#!/usr/bin/env python3
"""
DSSO MS2–MS3 Cross-link Pipeline (multithreaded, optimized)
----------------------------------------------------------
* Pairs the two best MS³ PSMs for each parent MS² scan
* Strict +54 / +86 stub complementarity ("true XL")
* Reporter-doublet QC: 31.99251 Da / z (+54 vs +86)
* Adds RT, DeltaCn, Conf %, gene names, cleaned protein lists
* Target–decoy link-level FDR
* ThreadPoolExecutor for speed (--threads N, default = CPU cores)
"""

import os, re, sys, argparse
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

# ── DSSO-specific constants (pre-compiled) ──────────────────────
STUB54_RE  = re.compile(r'\(54\.010')
STUB86_RE  = re.compile(r'\(86\.003')
DELTA_NEUT = 31.99251  # Pre-calculated
DECOY_RE   = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── Optimized .ms2 helpers ──────────────────────────────────────
def parse_ms3_child(path):
    """Return mapping {ms3_scan: (parent_scan, parent_charge)}."""
    mapping = {}
    scan = parent = charge = None
    
    with open(path, 'r', buffering=16384) as fh:  # Increased buffer
        for ln in fh:
            if ln[0] == 'S':  # Faster than startswith
                if scan and parent:
                    mapping[scan] = (parent, charge)
                parts = ln.split()
                scan = int(parts[1])
                parent = charge = None
            elif ln[0] == 'I' and 'PrecursorScan' in ln:
                parent = int(ln.split('\t')[2])
            elif ln[0] == 'Z':
                charge = int(float(ln.split()[1]))
    
    if scan and parent:
        mapping[scan] = (parent, charge)
    
    return os.path.basename(path)[:-4], mapping

@lru_cache(maxsize=2048)
def read_ms2_scan_cached(path, wanted):
    """Cached version to avoid repeated reads of same scan."""
    return read_ms2_scan(path, wanted)

def read_ms2_scan(path, wanted):
    """Return (charge, peaklist) for one parent MS² scan."""
    charge = None
    mz_list, int_list = [], []
    reading = False
    
    with open(path, 'r', buffering=16384) as fh:
        for ln in fh:
            if ln[0] == 'S':
                scan_num = int(ln.split()[1])
                if reading and scan_num != wanted:
                    break  # Exit early
                reading = (scan_num == wanted)
                if reading:
                    mz_list, int_list = [], []
                    charge = None
            elif reading:
                if ln[0] == 'Z':
                    charge = int(float(ln.split()[1]))
                elif ln and ln[0].isdigit():
                    parts = ln.split()
                    mz_list.append(float(parts[0]))
                    int_list.append(float(parts[1]))
                elif ln == '\n':
                    break
    
    peaks = list(zip(mz_list, int_list)) if mz_list else []
    return charge, peaks

def reporter_doublet(peaks, z, ppm=10, min_ratio=0.10):
    """Optimized search for 31.993 Da / z pair using numpy."""
    if not peaks or not z:
        return False, 0.0
    
    # Convert to numpy for vectorized operations
    peaks_arr = np.array(peaks)
    if len(peaks_arr) == 0:
        return False, 0.0
    
    # Sort by m/z
    sort_idx = np.argsort(peaks_arr[:, 0])
    mz_sorted = peaks_arr[sort_idx, 0]
    int_sorted = peaks_arr[sort_idx, 1]
    
    dmz = DELTA_NEUT / z
    topI = np.max(int_sorted)
    if topI == 0:
        topI = 1.0
    
    # Vectorized search
    for i in range(len(mz_sorted) - 1):
        tgt = mz_sorted[i] + dmz
        tol = ppm * tgt / 1e6
        
        # Find peaks within tolerance window
        candidates = np.where((mz_sorted[i+1:] >= tgt - tol) & 
                             (mz_sorted[i+1:] <= tgt + tol))[0]
        
        if len(candidates) > 0:
            j = i + 1 + candidates[0]
            ratio = min(int_sorted[i], int_sorted[j]) / topI
            if ratio >= min_ratio:
                return True, round(ratio, 3)
    
    return False, 0.0

@lru_cache(maxsize=10000)
def stub_type(seq):
    """Cached stub type detection."""
    if STUB54_RE.search(seq):
        return '54'
    if STUB86_RE.search(seq):
        return '86'
    return '?'

# ── Optimized main builder ──────────────────────────────────────
def build_links(pep_csv, spectra_dir, threads):
    # Read with optimized settings
    pep = pd.read_csv(pep_csv, low_memory=False)
    
    # Vectorized operations for metadata
    proteins_has_decoy = pep['Proteins'].fillna('').str.contains(DECOY_RE, regex=True)
    desc_has_decoy = pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE, regex=True)
    
    pep['is_decoy'] = proteins_has_decoy | desc_has_decoy
    pep['gene'] = pep['Protein Descriptions'].str.extract(GN_RE, expand=False)
    pep['conf'] = pep.get('Conf%', pd.Series(dtype=float))
    
    # Create sequence metadata more efficiently
    seq_cols = ['is_decoy', 'gene', 'conf', 'DeltCN', 'RetTime']
    seq_meta = pep.groupby('sequence').first()[seq_cols].rename(
        columns={'DeltCN': 'deltcn', 'RetTime': 'rt'}
    ).to_dict('index')
    
    # Filter MS3 data
    ms3 = pep[pep['File Name'].str.contains('_ms3', regex=False, na=False)].copy()
    
    if ms3.empty:
        sys.exit('No *_ms3.ms2 rows in pep_list.csv')
    
    # ── 1. Build MS3→parent maps (parallel) ─────────────────────
    child_files = ms3['File Name'].unique()
    parent_map = {}
    
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futures = {
            ex.submit(parse_ms3_child, os.path.join(spectra_dir, f'{c}.ms2')): c 
            for c in child_files
        }
        for future in as_completed(futures):
            base, mp = future.result()
            parent_map[base] = mp
    
    # Vectorized parent info assignment
    def get_parent_info(row):
        mapping = parent_map.get(row['File Name'], {})
        return mapping.get(int(row['Scan Num']), (None, None))
    
    parent_info = ms3.apply(get_parent_info, axis=1, result_type='expand')
    ms3[['parent_scan', 'parent_z']] = parent_info
    ms3 = ms3.dropna(subset=['parent_scan'])
    
    groups = list(ms3.groupby(['File Name', 'parent_scan']))
    
    # ── 2. Optimized worker function ────────────────────────────
    def process_group(item):
        (child, pscan), grp = item
        grp = grp.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        
        if len(grp) < 2:
            return None
        
        a, b = grp.iloc[0], grp.iloc[1]
        
        parent_path = os.path.join(spectra_dir, child.replace('_ms3', '') + '.ms2')
        if not os.path.isfile(parent_path):
            return None
        
        # Use cached read
        z, peaks = read_ms2_scan_cached(parent_path, int(pscan))
        dbl, ratio = reporter_doublet(peaks, z)
        
        # Use cached stub detection
        stA, stB = stub_type(a['sequence']), stub_type(b['sequence'])
        
        if stA in ('54', '86') and stB in ('54', '86'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('54', '86')) ^ (stB in ('54', '86')):
            note = 'Mono'
        else:
            note = 'Not sure'
        
        meta = lambda s, k: seq_meta.get(s, {}).get(k)
        clean = lambda p: p.strip('[] ')
        
        return {
            'MS2_file': os.path.basename(parent_path),
            'MS2_scan': int(pscan),
            'parent_charge': z,
            'MS3_file_A': f'{child}.ms2',
            'MS3_scan_A': int(a['Scan Num']),
            'sequence_A': a['sequence'],
            'XCorr_A': a['XCorr'],
            'DeltaCn_A': meta(a['sequence'], 'deltcn'),
            'RT_A': meta(a['sequence'], 'rt'),
            'Conf_A': meta(a['sequence'], 'conf'),
            'Proteins_A': clean(a['Proteins']),
            'Gene_A': meta(a['sequence'], 'gene'),
            'MS3_file_B': f'{child}.ms2',
            'MS3_scan_B': int(b['Scan Num']),
            'sequence_B': b['sequence'],
            'XCorr_B': b['XCorr'],
            'DeltaCn_B': meta(b['sequence'], 'deltcn'),
            'RT_B': meta(b['sequence'], 'rt'),
            'Conf_B': meta(b['sequence'], 'conf'),
            'Proteins_B': clean(b['Proteins']),
            'Gene_B': meta(b['sequence'], 'gene'),
            'doublet_present': dbl,
            'doublet_int_ratio': ratio,
            'Note': note
        }
    
    # ── 3. Parallel processing with better scheduling ────────────
    rows = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(process_group, group) for group in groups]
        for future in as_completed(futures):
            result = future.result()
            if result:
                rows.append(result)
    
    if not rows:
        sys.exit('No linkable parent scans found.')
    
    df = pd.DataFrame(rows)
    
    # ── 4. Vectorized decoy flags + FDR ─────────────────────────
    df['decoy_A'] = df['sequence_A'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']
    
    # Optimized sorting and FDR calculation
    df = df.sort_values('link_score', ascending=False, ignore_index=True)
    df['rank'] = np.arange(1, len(df) + 1)
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = np.round(2 * df['cum_decoys'] / df['cum_links'], 4)
    
    cols = [c for c in df.columns if c != 'Note'] + ['Note']
    return df[cols]

# ── Command-line driver ─────────────────────────────────────────
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='DSSO MS2-MS3 multithreaded XL pipeline')
    ap.add_argument('pep_list', help='IP2 pep_list.csv')
    ap.add_argument('spectra_dir', help='directory with *.ms2 and *_ms3.ms2')
    ap.add_argument('out_csv', help='output CSV file')
    ap.add_argument('--threads', type=int, default=os.cpu_count(),
                    help='worker threads (default = CPU cores)')
    args = ap.parse_args()
    
    print(f'Processing with {args.threads} threads...')
    tbl = build_links(args.pep_list, args.spectra_dir, args.threads)
    tbl.to_csv(args.out_csv, index=False)
    print(f'✓ {len(tbl)} links → {args.out_csv}')
    
    for lim in (0.01, 0.05):
        subset = tbl[tbl['FDR'] <= lim]
        if subset.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
        else:
            best = subset.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {best["link_score"]:.3f} '
                  f'({int(best["cum_links"])} links, {int(best["cum_decoys"])} decoys)')

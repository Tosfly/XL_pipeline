#!/usr/bin/env python3
"""
dsso_link_fdr_allin1_multi_optimized.py – Optimized DSSO MS2-MS3 XL analysis
"""
import os, re, sys, argparse
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from collections import defaultdict
import logging

# Pre-compiled regexes (compiled once)
STUB54_RE = re.compile(r'\(54\.010')
STUB86_RE = re.compile(r'\(86\.003')
DECOY_RE = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# Constants
DELTA_NEUT = 31.99251  # Pre-calculated
STUB_CACHE = {}  # Cache for stub detection

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# ── Optimized MS2 helpers ────────────────────────────────────────
@lru_cache(maxsize=1024)
def read_ms2_cached(path, scan):
    """Cached MS2 reading to avoid repeated file I/O"""
    return read_ms2(path, scan)

def map_ms3_to_parent_batch(path):
    """Optimized batch reading with single pass"""
    mapping = {}
    scan = parent = charge = None
    
    with open(path, 'r', buffering=8192) as fh:  # Larger buffer
        for ln in fh:
            if ln[0] == 'S':  # Faster than startswith
                if scan and parent:
                    mapping[scan] = (parent, charge)
                scan = int(ln.split()[1])
                parent = charge = None
            elif ln[0] == 'I' and 'PrecursorScan' in ln:
                parent = int(ln.split('\t')[2])
            elif ln[0] == 'Z':
                charge = int(float(ln.split()[1]))
    
    if scan and parent:
        mapping[scan] = (parent, charge)
    
    return os.path.basename(path)[:-4], mapping

def read_ms2(parent_path, scan):
    """Optimized MS2 reading with numpy arrays"""
    charge = None
    mz_list, int_list = [], []
    reading = False
    
    with open(parent_path, 'r', buffering=8192) as fh:
        for ln in fh:
            if ln[0] == 'S':
                if reading:  # Early exit
                    break
                reading = (int(ln.split()[1]) == scan)
                if reading:
                    mz_list, int_list = [], []
                    charge = None
            elif reading:
                if ln[0] == 'Z':
                    charge = int(float(ln.split()[1]))
                elif ln[0].isdigit():
                    parts = ln.split()
                    mz_list.append(float(parts[0]))
                    int_list.append(float(parts[1]))
                elif ln == '\n':
                    break
    
    if mz_list:
        # Use numpy for faster operations
        peaks = list(zip(mz_list, int_list))
        return charge, peaks
    return charge, []

def has_doublet_optimized(peaks, z, ppm=10, min_ratio=0.10):
    """Optimized doublet detection using numpy"""
    if not peaks or not z:
        return False, 0.0
    
    # Convert to numpy arrays for vectorized operations
    peaks_arr = np.array(peaks)
    if len(peaks_arr) == 0:
        return False, 0.0
    
    mz_arr = peaks_arr[:, 0]
    int_arr = peaks_arr[:, 1]
    
    # Sort by m/z
    sort_idx = np.argsort(mz_arr)
    mz_sorted = mz_arr[sort_idx]
    int_sorted = int_arr[sort_idx]
    
    dmz = DELTA_NEUT / z
    topI = np.max(int_sorted)
    if topI == 0:
        topI = 1.0
    
    # Vectorized search for doublets
    tgt_mz = mz_sorted + dmz
    tol = ppm * tgt_mz / 1e6
    
    for i in range(len(mz_sorted) - 1):
        # Find peaks within tolerance
        diff = np.abs(mz_sorted[i+1:] - tgt_mz[i])
        matches = np.where(diff <= tol[i])[0]
        
        if len(matches) > 0:
            # Get the first match
            j = i + 1 + matches[0]
            ratio = min(int_sorted[i], int_sorted[j]) / topI
            if ratio >= min_ratio:
                return True, round(ratio, 3)
    
    return False, 0.0

@lru_cache(maxsize=10000)
def stub_cached(seq):
    """Cached stub detection"""
    if STUB54_RE.search(seq):
        return '54'
    if STUB86_RE.search(seq):
        return '86'
    return '?'

# ── Optimized main build ─────────────────────────────────────────
def build_links_optimized(pep_csv, spectra_dir, n_threads):
    # Read peptide data
    pep = pd.read_csv(pep_csv, low_memory=False)
    
    # Vectorized operations for metadata
    pep['is_decoy'] = (pep['Proteins'].fillna('').str.contains(DECOY_RE, regex=True) |
                       pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE, regex=True))
    pep['gene'] = pep['Protein Descriptions'].str.extract(GN_RE, expand=False)
    pep['conf'] = pep.get('Conf%', pd.Series(dtype=float))
    
    # Create sequence metadata dict more efficiently
    seq_meta = pep.groupby('sequence').first()[
        ['is_decoy', 'gene', 'conf', 'DeltCN', 'RetTime']
    ].rename(columns={'DeltCN': 'deltcn', 'RetTime': 'rt'}).to_dict('index')
    
    # Filter MS3 data
    ms3_mask = pep['File Name'].str.contains('_ms3', regex=False)
    ms3 = pep[ms3_mask].copy()
    
    if ms3.empty:
        sys.exit('No *_ms3.ms2 rows.')
    
    # Parallel parent mapping
    child_files = ms3['File Name'].unique()
    parent_map = {}
    
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futures = {
            ex.submit(map_ms3_to_parent_batch, os.path.join(spectra_dir, f'{c}.ms2')): c
            for c in child_files
        }
        for future in as_completed(futures):
            child, mp = future.result()
            parent_map[child] = mp
    
    # Vectorized parent scan assignment
    def get_parent_info(row):
        mapping = parent_map.get(row['File Name'], {})
        scan_num = int(row['Scan Num'])
        return mapping.get(scan_num, (None, None))
    
    parent_info = ms3.apply(get_parent_info, axis=1, result_type='expand')
    ms3[['parent_scan', 'parent_z']] = parent_info
    ms3 = ms3.dropna(subset=['parent_scan'])
    
    # Group processing
    groups = list(ms3.groupby(['File Name', 'parent_scan']))
    
    def process_group(item):
        (child, pscan), g = item
        
        # Sort and deduplicate
        g = g.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        if len(g) < 2:
            return None
        
        a, b = g.iloc[0], g.iloc[1]
        parent = os.path.join(spectra_dir, child.replace('_ms3', '') + '.ms2')
        
        if not os.path.isfile(parent):
            return None
        
        # Use cached read
        z, peaks = read_ms2_cached(parent, int(pscan))
        dbl, ratio = has_doublet_optimized(peaks, z)
        
        # Use cached stub detection
        stA, stB = stub_cached(a['sequence']), stub_cached(b['sequence'])
        
        # Determine note
        if stA in ('54', '86') and stB in ('54', '86'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('54', '86')) ^ (stB in ('54', '86')):
            note = 'Mono'
        else:
            note = 'Not sure'
        
        # Build result dict
        get = lambda s, k: seq_meta.get(s, {}).get(k)
        clean = lambda p: p.strip('[] ')
        
        return {
            'MS2_file': os.path.basename(parent),
            'MS2_scan': int(pscan),
            'parent_charge': z,
            'MS3_file_A': f'{child}.ms2',
            'MS3_scan_A': int(a['Scan Num']),
            'sequence_A': a['sequence'],
            'XCorr_A': a['XCorr'],
            'DeltaCn_A': get(a['sequence'], 'deltcn'),
            'RT_A': get(a['sequence'], 'rt'),
            'Conf_A': get(a['sequence'], 'conf'),
            'Proteins_A': clean(a['Proteins']),
            'Gene_A': get(a['sequence'], 'gene'),
            'MS3_file_B': f'{child}.ms2',
            'MS3_scan_B': int(b['Scan Num']),
            'sequence_B': b['sequence'],
            'XCorr_B': b['XCorr'],
            'DeltaCn_B': get(b['sequence'], 'deltcn'),
            'RT_B': get(b['sequence'], 'rt'),
            'Conf_B': get(b['sequence'], 'conf'),
            'Proteins_B': clean(b['Proteins']),
            'Gene_B': get(b['sequence'], 'gene'),
            'doublet_present': dbl,
            'doublet_int_ratio': ratio,
            'Note': note
        }
    
    # Process groups in parallel
    rows = []
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futures = [ex.submit(process_group, group) for group in groups]
        for future in as_completed(futures):
            result = future.result()
            if result:
                rows.append(result)
    
    if not rows:
        sys.exit('No linkable groups.')
    
    # Build final dataframe
    df = pd.DataFrame(rows)
    
    # Vectorized operations for decoy detection and scoring
    df['decoy_A'] = df['sequence_A'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']
    
    # Sort and calculate FDR
    df = df.sort_values('link_score', ascending=False, ignore_index=True)
    df['rank'] = np.arange(1, len(df) + 1)
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = np.round(2 * df['cum_decoys'] / df['cum_links'], 4)
    
    # Reorder columns
    cols = [c for c in df.columns if c != 'Note'] + ['Note']
    return df[cols]

# ── Main entry point ─────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Optimized DSSO XL multithreaded pipeline'
    )
    parser.add_argument('pep_list', help='Input peptide CSV file')
    parser.add_argument('spectra_dir', help='Directory containing MS2 files')
    parser.add_argument('out_csv', help='Output CSV file')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads (default: 8)')
    
    args = parser.parse_args()
    
    # Process
    logging.info(f'Processing with {args.threads} threads...')
    table = build_links_optimized(args.pep_list, args.spectra_dir, args.threads)
    
    # Save results
    table.to_csv(args.out_csv, index=False)
    logging.info(f'✓ {len(table)} links → {args.out_csv}')
    
    # Report FDR thresholds
    for lim in (0.01, 0.05):
        sub = table[table['FDR'] <= lim]
        if sub.empty:
            logging.info(f'No links ≤ {lim:.0%} FDR.')
        else:
            top = sub.iloc[0]
            logging.info(
                f'{lim:.0%} FDR ⇒ link_score ≥ {top["link_score"]:.3f} '
                f'({int(top["cum_links"])} links, {int(top["cum_decoys"])} decoys)'
            )

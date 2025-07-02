#!/usr/bin/env python3
"""
DSSO MS2–MS3 XL pipeline (hybrid version - memory efficient)
Uses threading for I/O and limited multiprocessing for CPU-intensive work
"""
import os, re, sys, argparse
import pandas as pd, numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import lru_cache
import pickle

# ── DSSO stub masses -------------------------------------------------------
STUB_A_MASS = 54.01056                # thiol fragment
STUB_B_MASS = 86.00307                # alkene fragment
DELTA_NEUT  = STUB_B_MASS - STUB_A_MASS   # 31.99251 Da

STUB54_RE  = re.compile(r'\(54\.010')
STUB86_RE  = re.compile(r'\(86\.003')

DECOY_RE   = re.compile(r'(?:decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── Global cache for MS2 data (shared across threads) ---------------------
_MS2_CACHE = {}

def parse_ms3_child(path):
    """Parse single MS3 file"""
    m, scan, parent, charge = {}, None, None, None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                if scan and parent:
                    m[scan] = (parent, charge)
                scan, parent, charge = int(ln.split()[1]), None, None
            elif ln.startswith('I') and 'PrecursorScan' in ln:
                parent = int(ln.split('\t')[2])
            elif ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
        if scan and parent:
            m[scan] = (parent, charge)
    return os.path.basename(path)[:-4], m

@lru_cache(maxsize=10_000)
def read_ms2_scan_cached(path, wanted):
    """Read specific scan from MS2 file with caching"""
    charge, peaks, reading = None, [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                scan_num = int(ln.split()[1])
                if scan_num == wanted:
                    reading = True
                elif reading:
                    break  # We've passed our scan
                continue
            
            if not reading:
                continue
                
            if ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
            elif ln and ln[0].isdigit():
                parts = ln.split()
                if len(parts) >= 2:
                    mz, intensity = float(parts[0]), float(parts[1])
                    peaks.append((mz, intensity))
    
    return charge, peaks

def has_doublet_fast(peaks, z, ppm=10, min_ratio=0.10):
    """Fast doublet detection"""
    if not peaks or not z:
        return False, 0.0
    
    dmz = DELTA_NEUT / z
    peaks_sorted = sorted(peaks, key=lambda x: x[0])
    mzs = [p[0] for p in peaks_sorted]
    ints = [p[1] for p in peaks_sorted]
    top = max(ints) if ints else 1.0
    
    for i in range(len(mzs) - 1):
        target = mzs[i] + dmz
        tol = ppm * target / 1e6
        
        # Binary search for efficiency
        j = i + 1
        while j < len(mzs) and mzs[j] - target <= tol:
            if abs(mzs[j] - target) <= tol:
                ratio = min(ints[i], ints[j]) / top
                if ratio >= min_ratio:
                    return True, round(ratio, 3)
            j += 1
            
    return False, 0.0

def stub(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

def process_single_group(group_data, spectra_dir, meta):
    """Process single MS2-MS3 group (for threading)"""
    (child, pscan), g = group_data
    
    g = g.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
    if len(g) < 2:
        return None
        
    a, b = g.iloc[0], g.iloc[1]
    pfile_base = child.replace('_ms3', '') + '.ms2'
    pfile = os.path.join(spectra_dir, pfile_base)
    
    if not os.path.isfile(pfile):
        return None
        
    z, peaks = read_ms2_scan_cached(pfile, int(pscan))
    if z is None:
        return None
        
    dbl, ratio = has_doublet_fast(peaks, z)
    stA, stB = stub(a['sequence']), stub(b['sequence'])
    
    note = 'true XL' if {stA, stB} == {'54', '86'} else \
           ('Rare, may not be true' if stA == stB != '?' else
            'Mono' if '?' not in (stA, stB) else 'Not sure')
    
    get = lambda s, k: meta.get(s, {}).get(k)
    clean = lambda p: p.strip('[] ')
    
    return dict(
        MS2_file=os.path.basename(pfile),
        MS2_scan=int(pscan),
        parent_charge=z,
        MS3_file_A=f'{child}.ms2',
        MS3_scan_A=int(a['Scan Num']),
        sequence_A=a['sequence'],
        XCorr_A=a['XCorr'],
        DeltaCn_A=get(a['sequence'], 'deltcn'),
        RT_A=get(a['sequence'], 'rt'),
        Conf_A=get(a['sequence'], 'conf'),
        Proteins_A=clean(a['Proteins']),
        Gene_A=get(a['sequence'], 'gene'),
        MS3_file_B=f'{child}.ms2',
        MS3_scan_B=int(b['Scan Num']),
        sequence_B=b['sequence'],
        XCorr_B=b['XCorr'],
        DeltaCn_B=get(b['sequence'], 'deltcn'),
        RT_B=get(b['sequence'], 'rt'),
        Conf_B=get(b['sequence'], 'conf'),
        Proteins_B=clean(b['Proteins']),
        Gene_B=get(b['sequence'], 'gene'),
        doublet_present=dbl,
        doublet_int_ratio=ratio,
        Note=note
    )

def build_links(pep_csv, spectra, threads):
    """Build cross-links with hybrid threading/multiprocessing approach"""
    print(f"Loading peptide data from {pep_csv}...")
    pep = pd.read_csv(pep_csv, low_memory=False)
    pep['is_decoy'] = pep['Proteins'].fillna('').str.contains(DECOY_RE) | \
                      pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE)
    pep['gene'] = pep['Protein Descriptions'].str.extract(GN_RE, expand=False)
    pep['conf'] = pep.get('Conf%', pd.Series(dtype=float))
    
    meta = pep.groupby('sequence').first()\
             .rename(columns={'DeltCN': 'deltcn', 'RetTime': 'rt'})\
             .to_dict('index')
    
    ms3 = pep[pep['File Name'].str.contains('_ms3', na=False)].copy()
    if ms3.empty:
        sys.exit('No *_ms3 rows.')
    
    print(f"Found {len(ms3)} MS3 entries")
    
    # Parse MS3 files using ThreadPoolExecutor (I/O bound)
    unique_ms3_files = ms3['File Name'].unique()
    print(f"Parsing {len(unique_ms3_files)} MS3 files...")
    
    maps = {}
    with ThreadPoolExecutor(max_workers=min(threads, 16)) as executor:
        futures = {
            executor.submit(parse_ms3_child, os.path.join(spectra, f'{c}.ms2')): c
            for c in unique_ms3_files
        }
        
        for future in as_completed(futures):
            try:
                child, m = future.result()
                maps[child] = m
            except Exception as e:
                print(f"Error parsing {futures[future]}: {e}")
    
    ms3[['parent_scan', 'parent_z']] = ms3.apply(
        lambda r: pd.Series(maps.get(r['File Name'], {}).get(int(r['Scan Num']), (None, None))),
        axis=1
    )
    ms3.dropna(subset=['parent_scan'], inplace=True)
    
    groups = list(ms3.groupby(['File Name', 'parent_scan']))
    print(f"Processing {len(groups)} MS2-MS3 pairs...")
    
    # Process groups using ThreadPoolExecutor (mostly I/O bound with some CPU)
    all_results = []
    batch_size = 100  # Process in smaller batches to show progress
    
    for i in range(0, len(groups), batch_size):
        batch = groups[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(process_single_group, group, spectra, meta)
                for group in batch
            ]
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result:
                        all_results.append(result)
                except Exception as e:
                    print(f"Error processing group: {e}")
        
        # Progress update
        processed = min(i + batch_size, len(groups))
        print(f"Processed {processed}/{len(groups)} groups ({processed/len(groups)*100:.1f}%)")
    
    if not all_results:
        sys.exit('No links found.')
    
    print(f"Found {len(all_results)} candidate links")
    
    # Build final dataframe
    df = pd.DataFrame(all_results)
    df['decoy_A'] = df['sequence_A'].map(lambda s: meta.get(s, {}).get('is_decoy', False))
    df['decoy_B'] = df['sequence_B'].map(lambda s: meta.get(s, {}).get('is_decoy', False))
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']
    
    # Sort and calculate FDR
    df.sort_values('link_score', ascending=False, inplace=True, ignore_index=True)
    df['rank'] = np.arange(1, len(df) + 1)
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = np.round(2 * df['cum_decoys'] / df['cum_links'], 4)
    
    # Reorder columns with Note at the end
    cols = [c for c in df.columns if c != 'Note'] + ['Note']
    return df[cols]

# ── CLI --------------------------------------------------------------------
if __name__ == '__main__':
    ar = argparse.ArgumentParser(description='DSSO hybrid XL pipeline')
    ar.add_argument('pep_list')
    ar.add_argument('spectra_dir')
    ar.add_argument('out_csv')
    ar.add_argument('--threads', type=int, default=min(8, os.cpu_count()))
    args = ar.parse_args()
    
    print(f'Processing with {args.threads} threads...')
    print(f'Memory-efficient hybrid mode')
    
    out = build_links(args.pep_list, args.spectra_dir, args.threads)
    out.to_csv(args.out_csv, index=False)
    print(f'✓ {len(out)} links → {args.out_csv}')
    
    # Print FDR summary
    for lim in (0.01, 0.05):
        sub = out[out['FDR'] <= lim]
        if sub.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
        else:
            best = sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {best["link_score"]:.3f} '
                  f'({int(best["cum_links"])} links, {int(best["cum_decoys"])} decoys)')

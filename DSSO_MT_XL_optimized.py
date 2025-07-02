#!/usr/bin/env python3
"""
DSSO MS2–MS3 XL pipeline (optimized version with better I/O and multiprocessing)
"""
import os, re, sys, argparse
import pandas as pd, numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import lru_cache
from collections import defaultdict
import pickle

# ── DSSO stub masses -------------------------------------------------------
STUB_A_MASS = 54.01056                # thiol fragment
STUB_B_MASS = 86.00307                # alkene fragment
DELTA_NEUT  = STUB_B_MASS - STUB_A_MASS   # 31.99251 Da

STUB54_RE  = re.compile(r'\(54\.010')
STUB86_RE  = re.compile(r'\(86\.003')

DECOY_RE   = re.compile(r'(?:decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── Optimized MS2 reading with batch loading ------------------------------
class MS2Cache:
    def __init__(self, spectra_dir):
        self.spectra_dir = spectra_dir
        self.cache = {}
        
    def preload_files(self, filenames):
        """Preload multiple MS2 files into memory"""
        for fname in filenames:
            if fname not in self.cache:
                path = os.path.join(self.spectra_dir, fname)
                if os.path.exists(path):
                    self.cache[fname] = self._parse_ms2_file(path)
    
    def _parse_ms2_file(self, path):
        """Parse entire MS2 file into memory-efficient structure"""
        scans = {}
        current_scan = None
        charge = None
        peaks = []
        
        with open(path) as fh:
            for ln in fh:
                if ln.startswith('S'):
                    if current_scan is not None:
                        scans[current_scan] = (charge, peaks)
                    current_scan = int(ln.split()[1])
                    charge = None
                    peaks = []
                elif ln.startswith('Z') and current_scan is not None:
                    charge = int(float(ln.split()[1]))
                elif ln and ln[0].isdigit() and current_scan is not None:
                    parts = ln.split()
                    if len(parts) >= 2:
                        mz, intensity = float(parts[0]), float(parts[1])
                        peaks.append((mz, intensity))
                        
        if current_scan is not None:
            scans[current_scan] = (charge, peaks)
            
        return scans
    
    def get_scan(self, filename, scan_num):
        """Get specific scan from cache"""
        if filename in self.cache and scan_num in self.cache[filename]:
            return self.cache[filename][scan_num]
        return None, []

def parse_ms3_child_batch(paths):
    """Parse multiple MS3 files in parallel"""
    results = {}
    for path in paths:
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
        results[os.path.basename(path)[:-4]] = m
    return results

def has_doublet_vectorized(peaks, z, ppm=10, min_ratio=0.10):
    """Vectorized doublet detection for better performance"""
    if not peaks or not z:
        return False, 0.0
    
    peaks = np.array(peaks)
    if len(peaks) == 0:
        return False, 0.0
        
    mzs, ints = peaks[:, 0], peaks[:, 1]
    dmz = DELTA_NEUT / z
    top = np.max(ints) if len(ints) > 0 else 1.0
    
    # Vectorized distance calculation
    for i in range(len(mzs) - 1):
        target = mzs[i] + dmz
        tol = ppm * target / 1e6
        
        # Find peaks within tolerance
        mask = (np.abs(mzs[i+1:] - target) <= tol)
        if np.any(mask):
            matched_ints = ints[i+1:][mask]
            ratio = min(ints[i], np.max(matched_ints)) / top
            if ratio >= min_ratio:
                return True, round(ratio, 3)
    
    return False, 0.0

def stub(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

# ── Process chunk of work items --------------------------------------------
def process_work_chunk(args):
    """Process a chunk of work items (for multiprocessing)"""
    items, spectra_dir, meta_pickle = args
    
    # Deserialize metadata
    with open(meta_pickle, 'rb') as f:
        meta = pickle.load(f)
    
    # Initialize MS2 cache for this process
    ms2_cache = MS2Cache(spectra_dir)
    
    # Collect all needed files
    needed_files = set()
    for (child, pscan), g in items:
        pfile = child.replace('_ms3', '') + '.ms2'
        needed_files.add(pfile)
    
    # Preload files
    ms2_cache.preload_files(needed_files)
    
    results = []
    for (child, pscan), g in items:
        g = g.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        if len(g) < 2:
            continue
            
        a, b = g.iloc[0], g.iloc[1]
        pfile = child.replace('_ms3', '') + '.ms2'
        
        z, peaks = ms2_cache.get_scan(pfile, int(pscan))
        if z is None:
            continue
            
        dbl, ratio = has_doublet_vectorized(peaks, z)
        stA, stB = stub(a['sequence']), stub(b['sequence'])
        
        note = 'true XL' if {stA, stB} == {'54', '86'} else \
               ('Rare, may not be true' if stA == stB != '?' else
                'Mono' if '?' not in (stA, stB) else 'Not sure')
        
        get = lambda s, k: meta.get(s, {}).get(k)
        clean = lambda p: p.strip('[] ')
        
        results.append(dict(
            MS2_file=os.path.basename(pfile[:-4] + '.ms2'),
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
        ))
    
    return results

# ── Build links with optimized I/O and multiprocessing --------------------
def build_links(pep_csv, spectra, threads):
    print(f"Loading peptide data from {pep_csv}...")
    pep = pd.read_csv(pep_csv, low_memory=False)
    pep['is_decoy'] = pep['Proteins'].fillna('').str.contains(DECOY_RE) | \
                      pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE)
    pep['gene'] = pep['Protein Descriptions'].str.extract(GN_RE, expand=False)
    pep['conf'] = pep.get('Conf%', pd.Series(dtype=float))
    
    meta = pep.groupby('sequence').first()\
             .rename(columns={'DeltCN': 'deltcn', 'RetTime': 'rt'})\
             .to_dict('index')
    
    # Serialize metadata for multiprocessing
    meta_pickle = 'meta_temp.pkl'
    with open(meta_pickle, 'wb') as f:
        pickle.dump(meta, f)
    
    ms3 = pep[pep['File Name'].str.contains('_ms3', na=False)].copy()
    if ms3.empty:
        os.remove(meta_pickle)
        sys.exit('No *_ms3 rows.')
    
    print(f"Found {len(ms3)} MS3 entries")
    
    # Parse MS3 files in batches
    unique_ms3_files = ms3['File Name'].unique()
    ms3_paths = [os.path.join(spectra, f'{c}.ms2') for c in unique_ms3_files]
    
    print(f"Parsing {len(unique_ms3_files)} MS3 files...")
    
    # Use ThreadPoolExecutor for I/O-bound MS3 parsing
    chunk_size = max(1, len(ms3_paths) // threads)
    chunks = [ms3_paths[i:i + chunk_size] for i in range(0, len(ms3_paths), chunk_size)]
    
    maps = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(parse_ms3_child_batch, chunk) for chunk in chunks]
        for future in as_completed(futures):
            maps.update(future.result())
    
    ms3[['parent_scan', 'parent_z']] = ms3.apply(
        lambda r: pd.Series(maps.get(r['File Name'], {}).get(int(r['Scan Num']), (None, None))),
        axis=1
    )
    ms3.dropna(subset=['parent_scan'], inplace=True)
    
    groups = list(ms3.groupby(['File Name', 'parent_scan']))
    print(f"Processing {len(groups)} MS2-MS3 pairs...")
    
    # Split work into chunks for multiprocessing
    chunk_size = max(1, len(groups) // threads)
    work_chunks = [groups[i:i + chunk_size] for i in range(0, len(groups), chunk_size)]
    
    # Process chunks in parallel using ProcessPoolExecutor
    all_results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        chunk_args = [(chunk, spectra, meta_pickle) for chunk in work_chunks]
        futures = [executor.submit(process_work_chunk, args) for args in chunk_args]
        
        for future in as_completed(futures):
            all_results.extend(future.result())
    
    # Clean up temporary file
    os.remove(meta_pickle)
    
    if not all_results:
        sys.exit('No links found.')
    
    print(f"Found {len(all_results)} candidate links")
    
    # Build final dataframe
    df = pd.DataFrame(all_results)
    df['decoy_A'] = df['sequence_A'].map(lambda s: meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: meta[s]['is_decoy'])
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
    ar = argparse.ArgumentParser(description='DSSO optimized XL pipeline')
    ar.add_argument('pep_list')
    ar.add_argument('spectra_dir')
    ar.add_argument('out_csv')
    ar.add_argument('--threads', type=int, default=min(8, os.cpu_count()))
    args = ar.parse_args()
    
    print(f'Processing with {args.threads} processes...')
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

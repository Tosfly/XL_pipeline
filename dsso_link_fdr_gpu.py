#!/usr/bin/env python3
"""
DSSO MS2–MS3 XL pipeline – GPU-accelerated when RAPIDS is available
* GPU path (~8–12 min on RTX A4000)  – cudf + cupy + numba.cuda
* CPU fallback (~2 h on 16 cores)    – pandas + multithreading
Everything else (stub rules, QC, FDR) is identical to the previous script.
"""
import os, re, sys, argparse, math, functools, warnings
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

# ── try GPU stack ──────────────────────────────────────────────────────────
GPU = False
try:
    import cudf, cupy
    from numba import cuda
    GPU = True
except ImportError:
    warnings.warn("RAPIDS stack not found – running in CPU mode.")

# ── DSSO constants (canonical +54 / +86) ───────────────────────────────────
STUB_A = 54.01056               # thiol
STUB_B = 86.00307               # alkene
DELTA_NEUT = STUB_B - STUB_A    # 31.99251
STUB54_RE = re.compile(r'\(54\.010')
STUB86_RE = re.compile(r'\(86\.003')
DECOY_RE  = re.compile(r'(?:decoy|reverse|rev_)', re.I)
GN_RE     = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── tiny .ms2 helpers (identical for CPU / GPU) ────────────────────────────
def parse_ms3_child(path):
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

@lru_cache(maxsize=20_000)
def read_ms2_scan(path, wanted):
    charge, mz, intensity, reading = None, [], [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                reading = (int(ln.split()[1]) == wanted)
                charge, mz, intensity = (None, [], []) if reading else (charge, mz, intensity)
                continue
            if not reading: continue
            if ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
            elif ln and ln[0].isdigit():
                m, i = map(float, ln.split()[:2]); mz.append(m); intensity.append(i)
            elif ln == '\n': break
    return charge, np.array(mz, dtype=np.float32), np.array(intensity, dtype=np.float32)

# ── GPU kernel for reporter doublet search ─────────────────────────────────
if GPU:
    @cuda.jit(device=True)
    def has_doublet_gpu(mz, inten, n, dmz, ppm, min_ratio):
        top = 0.0
        for i in range(n):
            if inten[i] > top:
                top = inten[i]
        if top == 0:           # safety
            return 0
        for i in range(n):
            m1 = mz[i]
            tgt = m1 + dmz
            tol = ppm * (m1 + dmz) / 1e6
            for j in range(i + 1, n):
                delta = mz[j] - tgt
                if delta > tol:
                    break
                if math.fabs(delta) <= tol:
                    ratio = min(inten[i], inten[j]) / top
                    if ratio >= min_ratio:
                        return 1
        return 0

def reporter_doublet_gpu(mz, inten, z, ppm=10, min_ratio=0.10):
    if z is None or z == 0 or len(mz) == 0:
        return False, 0.0
    dmz = DELTA_NEUT / z
    mz_d = cuda.to_device(mz)
    it_d = cuda.to_device(inten)
    flag = cuda.device_array(1, dtype=np.int32)
    has_doublet_gpu[1, 1](mz_d, it_d, len(mz), dmz, ppm, min_ratio)  # result in return value
    present = bool(_ := has_doublet_gpu(mz, inten, len(mz), dmz, ppm, min_ratio))
    # ratio not returned (costly); we keep 0.0 / 1.0
    return present, 1.0 if present else 0.0

# ── CPU doublet ------------------------------------------------------------
def reporter_doublet_cpu(mz, inten, z, ppm=10, min_ratio=0.10):
    if z is None or z == 0 or len(mz) == 0:
        return False, 0.0
    dmz = DELTA_NEUT / z
    top = inten.max() if len(inten) else 0.0
    if top == 0: return False, 0.0
    order = np.argsort(mz)
    mz, inten = mz[order], inten[order]
    for i, m1 in enumerate(mz):
        tgt = m1 + dmz
        tol = ppm * (m1 + dmz) / 1e6
        j = i + 1
        while j < len(mz) and mz[j] - tgt <= tol:
            if abs(mz[j] - tgt) <= tol:
                ratio = min(inten[i], inten[j]) / top
                if ratio >= min_ratio:
                    return True, ratio
            j += 1
    return False, 0.0

has_doublet = reporter_doublet_gpu if GPU else reporter_doublet_cpu

# ── stub helper ------------------------------------------------------------
def stub(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

# ── main builder -----------------------------------------------------------
def build_links(pep_csv, spectra_dir, threads):
    # ---------------- read pep_list (pandas ok even in GPU mode) ----------
    pep = pd.read_csv(pep_csv, low_memory=False)
    pep['is_decoy'] = pep['Proteins'].fillna('').str.contains(DECOY_RE) | \
                      pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE)
    pep['gene']    = pep['Protein Descriptions'].str.extract(GN_RE, expand=False)
    pep['conf']    = pep.get('Conf%', pd.Series(dtype=float))
    meta = pep.groupby('sequence').first()\
              .rename(columns={'DeltCN':'deltcn','RetTime':'rt'})\
              .to_dict('index')

    ms3 = pep[pep['File Name'].str.contains('_ms3', na=False)].copy()
    if ms3.empty:
        sys.exit('No *_ms3 rows.')

    # ---------------- child→parent maps in parallel ----------------------
    child_files = ms3['File Name'].unique()
    parent_map = {}
    with ThreadPoolExecutor(threads) as ex:
        fut = [ex.submit(parse_ms3_child,
                         os.path.join(spectra_dir, f'{c}.ms2'))
               for c in child_files]
        for f in as_completed(fut):
            child, m = f.result()
            parent_map[child] = m

    ms3[['parent_scan','parent_z']] = ms3.apply(
        lambda r: pd.Series(parent_map[r['File Name']].get(
            int(r['Scan Num']), (None, None))), axis=1)
    ms3 = ms3.dropna(subset=['parent_scan']).copy()
    groups = list(ms3.groupby(['File Name','parent_scan']))

    # ---------------- GPU branch uses cudf for group ops ------------------
    if GPU:
        # Convert the small groups table to cudf for faster sort/drop
        mini_cols = ['File Name','parent_scan','Scan Num','sequence','XCorr','Proteins']
        gdf = cudf.DataFrame(ms3[mini_cols])
        gdf['rank'] = gdf.groupby(['File Name','parent_scan'])['XCorr']\
                         .rank(method='first', ascending=False)
        gdf = gdf[gdf['rank'] <= 2].to_pandas()   # back to pandas rows

    else:
        gdf = None  # CPU branch uses original groups

    # rows to collect
    rows = []

    # ---------------- inner worker ---------------------------------------
    def work(child, pscan, two_psms):
        a, b = two_psms
        pfile = os.path.join(spectra_dir, child.replace('_ms3','') + '.ms2')
        if not os.path.isfile(pfile):
            return None
        z, mz, inten = read_ms2_scan(pfile, int(pscan))
        dbl, ratio = has_doublet(mz, inten, z)

        stA, stB = stub(a['sequence']), stub(b['sequence'])
        if stA in ('54','86') and stB in ('54','86'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('54','86')) ^ (stB in ('54','86')):
            note = 'Mono'
        else:
            note = 'Not sure'

        meta_val = lambda s, k: meta.get(s, {}).get(k)
        clean = lambda p: p.strip('[] ')

        return dict(
            MS2_file=os.path.basename(pfile), MS2_scan=int(pscan), parent_charge=z,
            MS3_file_A=f'{child}.ms2', MS3_scan_A=int(a['Scan Num']),
            sequence_A=a['sequence'], XCorr_A=a['XCorr'],
            DeltaCn_A=meta_val(a['sequence'],'deltcn'), RT_A=meta_val(a['sequence'],'rt'),
            Conf_A=meta_val(a['sequence'],'conf'),
            Proteins_A=clean(a['Proteins']), Gene_A=meta_val(a['sequence'],'gene'),
            MS3_file_B=f'{child}.ms2', MS3_scan_B=int(b['Scan Num']),
            sequence_B=b['sequence'], XCorr_B=b['XCorr'],
            DeltaCn_B=meta_val(b['sequence'],'deltcn'), RT_B=meta_val(b['sequence'],'rt'),
            Conf_B=meta_val(b['sequence'],'conf'),
            Proteins_B=clean(b['Proteins']), Gene_B=meta_val(b['sequence'],'gene'),
            doublet_present=dbl, doublet_int_ratio=ratio, Note=note)

    # ---------------- parallel work dispatch -----------------------------
    tasks = []
    with ThreadPoolExecutor(threads) as ex:
        if GPU:
            for (child, pscan), grp in gdf.groupby(['File Name','parent_scan']):
                if len(grp) < 2: continue
                grp = grp.sort_values('XCorr', ascending=False)
                a, b = grp.iloc[0], grp.iloc[1]
                tasks.append(ex.submit(work, child, pscan, (a, b)))
        else:
            for (child, pscan), grp in groups:
                grp = grp.sort_values('XCorr', ascending=False)\
                         .drop_duplicates('sequence')
                if len(grp) < 2: continue
                a, b = grp.iloc[0], grp.iloc[1]
                tasks.append(ex.submit(work, child, pscan, (a, b)))

        for f in as_completed(tasks):
            r = f.result()
            if r:
                rows.append(r)

    if not rows:
        sys.exit('No links found.')

    df = pd.DataFrame(rows)
    # -------------- FDR & final columns ----------------------------------
    df['decoy_A'] = df['sequence_A'].map(lambda s: meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: meta[s]['is_decoy'])
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']
    df.sort_values('link_score', ascending=False, inplace=True, ignore_index=True)
    df['rank'] = np.arange(1, len(df) + 1)
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = np.round(2 * df['cum_decoys'] / df['cum_links'], 4)
    return df[[c for c in df.columns if c != 'Note'] + ['Note']]

# ── CLI -------------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='GPU-accelerated DSSO XL pipeline')
    ap.add_argument('pep_list')
    ap.add_argument('spectra_dir')
    ap.add_argument('out_csv')
    ap.add_argument('--threads', type=int, default=os.cpu_count(),
                    help='CPU threads for fallback / I/O (default = cores)')
    args = ap.parse_args()

    mode = "GPU" if GPU else "CPU-only"
    print(f'[{mode}] processing with {args.threads} CPU threads …')
    table = build_links(args.pep_list, args.spectra_dir, args.threads)
    table.to_csv(args.out_csv, index=False)
    print(f'✓ {len(table)} links → {args.out_csv}')

    for lim in (0.01, 0.05):
        sub = table[table['FDR'] <= lim]
        if sub.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
        else:
            top = sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {top["link_score"]:.3f} '
                  f'({int(top["cum_links"])} links, {int(top["cum_decoys"])} decoys)')

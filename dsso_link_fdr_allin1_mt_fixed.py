#!/usr/bin/env python3
"""
DSSO MS2–MS3 Cross-link Pipeline (multithreaded, fixed)
-------------------------------------------------------
* Pairs the two best MS³ PSMs for each parent MS² scan
* Strict +54 / +86 stub complementarity (“true XL”)
* Reporter-doublet QC: 31.99251 Da / z (+54 vs +86)
* Adds RT, DeltaCn, Conf %, gene names, cleaned protein lists
* Target–decoy link-level FDR
* ThreadPoolExecutor for speed (--threads N, default = CPU cores)
"""

import os, re, sys, argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── DSSO-specific constants ─────────────────────────────────────
STUB54_RE  = re.compile(r'\(54\.010')
STUB86_RE  = re.compile(r'\(86\.003')
DELTA_NEUT = 86.00307 - 54.01056        # 31.99251 Da
DECOY_RE   = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── tiny .ms2 helpers ───────────────────────────────────────────
def parse_ms3_child(path):
    """Return mapping {ms3_scan: (parent_scan, parent_charge)}."""
    mapping, scan, parent, charge = {}, None, None, None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                if scan and parent:
                    mapping[scan] = (parent, charge)
                scan, parent, charge = int(ln.split()[1]), None, None
            elif ln.startswith('I') and 'PrecursorScan' in ln:
                parent = int(ln.split('\t')[2])
            elif ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
        if scan and parent:
            mapping[scan] = (parent, charge)
    return os.path.basename(path)[:-4], mapping   # strip ".ms2"

def read_ms2_scan(path, wanted):
    """Return (charge, peaklist) for one parent MS² scan."""
    charge, peaks, reading = None, [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                reading = (int(ln.split()[1]) == wanted)
                peaks, charge = ([], None) if reading else (peaks, charge)
                continue
            if not reading:
                continue
            if ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
            elif ln and ln[0].isdigit():
                mz, inten = map(float, ln.split()[:2])
                peaks.append((mz, inten))
            elif ln == '\n':
                break
    return charge, peaks

def reporter_doublet(peaks, z, ppm=10, min_ratio=0.10):
    """Search for 31.993 Da / z pair; return (present?, ratio)."""
    if not peaks or not z:
        return False, 0.0
    peaks = sorted(peaks)
    dmz   = DELTA_NEUT / z
    topI  = max(i for _, i in peaks) or 1.0
    for i, (mz1, i1) in enumerate(peaks):
        tgt, tol = mz1 + dmz, ppm * (mz1 + dmz) / 1e6
        for mz2, i2 in peaks[i + 1:]:
            if mz2 - tgt > tol:
                break
            if abs(mz2 - tgt) <= tol:
                return (min(i1, i2) / topI) >= min_ratio, round(min(i1, i2) / topI, 3)
    return False, 0.0

def stub_type(seq):
    if STUB54_RE.search(seq):
        return '54'
    if STUB86_RE.search(seq):
        return '86'
    return '?'

# ── main builder ────────────────────────────────────────────────
def build_links(pep_csv, spectra_dir, threads):
    pep = pd.read_csv(pep_csv)

    # consolidate peptide-level metadata
    seq_meta = (
        pep.assign(
            is_decoy = pep['Proteins'].fillna('').str.contains(DECOY_RE) |
                       pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE),
            gene     = pep['Protein Descriptions'].str.extract(GN_RE, expand=False),
            conf     = pep.get('Conf%', None),
            deltcn   = pep['DeltCN'],
            rt       = pep['RetTime'])
        .groupby('sequence').first().to_dict('index'))

    ms3 = pep[pep['File Name'].str.contains('_ms3', regex=False)]
    if ms3.empty:
        sys.exit('No *_ms3.ms2 rows in pep_list.csv')

    # ── 1. build MS3→parent maps (parallel) ─────────────────────
    child_files = ms3['File Name'].unique()
    parent_map = {}
    with ThreadPoolExecutor(max_workers=threads) as ex:
        fut = {ex.submit(parse_ms3_child,
                         os.path.join(spectra_dir, f'{c}.ms2')): c for c in child_files}
        for f in as_completed(fut):
            base, mp = f.result()
            parent_map[base] = mp

    ms3[['parent_scan', 'parent_z']] = ms3.apply(
        lambda r: pd.Series(parent_map[r['File Name']].get(int(r['Scan Num']), (None, None))),
        axis=1)
    ms3.dropna(subset=['parent_scan'], inplace=True)

    groups = list(ms3.groupby(['File Name', 'parent_scan']))

    # ── 2. worker function for each (child,parent) group ────────
    def process_group(item):
        (child, pscan), grp = item
        grp = grp.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        if len(grp) < 2:
            return None
        a, b = grp.iloc[0], grp.iloc[1]

        parent_path = os.path.join(spectra_dir, child.replace('_ms3', '') + '.ms2')
        if not os.path.isfile(parent_path):
            return None
        z, peaks = read_ms2_scan(parent_path, int(pscan))
        dbl, ratio = reporter_doublet(peaks, z)

        stA, stB = stub_type(a['sequence']), stub_type(b['sequence'])
        if stA in ('54', '86') and stB in ('54', '86'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('54', '86')) ^ (stB in ('54', '86')):
            note = 'Mono'
        else:
            note = 'Not sure'

        meta = lambda s, k: seq_meta.get(s, {}).get(k)
        clean = lambda p: p.strip('[] ')

        return dict(
            MS2_file=os.path.basename(parent_path), MS2_scan=int(pscan), parent_charge=z,
            MS3_file_A=f'{child}.ms2', MS3_scan_A=int(a['Scan Num']),
            sequence_A=a['sequence'], XCorr_A=a['XCorr'],
            DeltaCn_A=meta(a['sequence'], 'deltcn'), RT_A=meta(a['sequence'], 'rt'),
            Conf_A=meta(a['sequence'], 'conf'), Proteins_A=clean(a['Proteins']),
            Gene_A=meta(a['sequence'], 'gene'),
            MS3_file_B=f'{child}.ms2', MS3_scan_B=int(b['Scan Num']),
            sequence_B=b['sequence'], XCorr_B=b['XCorr'],
            DeltaCn_B=meta(b['sequence'], 'deltcn'), RT_B=meta(b['sequence'], 'rt'),
            Conf_B=meta(b['sequence'], 'conf'), Proteins_B=clean(b['Proteins']),
            Gene_B=meta(b['sequence'], 'gene'),
            doublet_present=dbl, doublet_int_ratio=ratio, Note=note)

    # ── 3. parallel processing of groups ────────────────────────
    rows = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for row in ex.map(process_group, groups):
            if row:
                rows.append(row)

    if not rows:
        sys.exit('No linkable parent scans found.')

    df = pd.DataFrame(rows)

    # ── 4. decoy flags + FDR ────────────────────────────────────
    df['decoy_A'] = df['sequence_A'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: seq_meta[s]['is_decoy'])
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']

    df.sort_values('link_score', ascending=False, inplace=True, ignore_index=True)
    df['rank'] = df.index + 1
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = (2 * df['cum_decoys'] / df['cum_links']).round(4)

    cols = [c for c in df.columns if c != 'Note'] + ['Note']
    return df[cols]

# ── command-line driver ─────────────────────────────────────────
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='DSSO MS2-MS3 multithreaded XL pipeline')
    ap.add_argument('pep_list', help='IP2 pep_list.csv')
    ap.add_argument('spectra_dir', help='directory with *.ms2 and *_ms3.ms2')
    ap.add_argument('out_csv', help='output CSV file')
    ap.add_argument('--threads', type=int, default=os.cpu_count(),
                    help='worker threads (default = CPU cores)')
    args = ap.parse_args()

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

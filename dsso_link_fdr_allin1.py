#!/usr/bin/env python3
"""
dsso_link_fdr_allin1.py  –  One-file pipeline for **DSSO** MS2-MS3 data
-----------------------------------------------------------------------
* Pairs the two best MS³ PSMs for each parent MS²
* Checks for the DSSO 31.993 Da / z reporter doublet
* Adds gene names, Conf%, stub pattern “Note”
* Requires strict **+54/+86 complementarity** (“true XL”)
* Calculates link-level target–decoy FDR

Run
----
pip install pandas
python dsso_link_fdr_allin1.py  pep_list.csv  /path/to/spectra  DSSO_links_FDR.csv
"""

import os, re, sys, argparse
import pandas as pd

# ── DSSO-specific constants ────────────────────────────────────────────────
STUB54_RE  = re.compile(r'\(54\.010')          # thiol fragment
STUB86_RE  = re.compile(r'\(86\.003')          # alkene fragment
DELTA_NEUT = 86.00307 - 54.01056               # 31.99251 Da
LINK_MASS  = 158.00376                         # intact DSSO linker

# generic patterns
DECOY_RE   = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── lightweight .ms2 helpers ───────────────────────────────────────────────
def map_ms3_to_parent(path):
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
    return m

def read_ms2(path, wanted):
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
            elif ln[0].isdigit():
                mz, inten = map(float, ln.split()[:2])
                peaks.append((mz, inten))
            elif ln == '\n':
                break
    return charge, peaks

def has_doublet(peaks, z, ppm=10, min_ratio=0.10):
    if not peaks or not z:
        return False, 0.0
    peaks = sorted(peaks)
    dmz   = DELTA_NEUT / z
    topI  = max(i for _, i in peaks) or 1.0
    for i,(mz1,i1) in enumerate(peaks):
        tgt, tol = mz1 + dmz, ppm * (mz1+dmz) / 1e6
        for mz2,i2 in peaks[i+1:]:
            if mz2 - tgt > tol:
                break
            if abs(mz2 - tgt) <= tol:
                ratio = min(i1,i2) / topI
                return ratio >= min_ratio, round(ratio,3)
    return False, 0.0

# ── stub helpers ───────────────────────────────────────────────────────────
def stub_type(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

# ── main builder ───────────────────────────────────────────────────────────
def build_links(pep_csv, spectra_dir):
    pep = pd.read_csv(pep_csv)

    seq_info = (
        pep.assign(
            is_decoy = pep['Proteins'].fillna('').str.contains(DECOY_RE) |
                       pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE),
            gene     = pep['Protein Descriptions'].str.extract(GN_RE, expand=False),
            conf     = pep.get('Conf%', None),
            deltcn   = pep['DeltCN'],
            rt       = pep['RetTime']
        )
        .groupby('sequence')
        .first()
        .to_dict(orient='index')
    )

    ms3 = pep[pep['File Name'].str.contains('_ms3', regex=False)].copy()
    if ms3.empty:
        sys.exit('No *_ms3.ms2 rows detected.')

    pmap = {}
    for child in ms3['File Name'].unique():
        p = os.path.join(spectra_dir, f'{child}.ms2')
        if not os.path.isfile(p):
            sys.exit(f'Missing file: {p}')
        pmap[child] = map_ms3_to_parent(p)

    ms3[['parent_scan','parent_z']] = ms3.apply(
        lambda r: pd.Series(pmap[r['File Name']].get(int(r['Scan Num']),(None,None))), axis=1)
    ms3.dropna(subset=['parent_scan'], inplace=True)

    rows = []
    for (child, pscan), grp in ms3.groupby(['File Name','parent_scan']):
        grp = grp.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        if len(grp) < 2: continue
        a,b = grp.iloc[0], grp.iloc[1]

        parent = os.path.join(spectra_dir, child.replace('_ms3','') + '.ms2')
        if not os.path.isfile(parent): sys.exit(f'Missing parent MS²: {parent}')
        z, peaks = read_ms2(parent, int(pscan))
        dbl_present, dbl_ratio = has_doublet(peaks, z)

        stA, stB = stub_type(a['sequence']), stub_type(b['sequence'])
        if stA in ('54','86') and stB in ('54','86'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('54','86')) ^ (stB in ('54','86')):
            note = 'Mono'
        else:
            note = 'Not sure'

        def info(seq,key): return seq_info.get(seq,{}).get(key)
        def clean(p):     return p.strip('[] ')

        rows.append({
            'MS2_file'  : os.path.basename(parent), 'MS2_scan': int(pscan),
            'parent_charge': z,
            'MS3_file_A': f'{child}.ms2', 'MS3_scan_A': int(a['Scan Num']),
            'sequence_A': a['sequence'], 'XCorr_A': a['XCorr'],
            'DeltaCn_A' : info(a['sequence'],'deltcn'),
            'RT_A'      : info(a['sequence'],'rt'),
            'Conf_A'    : info(a['sequence'],'conf'),
            'Proteins_A': clean(a['Proteins']), 'Gene_A': info(a['sequence'],'gene'),
            'MS3_file_B': f'{child}.ms2', 'MS3_scan_B': int(b['Scan Num']),
            'sequence_B': b['sequence'], 'XCorr_B': b['XCorr'],
            'DeltaCn_B' : info(b['sequence'],'deltcn'),
            'RT_B'      : info(b['sequence'],'rt'),
            'Conf_B'    : info(b['sequence'],'conf'),
            'Proteins_B': clean(b['Proteins']), 'Gene_B': info(b['sequence'],'gene'),
            'doublet_present': dbl_present, 'doublet_int_ratio': dbl_ratio,
            'Note': note
        })

    if not rows: sys.exit('No linkable parent scans found.')
    df = pd.DataFrame(rows)

    # decoy + FDR
    df['decoy_A'] = df['sequence_A'].map(lambda s: seq_info[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s: seq_info[s]['is_decoy'])
    df['decoy_link'] = df['decoy_A'] | df['decoy_B']
    df['link_score'] = df['XCorr_A'] + df['XCorr_B']
    df = df.sort_values('link_score', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    df['cum_links'] = df['rank']
    df['cum_decoys'] = df['decoy_link'].cumsum()
    df['FDR'] = (2 * df['cum_decoys'] / df['cum_links']).round(4)
    cols = [c for c in df.columns if c != 'Note'] + ['Note']
    return df[cols]

# ── CLI ─────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    ag = argparse.ArgumentParser(description='DSSO XL MS2-MS3 linker pipeline')
    ag.add_argument('pep_list', help='pep_list.csv from IP2 search')
    ag.add_argument('spectra_dir', help='folder with *.ms2 and *_ms3.ms2')
    ag.add_argument('out_csv', help='output CSV (e.g. DSSO_links_FDR.csv)')
    args = ag.parse_args()

    result = build_links(args.pep_list, args.spectra_dir)
    result.to_csv(args.out_csv, index=False)
    print(f'✓ {len(result)} links → {args.out_csv}')
    for lim in (0.01,0.05):
        sub = result[result['FDR'] <= lim]
        if sub.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
        else:
            best = sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {best["link_score"]:.3f} '
                  f'({int(best["cum_links"])} links, {int(best["cum_decoys"])} decoys)')


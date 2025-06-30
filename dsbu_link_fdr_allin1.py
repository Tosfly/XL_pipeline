#!/usr/bin/env python3
"""
dsbu_link_fdr_allin1.py  –  Build DSBU cross-link table with QC + FDR  
rev 4 • adds stricter stub-complementarity logic  
      • +85/+85 or +111/+111 pairs are flagged “Rare, may not be true”  
      • keeps all previous columns (Conf, Gene, cleaned Proteins, etc.) and
        appends Note at the far right
"""

import os, re, sys, argparse
import pandas as pd

# ── patterns & constants ────────────────────────────────────────────────────
DECOY_RE   = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')
STUB85_RE  = re.compile(r'\(85\.052')
STUB111_RE = re.compile(r'\(111\.032')
DSBU_DELTA = 26.016      # 85- vs 111-stub gap (Da)

# ── lightweight *.ms2 helpers ───────────────────────────────────────────────
def map_ms3_to_parent(path: str):
    """Return {ms3_scan:int → (parent_scan:int, parent_charge:int)}."""
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
    return mapping


def read_ms2(path: str, wanted: int):
    """Return (charge:int, list[(mz,float,inten,float)]) for one MS² scan."""
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
    """Check for 26.016 Da reporter doublet; return (present?, ratio)."""
    if not peaks or not z:
        return False, 0.0
    peaks = sorted(peaks)
    dmz   = DSBU_DELTA / z
    topI  = max(i for _, i in peaks) or 1.0
    for i, (mz1, i1) in enumerate(peaks):
        tgt, tol = mz1 + dmz, ppm * (mz1 + dmz) / 1e6
        for mz2, i2 in peaks[i + 1:]:
            if mz2 - tgt > tol:
                break
            if abs(mz2 - tgt) <= tol:
                ratio = min(i1, i2) / topI
                return ratio >= min_ratio, round(ratio, 3)
    return False, 0.0

# ── core helpers ────────────────────────────────────────────────────────────
def stub_type(seq: str):
    if STUB85_RE.search(seq):
        return '85'
    if STUB111_RE.search(seq):
        return '111'
    return '?'


def build_links(pep_list: str, spectra_dir: str) -> pd.DataFrame:
    pep = pd.read_csv(pep_list)

    # dictionary: sequence → consolidated info
    seq_info = (
        pep.assign(
            is_decoy=pep['Proteins'].fillna('').str.contains(DECOY_RE) |
                     pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE),
            gene=pep['Protein Descriptions'].str.extract(GN_RE, expand=False),
            conf=pep.get('Conf%', None),
            deltcn=pep['DeltCN'],
            rt=pep['RetTime'],
        )
        .groupby('sequence')
        .first()
        .to_dict(orient='index')
    )

    # keep only MS³ rows
    ms3 = pep[pep['File Name'].str.contains('_ms3', regex=False)].copy()
    if ms3.empty:
        sys.exit('No “_ms3” rows in pep_list.csv.')

    # map child → parent scans
    pmap = {}
    for child in ms3['File Name'].unique():
        child_path = os.path.join(spectra_dir, f'{child}.ms2')
        if not os.path.isfile(child_path):
            sys.exit(f'Missing MS³ file: {child_path}')
        pmap[child] = map_ms3_to_parent(child_path)

    # add parent_scan & parent_z to each MS³ PSM
    ms3[['parent_scan', 'parent_z']] = ms3.apply(
        lambda r: pd.Series(pmap[r['File Name']].get(int(r['Scan Num']), (None, None))),
        axis=1,
    )
    ms3.dropna(subset=['parent_scan'], inplace=True)

    # ---------------- pair peptides per parent MS² --------------------------
    rows = []
    for (child, pscan), grp in ms3.groupby(['File Name', 'parent_scan']):
        grp = (
            grp.sort_values('XCorr', ascending=False)
               .drop_duplicates(subset='sequence', keep='first')
        )
        if len(grp) < 2:
            continue
        a, b = grp.iloc[0], grp.iloc[1]

        parent_path = os.path.join(spectra_dir, child.replace('_ms3', '') + '.ms2')
        if not os.path.isfile(parent_path):
            sys.exit(f'Missing parent MS²: {parent_path}')
        z, peaks = read_ms2(parent_path, int(pscan))
        dbl_present, dbl_ratio = has_doublet(peaks, z)

        # classify stubs
        stA, stB = stub_type(a['sequence']), stub_type(b['sequence'])
        if stA in ('85', '111') and stB in ('85', '111'):
            note = 'true XL' if stA != stB else 'Rare, may not be true'
        elif (stA in ('85', '111')) ^ (stB in ('85', '111')):
            note = 'Mono'
        else:
            note = 'Not sure'

        def get(seq, key):
            return seq_info.get(seq, {}).get(key)

        def clean(p):
            return p.strip('[] ')

        rows.append(
            {
                # parent
                'MS2_file': os.path.basename(parent_path),
                'MS2_scan': int(pscan),
                'parent_charge': z,
                # peptide A
                'MS3_file_A': f'{child}.ms2',
                'MS3_scan_A': int(a['Scan Num']),
                'sequence_A': a['sequence'],
                'XCorr_A': a['XCorr'],
                'DeltaCn_A': get(a['sequence'], 'deltcn'),
                'RT_A': get(a['sequence'], 'rt'),
                'Conf_A': get(a['sequence'], 'conf'),
                'Proteins_A': clean(a['Proteins']),
                'Gene_A': get(a['sequence'], 'gene'),
                # peptide B
                'MS3_file_B': f'{child}.ms2',
                'MS3_scan_B': int(b['Scan Num']),
                'sequence_B': b['sequence'],
                'XCorr_B': b['XCorr'],
                'DeltaCn_B': get(b['sequence'], 'deltcn'),
                'RT_B': get(b['sequence'], 'rt'),
                'Conf_B': get(b['sequence'], 'conf'),
                'Proteins_B': clean(b['Proteins']),
                'Gene_B': get(b['sequence'], 'gene'),
                # parent QC
                'doublet_present': dbl_present,
                'doublet_int_ratio': dbl_ratio,
                # stub note
                'Note': note,
            }
        )

    if not rows:
        sys.exit('No linkable parent scans found.')

    links = pd.DataFrame(rows)

    # -------------- decoy flags, score, FDR ---------------------------------
    links['decoy_A'] = links['sequence_A'].map(lambda s: seq_info[s]['is_decoy'])
    links['decoy_B'] = links['sequence_B'].map(lambda s: seq_info[s]['is_decoy'])
    links['decoy_link'] = links['decoy_A'] | links['decoy_B']
    links['link_score'] = links['XCorr_A'] + links['XCorr_B']

    links = links.sort_values('link_score', ascending=False).reset_index(drop=True)
    links['rank'] = links.index + 1
    links['cum_links'] = links['rank']
    links['cum_decoys'] = links['decoy_link'].cumsum()
    links['FDR'] = (2 * links['cum_decoys'] / links['cum_links']).round(4)

    # ensure Note is the right-most column
    cols = [c for c in links.columns if c != 'Note'] + ['Note']
    return links[cols]

# ── command-line entry point ────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description='DSBU cross-link builder with strict +85/+111 complementarity'
    )
    ap.add_argument('pep_list', help='pep_list.csv (IP2 export)')
    ap.add_argument('spectra_dir', help='directory containing *.ms2 and *_ms3.ms2')
    ap.add_argument('out_csv', help='output CSV (e.g. DSBU_links_FDR.csv)')
    args = ap.parse_args()

    df = build_links(args.pep_list, args.spectra_dir)
    df.to_csv(args.out_csv, index=False)
    print(f'✓ wrote {len(df)} links → {args.out_csv}')

    for lim in (0.01, 0.05):
        sub = df[df['FDR'] <= lim]
        if sub.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
            continue
        top = sub.iloc[0]
        print(
            f'{lim:.0%} FDR ⇒ link_score ≥ {top["link_score"]:.3f} '
            f'({int(top["cum_links"])} links, {int(top["cum_decoys"])} decoys)'
        )


if __name__ == '__main__':
    main()


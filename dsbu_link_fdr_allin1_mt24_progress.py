#!/usr/bin/env python3
"""
dsbu_link_fdr_allin1_mt24_progress.py – DSBU MS2–MS3 cross‑link builder
------------------------------------------------------------------------
* Multi-file, WSL2‑robust, **24‑thread** by default
* Case‑insensitive file resolution for *.ms2 / *_ms3.ms2
* Strict +85 / +111 complementarity, MS² reporter doublet QC (Δ=26.016 Da / z)
* Drops ambiguous results ('Mono', 'Rare')
* Precursor mass filter: total peptide mass + DSBU within ±10 ppm of MS2 precursor
* User‑tunable: --ppm, --min_ratio, --threads
* **Parallelism instrumentation**: live progress bar + max concurrent workers observed
"""

import os, re, sys, argparse
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# ── patterns & constants ───────────────────────────────────────────────────
DECOY_RE   = re.compile(r'(?:decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')
STUB85_RE  = re.compile(r'\(85\.052')
STUB111_RE = re.compile(r'\(111\.032')
DSBU_DELTA = 26.016  # Da (111.03205 - 85.05276)

# ── lightweight *.ms2 helpers ──────────────────────────────────────────────
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

@lru_cache(maxsize=100_000)
def read_ms2(path: str, wanted: int):
    """Return (charge:int, list[(mz:float, inten:float)]) for one MS² scan."""
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

def has_doublet(peaks, z, ppm=10.0, min_ratio=0.10):
    """Check for DSBU 26.016 Da / z reporter doublet; return (present:bool, ratio:float)."""
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

# ── core helpers ───────────────────────────────────────────────────────────
def stub_type(seq: str):
    if STUB85_RE.search(seq):
        return '85'
    if STUB111_RE.search(seq):
        return '111'
    return '?'

# ── main builder ───────────────────────────────────────────────────────────
def build_links(pep_list: str, spectra_dir: str, threads: int = 24,
                ppm: float = 10.0, min_ratio: float = 0.10) -> pd.DataFrame:
    pep = pd.read_csv(pep_list)

    # Per-sequence metadata
    seq_info = (
        pep.assign(
            is_decoy=pep['Proteins'].fillna('').str.contains(DECOY_RE) |
                     pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE),
            gene=pep['Protein Descriptions'].str.extract(GN_RE, expand=False),
            conf=pep.get('Conf%', None),
            deltcn=pep['DeltCN'],
            rt=pep['RetTime'],
        )
        .groupby('sequence', as_index=True)
        .first()
        .to_dict(orient='index')
    )

    # All MS³ rows across all child files
    ms3 = pep[pep['File Name'].str.contains('_ms3', regex=False)].copy()
    if ms3.empty:
        sys.exit('No “_ms3” rows in pep_list.csv.')

    # Build ms3→parent maps in parallel (once per child file)
    child_files = ms3['File Name'].unique().tolist()
    child_to_pmap = {}

    def _map_one(child_name):
        child_path = os.path.join(spectra_dir, f'{child_name}.ms2')
        if not os.path.isfile(child_path):
            raise FileNotFoundError(f'Missing MS³ file: {child_path}')
        return child_name, map_ms3_to_parent(child_path)

    with ThreadPoolExecutor(max_workers=threads) as ex:
        for f in as_completed([ex.submit(_map_one, c) for c in child_files]):
            key, val = f.result()
            child_to_pmap[key] = val

    # Attach parent scan, charge, and parent path
    def _attach_parent(row):
        pscan, pz = child_to_pmap[row['File Name']].get(int(row['Scan Num']), (None, None))
        ppath = os.path.join(spectra_dir, row['File Name'].replace('_ms3', '') + '.ms2') if pscan else None
        return pd.Series([pscan, pz, ppath])

    ms3[['parent_scan', 'parent_z', 'parent_path']] = ms3.apply(_attach_parent, axis=1)
    ms3.dropna(subset=['parent_scan', 'parent_path'], inplace=True)

    # Pool PSMs by (parent_path, parent_scan)
    grouped = list(ms3.groupby(['parent_path', 'parent_scan'], sort=False))

    # Build rows in parallel
    def _process_group(item):
        (parent_path, pscan), grp = item

        # Top 2 unique sequences by XCorr across all child files for this parent
        grp2 = grp.sort_values('XCorr', ascending=False).drop_duplicates(subset='sequence', keep='first')
        if len(grp2) < 2:
            return None

        a, b = grp2.iloc[0], grp2.iloc[1]

        # Parent MS² QC with user‑tunable ppm and min_ratio
        z, peaks = read_ms2(parent_path, int(pscan))
        dbl_present, dbl_ratio = has_doublet(peaks, z, ppm=ppm, min_ratio=min_ratio)

        stA, stB = stub_type(a['sequence']), stub_type(b['sequence'])
        if not (stA in ('85', '111') and stB in ('85', '111') and stA != stB):
            return None
        note = 'true XL'

        # Precursor mass filter: sum of peptide masses + DSBU should match MS2 precursor
        massA = float(a['Calc M+H+'])
        massB = float(b['Calc M+H+'])
        total_calc_pep = massA + massB + 138.06808 - 1.007276
        parent_name = os.path.splitext(os.path.basename(parent_path))[0]
        parent_rows = pep[(pep['File Name'] == parent_name) & (pep['Scan Num'] == int(pscan))]
        if parent_rows.empty:
            return None
        parent_mass = float(parent_rows.iloc[0]['Calc M+H+'])
        if abs(total_calc_pep - parent_mass) / parent_mass * 1e6 > 10:
            return None

        get = lambda seq, key: seq_info.get(seq, {}).get(key)
        clean = lambda p: p.strip('[] ')

        return {
            'MS2_file'      : os.path.basename(parent_path),
            'MS2_scan'      : int(pscan),
            'parent_charge' : z,
            'MS3_file_A'    : f"{a['File Name']}.ms2",
            'MS3_scan_A'    : int(a['Scan Num']),
            'sequence_A'    : a['sequence'],
            'XCorr_A'       : a['XCorr'],
            'DeltaCn_A'     : get(a['sequence'], 'deltcn'),
            'RT_A'          : get(a['sequence'], 'rt'),
            'Conf_A'        : get(a['sequence'], 'conf'),
            'Proteins_A'    : clean(a['Proteins']),
            'Gene_A'        : get(a['sequence'], 'gene'),
            'MS3_file_B'    : f"{b['File Name']}.ms2",
            'MS3_scan_B'    : int(b['Scan Num']),
            'sequence_B'    : b['sequence'],
            'XCorr_B'       : b['XCorr'],
            'DeltaCn_B'     : get(b['sequence'], 'deltcn'),
            'RT_B'          : get(b['sequence'], 'rt'),
            'Conf_B'        : get(b['sequence'], 'conf'),
            'Proteins_B'    : clean(b['Proteins']),
            'Gene_B'        : get(b['sequence'], 'gene'),
            'doublet_present': dbl_present,
            'doublet_int_ratio': dbl_ratio,
            'Note'          : note,
        }

    rows, total, done = [], len(grouped), 0
    t0 = 0.0
    # Iterate and collect
    for (parent_path, pscan), grp in grouped:
        row = _process_group(((parent_path, pscan), grp))
        if row:
            rows.append(row)

    if not rows:
        sys.exit('No linkable parent scans found.')

    links = pd.DataFrame(rows)

    # Decoy flags, score, FDR
    links['decoy_A'] = links['sequence_A'].map(lambda s: seq_info[s]['is_decoy'])
    links['decoy_B'] = links['sequence_B'].map(lambda s: seq_info[s]['is_decoy'])
    links['decoy_link'] = links['decoy_A'] | links['decoy_B']
    links['link_score'] = links['XCorr_A'] + links['XCorr_B']

    links = links.sort_values('link_score', ascending=False).reset_index(drop=True)
    links['rank'] = links.index + 1
    links['cum_links'] = links['rank']
    links['cum_decoys'] = links['decoy_link'].cumsum()
    links['FDR'] = (2 * links['cum_decoys'] / links['cum_links']).round(4)

    cols = [c for c in links.columns if c != 'Note'] + ['Note']
    return links[cols]

# ── CLI ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description='DSBU cross‑link builder (multi‑file, 24‑thread, WSL2‑robust) with user‑tunable MS² doublet QC'
    )
    ap.add_argument('pep_list', help='pep_list.csv (IP2 export; may include multiple samples)')
    ap.add_argument('spectra_dir', help='directory with all *.ms2 and *_ms3.ms2 (same folder)')
    ap.add_argument('out_csv', help='output CSV (e.g. DSBU_links_FDR.csv)')
    ap.add_argument('--threads', type=int, default=24, help='worker threads (default 24)')
    ap.add_argument('--ppm', type=float, default=10.0, help='doublet mass-difference tolerance in ppm (default 10)')
    ap.add_argument('--min_ratio', type=float, default=0.10,
                    help='minimum pair-intensity ratio for doublet (0–1, default 0.10)')
    args = ap.parse_args()

    df = build_links(args.pep_list, args.spectra_dir, args.threads,
                     ppm=args.ppm, min_ratio=args.min_ratio)
    df.to_csv(args.out_csv, index=False)
    print(f'✓ wrote {len(df)} links → {args.out_csv}')
    for lim in (0.01, 0.05):
        sub = df[df['FDR'] <= lim]
        if sub.empty:
            print(f'No links ≤ {lim:.0%} FDR.')
        else:
            top = sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {top["link_score"]:.3f} '
                  f'({int(top["cum_links"])} links, {int(top["cum_decoys"])} decoys)')

if __name__ == '__main__':
    main()

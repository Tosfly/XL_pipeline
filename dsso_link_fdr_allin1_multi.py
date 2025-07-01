#!/usr/bin/env python3
"""
dsso_link_fdr_allin1_multi.py  –  DSSO MS2-MS3 XL analysis (multithreaded)
--------------------------------------------------------------------
* Pairs the two best MS³ PSMs for each parent MS²
* Checks DSSO reporter doublet (31.993 Da / z)
* Enforces strict +54/+86 complementarity
* Adds Conf%, gene names, reporter-QC, Note
* Calculates target-decoy FDR
* Uses ThreadPoolExecutor (–-threads N) for speed-up
--------------------------------------------------------------------
CLI
  pip install pandas
  python dsso_link_fdr_allin1.py pep_list.csv /path/to/spectra DSSO_links_FDR.csv --threads 8
"""
import os, re, sys, argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# DSSO constants
STUB54_RE = re.compile(r'\(54\.010')
STUB86_RE = re.compile(r'\(86\.003')
DELTA_NEUT = 86.00307 - 54.01056          # 31.99251 Da
DECOY_RE = re.compile(r'(decoy|reverse|rev_)', re.I)
GN_RE = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ──  mini .ms2 helpers ───────────────────────────────────────────
def map_ms3_to_parent(path):
    """Return (child_basename, {ms3_scan:(parent_scan,charge)})"""
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
    return os.path.basename(path)[:-4], mapping  # strip ".ms2"

def read_ms2(parent_path, scan):
    charge, peaks, reading = None, [], False
    with open(parent_path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                reading = (int(ln.split()[1]) == scan)
                peaks, charge = ([], None) if reading else (peaks, charge)
                continue
            if not reading: continue
            if ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
            elif ln[0].isdigit():
                mz, inten = map(float, ln.split()[:2]); peaks.append((mz,inten))
            elif ln == '\n':
                break
    return charge, peaks

def has_doublet(peaks, z, ppm=10, min_ratio=0.10):
    if not peaks or not z: return False, 0.0
    dmz = DELTA_NEUT / z
    peaks = sorted(peaks); topI = max(i for _,i in peaks) or 1.0
    for i,(mz1,i1) in enumerate(peaks):
        tgt, tol = mz1+dmz, ppm*(mz1+dmz)/1e6
        for mz2,i2 in peaks[i+1:]:
            if mz2 - tgt > tol: break
            if abs(mz2 - tgt) <= tol:
                return (min(i1,i2)/topI)>=min_ratio, round(min(i1,i2)/topI,3)
    return False, 0.0

def stub(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

# ── main build ──────────────────────────────────────────────────
def build_links(pep_csv, spectra_dir, n_threads):
    pep = pd.read_csv(pep_csv)
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
    if ms3.empty: sys.exit('No *_ms3.ms2 rows.')

    # 1) parent maps in parallel
    child_files = ms3['File Name'].unique()
    parent_map = {}
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        fut = {ex.submit(map_ms3_to_parent, os.path.join(spectra_dir,f'{c}.ms2')): c for c in child_files}
        for f in as_completed(fut):
            child, mp = f.result(); parent_map[child] = mp

    ms3[['parent_scan','parent_z']] = ms3.apply(
        lambda r: pd.Series(parent_map[r['File Name']].get(int(r['Scan Num']),(None,None))), axis=1)
    ms3.dropna(subset=['parent_scan'], inplace=True)

    groups = list(ms3.groupby(['File Name','parent_scan']))

    def do_group(item):
        (child, pscan), g = item
        g = g.sort_values('XCorr', ascending=False).drop_duplicates('sequence')
        if len(g) < 2: return None
        a,b = g.iloc[0], g.iloc[1]
        parent = os.path.join(spectra_dir, child.replace('_ms3','')+'.ms2')
        if not os.path.isfile(parent): return None
        z, peaks = read_ms2(parent, int(pscan))
        dbl, ratio = has_doublet(peaks,z)

        stA, stB = stub(a['sequence']), stub(b['sequence'])
        note = ('true XL' if stA!=stB else 'Rare, may not be true') if stA in ('54','86') and stB in ('54','86') \
               else ('Mono' if (stA in ('54','86')) ^ (stB in ('54','86')) else 'Not sure')
        get = lambda s,k: seq_meta.get(s,{}).get(k)
        clean = lambda p: p.strip('[] ')
        return dict(
            MS2_file=os.path.basename(parent), MS2_scan=int(pscan), parent_charge=z,
            MS3_file_A=f'{child}.ms2', MS3_scan_A=int(a['Scan Num']),
            sequence_A=a['sequence'], XCorr_A=a['XCorr'],
            DeltaCn_A=get(a['sequence'],'deltcn'), RT_A=get(a['sequence'],'rt'),
            Conf_A=get(a['sequence'],'conf'), Proteins_A=clean(a['Proteins']),
            Gene_A=get(a['sequence'],'gene'),
            MS3_file_B=f'{child}.ms2', MS3_scan_B=int(b['Scan Num']),
            sequence_B=b['sequence'], XCorr_B=b['XCorr'],
            DeltaCn_B=get(b['sequence'],'deltcn'), RT_B=get(b['sequence'],'rt'),
            Conf_B=get(b['sequence'],'conf'), Proteins_B=clean(b['Proteins']),
            Gene_B=get(b['sequence'],'gene'),
            doublet_present=dbl, doublet_int_ratio=ratio, Note=note)

    rows=[]
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        for row in ex.map(do_group, groups):
            if row: rows.append(row)
    if not rows: sys.exit('No linkable groups.')

    df = pd.DataFrame(rows)
    df['decoy_A'] = df['sequence_A'].map(lambda s:seq_meta[s]['is_decoy'])
    df['decoy_B'] = df['sequence_B'].map(lambda s:seq_meta[s]['is_decoy'])
    df['decoy_link']= df['decoy_A'] | df['decoy_B']
    df['link_score']= df['XCorr_A'] + df['XCorr_B']
    df.sort_values('link_score',ascending=False,inplace=True,ignore_index=True)
    df['rank']=df.index+1; df['cum_links']=df['rank']; df['cum_decoys']=df['decoy_link'].cumsum()
    df['FDR']=(2*df['cum_decoys']/df['cum_links']).round(4)
    cols=[c for c in df.columns if c!='Note']+['Note']; return df[cols]

# ── main ─────────────────────────────────────────────────────────
if __name__=='__main__':
    p=argparse.ArgumentParser(description='DSSO XL multithreaded pipeline')
    p.add_argument('pep_list'); p.add_argument('spectra_dir'); p.add_argument('out_csv')
    p.add_argument('--threads',type=int,default=8)
    a=p.parse_args()
    table=build_links(a.pep_list,a.spectra_dir,a.threads)
    table.to_csv(a.out_csv,index=False)
    print(f'✓ {len(table)} links → {a.out_csv}')
    for lim in (0.01,0.05):
        sub=table[table['FDR']<=lim]
        if sub.empty: print(f'No links ≤ {lim:.0%} FDR.')
        else:
            top=sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {top["link_score"]:.3f} '
                  f'({int(top["cum_links"])} links, {int(top["cum_decoys"])} decoys)')

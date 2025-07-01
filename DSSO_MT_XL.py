#!/usr/bin/env python3
"""
DSSO MS2–MS3 XL pipeline (multithreaded, stub-mass fixed)
"""
import os, re, sys, argparse
import pandas as pd, numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

# ── DSSO stub masses -------------------------------------------------------
STUB_A_MASS = 54.01056                # thiol fragment
STUB_B_MASS = 85.98264                # alkene fragment
DELTA_NEUT  = STUB_B_MASS - STUB_A_MASS   # 31.97208 Da

STUB54_RE  = re.compile(r'\(54\.010')
STUB86_RE  = re.compile(r'\(85\.982')

DECOY_RE   = re.compile(r'(?:decoy|reverse|rev_)', re.I)
GN_RE      = re.compile(r'GN=([A-Za-z0-9_.-]+)')

# ── tiny .ms2 helpers ------------------------------------------------------
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

@lru_cache(maxsize=10_000)
def read_ms2_scan(path, wanted):
    charge, peaks, reading = None, [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('S'):
                reading = int(ln.split()[1]) == wanted
                charge = None; peaks = [] if reading else peaks
                continue
            if not reading: continue
            if ln.startswith('Z'):
                charge = int(float(ln.split()[1]))
            elif ln and ln[0].isdigit():
                mz,i = map(float, ln.split()[:2]); peaks.append((mz,i))
            elif ln=='\n': break
    return charge, peaks

def has_doublet(peaks,z,ppm=10,min_ratio=0.10):
    if not peaks or not z: return False,0.0
    dmz = DELTA_NEUT/z
    mzs, ints = zip(*sorted(peaks))
    top = max(ints) or 1.0
    for i,m1 in enumerate(mzs):
        tgt, tol = m1+dmz, ppm*(m1+dmz)/1e6
        for m2, i2 in zip(mzs[i+1:],ints[i+1:]):
            if m2-tgt > tol: break
            if abs(m2-tgt)<=tol:
                ratio = min(ints[i],i2)/top
                return ratio>=min_ratio, round(ratio,3)
    return False,0.0

def stub(seq):
    if STUB54_RE.search(seq): return '54'
    if STUB86_RE.search(seq): return '86'
    return '?'

# ── build table -----------------------------------------------------------
def build_links(pep_csv, spectra, threads):
    pep = pd.read_csv(pep_csv, low_memory=False)
    pep['is_decoy'] = pep['Proteins'].fillna('').str.contains(DECOY_RE) | \
                      pep['Protein Descriptions'].fillna('').str.contains(DECOY_RE)
    pep['gene']  = pep['Protein Descriptions'].str.extract(GN_RE,expand=False)
    pep['conf']  = pep.get('Conf%', pd.Series(dtype=float))
    meta = pep.groupby('sequence').first()\
             .rename(columns={'DeltCN':'deltcn','RetTime':'rt'})\
             .to_dict('index')

    ms3 = pep[pep['File Name'].str.contains('_ms3',na=False)].copy()
    if ms3.empty: sys.exit('No *_ms3 rows.')

    # build child→parent map
    maps={}
    with ThreadPoolExecutor(threads) as ex:
        fut = [ex.submit(parse_ms3_child,os.path.join(spectra,f'{c}.ms2'))
               for c in ms3['File Name'].unique()]
        for f in as_completed(fut):
            child,m = f.result(); maps[child]=m

    ms3[['parent_scan','parent_z']] = ms3.apply(
        lambda r: pd.Series(maps[r['File Name']].get(int(r['Scan Num']),(None,None))),axis=1)
    ms3.dropna(subset=['parent_scan'], inplace=True)

    groups = list(ms3.groupby(['File Name','parent_scan']))

    def work(item):
        (child,pscan),g=item
        g=g.sort_values('XCorr',ascending=False).drop_duplicates('sequence')
        if len(g)<2: return
        a,b=g.iloc[0],g.iloc[1]
        pfile=os.path.join(spectra,child.replace('_ms3','')+'.ms2')
        if not os.path.isfile(pfile): return
        z,peaks=read_ms2_scan(pfile,int(pscan))
        dbl,ratio=has_doublet(peaks,z)
        stA,stB=stub(a['sequence']),stub(b['sequence'])
        note='true XL' if {stA,stB}=={'54','86'} else \
             ('Rare, may not be true' if stA==stB!='?' else
              'Mono' if '?' not in (stA,stB) else 'Not sure')
        get=lambda s,k:meta.get(s,{}).get(k)
        clean=lambda p:p.strip('[] ')
        return dict(
            MS2_file=os.path.basename(pfile),MS2_scan=int(pscan),parent_charge=z,
            MS3_file_A=f'{child}.ms2',MS3_scan_A=int(a['Scan Num']),
            sequence_A=a['sequence'],XCorr_A=a['XCorr'],
            DeltaCn_A=get(a['sequence'],'deltcn'),RT_A=get(a['sequence'],'rt'),
            Conf_A=get(a['sequence'],'conf'),Proteins_A=clean(a['Proteins']),
            Gene_A=get(a['sequence'],'gene'),
            MS3_file_B=f'{child}.ms2',MS3_scan_B=int(b['Scan Num']),
            sequence_B=b['sequence'],XCorr_B=b['XCorr'],
            DeltaCn_B=get(b['sequence'],'deltcn'),RT_B=get(b['sequence'],'rt'),
            Conf_B=get(b['sequence'],'conf'),Proteins_B=clean(b['Proteins']),
            Gene_B=get(b['sequence'],'gene'),
            doublet_present=dbl,doublet_int_ratio=ratio,Note=note)

    with ThreadPoolExecutor(threads) as ex:
        rows=[r for r in ex.map(work,groups) if r]

    if not rows: sys.exit('No links.')
    df=pd.DataFrame(rows)
    df['decoy_A']=df['sequence_A'].map(lambda s:meta[s]['is_decoy'])
    df['decoy_B']=df['sequence_B'].map(lambda s:meta[s]['is_decoy'])
    df['decoy_link']=df['decoy_A']|df['decoy_B']
    df['link_score']=df['XCorr_A']+df['XCorr_B']
    df.sort_values('link_score',ascending=False,inplace=True,ignore_index=True)
    df['rank']=np.arange(1,len(df)+1)
    df['cum_links']=df['rank']; df['cum_decoys']=df['decoy_link'].cumsum()
    df['FDR']=np.round(2*df['cum_decoys']/df['cum_links'],4)
    return df[[c for c in df.columns if c!='Note']+['Note']]

# ── CLI --------------------------------------------------------------------
if __name__=='__main__':
    ar=argparse.ArgumentParser(description='DSSO multithread XL pipeline')
    ar.add_argument('pep_list'); ar.add_argument('spectra_dir'); ar.add_argument('out_csv')
    ar.add_argument('--threads',type=int,default=os.cpu_count())
    args=ar.parse_args()

    print(f'Processing with {args.threads} threads …')
    out=build_links(args.pep_list,args.spectra_dir,args.threads)
    out.to_csv(args.out_csv,index=False); print(f'✓ {len(out)} links → {args.out_csv}')
    for lim in (0.01,0.05):
        sub=out[out['FDR']<=lim]
        if sub.empty: print(f'No links ≤ {lim:.0%} FDR.')
        else:
            best=sub.iloc[0]
            print(f'{lim:.0%} FDR ⇒ link_score ≥ {best["link_score"]:.3f} '
                  f'({int(best["cum_links"])} links, {int(best["cum_decoys"])} decoys)')

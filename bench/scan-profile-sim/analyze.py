"""Group-stratified TDC on sage_discriminant_score; combine arms per scan group.
Usage: analyze.py  -> prints tables, writes results.tsv"""
import pandas as pd, numpy as np, os
O = '/mnt/data1/scan-profile-bench/sim'
GMAP = {'ot2': 'SF_200217_U2OS_TiO2_HCD_nlEThcD_OT_rep2', 'quad1': 'SF_200217_U2OS_TiO2_HCD_nlETcaD_quad_OT_rep1',
        'hyb': '01625b_GE2-TUM_first_pool_13_01_01-2xIT_2xHCD-1h-R1',
        'ethcd': 'B191108_10_Lumos_AB_DE_165_Ec_dualHCD_EThcD_1', 'etd': 'B191108_09_Lumos_AB_DE_165_Ec_dualHCD_ETD_1'}
_cache = {}
def load(run, f):
    if run in _cache: return _cache[run]
    p = f'{O}/out/{run}/results.sage.parquet'
    if not os.path.exists(p): return None
    d = pd.read_parquet(p, columns=['scannr', 'peptide', 'rank', 'is_decoy', 'sage_discriminant_score', 'spectrum_q', 'peptide_q'])
    d = d[d['rank'] == 1].copy()
    d['peptide'] = d.peptide.str.replace('-[+1.007825]', '', regex=False)  # z-dot hack arms
    d['scan'] = d.scannr.str.extract(r'scan=(\d+)')[0]
    g = pd.read_csv(f'{O}/{GMAP[f]}.groups.tsv', sep='\t', names=['scan', 'group'], dtype=str)
    d = d.merge(g, on='scan', how='left')
    d['gq'] = np.nan
    for _, idx in d.groupby('group').groups.items():
        s = d.loc[idx].sort_values('sage_discriminant_score', ascending=False)
        dec = s.is_decoy.cumsum(); tgt = (~s.is_decoy).cumsum()
        q = (dec / tgt.clip(lower=1)).iloc[::-1].cummin().iloc[::-1]
        d.loc[s.index, 'gq'] = q.values
    _cache[run] = d
    return d
rows = []
def arm(f, label, choice, default=None):
    """choice: dict group->run; default run for other groups"""
    parts = []
    groups = pd.read_csv(f'{O}/{GMAP[f]}.groups.tsv', sep='\t', names=['s', 'g']).g.unique()
    for g in sorted(groups):
        run = choice.get(g, default)
        d = load(run, f)
        if d is None: return
        parts.append(d[d.group == g])
    d = pd.concat(parts)
    ok = d[(d.gq <= 0.01) & ~d.is_decoy]
    # peptide-level TDC within each group: best score per (peptide, decoy) in group, union of passing targets
    peps = set()
    for g, x in d.groupby('group'):
        b = x.sort_values('sage_discriminant_score', ascending=False).drop_duplicates(['peptide', 'is_decoy'])
        dec = b.is_decoy.cumsum(); tgt = (~b.is_decoy).cumsum()
        q = (dec / tgt.clip(lower=1)).iloc[::-1].cummin().iloc[::-1]
        peps |= set(b.peptide[(q.values <= 0.01) & ~b.is_decoy.values])
    per = ok.groupby('group').size().to_dict()
    r = dict(file=f, arm=label, psms=len(ok), peps_from_psms=ok.peptide.nunique(), peptides=len(peps),
             per_group=' '.join(f'{k}={v}' for k, v in sorted(per.items())))
    if not choice:
        a = load(default, f)
        r['native_psms'] = int(((a.spectrum_q <= 0.01) & ~a.is_decoy).sum())
        r['native_peps'] = a[(a.peptide_q <= 0.01) & ~a.is_decoy].peptide.nunique()
    rows.append(r)
IT = ['ion_trap/cid', 'ion_trap/hcd', 'ion_trap/ethcd']
# ---- A: fragment tolerance
for f in ('ot2', 'quad1', 'hyb'):
    arm(f, 'ppm-all', {}, f'A-{f}-ppm'); arm(f, 'da-all', {}, f'A-{f}-da')
    arm(f, 'combined(OT ppm, IT Da)', {g: f'A-{f}-da' for g in IT}, f'A-{f}-ppm')
# ---- B: ion kinds, Orbitrap-only files
for f, act in (('ethcd', 'orbitrap/ethcd'), ('etd', 'orbitrap/etd')):
    for ions in ('by', 'cz', 'bycz'): arm(f, f'{ions}-all', {}, f'B-{f}-{ions}')
    arm(f, 'combined(HCD by, ETx cz)', {act: f'B-{f}-cz'}, f'B-{f}-by')
    arm(f, 'combined(HCD by, ETx bycz)', {act: f'B-{f}-bycz'}, f'B-{f}-by')
for f, act in (('ethcd', 'orbitrap/ethcd'), ('etd', 'orbitrap/etd')):
    arm(f, 'czdot-all (hack)', {}, f'B-{f}-czdot')
    arm(f, 'combined(HCD by, ETx czdot)', {act: f'B-{f}-czdot'}, f'B-{f}-by')
    arm(f, 'combined(HCD by, ETx b+c+zdot+y+1)', {act: f'B-{f}-byczdot'}, f'B-{f}-by')
# ---- B on PXD007058 (IT EThcD); HCD at ppm b/y
for f in ('ot2', 'quad1'):
    for ions in ('cz', 'bycz'):
        for t in ('ppm', 'da'): arm(f, f'{ions}-{t}-all', {}, f'B-{f}-{ions}-{t}')
    arm(f, 'combined(HCD ppm by, IT Da cz)', {'ion_trap/ethcd': f'B-{f}-cz-da'}, f'A-{f}-ppm')
    arm(f, 'combined(HCD ppm by, IT Da bycz)', {'ion_trap/ethcd': f'B-{f}-bycz-da'}, f'A-{f}-ppm')
r = pd.DataFrame(rows)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
print(r.to_string(index=False))
r.to_csv(f'{O}/results.tsv', sep='\t', index=False)

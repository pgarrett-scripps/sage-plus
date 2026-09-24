import json
S = '/mnt/data1/scan-profile-bench'; O = f'{S}/sim'
base = json.load(open('/mnt/data1/hybrid-recal-bench/it-ppm-off.json'))
for k in ('output_directory', 'mzml_paths'): base.pop(k, None)
base['mass_recalibration'] = 'off'
base.pop('tolerance_mode', None); base['report_psms'] = 1; base['output_filter'] = {'psm_q_value': 1.0}
PPM = {'ppm': [-20, 20]}; DA = {'da': [-0.5, 0.5]}
phos = {'M': [15.9949], 'S': [79.966331], 'T': [79.966331], 'Y': [79.966331]}
files = {
 'ot2':   (f'{S}/PXD007058/mzml/SF_200217_U2OS_TiO2_HCD_nlEThcD_OT_rep2.mzML', 'human', phos),
 'quad1': (f'{S}/PXD007058/mzml/SF_200217_U2OS_TiO2_HCD_nlETcaD_quad_OT_rep1.mzML', 'human', phos),
 'hyb':   ('/mnt/data1/hybrid-recal-bench/mzml/01625b_GE2-TUM_first_pool_13_01_01-2xIT_2xHCD-1h-R1.mzML', 'human', {'M': [15.9949]}),
 'ethcd': (f'{S}/PXD018176/mzml/B191108_10_Lumos_AB_DE_165_Ec_dualHCD_EThcD_1.mzML', 'ecoli', {'M': [15.9949]}),
 'etd':   (f'{S}/PXD018176/mzml/B191108_09_Lumos_AB_DE_165_Ec_dualHCD_ETD_1.mzML', 'ecoli', {'M': [15.9949]}),
}
fasta = {'human': f'{S}/fasta/human_sp.fasta', 'ecoli': f'{S}/fasta/ecoli_k12_sp.fasta'}
runs = {}
# Experiment A: fragment tolerance (b,y)
for f in ('ot2', 'quad1', 'hyb'):
    runs[f'A-{f}-ppm'] = (f, PPM, ['b', 'y'])
    runs[f'A-{f}-da'] = (f, DA, ['b', 'y'])
# Experiment B: ion kinds. Orbitrap files at 20 ppm; PXD007058 at 0.5 Da too (ion-trap EThcD)
for f in ('ethcd', 'etd'):
    for ions in ('by', 'cz', 'bycz'):
        runs[f'B-{f}-{ions}'] = (f, PPM, list(ions))
for f in ('ot2', 'quad1'):
    for ions in ('cz', 'bycz'):
        runs[f'B-{f}-{ions}-da'] = (f, DA, list(ions))
        runs[f'B-{f}-{ions}-ppm'] = (f, PPM, list(ions))
for name, (f, tol, ions) in runs.items():
    c = json.loads(json.dumps(base))
    path, sp, vm = files[f]
    c['mzml_paths'] = [path]; c['fragment_tol'] = tol
    c['database']['fasta'] = fasta[sp]; c['database']['ion_kinds'] = ions
    c['database']['variable_mods'] = vm
    if f in ('ot2', 'quad1'):
        c['database']['prefilter'] = True; c['database']['enzyme']['missed_cleavages'] = 1; c['database']['enzyme']['max_len'] = 35
    c['output_directory'] = f'{O}/out/{name}'
    json.dump(c, open(f'{O}/cfg/{name}.json', 'w'), indent=1)
    print(name)

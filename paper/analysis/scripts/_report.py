"""Shared report statistics from traceable result snapshots."""
import json
from statistics import median
from _scientific import PAPER, load, timing

REPORT_INPUTS = [f'analysis/data/report-extension/{name}.json' for name in ('public','ptm','lfq','scaling','verification')]


def report(name):
    return json.loads((PAPER/f'analysis/data/report-extension/{name}.json').read_text())['data']


def scale_summary():
    rows = report('scaling')
    return [dict(engine=e, threads=t, seconds=median(r['seconds'] for r in rows if r['engine']==e and r['threads']==t and not r['warmup']), rss=median(r['rss'] for r in rows if r['engine']==e and r['threads']==t and not r['warmup'])) for e in ('upstream','plus') for t in (1,2,4,8)]


def workloads():
    pilot,_ = load()
    output = []
    for study,label in (('PXD001468','HEK public'),('PXD028735','Mixture public')):
        output.append(dict(label=label,kind='Repeated fixed input',engines={r['engine']:dict(seconds=r['seconds'],rss=r['rss']) for r in timing(pilot) if r['study']==study}))
    for study,label in (('human','HEK entrapment'),('hye','Mixture entrapment')):
        output.append(dict(label=label,kind='Across files and seeds',engines={e:{k:median(r[source] for r in pilot['jobs'] if r['suite'].startswith(f'entrapment-{study}-') and r['engine']==e) for k,source in (('seconds','wall_seconds'),('rss','peak_rss_mib'))} for e in ('upstream','plus')}))
    for work,label in (('standard','Local standard'),('common-mods','Local modified')):
        output.append(dict(label=label,kind='Contextual local timing',engines={e:{k:median(r[source] for r in pilot['jobs'] if r['suite']=='local-paired' and r['workload']==work and r['engine']==e and not r['warmup']) for k,source in (('seconds','wall_seconds'),('rss','peak_rss_mib'))} for e in ('upstream','plus')}))
    return output


def add_report_stats(st):
    audit = json.loads((PAPER/'analysis/data/report-extension/verification.json').read_text())
    for key in ('completed_jobs','source_files_verified'):
        st.add(f'report.extension.{key}',audit[key],fmt=',',desc=f'Report extension {key}',sign='+')
    public = report('public')
    for study in ('PXD001468','PXD028735'):
        group = [p for p in public if study in p['pair']]
        for measure in ('target_psms','target_peptidoforms'):
            delta = [100*(p['engines']['plus']['0.01'][measure]/p['engines']['upstream']['0.01'][measure]-1) for p in group]
            for bound,value in (('min',min(delta)),('max',max(delta))):
                st.add(f'report.{study}.{measure}.{bound}',value,fmt='+.2f',desc=f'{study} {bound} per-file change in {measure} percent',between=(-100,100))
    for engine in ('upstream','plus'):
        counts = [p['disagreement'][engine] for p in public]
        total = sum(sum(c.values()) for c in counts)
        same = sum(c.get('same_assignment_above_threshold',0) for c in counts)
        st.add(f'report.{engine}.threshold_disagreement',100*same/total,fmt='.1f',desc='Percent of engine-only accepted PSMs with the identical assignment above threshold in the other engine',between=(0,100))
    for row in report('ptm'):
        for key in ('spectrum_accepted','joint_accepted'):
            if row['engine'] != 'upstream' or (key == 'joint_accepted' and row['library'] != 1):
                continue
            st.add(f"report.ptm.{row['library']}.{row['engine']}.{key}",row[key],fmt=',',desc=f"Synthetic HCD {row['library']} {row['engine']} {key}",between=(0,1000000))
    for r in workloads()[2:4]:
        for key in ('seconds','rss'):
            value=100*(r['engines']['plus'][key]/r['engines']['upstream'][key]-1)
            st.add(f"report.{r['label'].split()[0].lower()}.entrapment.{key}",value,fmt='+.1f',desc=f"{r['label']} relative change in {key} percent",between=(-100,100))
    q = report('lfq')
    quant_engines = {r['engine']:r for r in q['engines']}
    for threshold in ('0.01','0.05'):
        change = 100*(quant_engines['plus']['thresholds'][threshold]['target_precursors']/quant_engines['upstream']['thresholds'][threshold]['target_precursors']-1)
        st.add(f'report.lfq.yield_change.{threshold}',change,fmt='+.2f',desc=f'Percent LFQ precursor yield change at {threshold}',between=(-100,100))
    differences = [r['median_absolute_ratio_difference'] for r in q['shared_ratio_pairs'].values()]
    for bound,value in (('min',min(differences)),('max',max(differences))):
        st.add(f'report.lfq.shared_difference.{bound}',value,fmt='.4f',desc=f'{bound} species median absolute inter-engine log2 ratio difference',sign='+')
    for row in q['engines']:
        e = row['engine']
        for key,value in row['thresholds']['0.01'].items():
            st.add(f'report.lfq.{e}.{key}',value,fmt=',',desc=f'{e} LFQ primary {key}',sign='+')
        for key,value in row['control'].items():
            if key not in ('quantified','foreign'):
                continue
            st.add(f'report.lfq.{e}.control.{key}',value,fmt=',',desc=f'{e} human-only control {key}',between=(0,1000000))
        fraction=100*row['control']['foreign_without_strict_ms2']/row['control']['without_strict_ms2']
        st.add(f'report.lfq.{e}.control.foreign_percent',fraction,fmt='.1f',desc=f'{e} percent foreign among control rows without strict MS2 evidence',between=(0,100))
        for species,metrics in row['species'].items():
            st.add(f'report.lfq.{e}.{species}.bias',metrics['median_log2_bias'],fmt='+.3f',desc=f'{e} {species} log2 B/A bias',between=(-10,10))
    for sp,row in q['shared_ratio_pairs'].items():
        st.add(f'report.lfq.shared.{sp}',row['pairs'],fmt=',',desc=f'{sp} paired ratios shared across releases',sign='+')
    for e in ('upstream','plus'):
        rows={r['threads']:r for r in scale_summary() if r['engine']==e}
        st.add(f'report.scaling.{e}.speedup',rows[1]['seconds']/rows[8]['seconds'],fmt='.2f',desc=f'{e} wall time speedup from one to eight workers',sign='+')
        st.add(f'report.scaling.{e}.eight_seconds',rows[8]['seconds'],fmt='.2f',desc=f'{e} eight worker median seconds in the extension',sign='+')

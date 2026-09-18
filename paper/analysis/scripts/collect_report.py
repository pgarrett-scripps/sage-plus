"""Derive auditable report extensions without changing frozen pilot analyses."""
import argparse
import csv
import hashlib
import json
import math
import statistics as stats
import sys
from collections import Counter, defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
REPO = PAPER.parent
sys.path.insert(0, str(REPO / 'benchmarks'))
from scientific_metrics import read_table, normalize_psms, canonical_peptide, boolean, stripped, quantification_metrics
from scientific_diagnostics import lfq_threshold_yields
from analyze_scientific import accepted_ms2_keys, reference_membership
from generate_ptm_library_benchmark import fasta_entries

ROOT = Path('/data/sage-plus-scientific/20260914')
EXT = Path('/data/sage-plus-scientific/report-extension-20260915')
OUT = PAPER / 'analysis/data/report-extension'
INPUTS = {}


def track(path):
    path = Path(path)
    INPUTS[str(path)] = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    return path


def load(path):
    return json.loads(track(path).read_text())


def psms(directory, engine):
    source = directory / ('results.sage.tsv' if engine == 'upstream' else 'results.sage.parquet')
    return normalize_psms(read_table(track(source)))


def public():
    result = []
    pilot = load(REPO / 'benchmarks/scientific-results/20260914/pilot-summary.json')
    for pair in pilot['public_identifications']:
        if pair['status'] != 'complete':
            continue
        directory = ROOT / 'runs' / pair['pair']
        entry = {'pair': pair['pair'], 'overlap': pair, 'engines': {}}
        rows = {}
        for engine in ('upstream', 'plus'):
            path = directory.with_name(directory.name + '-' + engine)
            rows[engine] = psms(path, engine)
            entry['engines'][engine] = load(path / 'identification-metrics.json')
        accepted = {e: {r['key']: r for r in rs if not r['is_decoy'] and r['spectrum_q'] <= .01} for e, rs in rows.items()}
        counts = {}
        for engine, other in (('upstream', 'plus'), ('plus', 'upstream')):
            all_other = {r['key']: r for r in rows[other]}
            other_spectra = {r['key'][:3]: r for r in rows[other]}
            c = Counter()
            for key in accepted[engine].keys() - accepted[other].keys():
                if key in all_other:
                    c['same_assignment_above_threshold'] += 1
                elif key[:3] in other_spectra:
                    c['different_assignment'] += 1
                else:
                    c['no_matching_rank_one_spectrum_charge'] += 1
            assert sum(c.values()) == pair['baseline_only' if engine == 'upstream' else 'candidate_only']
            counts[engine] = dict(c)
        entry['disagreement'] = counts
        result.append(entry)
    return result


def ptm():
    output = []
    for library in (1, 2):
        for engine in ('upstream', 'plus'):
            path = ROOT / f'runs/ptm/library-{library}-{engine}-all'
            rows = psms(path, engine)
            primary = [r for r in rows if not r['is_decoy'] and r['spectrum_q'] <= .01]
            joint = [r for r in primary if r['peptide_q'] <= .01]
            result = load(path / 'result.json')
            output.append(dict(library=library, engine=engine, spectrum_accepted=len(primary), joint_accepted=len(joint), peptide_q_one_fraction=sum(r['peptide_q'] == 1 for r in rows)/len(rows), thresholds=load(path/'identification-metrics.json'), seconds=result['wall_seconds'], rss=result['peak_rss_mib']))
    return output


def lfq():
    references = {name: '\x00'.join(sequence.replace('I', 'L') for _, sequence in fasta_entries(track(ROOT/f'references/{name}.fasta'))) for name in ('human', 'yeast', 'ecoli')}
    protein_species = load(ROOT/'references/hye-species.json')
    membership = {}
    result = []
    expected = {'human': 0, 'yeast': -1, 'ecoli': 2}
    for engine in ('upstream', 'plus'):
        path = EXT/f'lfq-{engine}'
        run = load(path/'result.json')
        assert run['status'] == 'complete'
        strict = accepted_ms2_keys(psms(path, engine))
        source = path/('lfq.tsv' if engine == 'upstream' else 'lfq.parquet')
        original = read_table(track(source))
        rows = []
        for r in original:
            if engine == 'upstream':
                samples = {k: v for k,v in r.items() if 'LFQ_Orbitrap_DDA_' in k}
                assert len(samples) == 5
                converted = [dict(peptide=r['peptide'],charge='',proteins=r['proteins'],q_value=r['q_value'],is_decoy=False,filename=Path(k).name,intensity=v) for k,v in samples.items()]
            else:
                converted = [dict(r, charge='' if r.get('charge') in ('', '-1', None) else str(r['charge']),filename=Path(r['filename']).name)]
            for item in converted:
                key = (item['filename'],canonical_peptide(item['peptide']),str(item['charge']))
                item['ms2_confirmed'] = key in strict
                rows.append(item)
        design = {r['filename']: {'condition': r['filename'].split('Condition_')[1][0], 'preparation': r['filename'].split('Sample_')[1].split('_')[0]} for r in rows if 'Condition_' in r['filename']}
        selected = []
        control = Counter()
        values = defaultdict(dict)
        for r in rows:
            if boolean(r['is_decoy']):
                continue
            seq = stripped(r['peptide']).replace('I','L')
            if seq not in membership:
                membership[seq] = reference_membership(seq, references)
            species = {protein_species.get(p) for p in r['proteins'].split(chr(59))}
            if len(species) != 1 or not species <= expected.keys():
                continue
            sp = next(iter(species))
            if r['filename'] in design and len(membership[seq]) <= 1:
                selected.append(r)
            accepted = float(r['q_value']) <= .01 and bool(r.get('intensity')) and math.isfinite(float(r['intensity'])) and float(r['intensity']) > 0
            if not accepted:
                continue
            if r['filename'] in design and len(membership[seq]) <= 1:
                key = (sp,canonical_peptide(r['peptide']),r['charge'])
                assert r['filename'] not in values[key]
                values[key][r['filename']] = float(r['intensity'])
            if 'DDA_Human_' in r['filename']:
                if sp != 'human' and 'human' in membership[seq]:
                    continue
                control['quantified'] += 1
                control['foreign'] += sp != 'human'
                if not r['ms2_confirmed']:
                    control['without_strict_ms2'] += 1
                    control['foreign_without_strict_ms2'] += sp != 'human'
        metrics = quantification_metrics(selected, design, protein_species, expected)
        ratios = defaultdict(dict)
        for (sp, peptide, charge), observed in values.items():
            for prep in ('Alpha','Beta'):
                pair = {d['condition']: observed[f] for f,d in design.items() if d['preparation'] == prep and f in observed}
                if set(pair) == {'A','B'}:
                    ratios[sp][f'{peptide}|{charge}|{prep}'] = math.log2(pair['B']/pair['A'])
        for sp in expected:
            assert len(ratios[sp]) == metrics['species'][sp]['ratio_pairs']
        result.append(dict(engine=engine, seconds=run['wall_seconds'], rss=run['peak_rss_mib'], thresholds=lfq_threshold_yields(rows), species=metrics['species'], control=dict(control), ratios=ratios))
    paired = {}
    for sp in expected:
        a,b = [r['ratios'][sp] for r in result]
        common = sorted(a.keys() & b.keys())
        paired[sp] = {'pairs':len(common),'median_absolute_error': {e: stats.median(abs(r[k]-expected[sp]) for k in common) if common else None for e,r in (('upstream',a),('plus',b))},'median_absolute_ratio_difference':stats.median(abs(a[k]-b[k]) for k in common) if common else None}
    return {'engines':result,'shared_ratio_pairs':paired}


def scaling():
    plan = load(EXT/'plan.json')
    result = []
    for job in plan['jobs']:
        if job['suite'] != 'report-threads':
            continue
        path = EXT/job['id']
        run = load(path/'result.json')
        assert run['status'] == 'complete'
        ids = load(path/'identification-metrics.json')
        result.append({k:job[k] for k in ('id','engine','threads','trial','warmup')} | dict(seconds=run['wall_seconds'],rss=run['peak_rss_mib'],target_psms=ids['0.01']['target_psms']))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('section', choices=('public','ptm','lfq','scaling'))
    args = p.parse_args()
    data = globals()[args.section]()
    sources = [Path(__file__)] + [REPO/'benchmarks'/name for name in ('scientific_metrics.py','scientific_diagnostics.py','analyze_scientific.py','generate_ptm_library_benchmark.py')]
    analysis = {str(path.resolve()):hashlib.file_digest(path.open('rb'),'sha256').hexdigest() for path in sources}
    record = {'data':data,'inputs_sha256':INPUTS,'analysis_sources_sha256':analysis}
    (OUT/f'{args.section}.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(args.section, len(INPUTS), 'input files hashed')


if __name__ == '__main__':
    main()

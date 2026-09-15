"""Verify retained evidence identities and close the report extension snapshot."""
import hashlib
import json
from pathlib import Path

PAPER=Path(__file__).resolve().parents[2]
DATA=PAPER/'analysis/data/report-extension'
ROOT=Path('/data/sage-plus-scientific/report-extension-20260915')


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    status=json.loads((ROOT/'matrix-status.json').read_text())
    assert not status['failures'] and not status['pending']
    assert len(status['completed'])==len(status['expected_jobs'])==34
    plan=json.loads((ROOT/'plan.json').read_text())
    assert plan==json.loads((DATA/'plan.json').read_text())
    checked={}
    for section in ('public','ptm','lfq','scaling'):
        snapshot=json.loads((DATA/f'{section}.json').read_text())
        for name,expected in (snapshot['inputs_sha256'] | snapshot['analysis_sources_sha256']).items():
            actual=checked[name] if name in checked else digest(name)
            checked[name]=actual
            assert actual==expected, f'Changed input: {name}'
    for job in plan['jobs']:
        record=json.loads((ROOT/job['id']/'result.json').read_text())
        assert record['status']=='complete'
        for name,expected in (record['outputs'] | record['signature']['inputs']).items():
            actual = checked[name] if name in checked else digest(name)
            checked[name] = actual
            assert actual==expected, f'Changed search evidence: {name}'
    final={'status':'verified','completed_jobs':len(status['completed']),'source_files_verified':len(checked),'snapshot_sha256':{name:digest(DATA/name) for name in ('public.json','ptm.json','lfq.json','scaling.json','plan.json')},'scope':'The original pilot and post hoc report extension use the frozen released pair. No claim of reserved validation. Search outputs and source identities verified.'}
    (DATA/'verification.json').write_text(json.dumps(final,indent=2,sort_keys=True)+'\n')
    print(json.dumps(final,indent=2))


if __name__=='__main__':
    main()

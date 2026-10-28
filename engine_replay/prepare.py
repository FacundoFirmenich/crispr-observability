"""Stage a new 144-job replay without requiring or running CRISPResso2.

The preparation preserves original frozen metadata separately. It authenticates
all input records and regenerates quality-filtered files in a new work directory.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parent
SCRIPTS = ['run_grid.py', 'run_profiles.py', 'checked_run.py', 'staged_run.py',
           'summarize_grid.py', 'build_findings.py', 'close_manifest.py']
COMMIT = '1fc4b629e5400a48d6569b0ab4b8d1766f33636f'


def sha(path):
    return hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()


def copy(source, destination, records):
    assert source.is_file(), 'Missing dependency: ' + str(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    assert sha(source) == sha(destination)
    records.append({'source': str(source), 'destination': str(destination),
                    'sha256': sha(destination)})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=BASE.parent,
                   help='Repository containing raw_replay/inputs/fastq/')
    p.add_argument('--work', type=Path, required=True,
                   help='New nonexistent output directory, outside the repository')
    a = p.parse_args()
    repo, work = a.repo.resolve(), a.work.resolve()
    assert not work.exists(), '--work must not exist; previous work is never overwritten'
    assert not work.is_relative_to(repo), '--work must be outside the source repository'
    for r in json.loads((BASE/'ORIGINAL_SOURCE_MANIFEST.json').read_text()):
        file = BASE/'original'/Path(r['archive_member']).name
        assert file.stat().st_size == r['size'] and sha(file) == r['sha256']
    for r in json.loads((BASE/'DEPENDENCY_MANIFEST.json').read_text()):
        file = BASE/r['path']
        assert file.stat().st_size == r['bytes'] and sha(file) == r['sha256']
    records, identities = [], []
    work.mkdir(parents=True)
    old = work/'restored/m2/advance_m2'
    new = work/'publication_02/sensitivity'
    for r in json.loads((BASE/'DEPENDENCY_MANIFEST.json').read_text()):
        name = Path(r['path']).relative_to('dependencies/m2')
        copy(BASE/r['path'], old/name, records)
    design = json.loads((old/'DESIGN.json').read_text())
    for run in design['runs']:
        rid = run['run_accession']
        source = repo/'raw_replay/inputs/fastq'/(rid+'.fastq.gz')
        dest = work/'publication_verification/P1/advance_m7/replay/project/evidence/v22/v21/sources/fastq'/source.name
        copy(source, dest, records)
        digest, count = hashlib.sha256(), 0
        with gzip.open(dest, 'rt') as f:
            while header := f.readline():
                seq, plus, quality = f.readline().strip(), f.readline(), f.readline().strip()
                assert header.startswith('@') and plus.startswith('+') and len(seq) == len(quality)
                count += 1
                digest.update((seq+'\n'+quality+'\n').encode())
        assert count == run['raw_reads'] and digest.hexdigest() == run['sequence_quality_stream_sha256']
        identities.append({'run': rid, 'N_raw': count,
                           'sequence_quality_stream_sha256': digest.hexdigest()})
    for source in sorted((BASE/'original').iterdir()):
        copy(source, work/'antecedent/original_sensitivity_metadata'/source.name, records)
    for name in SCRIPTS:
        copy(BASE/'original'/name, new/name, records)
    stage = new/'staged_run.py'
    before = sha(stage)
    text = stage.read_text()
    needle = "'/tmp/mp02_engine_jobs'"
    assert text.count(needle) == 1
    scratch = Path(tempfile.mkdtemp(prefix='mp03_engine_replay_'))/'jobs'
    stage.write_text(text.replace(needle, repr(str(scratch))))
    subprocess.run([sys.executable, str(new/'run_grid.py'), 'prepare'], check=True)
    expected = {x['run']: x for x in json.loads((BASE/'original/INPUT_MANIFEST.json').read_text())}
    actual = json.loads((new/'INPUT_MANIFEST.json').read_text())
    for x in actual:
        original = expected[x['run']]
        for key in ['N_raw', 'Q30_pass', 'raw_sequence_quality_stream_sha256', 'Q30_sequence_quality_stream_sha256']:
            assert x[key] == original[key], (x['run'], key)
    assert len(json.loads((new/'JOBS.json').read_text())) == 126
    result = {'state': 'PREPARATION_PASS_NO_ENGINE_ATTACHED_NO_ENGINE_JOBS_RUN',
              'raw_records': sum(x['N_raw'] for x in identities),
              'quality_pass_records': sum(x['Q30_pass'] for x in actual),
              'OAT_grid_jobs_prepared': 126, 'stock_jobs_awaiting_engine_attachment': 18,
              'identities': identities, 'copied_dependencies': records,
              'engine_commit_required': COMMIT,
              'engine_redistributed': False,
              'stock_defaults_verified_against_engine': False,
              'staging_patch': {'before_sha256': before, 'after_sha256': sha(stage),
                                'reason': 'Unique scratch root; scientific parameters unchanged'},
              'retrospective_scope': 'Historical outcomes are already known. New freeze timestamps describe only this fresh computational replay; they are not prospective validation.',
              'excluded_original_outputs': 'Original 144 full alignment directories and failed attempts are not included in this compact source package. The original closed manifest is retained as provenance, not as a claim those bytes are present.'}
    (work/'PREPARATION_VALIDATION.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['state', 'raw_records', 'quality_pass_records', 'OAT_grid_jobs_prepared']}))


if __name__ == '__main__':
    main()

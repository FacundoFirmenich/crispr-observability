"""Verify and attach a separately installed official engine to prepared work.

This runs only --version, --help and stock-profile preparation. It starts no
alignment jobs. All scientific run options remain in the original scripts.
"""
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys

COMMIT = '1fc4b629e5400a48d6569b0ab4b8d1766f33636f'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--engine-source', type=Path, required=True)
    p.add_argument('--engine-executable', type=Path, required=True)
    a = p.parse_args()
    work, source, executable = a.work.resolve(), a.engine_source.resolve(), a.engine_executable.resolve()
    preparation = json.loads((work/'PREPARATION_VALIDATION.json').read_text())
    assert preparation['state'] == 'PREPARATION_PASS_NO_ENGINE_ATTACHED_NO_ENGINE_JOBS_RUN'
    assert subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip() == COMMIT
    subprocess.run(['git', '-C', str(source), 'diff', '--exit-code', 'HEAD'], check=True)
    version = subprocess.check_output([str(executable), '--version'], text=True).strip()
    assert re.search(r'(?<![\d.])2\.3\.4(?![\d.])', version)
    tools = work/'_mp02_engine'
    assert not tools.exists(), 'An engine is already attached; preserve that work directory'
    (tools/'env/bin').mkdir(parents=True)
    (tools/'CRISPResso2').symlink_to(source, target_is_directory=True)
    (tools/'env/bin/CRISPResso').symlink_to(executable)
    scripts = work/'publication_02/sensitivity'
    for name in ['run_grid.py', 'run_profiles.py', 'checked_run.py', 'staged_run.py', 'summarize_grid.py']:
        subprocess.run([sys.executable, str(scripts/name), '--help'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(scripts/'run_profiles.py'), 'prepare'], check=True)
    historical = json.loads((work/'antecedent/original_sensitivity_metadata/STOCK_PROFILE.json').read_text())
    fresh = json.loads((scripts/'STOCK_PROFILE.json').read_text())
    assert fresh['upstream_default_definitions'] == historical['upstream_default_definitions']
    subprocess.run([str(executable), '--help'], check=True, stdout=subprocess.DEVNULL)
    result = {'state': 'ENGINE_ATTACHED_VERSION_COMMIT_HELP_STOCK_DEFAULTS_PASS_NO_ALIGNMENT_JOBS_RUN',
              'engine_version': version, 'engine_commit': COMMIT,
              'source': str(source), 'executable': str(executable),
              'stock_default_definitions_equal_original': True,
              'engine_redistributed': False, 'alignment_jobs_run_by_attachment': 0}
    (work/'ENGINE_ATTACHMENT_VALIDATION.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

"""Check the declared public inventory and its text-withholding boundary."""
from pathlib import Path
import hashlib
import json
import re
import sys
import zipfile
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    inventory = json.loads((ROOT / 'metadata/REPOSITORY_MANIFEST.json').read_text())
    policy = json.loads((ROOT / 'metadata/PUBLIC_DISTRIBUTION_POLICY.json').read_text())
    assert policy['manuscript_public'] is False and policy['supplement_public'] is False
    assert policy['current_acknowledgment_person_count'] == 5
    manifest = inventory['files']
    tracked_scope = 'Declared source snapshot; no matching Git checkout available.'
    try:
        checkout = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', '--show-toplevel'],
                                  capture_output=True, text=True, check=False)
    except FileNotFoundError:
        checkout = None
    if checkout is not None and checkout.returncode == 0 and Path(checkout.stdout.strip()).resolve() == ROOT:
        tracked = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '-z'],
                                 capture_output=True, check=True)
        tracked_paths = set(tracked.stdout.decode('utf-8').rstrip('\0').split('\0'))
        declared = set(manifest) | {'metadata/REPOSITORY_MANIFEST.json'}
        assert tracked_paths == declared, {'unlisted_tracked_files': sorted(tracked_paths - declared),
                                           'declared_but_untracked_files': sorted(declared - tracked_paths)}
        tracked_scope = 'Exact matching Git checkout: all tracked files equal the declared inventory.'
    for name, expected in manifest.items():
        relative = Path(name)
        assert not relative.is_absolute() and '..' not in relative.parts
        assert not any(part in {'manuscript', 'private_current', 'private_prechange',
                                'source_snapshots', 'source_baseline', 'authoring_MP04',
                                'node_modules'} for part in relative.parts)
        assert relative.suffix.lower() not in {'.tex', '.docx', '.bundle'}
        if relative.suffix.lower() == '.zip':
            assert re.fullmatch(r'raw_replay/inputs/alleles/(?:mutant|programmed)/SRR[0-9]+\.zip', name)
            with zipfile.ZipFile(ROOT / name) as archive:
                assert archive.namelist() == ['Alleles_frequency_table.txt']
        assert not (relative.suffix.lower() == '.pdf' and not name.startswith('figures/'))
        data = (ROOT / name).read_bytes()
        assert len(data) == expected['bytes'], name
        assert hashlib.sha256(data).hexdigest() == expected['sha256'], name
        if relative.suffix in {'.json', '.html', '.md', '.txt', '.csv', '.yml', '.svg', '.js'}:
            text = data.decode('utf-8')
            assert '\\begin{document}' not in text, name
            assert 'P1_Revised_Integrated_EN.txt' not in text, name
    release = json.loads((ROOT / 'metadata/RELEASE_STATE.json').read_text())
    assert release['version'] == inventory['version'] == '0.4.0'
    assert release['edition'] == inventory['edition'] == policy['edition'] == 'MP04'
    assert release['manuscript_public'] is False and release['supplement_public'] is False
    # Existing tools and scientific observations are frozen; PaperGraph is a
    # separate communication artifact, not a rewrite of historical analyses.
    preserved = {
        'interactive/epistemic_graph.html':
            '29fc7fce6805bb50b4a6b2c8340d745d07552bad4b8fa7d5b359a7d35bad6445',
        'interactive/supergraph/crispr_supergraph.html':
            '19ae4eac0eb949faeea90027be2151a0257613766ea3370b2390d6f9dd75c31c'}
    for name, expected in preserved.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    model = json.loads((ROOT / 'interactive/papergraph/papergraph.json').read_text())
    def inspect(value):
        if isinstance(value, dict):
            assert not set(value) & {'literalAnchor', 'originalPath', 'sectionCoverage',
                                     'integrationRecommendation'}
            for child in value.values():
                inspect(child)
        elif isinstance(value, list):
            for child in value:
                inspect(child)
    inspect(model)
    for source in model['sources']:
        if source['id'] in {'article', 'supplement'}:
            assert source['access'] == 'WITHHELD_PRIVATE_SOURCE_IDENTITY'
            assert 'path' not in source
    figure_inventory = json.loads((ROOT / 'figures/FIGURE_PACKAGE_MANIFEST.json').read_text())
    for name, expected in figure_inventory['files'].items():
        data = (ROOT / 'figures' / name).read_bytes()
        assert len(data) == expected['size_bytes'], name
        assert hashlib.sha256(data).hexdigest() == expected['sha256'], name
    print(json.dumps({'status': 'PASS', 'inventoried_files_verified': len(manifest),
                      'manuscript_and_formal_supplement_withheld': True,
                      'papergraph_private_text_fields_withheld': True,
                      'existing_interactive_tools_preserved': len(preserved),
                      'tracked_file_scope': tracked_scope,
                      'scope': 'Current declared source distribution; does not erase prior publication or caches.'}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

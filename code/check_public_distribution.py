"""Check the declared public inventory and its text-withholding boundary."""
from pathlib import Path
import hashlib
import os
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
    # Every declared frozen replay source must actually be distributed; a clean
    # checkout cannot rely on ignored local files with historical results/ names.
    reproduction = json.loads((ROOT / 'extensions/reproduction/SOURCE_MANIFEST.json').read_text())
    for name, source in reproduction['files'].items():
        key = 'extensions/reproduction/' + name
        assert key in manifest, 'Missing distributed scientific input: ' + key
        assert manifest[key]['bytes'] == source['bytes'] and manifest[key]['sha256'] == source['sha256'], key
    frontier = json.loads((ROOT / 'extensions/correspondence_frontier/INPUT_MANIFEST.json').read_text())
    for name, source in frontier['files'].items():
        key = 'extensions/correspondence_frontier/' + name
        assert key in manifest and manifest[key] == source, key
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
        if relative.suffix.lower() == '.zip' and name == 'validation/rol-v03-inputs.zip':
            with zipfile.ZipFile(ROOT / name) as archive:
                for member in archive.namelist():
                    part = Path(member)
                    assert not part.is_absolute() and '..' not in part.parts and '\\' not in member, member
                    assert part.suffix.lower() in {'.bam', '.bai', '.sam', '.json', '.html', '.mjs', '.js', '.cif', '.sha256'}, member
        elif relative.suffix.lower() == '.zip':
            assert (re.fullmatch(r'raw_replay/inputs/alleles/(?:mutant|programmed)/SRR[0-9]+\.zip', name) or re.fullmatch(r'extensions/reproduction/external/inputs/(?:nhej|allele_specific|base_editor)_Alleles_frequency_table\.zip', name) or re.fullmatch(r'extensions/correspondence_frontier/inputs/(?:cftr/runs/(?:M_rc|P_rc|MP_fw|MP_rc|PE_zero)|PE_default|hek3/runs/(?:M_fw|P_fw|M_rc|P_rc))/CRISPResso_on_SRR[0-9]+/Alleles_frequency_table\.zip', name))
            with zipfile.ZipFile(ROOT / name) as archive:
                assert archive.namelist() == ['Alleles_frequency_table.txt']
        assert not (relative.suffix.lower() == '.pdf' and not name.startswith('figures/'))
        path = ROOT / name
        if os.name == "nt" and len(str(path)) >= 260:
            path = Path(chr(92)*2 + "?" + chr(92) + str(path))
        data = path.read_bytes()
        assert len(data) == expected['bytes'], name
        assert hashlib.sha256(data).hexdigest() == expected['sha256'], name
        if relative.suffix in {'.json', '.html', '.md', '.txt', '.csv', '.yml', '.svg', '.js'}:
            # The text markers are ASCII. Scan exact bytes so archived Windows
            # code-page files can be checked without rewriting their evidence.
            assert b'\\begin{document}' not in data, name
            assert b'P1_Revised_Integrated_EN.txt' not in data, name
    release = json.loads((ROOT / 'metadata/RELEASE_STATE.json').read_text())
    assert release['version'] == inventory['version'] == '0.9.0'
    assert release['edition'] == inventory['edition'] == policy['edition'] == 'NATIVE_REPRODUCTION'
    assert release['manuscript_public'] is False and release['supplement_public'] is False
    license_scope = json.loads((ROOT / 'metadata/LICENSE_SCOPE.json').read_text())
    assert license_scope['original_code_license'] == 'GPL-3.0-only'
    zenodo = json.loads((ROOT / '.zenodo.json').read_text())
    assert zenodo['license'] == 'gpl-3.0-only' and zenodo['version'] == release['version']
    assert not zenodo.get('doi')
    # Verify the current navigation files without changing frozen scientific observations.
    preserved = {
        'interactive/epistemic_graph.html':
            'd7630fb7f3ff4f74cca55a91dcf41e40c8825ffb394bea130946ebb162206a57',
        'interactive/supergraph/crispr_supergraph.html':
            (json.loads((ROOT / 'interactive/supergraph/BUILD_RECEIPT.json').read_bytes())['html_sha256'])}
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

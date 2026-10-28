"""Refresh scoped asset notices and the exact public source inventory."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = 'metadata/REPOSITORY_MANIFEST.json'


def public_files():
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if any(part in {'.git', '__pycache__', 'node_modules', 'results', '.qa'} for part in relative.parts):
            continue
        if path.suffix in {'.pyc', '.stdout', '.stderr'}:
            continue
        yield relative.as_posix(), path


def main():
    scope = json.loads((ROOT / 'metadata/LICENSE_SCOPE.json').read_text())
    code = set(scope['exact_original_code_allowlist'])
    paths = dict(public_files())
    assert code <= set(paths), sorted(code - set(paths))
    notice = ROOT / 'THIRD_PARTY_NOTICES.csv'
    fields = ['asset', 'bundled', 'asset_kind', 'attribution', 'rights_holder_status',
              'source_or_provenance', 'applicable_terms', 'license_assigned_by_candidate',
              'third_party_relicensing']
    with notice.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for name in sorted(paths):
            original = name in code
            writer.writerow({'asset': name, 'bundled': 'yes',
                'asset_kind': 'ORIGINAL_CODE' if original else 'SCIENTIFIC_DATA_FIGURE_OR_METADATA',
                'attribution': 'Facundo Firmenich programme; original CFTR experiments Sousa et al.; upstream sources retain attribution',
                'rights_holder_status': 'SCOPED_ORIGINAL_CODE' if original else 'SOURCE_AND_DATA_TERMS_RETAINED',
                'source_or_provenance': 'Current computational companion with source identities retained; independently implemented interactive interface, no KwanTube source code copied',
                'applicable_terms': 'MIT scoped by metadata/LICENSE_SCOPE.json' if original else 'Outside software MIT grant; source attribution and applicable data/text/figure terms retained',
                'license_assigned_by_candidate': 'MIT_SCOPED' if original else 'NONE',
                'third_party_relicensing': 'no'})
    files = {name: {'bytes': path.stat().st_size,
                    'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
             for name, path in public_files() if name != INVENTORY}
    data = {'edition': 'MP03D', 'version': '0.3.3',
            'scope': 'Exact current public computational companion source. Complete manuscript and formal supplement remain private. Inventory excludes itself and runtime results/screenshots/dependencies.',
            'files': files}
    (ROOT / INVENTORY).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'REFRESHED', 'public_files_including_inventory': len(files)+1}))


if __name__ == '__main__':
    main()

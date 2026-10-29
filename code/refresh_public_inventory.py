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
    gpl_code = set(scope.get('exact_gpl_code_allowlist', []))
    paths = dict(public_files())
    assert code <= set(paths), sorted(code - set(paths))
    assert gpl_code <= set(paths), sorted(gpl_code - set(paths))
    assert not (code & gpl_code), sorted(code & gpl_code)
    notice = ROOT / 'THIRD_PARTY_NOTICES.csv'
    fields = ['asset', 'bundled', 'asset_kind', 'attribution', 'rights_holder_status',
              'source_or_provenance', 'applicable_terms', 'license_assigned_by_candidate',
              'third_party_relicensing']
    with notice.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for name in sorted(paths):
            original = name in code
            gpl = name in gpl_code
            gpl_html = name == 'interactive/supergraph/crispr_supergraph.html'
            gpl_license = name == 'interactive/supergraph/LICENSE'
            mit_html = name in {'interactive/epistemic_graph.html',
                                'interactive/papergraph/PaperGraph_CRISPR.html'}
            if gpl:
                kind = 'GPL_COMPONENT_CODE'
                rights = 'SCOPED_GPL_COMPONENT'
                terms = 'GPL-3.0-only; complete license and modification notice in interactive/supergraph/'
                assigned = 'GPL-3.0-only'
            elif gpl_html:
                kind = 'MIXED_HTML_GPL_CODE_AND_SCIENTIFIC_DATA'
                rights = 'CODE_GPL_DATA_TERMS_RETAINED'
                terms = 'Interface code GPL-3.0-only; embedded scientific data and claims retain separate source terms'
                assigned = 'GPL-3.0-only_CODE_ONLY'
            elif gpl_license:
                kind = 'LICENSE_DOCUMENT'
                rights = 'FREE_SOFTWARE_FOUNDATION_LICENSE_TEXT'
                terms = 'Verbatim GNU GPL version 3 license text; copyright Free Software Foundation'
                assigned = 'NONE'
            elif original:
                kind, rights = 'ORIGINAL_CODE', 'SCOPED_ORIGINAL_CODE'
                terms, assigned = 'MIT scoped by metadata/LICENSE_SCOPE.json', 'MIT_SCOPED'
            elif mit_html:
                kind, rights = 'MIXED_HTML_MIT_CODE_AND_SCIENTIFIC_DATA', 'CODE_MIT_DATA_TERMS_RETAINED'
                terms, assigned = 'Original interface code MIT; embedded scientific data retain separate source terms', 'MIT_CODE_ONLY'
            else:
                kind, rights = 'SCIENTIFIC_DATA_FIGURE_OR_METADATA', 'SOURCE_AND_DATA_TERMS_RETAINED'
                terms = 'Outside software code grants; source attribution and applicable data/text/figure terms retained'
                assigned = 'NONE'
            spatial = name.startswith('interactive/supergraph/') or gpl
            papergraph = name.startswith(('interactive/papergraph/', 'figures/papergraph/')) or name in {
                'code/check_papergraph_browser.cjs', '.github/workflows/papergraph-checks.yml'}
            provenance = ('New separate spatial interface adapts the actual published KwanTube Beyond Redfield companion; exact upstream identity and modifications in interactive/supergraph/' if spatial else
                          'Independent PaperGraph interface and renderer; public authored claim summaries, typed relations and existing numerical source identities; formal manuscript and supplement withheld' if papergraph else
                          'Current computational companion; the existing white interactive interface was independently implemented and remains unchanged')
            writer.writerow({'asset': name, 'bundled': 'yes',
                'asset_kind': kind,
                'attribution': 'Facundo Firmenich programme; original CFTR experiments Sousa et al.; upstream sources retain attribution',
                'rights_holder_status': rights,
                'source_or_provenance': provenance,
                'applicable_terms': terms,
                'license_assigned_by_candidate': assigned,
                'third_party_relicensing': 'no'})
    # The figure inventory includes the original four figures and the separate
    # PaperGraph figures. Neither is a full manuscript distribution.
    figure_manifest = ROOT / 'figures/FIGURE_PACKAGE_MANIFEST.json'
    figure_files = {name.removeprefix('figures/'): {
                       'size_bytes': path.stat().st_size,
                       'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                    for name, path in public_files()
                    if name.startswith('figures/') and path != figure_manifest}
    figure_manifest.write_text(json.dumps({
        'edition': 'MP04', 'files': figure_files,
        'full_manuscript_or_supplement_included': False,
        'coordinate_contract': 'sources/anchor_geometry_contract.json',
        'papergraph_source': '../interactive/papergraph/papergraph.json'
    }, ensure_ascii=False, indent=2) + '\n')
    # Fingerprint after refreshing nested inventories.
    files = {name: {'bytes': path.stat().st_size,
                    'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
             for name, path in public_files() if name != INVENTORY}
    data = {'edition': 'MP04', 'version': '0.4.0',
            'scope': 'Exact current public computational companion source. Complete manuscript and formal supplement remain private. Inventory excludes itself and runtime results/screenshots/dependencies.',
            'files': files}
    (ROOT / INVENTORY).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'REFRESHED', 'public_files_including_inventory': len(files)+1}))


if __name__ == '__main__':
    main()

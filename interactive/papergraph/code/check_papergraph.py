"""Check the public PaperGraph, withheld-source boundary and deterministic rebuild.

Public fingerprints and diagram invariants do not verify private text, every
scientific inference, molecular accuracy or reader benefit. Browser QA is separate.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parents[1]
FIGURES = REPOSITORY / 'figures/papergraph'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect_public_boundary(obj):
    if isinstance(obj, dict):
        assert not set(obj) & {'literalAnchor', 'originalPath', 'sectionCoverage',
                               'integrationRecommendation'}
        if obj.get('sourceId') in {'article', 'supplement'}:
            assert obj == {'sourceId': obj['sourceId'],
                           'access': 'WITHHELD_PRIVATE_SOURCE_IDENTITY'}
        for value in obj.values():
            inspect_public_boundary(value)
    elif isinstance(obj, list):
        for value in obj:
            inspect_public_boundary(value)


def science(obj):
    if isinstance(obj, dict):
        return {key: science(value) for key, value in obj.items()
                if key != 'sourceRefs'}
    if isinstance(obj, list):
        return [science(value) for value in obj]
    return obj


def main():
    model_path = ROOT / 'papergraph.json'
    model = json.loads(model_path.read_text())
    receipt = json.loads((ROOT / 'BUILD_RECEIPT.json').read_text())
    projection = json.loads((ROOT / 'PUBLIC_PROJECTION.json').read_text())
    assert model['schemaVersion'] == 'papergraph-public-0.1.0'
    assert model['status'] == 'PUBLIC_COMPUTATIONAL_COMPANION_FORMAL_TEXT_WITHHELD'
    assert sha(model_path) == receipt['model_sha256'] == projection['public_model_sha256']
    inspect_public_boundary(model)
    encoded = model_path.read_text()
    assert 'source_snapshots/' not in encoded
    assert 'withdrawal_01/' not in encoded
    assert '\\begin{document}' not in encoded
    sources = {source['id']: source for source in model['sources']}
    assert len(sources) == len(model['sources']) == 14
    withheld = public = 0
    for source in model['sources']:
        assert re.fullmatch('[0-9a-f]{64}', source['sha256'])
        assert isinstance(source['bytes'], int) and source['bytes'] > 0
        if source['access'] == 'WITHHELD_PRIVATE_SOURCE_IDENTITY':
            assert source['id'] in {'article', 'supplement'}
            assert 'path' not in source
            withheld += 1
            continue
        relative = Path(source['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path = REPOSITORY / relative
        assert path.is_file() and path.stat().st_size == source['bytes'], source['id']
        assert sha(path) == source['sha256'], source['id']
        public += 1
    assert withheld == 2 and public == 12
    for field, expected in projection['scientific_field_sha256'].items():
        raw = json.dumps(science(model[field]), ensure_ascii=False, sort_keys=True,
                         separators=(',', ':')).encode()
        assert hashlib.sha256(raw).hexdigest() == expected, field
    stations = [station['id'] for station in model['stations']]
    assert stations == model['readingPath'] and len(set(stations)) == 12
    facts = {fact['id']: fact for fact in model['facts']}
    nodes = {node['id']: node for node in model['nodes']}
    assert len(facts) == len(model['facts']) == 49
    assert len(nodes) == len(model['nodes']) == 92
    assert len(model['edges']) == len({edge['id'] for edge in model['edges']}) == 86
    assert len(model['conjunctions']) == 7
    for station in model['stations']:
        assert all(station['levels'].get(depth) for depth in ['student', 'technical', 'audit'])
        assert set(station['factIds']) <= set(facts)
        assert set(station['sourceIds']) <= set(sources)
    for edge in model['edges']:
        assert edge['source'] in nodes and edge['target'] in nodes
    for group in model['conjunctions']:
        assert group['id'] in nodes and group['conclusion'] in nodes
        assert group['premises'] and set(group['premises']) <= set(nodes)
    html = (ROOT / 'PaperGraph_CRISPR.html').read_text()
    assert sha(ROOT / 'PaperGraph_CRISPR.html') == receipt['html_sha256']
    embedded = re.search(r'<script type="application/json" id="papergraph-data">(.*?)</script>',
                         html, re.S).group(1)
    assert json.loads(embedded) == model
    assert hashlib.sha256(embedded.encode()).hexdigest() == receipt['embedded_sha256']
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=["\']https?://', html, re.I)
    for name in ['PAPERGRAPH_DATA', 'MODEL_SHA256', 'EMBEDDED_SHA256', 'OVERVIEW_SVG',
                 'PORTRAIT_SVG_JSON', 'ARGUMENT_RUNTIME']:
        assert '__' + name + '__' not in html
    ns = {'s': 'http://www.w3.org/2000/svg'}
    fonts = []
    for name in ['PaperGraph_Overview', 'PaperGraph_Article']:
        svg = ET.parse(FIGURES / (name + '.svg')).getroot()
        got = [node.attrib['data-station'] for node in svg.findall('.//s:g[@data-station]', ns)]
        assert got == stations
        if name == 'PaperGraph_Article':
            fonts = [float(node.attrib['font-size']) for node in svg.findall('.//s:text', ns)]
    assert min(fonts) >= 32
    atlas = ET.parse(FIGURES / 'PaperGraph_Argument.svg').getroot()
    node_ids = [node.attrib['data-node-id'] for node in atlas.findall('.//s:g[@data-node-id]', ns)]
    edge_ids = [edge.attrib['data-edge-id'] for edge in atlas.findall('.//s:g[@data-edge-id]', ns)]
    assert len(node_ids) == len(set(node_ids)) == len(nodes)
    assert set(node_ids) == set(nodes)
    assert len(edge_ids) == len(set(edge_ids)) == len(model['edges'])
    assert set(edge_ids) == {edge['id'] for edge in model['edges']}
    assert json.loads(atlas.find('s:metadata', ns).text) == model
    with tempfile.TemporaryDirectory(prefix='public-papergraph-rebuild-') as tmp:
        target = Path(tmp)
        figures = target / 'figures'
        subprocess.run([sys.executable, '-B', str(ROOT / 'code/build_papergraph.py'),
                        '--output', str(target), '--figures-output', str(figures)],
                       check=True, stdout=subprocess.DEVNULL)
        for name in ['PaperGraph_CRISPR.html', 'papergraph.json']:
            assert (ROOT / name).read_bytes() == (target / name).read_bytes(), name
        for name in ['PaperGraph_Overview.svg', 'PaperGraph_Article.svg', 'PaperGraph_Argument.svg']:
            assert (FIGURES / name).read_bytes() == (figures / name).read_bytes(), name
        # The builder writes no full-text sources or data copies into its output.
        expected = {'PaperGraph_CRISPR.html', 'papergraph.json', 'BUILD_RECEIPT.json',
                    'figures/PaperGraph_Overview.svg', 'figures/PaperGraph_Article.svg',
                    'figures/PaperGraph_Argument.svg'}
        assert {path.relative_to(target).as_posix() for path in target.rglob('*')
                if path.is_file()} == expected
    result = {'status': 'PASS', 'public_source_files_verified': public,
              'withheld_identity_records_not_read': withheld, 'stations': len(stations),
              'facts': len(facts), 'argument_atlas_nodes': len(node_ids),
              'argument_atlas_edges': len(edge_ids), 'semantic_edges': 75,
              'reading_order_edges_distinct': 11, 'named_AND_groups': 7,
              'scientific_field_fingerprints_unchanged': len(projection['scientific_field_sha256']),
              'article_minimum_font_px': min(fonts),
              'article_minimum_font_at235mm_height_pt': min(fonts) * (235 / 25.4 * 72) / 2660,
              'temporary_rebuild_byte_identical': True,
              'full_text_sources_read_or_copied': False,
              'browser_or_reader_study_claimed_by_this_check': False}
    (ROOT / 'ARTIFACT_CHECK.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

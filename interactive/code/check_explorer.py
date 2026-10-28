"""Verify frozen source identities, independent rebuilds and embedded JSON bytes."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    model_path = ROOT / 'evidence_model.json'
    model = json.loads(model_path.read_text())
    canonical_sources = []
    for source in model['sources']:
        relative = Path(source['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Source path must stay inside the public repository')
        canonical = REPOSITORY / relative
        frozen = ROOT / 'sources' / ('tables' if relative.parts[0] == 'tables' else 'numerical') / relative.name
        assert sha(canonical) == sha(frozen) == source['sha256'], source['path']
        canonical_sources.append(source['path'])
    html = (ROOT / 'epistemic_graph.html').read_text()
    embedded = re.search(r'<script type="application/json" id="evidence-data">(.*?)</script>', html, re.S)
    assert embedded, 'Embedded evidence-data script is missing'
    raw = embedded.group(1)
    assert json.loads(raw) == model
    receipt = json.loads((ROOT / 'BUILD_RECEIPT.json').read_text())
    assert receipt['input_json_sha256'] == sha(model_path)
    assert receipt['embedded_json_sha256'] == hashlib.sha256(raw.encode()).hexdigest()
    assert receipt['html_sha256'] == sha(ROOT / 'epistemic_graph.html')
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=["\']https?://', html, re.I)
    assert not re.search(r'/(?:blob|commit|raw)/[0-9a-f]{7,40}(?:/|\b)', html)
    assert '\\begin{document}' not in html
    with tempfile.TemporaryDirectory(prefix='crispr-explorer-check-') as folder:
        target = Path(folder)
        generated_model = target / 'evidence_model.json'
        subprocess.run([sys.executable, '-B', str(ROOT / 'code/build_model.py'),
                        '--output', str(generated_model)], check=True, stdout=subprocess.DEVNULL)
        assert generated_model.read_bytes() == model_path.read_bytes()
        generated_html = target / 'epistemic_graph.html'
        subprocess.run([sys.executable, '-B', str(ROOT / 'code/build_explorer.py'),
                        '--model', str(generated_model), '--output', str(generated_html)],
                       check=True, stdout=subprocess.DEVNULL)
        assert generated_html.read_bytes() == (ROOT / 'epistemic_graph.html').read_bytes()
    result = {'status': 'PASS', 'canonical_source_fingerprints': len(canonical_sources),
              'model_rebuild_byte_identical': True, 'html_rebuild_byte_identical': True,
              'embedded_model_object_equal': True, 'embedded_json_sha256': receipt['embedded_json_sha256'],
              'input_json_sha256': receipt['input_json_sha256'], 'html_sha256': receipt['html_sha256'],
              'browser_rendering_tested_by_this_check': False}
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

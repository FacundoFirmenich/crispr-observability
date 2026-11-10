"""Reconstruct exact large source files from digest-verified Git-tracked parts."""
from pathlib import Path
import argparse,hashlib,json
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
for descriptor in sorted(root.rglob('*.parts.json')):
    item=json.loads(descriptor.read_text());target=descriptor.parent/item['original_name']
    assert target.resolve().is_relative_to(root)
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest()==item['sha256'];continue
    digest=hashlib.sha256();total=0
    for row in item['parts']:
        p=descriptor.parent/row['name'];assert p.resolve().is_relative_to(root)
        b=p.read_bytes();assert len(b)==row['size'] and hashlib.sha256(b).hexdigest()==row['sha256']
        digest.update(b);total+=len(b)
    assert total==item['size'] and digest.hexdigest()==item['sha256']
    if not args.check:
        with target.open('xb') as out:
            for row in item['parts']:out.write((descriptor.parent/row['name']).read_bytes())
    print('VERIFIED' if args.check else 'RESTORED',target.relative_to(root),round(total/1e6,2),'MB')

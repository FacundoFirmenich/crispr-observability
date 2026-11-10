"""Known adverse panel, explicitly not held-out validation."""
import json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from crispr_extended import certify
cases=json.loads((R/'tests/fixtures/persistent_source_flanks.json').read_text())
for case in cases:
    q=case['request'];r=certify(q['read'],q['target'],q['haplotypes'],q['budgets'])
    for d in r['decisions']:
        if d['budget']>=case['literal_errors']:
            assert case['expected_origin_label'] in d['possible_observed_labels']
            assert case['source_allele'] in d['possible_source_alleles']
            assert d['observed_tract_certificate']=='U'
print('PASS',len(cases),'known adverse records')

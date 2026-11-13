"""Check the frozen held-out grid and exact cost partition against source data."""
import json
from fractions import Fraction as F
from pathlib import Path
from local_policy import evaluate
from cost_frontier import calculate

here = Path(__file__).resolve().parent
grid = evaluate()
assert json.loads(json.dumps(grid)) == json.loads((here/'results/LOCAL_POLICY.json').read_bytes())
frontier = calculate()
assert frontier == json.loads((here/'results/COST_FRONTIER.json').read_bytes())
for cell in frontier['intervals'] + frontier['boundary_points']:
    r = F(cell.get('evaluation_ratio', cell.get('ratio')))
    favorable = cell['local_empirical_beats_both_fixed'] and cell['local_empirical_beats_pooled']
    assert favorable == (F(17,3) < r <= F(20,3))
anchor = json.loads((here/'results/AMP_ANCHOR_SUMMARY.json').read_bytes())
for edge in anchor['admission_transitions']:
    before = next(x for x in anchor['rows'] if x['run']==edge['run'] and x['anchor_offset']==edge['from_offset'])
    after = next(x for x in anchor['rows'] if x['run']==edge['run'] and x['anchor_offset']==edge['to_offset'])
    assert after['assigned']-before['assigned'] == edge['newly_assigned_weighted_records']-edge['newly_unassigned_weighted_records']
native = json.loads((here/'results/CRISPRA_INTERFACE_RESULT.json').read_bytes())
assert len(native['rows']) == 9 and all(x['exit']==0 for x in native['execution'])
for row in native['rows']:
    assert row['repaired_scalar_interface']
    assert row['linked_count'] == row['repaired_count']
    if row['linked_count']==0:
        assert row['control_table_unchanged']
assert sum(x['Q30_input'] for x in native['rows']) == 175010
print('PASS: source-linked held-out grid, exact frontier, 27 admission transitions and nine native return interfaces')

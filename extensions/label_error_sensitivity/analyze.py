"""Exact single-target sensitivity to bounded accepted-label errors.

Uses frozen actual-record counts, not synthetic biological data. Budgets are
hypothetical deterministic allowances, never estimated error rates. Does not
modify MP06 data, calibrate errors, rerun an aligner or estimate accuracy.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import csv
import hashlib
import itertools
import json

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'extensions/correspondence_frontier/expected/mechanism/FAMILY_CERTIFICATE_COUNTS.csv'
RATES = tuple(map(F, ('0', '1/1000', '1/200', '1/100', '1/50', '1/20', '1/10', '1/5', '1')))


def bounds(e, ne, u, a, b):
    """At most a false accepted events and b false accepted non-events."""
    assert min(e, ne, u, a, b) >= 0 and a <= e and b <= ne
    n = e + ne + u
    assert n > 0
    return F(e-a, n), F(e+u+b, n)


def floor_product(rate, count):
    return rate.numerator * count // rate.denominator


def contrast(t, c, rate):
    tl, th = bounds(t['E'], t['NE'], t['U'], floor_product(rate, t['E']), floor_product(rate, t['NE']))
    cl, ch = bounds(c['E'], c['NE'], c['U'], floor_product(rate, c['E']), floor_product(rate, c['NE']))
    return tl-ch, th-cl


def critical_rate(t, c):
    # Lower endpoints change only at these exact breakpoints. This is a shared
    # cap fraction within each accepted class, not the observed error fraction.
    cuts = {F(0), F(1)}
    for cap in (t['E'], c['NE']):
        if cap:
            cuts.update(F(k, cap) for k in range(cap+1))
    previous = None
    for rate in sorted(cuts):
        low, _ = contrast(t, c, rate)
        if low <= 0:
            if rate > 0:
                assert previous is not None and contrast(t, c, previous)[0] > 0
            return rate, low, previous
        previous = rate
    raise AssertionError('A cap of one must permit a nonpositive contrast')


def check_sharpness():
    # Exhaustive logical label assignments. These are mathematical objects,
    # never biological observations. Check all count triples with 1 <= N <= 7.
    cases = 0
    for n in range(1, 8):
        for e in range(n+1):
            for ne in range(n-e+1):
                u = n-e-ne
                truth = tuple(itertools.product((0, 1), repeat=n))
                for a in range(e+1):
                    for b in range(ne+1):
                        attainable = {F(sum(x), n) for x in truth
                                      if e-sum(x[:e]) <= a and sum(x[e:e+ne]) <= b}
                        lo, hi = bounds(e, ne, u, a, b)
                        assert min(attainable) == lo and max(attainable) == hi
                        cases += 1
    return cases


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    out = args.output.resolve()
    assert not out.exists() and not out.is_relative_to(ROOT)
    # Authenticate this frozen source against the existing committed manifest.
    record = json.loads((ROOT/'metadata/REPOSITORY_MANIFEST.json').read_text())['files'][str(SOURCE.relative_to(ROOT))]
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert SOURCE.stat().st_size == record['bytes'] and digest == record['sha256']
    with SOURCE.open(newline='') as f:
        rows = list(csv.DictReader(f))
    groups = {}
    for row in rows:
        x = {**row, **{k: int(row[k]) for k in ('N', 'E', 'NE', 'U')}}
        assert x['E']+x['NE']+x['U'] == x['N']
        assert bounds(x['E'], x['NE'], x['U'], 0, 0) == (F(x['L_exact']), F(x['H_exact']))
        key = tuple(x[k] for k in ('preparation', 'family', 'certificate', 'event'))
        assert x['condition'] not in groups.setdefault(key, {})
        groups[key][x['condition']] = x
    curves, tips = [], []
    for key, arms in sorted(groups.items()):
        assert set(arms) == {'untreated', 'mock', 'PE treated'}, set(arms)
        t, c = arms['PE treated'], arms['mock']
        common = dict(zip(('preparation', 'family', 'certificate', 'event'), key))
        common.update(treated_run=t['run'], mock_run=c['run'])
        old = None
        for rate in RATES:
            lo, hi = contrast(t, c, rate)
            if old:
                assert lo <= old[0] and hi >= old[1]
            old = lo, hi
            curves.append({**common, 'class_cap': str(rate),
                           'a_t': floor_product(rate, t['E']), 'b_t': floor_product(rate, t['NE']),
                           'a_c': floor_product(rate, c['E']), 'b_c': floor_product(rate, c['NE']),
                           'lower_exact': str(lo), 'upper_exact': str(hi),
                           'lower_pp': f'{100*float(lo):.8f}', 'upper_pp': f'{100*float(hi):.8f}'})
        rate, low, prev = critical_rate(t, c)
        tips.append({**common, 'first_nonpositive_class_cap': str(rate),
                     'cap_percent': f'{100*float(rate):.8f}',
                     'a_t': floor_product(rate, t['E']), 'b_c': floor_product(rate, c['NE']),
                     'lower_at_cap_exact': str(low),
                     'previous_breakpoint': str(prev) if prev is not None else '',
                     'lower_at_previous_exact': str(contrast(t,c,prev)[0]) if prev is not None else '',
                     'baseline_lower_exact': str(contrast(t,c,F(0))[0])})
    logical_cases = check_sharpness()
    out.mkdir(parents=True)
    for name, data in [('SENSITIVITY.csv', curves), ('TIPPING_CAPS.csv', tips)]:
        with (out/name).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader(); writer.writerows(data)
    receipt = {'status': 'PASS_CONDITIONAL_SENSITIVITY_ONLY',
               'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': digest,
               'source_count_rows': len(rows), 'contrasts': len(groups),
               'grid_rows': len(curves), 'sharpness_logical_cases': logical_cases,
               'new_biological_data': False, 'alignment_engine_rerun': False,
               'error_rates_calibrated': False, 'joint_target_hulls_claimed': False,
               'assumptions': ['Each target is evaluated separately.',
                  'At most floor(r E) false events and floor(r NE) false non-events in each arm.',
                  'Unknown labels are unrestricted; no cross-record coupling.',
                  'Allowed budgets are not estimated probabilities or confidence levels.']}
    (out/'RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

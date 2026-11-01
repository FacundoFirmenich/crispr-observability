from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib

BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = dict(zip((''.join(c) for c in product(BASES, repeat=3)), AA))
NATIVE = 'GTIKENIIFGVSYD'
PROPERTIES = ('observed_coding', 'peptide_or_frame_status', 'coding_length', 'native_event', 'exact_programmed_event')
TARGETS = ('observed_coding', 'peptide_or_frame_status', 'native_event', 'exact_programmed_event')

class ScopeError(ValueError):
    pass

@lru_cache(maxsize=None)
def peptide(coding):
    if not coding or set(coding) - set('ACGT'):
        raise ValueError('A resolved coding string must be nonempty ACGT')
    if len(coding) % 3:
        return None
    return ''.join(CODE[coding[i:i+3]] for i in range(0, len(coding), 3))

def property_value(coding, name, exact_coding):
    if name == 'observed_coding':
        return coding
    if name == 'peptide_or_frame_status':
        p = peptide(coding)
        return 'PEPTIDE:' + p if p is not None else 'OUT_OF_FRAME:' + str(len(coding) % 3)
    if name == 'coding_length':
        return len(coding)
    if name == 'native_event':
        return peptide(coding) == NATIVE
    if name == 'exact_programmed_event':
        return coding == exact_coding
    raise ValueError('Unknown declared property: ' + name)

def factorization(catalog, representation, target, exact_coding):
    """Certificate on a supplied finite catalog, including actual counterexamples.

    Catalog entries carry coding, multiplicity and a real record locator.
    Counts describe this catalog only; homogeneous observed fibers do not prove
    sufficiency outside the observed support.
    """
    fibers = defaultdict(dict)
    total = 0
    for row in catalog:
        c, n = row['coding'], row['count']
        if type(n) is not int or n <= 0:
            raise ValueError('Positive integer multiplicity required')
        r = property_value(c, representation, exact_coding)
        y = property_value(c, target, exact_coding)
        cell = fibers[r].setdefault(y, {'count': 0, 'example': row['example']})
        cell['count'] += n
        total += n
    bad = [(r, ys) for r, ys in fibers.items() if len(ys) > 1]
    examples = []
    for r, ys in sorted(bad, key=lambda x: repr(x[0]))[:8]:
        entries = list(ys.items())
        examples.append({'representation_value': r, 'target_values': [entries[0][0], entries[1][0]], 'actual_records': [entries[0][1]['example'], entries[1][1]['example']]})
    return {'representation': representation, 'target': target, 'resolved_records': total, 'distinct_coding_strings': len(catalog), 'observed_fibers': len(fibers), 'heterogeneous_fibers': len(bad), 'records_in_heterogeneous_fibers': sum(sum(v['count'] for v in ys.values()) for _, ys in bad), 'factors_on_observed_catalog': not bad, 'counterexamples': examples, 'authority': 'Finite supplied observed catalog only; no out-of-support sufficiency or accuracy claim.'}

def require_scope(requested, certificate):
    if requested != 'finite_observed_catalog':
        raise ScopeError('This certificate does not authenticate full alleles, function, biological truth or future sufficiency')
    if not certificate['factors_on_observed_catalog']:
        raise ScopeError('The requested target is not constant in every observed representation fiber')
    return True

def compile_partition(catalog, targets, exact_coding):
    """Coarsest observed partition preserving every requested target.

    Distinct target signatures must remain separate in any sufficient encoding.
    This classical factorization construction is applied to actual coding rows;
    it does not authenticate biological truth or generalize beyond the catalog.
    """
    if not targets or len(set(targets)) != len(targets):
        raise ValueError('Declare a nonempty sequence of distinct targets')
    cells = {}
    for r in catalog:
        signature = tuple(property_value(r['coding'],t,exact_coding) for t in targets)
        g = cells.setdefault(signature, {'target_signature':list(signature),'records':0,'distinct_coding_strings':0,'actual_example':r['example']})
        g['records'] += r['count']
        g['distinct_coding_strings'] += 1
    return {'targets':list(targets),'minimum_observed_classes':len(cells),'resolved_records':sum(g['records'] for g in cells.values()),'classes':sorted(cells.values(),key=lambda g:repr(g['target_signature'])),'scope':'Coarsest partition on this observed catalog for these targets only; not a state ontology or full-sequence/function certificate.'}

def grouped_withdrawal(t_catalog, c_catalog, t_counts, c_counts, event, exact_coding):
    """Minimum number of whole identical-coding families to remove positivity.

    Families span both arms. All resolved records of a selected family are
    withdrawn, including labels neutral to the lower endpoint. This is a
    deterministic sensitivity convention, not a measured fault mechanism.
    """
    nt, nc = t_counts['N'], c_counts['N']
    lo = Fraction(t_counts['E'], nt) - Fraction(c_counts['E'] + c_counts['U'], nc)
    if lo <= 0:
        raise ValueError('A positive original certificate is required')
    ts = {r['coding']: r['count'] for r in t_catalog}
    cs = {r['coding']: r['count'] for r in c_catalog}
    groups = []
    for s in sorted(ts.keys() | cs.keys()):
        e = property_value(s, event, exact_coding)
        impact = Fraction(ts.get(s, 0), nt) if e else Fraction(cs.get(s, 0), nc)
        if impact:
            groups.append({'coding': s, 'coding_sha256': hashlib.sha256(s.encode()).hexdigest(), 'treated_resolved': ts.get(s, 0), 'control_resolved': cs.get(s, 0), 'event': e, 'lower_bound_loss_exact': str(impact)})
    groups.sort(key=lambda g: (-Fraction(g['lower_bound_loss_exact']), g['coding']))
    loss = Fraction(0)
    chosen = []
    for g in groups:
        chosen.append(g)
        loss += Fraction(g['lower_bound_loss_exact'])
        if loss >= lo:
            break
    if loss < lo:
        raise AssertionError('Available full withdrawal must remove the original certificate')
    k = len(chosen)
    # When one family suffices, select the least-record witness among all
    # sufficient singleton families. The primary minimum-family count is unchanged.
    if k == 1:
        candidates = [g for g in groups if Fraction(g['lower_bound_loss_exact']) >= lo]
        chosen = [min(candidates, key=lambda g: (g['treated_resolved']+g['control_resolved'], Fraction(g['lower_bound_loss_exact']), g['coding']))]
        loss = Fraction(chosen[0]['lower_bound_loss_exact'])
    before = sum((Fraction(g['lower_bound_loss_exact']) for g in groups[:k-1]), Fraction(0))
    assert before < lo <= loss
    after = []
    for arm, counts in [('treated', t_counts), ('control', c_counts)]:
        out = dict(counts)
        for g in chosen:
            n = g[arm+'_resolved']
            out['E' if g['event'] else 'NE'] -= n
            out['U'] += n
        assert out['E'] + out['NE'] + out['U'] == out['N']
        after.append(out)
    new_lo = Fraction(after[0]['E'], nt) - Fraction(after[1]['E']+after[1]['U'], nc)
    assert new_lo == lo-loss and new_lo <= 0
    return {'event': event, 'original_lower_exact': str(lo), 'minimum_sequence_families': k, 'secondary_record_minimization': 'exact_among_sufficient_singletons' if k == 1 else 'not_optimized', 'available_adverse_families': len(groups), 'max_loss_with_k_minus_one_exact': str(before), 'selected_loss_exact': str(loss), 'post_withdrawal_lower_exact': str(new_lo), 'selected_families': chosen, 'after_counts': after, 'total_resolved_records_withdrawn': sum(g['treated_resolved']+g['control_resolved'] for g in chosen), 'all_group_impacts': groups, 'scope': 'Hypothetical whole-family withdrawal on original observed coding strings; not observed errors or a probability model.'}

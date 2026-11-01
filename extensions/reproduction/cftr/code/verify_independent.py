"""Second grouping/translation implementation over the exported real records.
This is computational verification on the same sources, not an independent
biological cohort, analyst or molecular truth benchmark.
"""
from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
import gzip, hashlib, json

R = Path(__file__).resolve().parents[1]
GROUPS = {'F':'TTT TTC','L':'TTA TTG CTT CTC CTA CTG','I':'ATT ATC ATA','M':'ATG','V':'GTT GTC GTA GTG','S':'TCT TCC TCA TCG AGT AGC','P':'CCT CCC CCA CCG','T':'ACT ACC ACA ACG','A':'GCT GCC GCA GCG','Y':'TAT TAC','*':'TAA TAG TGA','H':'CAT CAC','Q':'CAA CAG','N':'AAT AAC','K':'AAA AAG','D':'GAT GAC','E':'GAA GAG','C':'TGT TGC','W':'TGG','R':'CGT CGC CGA CGG AGA AGG','G':'GGT GGC GGA GGG'}
TABLE = {c:a for a,cs in GROUPS.items() for c in cs.split()}
def val(s, prop, exact):
    p = None if len(s)%3 else ''.join(TABLE[s[i:i+3]] for i in range(0,len(s),3))
    return {'observed_coding':s,'peptide_or_frame_status':('OUT_OF_FRAME:'+str(len(s)%3) if p is None else 'PEPTIDE:'+p),'coding_length':len(s),'native_event':p=='GTIKENIIFGVSYD','exact_programmed_event':s==exact}[prop]

def main():
    data = json.loads((R/'results/CODING_CATALOGS.json').read_text())
    atlas = json.loads((R/'results/FACTORIZATION_ATLAS.json').read_text())
    needed = {(e['operator'],e['run'],e['read_id']) for a in atlas for c in a['counterexamples'] for e in c['actual_records']}
    locators = {}
    counts = defaultdict(Counter)
    seen = set()
    disagreements, equal_ne, shared_peptide = 0,0,0
    with gzip.open(R/'results/OBSERVED_RECORD_CODINGS.jsonl.gz','rt') as f:
        for line in f:
            r = json.loads(line)
            key = r['run'],r['read_id']
            assert key not in seen
            seen.add(key)
            codes = r['codings']
            assert r['common_Q30'] or all(c is None for c in codes.values())
            for op,s in codes.items():
                if s is not None:
                    counts[op,r['run']][s] += 1
                    counts[op,'ALL'][s] += 1
                if (op,*key) in needed: locators[op,*key] = s
            witness = all(s is not None for s in codes.values()) and len(set(codes.values())) == 1
            assert r['witness'] == witness
            if witness:
                counts['witness',r['run']][next(iter(codes.values()))] += 1
                counts['witness','ALL'][next(iter(codes.values()))] += 1
            if all(s is not None for s in codes.values()) and not witness:
                disagreements += 1
                labels = [val(s,'native_event',data['exact_programmed_coding']) for s in codes.values()]
                equal_ne += int(labels == [False,False,False])
                ps = [val(s,'peptide_or_frame_status',data['exact_programmed_coding']) for s in codes.values()]
                shared_peptide += int(all(p.startswith('PEPTIDE:') for p in ps) and len(set(ps)) == 1)
    assert len(seen) == 177006 and disagreements == 83 and equal_ne == 18 and shared_peptide == 0
    for op,byrun in data['catalogs'].items():
        for run,rows in byrun.items():
            assert counts[op,run] == Counter({r['coding']:r['count'] for r in rows})
    for a in atlas:
        rows = data['catalogs'][a['operator']][a['run']]
        partition = defaultdict(set)
        weights = Counter()
        for r in rows:
            x = val(r['coding'],a['representation'],data['exact_programmed_coding'])
            y = val(r['coding'],a['target'],data['exact_programmed_coding'])
            partition[x].add(y)
            weights[x] += r['count']
        bad = {k for k,v in partition.items() if len(v)>1}
        assert a['factors_on_observed_catalog'] == (len(bad)==0)
        assert a['heterogeneous_fibers'] == len(bad)
        assert a['records_in_heterogeneous_fibers'] == sum(weights[k] for k in bad)
        for c in a['counterexamples']:
            es = c['actual_records']
            ss = [locators[e['operator'],e['run'],e['read_id']] for e in es]
            assert all(s == e['coding'] for s,e in zip(ss,es))
            assert val(ss[0],a['representation'],data['exact_programmed_coding']) == val(ss[1],a['representation'],data['exact_programmed_coding'])
            assert val(ss[0],a['target'],data['exact_programmed_coding']) != val(ss[1],a['target'],data['exact_programmed_coding'])
    blocks = json.loads((R/'results/SEQUENCE_FAMILY_WITHDRAWAL.json').read_text())
    for b in blocks:
        assert b['minimum_sequence_families'] == 1 and F(b['original_lower_exact']) > 0
        sufficient = [g for g in b['all_group_impacts'] if F(g['lower_bound_loss_exact']) >= F(b['original_lower_exact'])]
        assert sufficient and b['selected_families'][0] in sufficient
        assert b['total_resolved_records_withdrawn'] == min(g['treated_resolved']+g['control_resolved'] for g in sufficient)
        t,c = b['after_counts']
        assert F(t['E'],t['N'])-F(c['E']+c['U'],c['N']) == F(b['post_withdrawal_lower_exact']) <= 0
    joint=json.loads((R/'results/JOINT_TARGET_PARTITIONS.json').read_text())
    for row in joint['rows']:
        run=row['run'];expected=Counter()
        for coding,n in counts['witness',run].items():
            expected[val(coding,'native_event',data['exact_programmed_coding']),val(coding,'exact_programmed_event',data['exact_programmed_coding'])]+=n
        assert {tuple(g['target_signature']):g['records'] for g in row['partition']['classes']}==dict(expected)
        assert row['partition']['minimum_observed_classes']==len(expected)
        if run!='ALL':
            q=row['raw_reporting_counts'];assert sum(q[k] for k in ['exact_programmed','native_other_observed_coding','non_native','unknown'])==q['N']
            hull=row['native_minus_exact_within_run'];naive=list(map(F,hull['naive_separate_marginal_hull_exact']));new=list(map(F,hull['joint_logically_compatible_hull_exact']))
            assert new[0]-naive[0]==F(q['unknown'],q['N']) and new[1]==naive[1]
            assert new==[F(q['native_other_observed_coding'],q['N']),F(q['native_other_observed_coding']+q['unknown'],q['N'])]
    result = {'status':'PASS','real_record_identities_checked':len(seen),'catalog_partitions_checked':len(atlas),'counterexample_locators_checked':len(locators),'whole_family_minima_checked':len(blocks),'all_minima_one_with_positive_original_lower_bound':True,'least_record_singleton_witnesses_checked':True,'joint_target_partitions_checked':len(joint['rows']),'joint_hulls_checked':len(joint['rows'])-1,'same_sources_not_independent_biological_truth':True}
    (R/'results/INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__ == '__main__': main()

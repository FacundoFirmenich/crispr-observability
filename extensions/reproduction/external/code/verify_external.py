"""Second implementation: relational sorting on independently parsed source rows.

Does not import the production adapter or factorization/partition functions.
Uses real tables only, and adjudicates finite computational claims, not biology.
"""
from pathlib import Path
from collections import Counter
from itertools import groupby
import csv,hashlib,io,json,zipfile
ROOT=Path(__file__).resolve().parents[1]
def enc(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    result=json.loads((ROOT/'results/EXTERNAL_RESULTS.json').read_text());checked=0;witnesses=0;rows_checked=0
    report=[]
    for case in result['cases']:
        name=case['case'];inputs=ROOT/'inputs'
        info=json.loads((inputs/f'{name}_CRISPResso2_info.json').read_text())
        a=info['running_info']['args']['value'];templates=dict(zip(a['amplicon_name'].split(','),a['amplicon_seq'].split(',')))
        p=inputs/f'{name}_Alleles_frequency_table.zip'
        with zipfile.ZipFile(p) as z:
            member=[m for m in z.namelist() if m.endswith('.txt')][0]
            raw=list(csv.DictReader(io.StringIO(z.read(member).decode()),delimiter='\t'))
        data=[]
        for line,r in enumerate(raw,2):
            s=r['Aligned_Sequence'].upper();ref=r['Reference_Sequence'].upper();seq=''.join(x for x in s if x!='-');rn=r['Reference_Name'];n=int(r['#Reads'])
            features={'source_status':r['Read_Status'],'source_reference':rn,'source_status_and_reference':[r['Read_Status'],rn],
              'observed_length':len(seq),'observed_sequence':seq,'full_alignment':[s,ref,rn],
              'exact_declared_template':None if rn not in templates else seq==templates[rn],
              'any_alignment_indel':any(c=='-' for c in s+ref)}
            data.append((line,n,features));rows_checked+=1
        lookup={line:(n,f) for line,n,f in data}
        assert sum(n for _,n,_ in data)==case['assignment_weight']
        supplied=json.loads((ROOT/'results'/f'{name}_CATALOG.json').read_text())['rows']
        assert len(supplied)==len(data)
        for row,(line,n,f) in zip(supplied,data):
            assert row['features']==f and row['weight']==n
            assert row['locator']=={'case':name,'zip_file':p.name,'zip_sha256':digest(p),'member':member,'line_including_header':line}
            assert row['original']==raw[line-2]
        for cert in case['atlas']:
            r,t=cert['representation'],cert['target'];tuples=sorted((enc(f[r]),enc(f[t]),n,line) for line,n,f in data if f[t] is not None)
            bad=0;mass=0;total_fibers=0
            for _,g in groupby(tuples,key=lambda v:v[0]):
                group=list(g);total_fibers+=1
                if len(set(v[1] for v in group))>1:bad+=1;mass+=sum(v[2] for v in group)
            unknown=sum(n for _,n,f in data if f[t] is None)
            support=len({enc(f[t]) for _,_,f in data if f[t] is not None})
            assert cert['distinct_resolved_target_values']==support
            assert cert['target_constant_on_resolved_domain']==(support==1)
            assert (bad,mass,unknown,total_fibers)==(cert['heterogeneous_fibers'],cert['assignment_weight_in_heterogeneous_fibers'],cert['target_unresolved_weight'],cert['observed_resolved_fibers'])
            assert cert['factors_on_resolved_subcatalog']==(bad==0)
            assert cert['complete_catalog_authorized']==(bad==0 and unknown==0)
            for witness in cert['counterexamples']:
                fs=[]
                for locator in witness['actual_rows']:
                    assert locator['zip_sha256']==digest(p) and locator['member']==member and locator['case']==name
                    line=locator['line_including_header'];n,f=lookup[line];assert n>0;fs.append(f)
                assert fs[0][r]==fs[1][r]==witness['representation_value']
                assert fs[0][t]!=fs[1][t]
                assert [fs[0][t],fs[1][t]]==witness['target_values'];witnesses+=1
            checked+=1
        targets=case['joint']['targets'];groups={};joint_unknown=0
        for line,n,f in data:
            sig=[f[t] for t in targets]
            if None in sig:joint_unknown+=n;continue
            groups[enc(sig)]=groups.get(enc(sig),0)+n
        assert len(groups)==case['joint']['minimum_classes_on_joint_resolved_domain']
        assert joint_unknown==case['joint']['joint_unresolved_weight']
        assert {enc(c['signature']):c['weight'] for c in case['joint']['classes']}==groups
        mapping=list(csv.DictReader(io.StringIO(case['source_mapping_text']),delimiter='\t'))[0]
        assert int(mapping['READS ALIGNED'])==case['assignment_weight']
        assert int(mapping['READS IN INPUTS'])==case['raw'][0]['records']
        quant=list(csv.DictReader(io.StringIO(case['source_quantification_text']),delimiter='\t'))
        qweight=0
        for q in quant:
            subset=[(n,f) for _,n,f in data if f['source_reference']==q['Amplicon']]
            assert sum(n for n,_ in subset)==int(q['Reads_aligned'])
            for status,col in [('MODIFIED','Modified'),('UNMODIFIED','Unmodified')]:
                assert sum(n for n,f in subset if f['source_status']==status)==int(q[col])
            qweight+=int(q['Reads_aligned'])
        ambiguous=case['assignment_weight']-qweight
        assert ambiguous==sum(n for _,n,f in data if f['source_reference'] not in templates)
        feature_counts={}
        for t in ['exact_declared_template','any_alignment_indel']:
            c=Counter()
            for _,n,f in data:c[enc(f[t])]+=n
            feature_counts[t]=dict(c)
        report.append({'case':name,'table_rows':len(data),'assignment_weight':case['assignment_weight'],
          'quantified_unambiguous_weight':qweight,'ambiguous_assignment_weight':ambiguous,
          'unrepresented_raw_count':int(mapping['READS IN INPUTS'])-case['assignment_weight'],
          'source_quantification_counts':'PASS','source_mapping_total':'PASS','feature_weights':feature_counts})
    output={'status':'PASS','parsed_actual_rows':rows_checked,'factorization_certificates_checked':checked,
      'actual_counterexample_pairs_checked':witnesses,'joint_partitions_checked':len(report),'cases':report,
      'limits':'Second computational implementation on identical public sources; no independent analyst or molecular truth. Raw sequence linkage in main runner is not an engine replay.'}
    (ROOT/'results/INDEPENDENT_VERIFICATION.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()

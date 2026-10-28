"""Raw-identity, allele-count, baseline, paired-label and interval verification."""
import argparse, csv, gzip, hashlib, io, json, os, zipfile
from collections import Counter
from fractions import Fraction
import run_grid as g

EVENTS=['CTT_plus_full_SE2','native_local_peptide_all_observed_sequences']
CODES={'U':0,'E':1,'N':2}; LETTERS='UEN'
AA={'F':'TTT TTC','L':'TTA TTG CTT CTC CTA CTG','I':'ATT ATC ATA','M':'ATG','V':'GTT GTC GTA GTG','S':'TCT TCC TCA TCG AGT AGC','P':'CCT CCC CCA CCG','T':'ACT ACC ACA ACG','A':'GCT GCC GCA GCG','Y':'TAT TAC','*':'TAA TAG TGA','H':'CAT CAC','Q':'CAA CAG','N':'AAT AAC','K':'AAA AAG','D':'GAT GAC','E':'GAA GAG','C':'TGT TGC','W':'TGG','R':'CGT CGC CGA CGG AGA AGG','G':'GGT GGC GGA GGG'}
CODE={c:a for a,cs in AA.items() for c in cs.split()}
def project_independent(ar,aq,r):
    # Distinct coordinate-slicing implementation from frozen verify_m2.py.
    assert ar.replace('-','')==r['sequence'] and len(ar)==len(aq)
    columns=[i for i,v in enumerate(ar) if v!='-'];lo=r['core_start_0'];hi=r['coding_end_exclusive_0']
    if any(aq[columns[i]] not in 'ACGT' for i in [lo-1,hi]):return ('U','U')
    s=aq[columns[lo-1]+1:columns[hi-1]+1].replace('-','')
    if not s or set(s)-set('ACGT'):return ('U','U')
    peptide=''.join(CODE[s[i:i+3]] for i in range(0,len(s),3)) if len(s)%3==0 else None
    return ('E' if s==r['exact_programmed_core']+'AT' else 'N','E' if peptide==r['native_peptide'] else 'N')
def plain(x):
    if isinstance(x,dict) and '_type' in x and 'value' in x:return plain(x['value'])
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,list):return [plain(v) for v in x]
    return x
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def bounds(c):
    lo=Fraction(c['E'],c['N']);hi=Fraction(c['E']+c['U'],c['N'])
    return {**c,'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi)}
def one_job(job,spec,raw,idx,manifest,quality_pass_ids):
    rid=job['run'];ref=spec['references'][job['reference']]
    f=g.BASE/'runs'/job['reference']/job['setting']/('CRISPResso_on_'+rid)
    c=[Counter(),Counter()]; labels=bytearray(2*manifest['N_raw']); seen=set(); cache={};nc=Counter()
    for name,s,plus,q in g.rows(f/'CRISPResso_output.fastq.gz'):
        assert name not in seen;seen.add(name);assert (s,q)==raw[name]
        assert name in quality_pass_ids
        a=dict(x.split('=',1) for x in plus[1:].split() if '=' in x)
        if a['ALN']=='NA':ll=('U','U');nc['not_aligned']+=1
        else:
            assert a['ALN']=='observed_mutant';nc[a['CLASS']]+=1
            k=a['ALN_REF'],a['ALN_SEQ']
            if k not in cache:
                cache[k]=g.project(*k,ref)
                assert cache[k]==project_independent(*k,ref)
            ll=cache[k]
        for i,l in enumerate(ll):c[i][l]+=1;labels[2*idx[name]+i]=CODES[l]
    assert len(seen)==manifest['Q30_pass'] and seen==quality_pass_ids
    for cc in c:cc['U']+=manifest['N_raw']-len(seen)
    allele=[Counter(),Counter()];aligned=0
    with zipfile.ZipFile(f/'Alleles_frequency_table.zip') as z:
        with z.open(z.namelist()[0]) as b:
            for a in csv.DictReader(io.TextIOWrapper(b),delimiter='\t'):
                n=int(a['#Reads']);aligned+=n;ll=g.project(a['Reference_Sequence'],a['Aligned_Sequence'],ref)
                for i,l in enumerate(ll):allele[i][l]+=n
    assert aligned==len(seen)-nc['not_aligned']
    for i in range(2):
        allele[i]['U']+=manifest['N_raw']-aligned
        assert all(allele[i][k]==c[i][k] for k in 'ENU')
    quant=table(f/'CRISPResso_quantification_of_editing_frequency.txt');assert len(quant)==1;qr=quant[0]
    assert int(qr['Reads_aligned'])==aligned
    assert int(qr['Modified'])==nc['observed_mutant_MODIFIED'] and int(qr['Unmodified'])==nc['observed_mutant_UNMODIFIED']
    info=plain(json.loads((f/'CRISPResso2_info.json').read_text()));args=info['running_info']['args'];rr=info['results']['refs']['observed_mutant']
    effective={'options':{k:args[k] for k in ['default_min_aln_score','needleman_wunsch_gap_open','needleman_wunsch_gap_extend','needleman_wunsch_gap_incentive','quantification_window_coordinates','exclude_bp_from_left','exclude_bp_from_right','trim_sequences']},
      'include_idxs':rr['include_idxs'],'exclude_idxs':rr['exclude_idxs'],'contains_guide':rr['contains_guide'],'effective_gap_incentive_vector':rr['gap_incentive']}
    summary={**job,'state':'VERIFIED','N_raw':manifest['N_raw'],'Q30_pass':len(seen),'aligned':aligned,'projected':{e:bounds(g.summarize(c[i],manifest['N_raw'])) for i,e in enumerate(EVENTS)},
      'native_software_quantification':qr,'native_annotation_classes':dict(nc),'effective_configuration':effective,
      'validation':{'original_read_identity_and_quality':'PASS','all_Q30_records_represented':'PASS','independent_coordinate_projection':'PASS','allele_table_vs_per_read_projection':'PASS','native_table_vs_per_read_classes':'PASS'}}
    return summary,labels
def baseline_check(summaries):
    old0=json.loads((g.OLD/'results/COUNTS.json').read_text());old1=json.loads((g.OLD/'results/REFERENCE_COUNTS.json').read_text());checks=[]
    for s in summaries:
        if s['setting']!='baseline':continue
        old=old0 if s['reference']=='mutant' else old1;op='crispresso_projected' if s['reference']=='mutant' else 'programmed_reference'
        for e in EVENTS:
            want=next(x for x in old if x['run_accession']==s['run'] and x['observer']==op and x['event']==e)
            match=all(s['projected'][e][k]==want[k] for k in ['N','E','NE','U'])
            checks.append({'reference':s['reference'],'run':s['run'],'event':e,'counts_exact_match':match})
        folder='runs' if s['reference']=='mutant' else 'runs_programmed_ref'
        oldtable=table(g.OLD/folder/('CRISPResso_on_'+s['run'])/'CRISPResso_quantification_of_editing_frequency.txt')
        checks.append({'reference':s['reference'],'run':s['run'],'native_quantification_table_exact_match':oldtable==[s['native_software_quantification']]})
    passed=all(x.get('counts_exact_match',x.get('native_quantification_table_exact_match')) for x in checks)
    result={'state':'PASS' if passed else 'FAIL','baseline_jobs':len([s for s in summaries if s['setting']=='baseline']),'checks':checks}
    g.dump(g.BASE/'BASELINE_VALIDATION.json',result);assert passed,'Baseline differs from frozen outputs'
def run(stage):
    spec=json.loads((g.BASE/'GRID.json').read_text());man={r['run']:r for r in json.loads((g.BASE/'INPUT_MANIFEST.json').read_text())}
    jobs=[]
    for p in sorted((g.BASE/'executions').rglob('*.json')):
        j=json.loads(p.read_text())
        if j['state']=='COMPLETED' and (stage=='all' or j['setting']=='baseline'):jobs.append(j)
    if stage=='all':assert len(jobs)==144,len(jobs)
    elif stage=='baseline':assert len(jobs)==18,len(jobs)
    summaries=[];transitions=[];joint=[];ranges=[];examples=[]
    joint_tmp=g.BASE/'JOINT_LABELS.jsonl.gz.tmp'
    joint_out=gzip.open(joint_tmp,'wt') if stage=='all' else None
    for run in spec['runs']:
        rid=run['run_accession'];raw={};ids=[];digest=hashlib.sha256()
        for name,s,plus,q in g.rows(g.RAW/(rid+'.fastq.gz')):
            assert name not in raw;raw[name]=(s,q);ids.append(name);digest.update((s+'\n'+q+'\n').encode())
        assert len(ids)==man[rid]['N_raw'] and digest.hexdigest()==man[rid]['raw_sequence_quality_stream_sha256'];idx={x:i for i,x in enumerate(ids)};per={}
        quality_pass_ids={name for name,(seq,qual) in raw.items() if sum(ord(z)-33 for z in qual)>=30*len(qual)}
        assert len(quality_pass_ids)==man[rid]['Q30_pass']
        local=[j for j in jobs if j['run']==rid]
        for j in local:
            summary,ll=one_job(j,spec,raw,idx,man[rid],quality_pass_ids);summaries.append(summary);per[j['reference']+'/'+j['setting']]=ll
        for ref in spec['references']:
            key=ref+'/baseline'
            if key not in per:continue
            b=per[key]
            for k,ll in per.items():
                if not k.startswith(ref+'/') or k==key:continue
                for ei,e in enumerate(EVENTS):
                    c=Counter(LETTERS[b[2*i+ei]]+'>'+LETTERS[ll[2*i+ei]] for i in range(len(ids)))
                    transitions.append({'run':rid,'reference':ref,'setting':k.split('/')[1],'event':e,'N':len(ids),'transitions':{a+'>'+z:c[a+'>'+z] for a in 'ENU' for z in 'ENU'}})
        if stage=='all':
            scopes={'mutant_OAT':[k for k in per if k.startswith('mutant/') and 'stock' not in k],
                    'programmed_OAT':[k for k in per if k.startswith('programmed/') and 'stock' not in k],
                    'both_references_OAT':[k for k in per if 'stock' not in k],
                    'both_references_all_profiles':list(per)}
            for scope,keys in scopes.items():
                for ei,e in enumerate(EVENTS):
                    cc=Counter()
                    for i,name in enumerate(ids):
                        vals={per[k][2*i+ei] for k in keys}
                        kind='all_'+LETTERS[next(iter(vals))] if len(vals)==1 else ('E_NE_conflict' if 1 in vals and 2 in vals else 'resolution_disagreement')
                        cc[kind]+=1
                        if kind=='E_NE_conflict' and len([x for x in examples if x['run']==rid and x['scope']==scope and x['event']==e])<10:
                            examples.append({'run':rid,'read_id':name,'scope':scope,'event':e,'labels':{k:LETTERS[per[k][2*i+ei]] for k in keys},'raw_sequence_sha256':__import__('hashlib').sha256(raw[name][0].encode()).hexdigest()})
                    joint.append({'run':rid,'scope':scope,'event':e,'N':len(ids),'profiles':keys,'counts':dict(cc)})
            for i,name in enumerate(ids):joint_out.write(json.dumps({'run':rid,'read_id':name,'labels':{k:LETTERS[v[2*i]]+LETTERS[v[2*i+1]] for k,v in per.items()}},separators=(',',':'))+'\n')
        print(json.dumps({'run':rid,'verified_jobs':len(local)}),flush=True)
    if joint_out:
        joint_out.close()
        with joint_tmp.open('rb') as f:os.fsync(f.fileno())
        gzip.open(joint_tmp,'rb').read()
        os.replace(joint_tmp,g.BASE/'JOINT_LABELS.jsonl.gz')
    baseline_check(summaries)
    for r in spec['runs']:
        for ref in spec['references']:
            for e in EVENTS:
                ss=[x for x in summaries if x['run']==r['run_accession'] and x['reference']==ref and x['setting']!='stock_default_no_guide']
                if not ss:continue
                ranges.append({'run':r['run_accession'],'reference':ref,'event':e,'scope':'OAT_only','settings':[x['setting'] for x in ss],
                  'ranges':{k:{'min':min(x['projected'][e][k] for x in ss),'max':max(x['projected'][e][k] for x in ss)} for k in ['E','NE','U','L_pp','H_pp']}})
    counts={(x['reference'],x['setting'],x['run']):x for x in summaries};contrasts=[]
    for x in summaries:
        r=next(z for z in spec['runs'] if z['run_accession']==x['run'])
        if r['condition']!='PE treated':continue
        for ctl in ['mock','untreated']:
            cr=next(z for z in spec['runs'] if z['condition']==ctl and z['source_replica']==r['source_replica']);y=counts.get((x['reference'],x['setting'],cr['run_accession']))
            if not y:continue
            for e in EVENTS:
                a=x['projected'][e];b=y['projected'][e];lo=Fraction(a['L_exact'])-Fraction(b['H_exact']);hi=Fraction(a['H_exact'])-Fraction(b['L_exact'])
                contrasts.append({'reference':x['reference'],'setting':x['setting'],'preparation':r['source_replica'],'PE':x['run'],'control':ctl,'control_run':cr['run_accession'],'event':e,'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi)})
    for name,data in [('SUMMARIES',summaries),('TRANSITIONS',transitions),('JOINT_STABILITY',joint),('SENSITIVITY_RANGES',ranges),('CONTRASTS',contrasts),('CONFLICT_EXAMPLES',examples)]:g.dump(g.BASE/(name+'.json'),data)
    native_conflicts=sum(x['transitions']['E>N']+x['transitions']['N>E'] for x in transitions if x['event']==EVENTS[1] and x['setting']!='stock_default_no_guide')
    g.dump(g.BASE/'VERIFICATION.json',{'state':'PASS','stage':stage,'jobs':len(summaries),'raw_records':177006,'baseline_validation':'PASS','all_jobs_identity_and_two_count_routes':'PASS','native_OAT_resolved_label_conflicts_against_baseline':native_conflicts,'limits':'Finite parameter sensitivity; neither independent biology nor molecular truth calibration. Stock multiple-change profile excluded from OAT ranges.'})
    print(json.dumps({'state':'PASS','stage':stage,'verified_jobs':len(summaries),'native_OAT_conflicts_vs_baseline':native_conflicts}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['pilot','baseline','all']);a=p.parse_args();run(a.stage)

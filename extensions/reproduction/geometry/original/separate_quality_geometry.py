"""Retrospective paired policies on unchanged real v22 FASTQ streams.

No biological simulation, alignment imputation or independently sampled reads.
Structural rules are reimplemented and regression checked against frozen v22.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter, defaultdict
import argparse, gzip, hashlib, itertools, json, time

OBS = {'target': ('TGGATTATGCCT', 'ATGAATATAGATACAGAAGC'),
       'broad_current': ('GATGATTATGGG', 'AAACTATGTGAA'),
       'broad_edge': ('GATGATTATGGG', 'GTGAAAACTTTT')}
EVENTS = ['CTT_plus_full_SE2', 'native_local_peptide_all_observed_sequences']
MODES = ['structural_only','own_meanQ30','common_raw_meanQ30','common_and_own_meanQ30']
PAIRS = [('broad_edge','target'),('broad_current','target'),('broad_edge','broad_current')]
CODON = dict(zip([''.join(c) for c in itertools.product('TCAG',repeat=3)],
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMP = str.maketrans('ACGTN','TGCAN')

def require(ok,message):
    if not ok: raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def dump(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def positions(s,motif):
    result=[];start=0
    while True:
        i=s.find(motif,start)
        if i<0:return result
        result.append(i);start=i+1

def inspect(oriented,name,exact):
    left,right=OBS[name];hits=[]
    for s,q,orientation in oriented:
        for a in positions(s,left):
            for b in positions(s,right):
                if a+len(left)<=b<=a+len(left)+350:
                    hits.append((s,q,a,b,orientation))
    result={'labels':('U','U'),'qpass':False,'reason':'no_outer_pair' if not hits else 'multiple_outer_pairs','pair_count':len(hits)}
    if len(hits)!=1:return result
    s,q,a,b,orientation=hits[0];payload=s[a+len(left):b];span=q[a:b+len(right)]
    require(bool(span),'empty quality span')
    result.update(qpass=sum(ord(x)-33 for x in span)>=30*len(span),orientation=orientation,
                  outer_start=a,outer_end=b+len(right),span_length=len(span))
    if any(x not in 'ACGT' for x in payload):
        result['reason']='non_ACGT_outer_payload';return result
    if name=='target':core=payload
    else:
        tl,tr=OBS['target']
        if payload.count(tl)!=1 or payload.count(tr)!=1:
            result['reason']='nested_target_motif_nonunique_or_missing';return result
        ta=payload.find(tl)+len(tl);tb=payload.find(tr,ta)
        if tb<0:
            result['reason']='nested_target_order';return result
        core=payload[ta:tb]
    coding=core+OBS['target'][1][:2]
    native=len(coding)%3==0 and ''.join(CODON[coding[i:i+3]] for i in range(0,len(coding),3))=='GTIKENIIFGVSYD'
    result.update(labels=('E' if core==exact else 'N','E' if native else 'N'),reason='resolved_structurally')
    return result

def geometry(oriented):
    """Literal expected edge location, only under unique observed target/current motifs."""
    candidates=[];left,right=OBS['target'];current=OBS['broad_current'][1];edge=OBS['broad_edge'][1]
    for s,q,orientation in oriented:
        if s.count(left)!=1 or s.count(right)!=1 or s.count(current)!=1:continue
        a=s.find(left);b=s.find(right);c=s.find(current)
        if not a+len(left)<=b<c:continue
        start=c+7;end=start+len(edge);available=s[start:min(end,len(s))]
        category=('full_edge_at_expected' if end<=len(s) and s[start:end]==edge else
                  'matching_prefix_but_read_ends' if end>len(s) and available and edge.startswith(available) else
                  'other_or_mismatch')
        candidates.append((category,len(s),c,start,end,orientation))
    return candidates[0] if len(candidates)==1 else ('not_uniquely_localizable',len(oriented[0][0]),None,None,None,None)

def mode_labels(info,mode,common):
    gate=(mode=='structural_only' or mode=='own_meanQ30' and info['qpass'] or
          mode=='common_raw_meanQ30' and common or mode=='common_and_own_meanQ30' and common and info['qpass'])
    return info['labels'] if gate else ('U','U')

def summarize(c):
    n=sum(c.values());require(n>0 and set(c)<=set('ENU'),'invalid label partition')
    e,ne,u=c['E'],c['N'],c['U'];require(e+ne+u==n,'partition failure')
    return {'N':n,'E':e,'NE':ne,'U':u,'L_exact':str(F(e,n)),'H_exact':str(F(e+u,n)),
            'L_pp':float(100*F(e,n)),'H_pp':float(100*F(e+u,n))}

def endpoints(r):return F(r['L_exact']),F(r['H_exact'])

def run(root,out,design_path):
    root=Path(root).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    require(not out.is_relative_to(root),'output must not modify frozen source')
    design=json.loads(Path(design_path).read_text(encoding='utf-8'))
    require(design['trims']==[0,8] and list(design['policies'])==MODES,'design mismatch')
    start=time.monotonic();prior=json.loads((root/'v21/results/READ_PHASE_RESULTS.json').read_text())
    expected=json.loads((root/'EXPECTED_PRODUCTS.json').read_text())
    refs=json.loads((root/'v21/results/REFERENCE_AND_TEMPLATES.json').read_text())
    exact=next(t['target_sequence'] for t in refs['templates'] if t['CTT_restored'] and t['full_SE2'])
    frozen=[]
    for rel in ['observer/results/OBSERVATION_COUNTS.json','observer/results/TRUNCATION_COUNTS.json']:
        require(sha(root/rel)==expected[rel],'frozen source identity mismatch')
        frozen+=json.loads((root/rel).read_text())
    lookup={(r['run_accession'],r['trim_last_raw_bases'],r['observer'],r['filter'],r['event']):r for r in frozen}
    runs=[];seen=set()
    for r in prior:
        if r['sequence_quality_stream_sha256'] not in seen:
            runs.append(r);seen.add(r['sequence_quality_stream_sha256'])
    require(len(runs)==15,'unique stream set changed')
    counts_rows=[];pair_rows=[];loss_rows=[];length_rows=[];checks=[];input_rows=[]
    for r in runs:
        rid=r['run_accession'];path=root/'v21/sources/fastq'/(rid+'.fastq.gz')
        n=0;digest=hashlib.sha256();common_pass=0;counters=defaultdict(Counter);paired=defaultdict(Counter)
        loss=Counter();lengths=defaultdict(Counter)
        with gzip.open(path,'rt',encoding='ascii',newline='') as f:
            while h:=f.readline():
                s=f.readline().rstrip('\n');plus=f.readline();q=f.readline().rstrip('\n')
                require(h.startswith('@') and plus.startswith('+') and len(s)==len(q) and bool(s),'invalid FASTQ')
                require(all(33<=ord(z)<=126 for z in q),'invalid Phred+33 encoding')
                n+=1;digest.update((s+'\n'+q+'\n').encode('ascii'))
                common=sum(ord(z)-33 for z in q)>=30*len(q);common_pass+=common
                rev=s.translate(COMP)[::-1];qr=q[::-1]
                for trim in design['trims']:
                    oriented=[(s,q,'forward'),(rev,qr,'reverse')] if not trim else [(s[:-trim],q[:-trim],'forward'),(rev[trim:],qr[trim:],'reverse')]
                    infos={name:inspect(oriented,name,exact) for name in OBS};geom=None
                    for mode in MODES:
                        labels={name:mode_labels(info,mode,common) for name,info in infos.items()}
                        union=[]
                        for ei,event in enumerate(EVENTS):
                            vals={v[ei] for v in labels.values()}-{'U'}
                            require(len(vals)<=1,f'CONFLICT {rid}:{n}:{trim}:{mode}:{event}')
                            union.append(next(iter(vals)) if vals else 'U')
                        labels['compatible_union']=tuple(union)
                        for name,v in labels.items():
                            for ei,event in enumerate(EVENTS):
                                counters[trim,mode,name,event][v[ei]]+=1
                                if name!='compatible_union':lengths[len(s),trim,mode,name,event][v[ei]]+=1
                        for a,b in PAIRS:
                            for ei,event in enumerate(EVENTS):
                                av,bv=labels[a][ei],labels[b][ei]
                                paired[trim,mode,a,b,event][av+'>'+bv]+=1
                                if bv!='U' and av=='U':
                                    if geom is None:geom=geometry(oriented)
                                    reason=infos[a]['reason']
                                    if reason=='resolved_structurally':reason='numerical_quality_gate'
                                    loss[trim,mode,a,b,event,bv,reason,common,*geom]+=1
        require(n==r['raw_reads'] and digest.hexdigest()==r['sequence_quality_stream_sha256'],'stream changed')
        fhash=sha(path);require(fhash==lookup[rid,0,'target','ACGT_only',EVENTS[0]]['fastq_sha256'],'compressed source changed')
        meta={k:r[k] for k in ['run_accession','figure','condition','source_replica','source_filename']}
        for (trim,mode,name,event),c in sorted(counters.items()):
            row={**meta,'trim':trim,'policy':mode,'observer':name,'event':event,**summarize(c)}
            require(row['N']==n,'denominator changed')
            if mode in ['structural_only','own_meanQ30']:
                flt='ACGT_only' if mode=='structural_only' else 'meanQ30';old=lookup[rid,trim,name,flt,event]
                require([row[k] for k in ['N','E','NE','U']]==[old[k] for k in ['N','E','non_E','U_E']],f'frozen mismatch {rid,trim,name,flt,event}')
                checks.append({'run':rid,'trim':trim,'observer':name,'filter':flt,'event':event,'status':'PASS'})
            counts_rows.append(row)
        for key,c in sorted(paired.items()):
            trim,mode,a,b,event=key;require(sum(c.values())==n and not c['E>N'] and not c['N>E'],'pair partition or conflict')
            pair_rows.append({**meta,'trim':trim,'policy':mode,'observer_A':a,'observer_B':b,'event':event,'N':n,
                              'transitions':{x+'>'+y:c[x+'>'+y] for x in 'ENU' for y in 'ENU'}})
        for key,count in sorted(loss.items(),key=lambda x:str(x[0])):
            trim,mode,a,b,event,label,reason,common,geo,length,c,gs,ge,orientation=key
            loss_rows.append({**meta,'trim':trim,'policy':mode,'observer_A':a,'observer_B':b,'event':event,'B_label':label,
                              'A_unresolved_reason':reason,'common_raw_Q30':common,'geometry':geo,'oriented_length':length,
                              'current_right_start':c,'expected_edge_start':gs,'expected_edge_end':ge,'orientation':orientation,'count':count})
        for (length,trim,mode,name,event),c in sorted(lengths.items()):
            length_rows.append({'run':rid,'raw_length':length,'trim':trim,'policy':mode,'observer':name,'event':event,**summarize(c)})
        input_rows.append({**meta,'N':n,'common_raw_Q30_pass':common_pass,'fastq_sha256':fhash,'sequence_quality_sha256':digest.hexdigest()})
        dump(out/'PROGRESS.json',{'completed_runs':len(input_rows),'runs':input_rows,'status':'RUNNING'})
        dump(out/'COUNTS.json',counts_rows)
        print(json.dumps({'run':rid,'N':n,'common_raw_Q30_pass':common_pass,'completed':len(input_rows)}),flush=True)
    index={(r['run_accession'],r['trim'],r['policy'],r['observer'],r['event']):r for r in counts_rows}
    contrasts=[]
    primary=[r for r in input_rows if r['figure']=='F4a' and r['condition']=='PE treated']
    for tr in primary:
        for control in ['mock','untreated']:
            cr=next(r for r in input_rows if r['figure']=='F4a' and r['source_replica']==tr['source_replica'] and r['condition']==control)
            for trim,mode,name,event in itertools.product([0,8],MODES,list(OBS)+['compatible_union'],EVENTS):
                t=index[tr['run_accession'],trim,mode,name,event];c=index[cr['run_accession'],trim,mode,name,event]
                tl,th=endpoints(t);cl,ch=endpoints(c);lo,hi=tl-ch,th-cl
                contrasts.append({'prep':{1:'CFB-26-16',2:'CFB-26-02',3:'CFB-06-04'}[tr['source_replica']],
                                  'treated':tr['run_accession'],'control':cr['run_accession'],'control_condition':control,'trim':trim,'policy':mode,
                                  'observer':name,'event':event,'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi),
                                  'sign':'positive' if lo>0 else 'negative' if hi<0 else 'includes_zero'})
    ci={(r['prep'],r['control_condition'],r['trim'],r['policy'],r['observer'],r['event']):r for r in contrasts};decompositions=[]
    for prep,control,trim,event,(a,b) in itertools.product(['CFB-26-16','CFB-26-02','CFB-06-04'],['mock','untreated'],[0,8],EVENTS,PAIRS):
        gaps={};arm={}
        for mode in MODES:
            ar=ci[prep,control,trim,mode,a,event];br=ci[prep,control,trim,mode,b,event];gap=F(br['L_exact'])-F(ar['L_exact'])
            ta=index[ar['treated'],trim,mode,a,event];tb=index[br['treated'],trim,mode,b,event]
            ca=index[ar['control'],trim,mode,a,event];cb=index[br['control'],trim,mode,b,event]
            gain_t=F(tb['E']-ta['E'],ta['N']);gain_c=F(cb['NE']-ca['NE'],ca['N'])
            require(gap==gain_t+gain_c,'lower bound decomposition failure')
            gaps[mode]=str(gap);arm[mode]={'treated_event_gain_pp':float(100*gain_t),'control_non_event_gain_pp':float(100*gain_c)}
        qnet=F(gaps['own_meanQ30'])-F(gaps['structural_only'])
        decompositions.append({'prep':prep,'control_condition':control,'trim':trim,'event':event,'A':a,'B':b,'gap_exact':gaps,
                               'gap_pp':{k:float(100*F(v)) for k,v in gaps.items()},'own_quality_net_component_pp':float(100*qnet),
                               'own_quality_net_exact':str(qnet),'arm_components':arm,'identity_residual_exact':'0'})
    for name,rows in [('CONTRASTS',contrasts),('PAIRS',pair_rows),('LOSS_INTERSECTIONS',loss_rows),('DECOMPOSITIONS',decompositions)]:dump(out/(name+'.json'),rows)
    with gzip.GzipFile(filename=str(out/'LENGTH_STRATA.json.gz'),mode='wb',mtime=0) as f:f.write((json.dumps(length_rows,sort_keys=True,separators=(',',':'))+'\n').encode())
    checks_doc={'status':'PASS','frozen_regression_checks':len(checks),'checks':checks,'raw_reads':sum(r['N'] for r in input_rows),
                'unique_streams':len(input_rows),'source_design_sha256':sha(design_path),'code_sha256':sha(__file__),'inputs':input_rows,
                'runtime_seconds':time.monotonic()-start,'counts_rows':len(counts_rows),'contrasts':len(contrasts),'decompositions':len(decompositions),
                'scope':'same historical streams; no label-error calibration or independent biological replication'}
    dump(out/'VALIDATION.json',checks_doc);dump(out/'PROGRESS.json',{'status':'COMPLETED','completed_runs':len(input_rows),'runs':input_rows})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);p.add_argument('--design',required=True)
    a=p.parse_args();run(a.root,a.out,a.design)

"""Post-hoc extension: retain original target ambiguity under all four M1 gates."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import argparse,gzip,hashlib,json
import separate_quality_geometry as m

def run(root,out):
    root=Path(root);out=Path(out);validation=json.loads((out/'VALIDATION.json').read_text())
    m.require(m.sha(Path(m.__file__))==validation['code_sha256'],'classification code changed')
    nominal=json.loads((out/'COUNTS.json').read_text());nominal={(r['run_accession'],r['policy'],r['event']):r for r in nominal if r['observer']=='target' and r['trim']==8}
    ref=json.loads((root/'v21/results/REFERENCE_AND_TEMPLATES.json').read_text());exact=next(t['target_sequence'] for t in ref['templates'] if t['CTT_restored'] and t['full_SE2'])
    old=json.loads((root/'observer/results/ORIGINAL_AMBIGUITY_VETO.json').read_text());old={r['run_accession']:r for r in old['samples'] if r['trim_last_raw_bases']==8}
    rows=[];checks=[]
    for r in validation['inputs']:
        path=root/'v21/sources/fastq'/(r['run_accession']+'.fastq.gz');m.require(m.sha(path)==r['fastq_sha256'],'input changed')
        c={p:[Counter(),Counter()] for p in m.MODES};nv={p:[Counter(),Counter()] for p in m.MODES};veto={p:Counter() for p in m.MODES};n=0
        with gzip.open(path,'rt',encoding='ascii') as f:
            while h:=f.readline():
                s=f.readline().rstrip('\n');plus=f.readline();q=f.readline().rstrip('\n');n+=1
                m.require(h.startswith('@') and plus.startswith('+') and len(s)==len(q),'invalid FASTQ')
                rev=s.translate(m.COMP)[::-1];qr=q[::-1];original=[(s,q,'forward'),(rev,qr,'reverse')]
                original_info=m.inspect(original,'target',exact)
                trimmed=m.inspect([(s[:-8],q[:-8],'forward'),(rev[8:],qr[8:],'reverse')],'target',exact)
                common=sum(ord(x)-33 for x in q)>=30*len(q)
                for policy in m.MODES:
                    labels=m.mode_labels(trimmed,policy,common)
                    for ei,label in enumerate(labels):
                        nv[policy][ei][label]+=1
                        if original_info['pair_count']>1 and label!='U':
                            veto[policy][ei,label]+=1;label='U'
                        c[policy][ei][label]+=1
        m.require(n==r['N'],'N mismatch')
        for policy in m.MODES:
            for ei,event in enumerate(m.EVENTS):
                x=m.summarize(nv[policy][ei]);base=nominal[r['run_accession'],policy,event]
                m.require(all(x[k]==base[k] for k in ['N','E','NE','U']),'nominal target changed')
                row={**{k:r[k] for k in ['run_accession','figure','condition','source_replica','source_filename']},'trim':8,'observer':'target',
                     'policy':policy,'event':event,**m.summarize(c[policy][ei]),'veto_E':veto[policy][ei,'E'],'veto_NE':veto[policy][ei,'N']}
                if policy=='own_meanQ30' and ei==1:
                    orig=old[r['run_accession']]
                    m.require([row[k] for k in ['N','E','NE','U']]==[orig[k] for k in ['N','E','non_E','U_E']],'frozen veto mismatch')
                    checks.append({'run':r['run_accession'],'status':'PASS'})
                rows.append(row)
    contrasts=[]
    for t in rows:
        if t['figure']!='F4a' or t['condition']!='PE treated':continue
        for baseline in ['mock','untreated']:
            c=next(r for r in rows if r['figure']=='F4a' and r['source_replica']==t['source_replica'] and r['condition']==baseline and r['policy']==t['policy'] and r['event']==t['event'])
            lo=F(t['E'],t['N'])-F(c['E']+c['U'],c['N']);hi=F(t['E']+t['U'],t['N'])-F(c['E'],c['N'])
            contrasts.append({'prep':{1:'CFB-26-16',2:'CFB-26-02',3:'CFB-06-04'}[t['source_replica']],'control_condition':baseline,
                              'policy':t['policy'],'event':t['event'],'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi)})
    m.dump(out/'COMMON_QUALITY_VETO.json',{'scope':'Post-hoc extension of original target ambiguity veto; original multiple candidate pairs cause trimmed resolved labels to return to U. Same original raw universes.',
             'extension_design':'Declared after main M1 results; not a prospective validation or error calibration','rows':rows,'contrasts':contrasts,
             'frozen_veto_checks':checks,'nominal_target_checks':len(rows),'code_sha256':m.sha(__file__)})
    print(json.dumps({'veto_rows':len(rows),'contrasts':len(contrasts),'original_veto_matches':len(checks),'status':'PASS'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.root,a.out)

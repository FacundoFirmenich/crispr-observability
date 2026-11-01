"""Execute the frozen S1-02 finite-row evaluation on public, real examples."""
from collections import Counter
from pathlib import Path
import csv,gzip,hashlib,io,json,sys,zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'software/src'))
from observed_contracts.tabular import audit,compile_joint,require_claim

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(name,value):
    (ROOT/'results'/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def canonical(s):return min(s,rc(s))

def raw_inventory(path):
    seqs=Counter();headers=[];examples={};logical=hashlib.sha256()
    with gzip.open(path,'rt',newline='') as f:
        while True:
            h=f.readline()
            if not h:break
            s,p,q=f.readline(),f.readline(),f.readline()
            if not h.startswith('@') or not p.startswith('+') or not q:
                raise ValueError('Malformed FASTQ')
            sequence=s.strip().upper();qual=q.rstrip('\r\n')
            if len(sequence)!=len(qual) or set(sequence)-set('ACGTN'):
                raise ValueError('Invalid FASTQ sequence/quality')
            if any(not 33<=ord(c)<=126 for c in qual):raise ValueError('Invalid quality ASCII')
            for v in [h,s,p,q]:logical.update(v.encode())
            ident=h.strip().split()[0][1:];headers.append(ident)
            key=canonical(sequence);seqs[key]+=1
            examples.setdefault(key,{'header':ident,'sequence':sequence})
    if len(set(headers))!=len(headers):raise ValueError('Duplicate FASTQ headers')
    return seqs,headers,examples,{'records':len(headers),'distinct_sequences_up_to_RC':len(seqs),'logical_fastq_sha256':logical.hexdigest(),'compressed_sha256':sha(path),'header_uniqueness':'PASS','structure_and_quality_lengths':'PASS'}

def main():
    protocol=json.loads((ROOT/'PROTOCOL_S1_02.json').read_text())
    freeze=json.loads((ROOT/'PROTOCOL_FREEZE.json').read_text())
    assert sha(ROOT/'PROTOCOL_S1_02.json')==freeze['sha256']
    results=[];catalogs={};inputs=ROOT/'inputs'
    for case in protocol['cases_fixed_order']:
        info=json.loads((inputs/f'{case}_CRISPResso2_info.json').read_text())
        args=info['running_info']['args']['value']
        refs=dict(zip(args['amplicon_name'].split(','),args['amplicon_seq'].split(',')))
        archive=inputs/f'{case}_Alleles_frequency_table.zip'
        with zipfile.ZipFile(archive) as z:
            member=[n for n in z.namelist() if n.endswith('.txt')][0]
            table=list(csv.DictReader(io.StringIO(z.read(member).decode()),delimiter='\t'))
        rows=[];assigned_sequences=Counter()
        for lineno,r in enumerate(table,2):
            s=r['Aligned_Sequence'].upper();ref=r['Reference_Sequence'].upper()
            if len(s)!=len(ref):raise ValueError('Unequal aligned lengths')
            if set(s+ref)-set('ACGTN-'):raise ValueError('Unsupported alignment symbol')
            n=int(r['#Reads']);observed=s.replace('-','');rn=r['Reference_Name']
            features={
              'source_status':r['Read_Status'],'source_reference':rn,
              'source_status_and_reference':[r['Read_Status'],rn],
              'observed_length':len(observed),'observed_sequence':observed,
              'full_alignment':[s,ref,rn],
              'exact_declared_template':observed==refs[rn] if rn in refs else None,
              'any_alignment_indel':'-' in s or '-' in ref}
            locator={'case':case,'zip_file':archive.name,'zip_sha256':sha(archive),'member':member,'line_including_header':lineno}
            rows.append({'locator':locator,'weight':n,'features':features,'original':r})
            assigned_sequences[canonical(observed)]+=n
        catalog={'case':case,'scope':'finite_external_assignment_catalog','rows':rows}
        catalogs[case]=catalog;dump(f'{case}_CATALOG.json',catalog)
        atlas=[audit(rows,a,b) for a in protocol['representations'] for b in protocol['targets']]
        joint=compile_joint(rows,['source_status','exact_declared_template','any_alignment_indel'])
        status=Counter();reference=Counter()
        for row in rows:
            status[row['features']['source_status']]+=row['weight']
            reference[row['features']['source_reference']]+=row['weight']
        rawfiles=['nhej_r1.fastq.gz','nhej_r2.fastq.gz'] if case=='nhej' else [case+'.fastq.gz']
        raw=[raw_inventory(inputs/p) for p in rawfiles]
        inventories=[{'file':p,**v[3]} for p,v in zip(rawfiles,raw)]
        if case=='nhej':
            h1=[h.removesuffix('/1') for h in raw[0][1]];h2=[h.removesuffix('/2') for h in raw[1][1]]
            paired=h1==h2
            linkage={'status':'NOT_AUTHENTICATED_PAIRED_TO_MERGED','paired_header_order_match':paired,
                     'raw_denominator_lift':'REJECTED','reason':'Published processed allele rows are not raw R1/R2 records. The merging, filtering and mapping engine was not rerun.'}
        else:
            inv=raw[0][0];overlap=sum(min(n,inv.get(s,0)) for s,n in assigned_sequences.items())
            deficits=sum(max(0,n-inv.get(s,0)) for s,n in assigned_sequences.items())
            leftovers=sum(max(0,n-assigned_sequences.get(s,0)) for s,n in inv.items())
            linkage={'status':'COMPLETE_SEQUENCE_MULTISET_MATCH' if assigned_sequences==inv else 'INCOMPLETE_SEQUENCE_MULTISET_MATCH',
                     'matching_assignment_weight':overlap,'assignment_weight_without_strict_raw_inventory':deficits,
                     'raw_records_without_table_inventory':leftovers,
                     'comparison':'Exact full observed sequence or its reverse complement; no trimming or alignment guessed.',
                     'raw_denominator_lift':'REJECTED',
                     'reason':'A sequence multiset check is not a processing replay or authenticated per-read assignment; unmatched mass cannot be silently discarded.'}
        claims=[]
        good=next(x for x in atlas if x['representation']=='full_alignment' and x['target']=='observed_sequence')
        for scope in ['finite_external_assignment_catalog','raw_record_catalog','biological_accuracy','future_sufficiency']:
            try:require_claim(good,scope);claims.append({'scope':scope,'decision':'ACCEPTED'})
            except ValueError as e:claims.append({'scope':scope,'decision':'REJECTED','reason':str(e)})
        result={'case':case,'table_rows':len(rows),'assignment_weight':sum(r['weight'] for r in rows),
          'source_status_counts':dict(status),'source_reference_counts':dict(reference),
          'source_version':info['running_info']['version'],'source_command':info['running_info']['command_used'],
          'source_arguments':{k:args[k] for k in ['amplicon_name','amplicon_seq','guide_seq','quantification_window_size','quantification_window_center','exclude_bp_from_left','exclude_bp_from_right','expand_ambiguous_alignments','assign_ambiguous_alignments_to_first_reference']},
          'atlas':atlas,'joint':joint,'raw':inventories,'raw_linkage':linkage,'scope_claims':claims,
          'source_mapping_text':(inputs/f'{case}_CRISPResso_mapping_statistics.txt').read_text(),
          'source_quantification_text':(inputs/f'{case}_CRISPResso_quantification_of_editing_frequency.txt').read_text(),
          'biological_genealogy_independence':'NOT_AUTHENTICATED_BY_TUTORIAL_DISTRIBUTION',
          'engine_rerun':False}
        results.append(result)
    dump('EXTERNAL_RESULTS.json',{'protocol_sha256':freeze['sha256'],'cases':results})
    for r in results:
        c=next(c for c in r['atlas'] if c['representation']=='source_status_and_reference' and c['target']=='observed_sequence')
        print(r['case'],'rows',r['table_rows'],'weight',r['assignment_weight'],'heterogeneous_status_ref_fibers',c['heterogeneous_fibers'],'raw',r['raw_linkage']['status'],flush=True)

if __name__=='__main__':main()

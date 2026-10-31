"""Independent raw FASTQ structure/quality summaries, without read filtering.

This is an original standard-library audit, not FastQC or MultiQC execution.
It retains original records and checks the historical whole-record Phred mask.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse,csv,gzip,hashlib,json
ROOT=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--raw-root',type=Path,default=ROOT/'frozen_sources/raw_replay')
 ap.add_argument('--output',type=Path,default=ROOT/'results/raw_qc')
 a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 manifest=json.loads((a.raw_root/'INPUT_MANIFEST.json').read_text())
 expected=json.loads((a.raw_root/'expected/WITNESS_BOUNDS.json').read_text())
 domains={x['run_accession']:x for x in expected['domains']}
 samples=[];positions=[]
 for item in sorted(manifest['files'],key=lambda x:x['path']):
  if not item['path'].startswith('inputs/fastq/'):continue
  p=a.raw_root/item['path'];assert p.stat().st_size==item['size']and sha(p)==item['sha256'],p
  rid=p.name.split('.')[0];lengths=Counter();sequences=Counter();headers=set();qpass=n=bases=qsum=ambiguous=0;sums=[];counts=[]
  with gzip.open(p,'rt',encoding='ascii',newline='')as f:
   while True:
    head=f.readline()
    if not head:break
    sequence=f.readline().rstrip('\r\n');plus=f.readline().rstrip('\r\n');quality=f.readline().rstrip('\r\n');head=head.rstrip('\r\n')
    assert head.startswith('@')and plus.startswith('+')and len(sequence)==len(quality)>0,(rid,n)
    assert head not in headers,(rid,'duplicate_header',head)
    headers.add(head);n+=1;L=len(sequence);lengths[L]+=1;sequences[sequence]+=1;bases+=L;ambiguous+=sum(c not in 'ACGT'for c in sequence)
    values=[ord(c)-33 for c in quality];assert min(values)>=0 and max(values)<=93
    total=sum(values);qsum+=total;qpass+=total>=30*L
    while len(sums)<L:sums.append(0);counts.append(0)
    for i,v in enumerate(values):sums[i]+=v;counts[i]+=1
  domain=domains[rid];assert n==domain['N_raw'];assert n-qpass==domain['counts']['common_quality_fail']
  samples.append({'run':rid,'records':n,'unique_record_headers':len(headers),'unique_recorded_sequences':len(sequences),'repeated_sequence_records':n-len(sequences),'shortest_read':min(lengths),'longest_read':max(lengths),'mean_read_length_exact':str(Fraction(bases,n)),'original_whole_record_meanQ30_pass':qpass,'original_whole_record_meanQ30_fail':n-qpass,'ambiguous_base_symbols':ambiguous,'mean_recorded_base_Phred_exact':str(Fraction(qsum,bases)),'source_sha256':item['sha256'],'read_length_histogram':dict(sorted(lengths.items()))})
  positions.extend({'run':rid,'position_1':i+1,'reads_covering_position':counts[i],'mean_recorded_Phred_exact':str(Fraction(sums[i],counts[i])),'mean_recorded_Phred':round(sums[i]/counts[i],8)}for i in range(len(sums)))
 assert len(samples)==9
 assert sum(x['records']for x in samples)==177006
 assert sum(x['original_whole_record_meanQ30_fail']for x in samples)==1996
 summary={'status':'PASS','analysis':'Independent standard-library FASTQ structure and recorded-quality audit','raw_records':sum(x['records']for x in samples),'quality_pass_records':sum(x['original_whole_record_meanQ30_pass']for x in samples),'quality_fail_records':sum(x['original_whole_record_meanQ30_fail']for x in samples),'sample_rows':samples,'input_manifest_sha256':sha(a.raw_root/'INPUT_MANIFEST.json'),'historical_witness_result_sha256':sha(a.raw_root/'expected/WITNESS_BOUNDS.json'),'unchanged_raw_source_files_verified':9,'reads_trimmed_or_removed':0,'new_alignment_jobs':0,'FastQC_or_MultiQC_executed':False,'limits':['Recorded Phred values are scores, not experimentally recalibrated error probabilities.','Sequence repeats in an amplicon are not automatically technical artifacts or independent molecules.','Q30 pass does not certify a correct editing label.','No altered observation rule, new sample, donor replication or molecular accuracy estimate is introduced.']}
 (a.output/'RAW_QC_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 for name,rows in [('RAW_QC_BY_STREAM.csv',samples),('RAW_QC_BY_POSITION.csv',positions)]:
  with (a.output/name).open('w',newline='')as f:
   fields=[x for x in rows[0]if x!='read_length_histogram'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:v for k,v in r.items()if k in fields}for r in rows)
 print(json.dumps({k:summary[k]for k in ['status','raw_records','quality_pass_records','quality_fail_records','reads_trimmed_or_removed','new_alignment_jobs','FastQC_or_MultiQC_executed']}))

if __name__=='__main__':main()

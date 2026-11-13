from pathlib import Path
import csv,gzip,hashlib,json,sys,datetime
ROOT=Path(__file__).resolve().parents[2]
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def prepare(out):
 out=Path(out).resolve();assert not out.exists();out.mkdir(parents=True)
 for f in ['inputs','refs','templates','evidence']:(out/f).mkdir()
 ref=json.loads((ROOT/'raw_replay/metadata/REFERENCE.json').read_bytes());edited=json.loads((ROOT/'raw_replay/metadata/REFERENCE_PROGRAMMED.json').read_bytes());construct=json.loads((ROOT/'extensions/reproduction/pe_mode/CONSTRUCT.json').read_bytes())
 runs=json.loads((ROOT/'engine_replay/dependencies/m2/DESIGN.json').read_bytes())['runs'];manifest=[];config=[];guide=construct['spacer_RNA_5to3'].replace('U','T');s=ref['sequence'];assert guide in rc(s)
 for row in runs:
  rid=row['run_accession'];src=ROOT/'raw_replay/inputs/fastq'/(rid+'.fastq.gz');dst=out/'inputs'/(rid+'_R1.fastq.gz');n=k=0;digest=hashlib.sha256()
  with gzip.open(src,'rt') as f,dst.open('wb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as g:
   while h:=f.readline():
    seq,plus,q=f.readline().strip(),f.readline(),f.readline().strip();assert len(seq)==len(q) and plus.startswith('+');n+=1;digest.update((seq+'\n'+q+'\n').encode())
    if sum(ord(c)-33 for c in q)>=30*len(q):k+=1;g.write((h.rstrip()+'\n'+seq+'\n+\n'+q+'\n').encode())
  assert n==row['raw_reads'] and digest.hexdigest()==row['sequence_quality_stream_sha256']
  manifest.append({**row,'Q30':k,'filtered_sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
  (out/'refs'/(rid+'.fasta')).write_text('>'+rid+'\n'+rc(s)+'\n');(out/'templates'/(rid+'.fasta')).write_text('>'+rid+'\n'+rc(edited['sequence'])+'\n')
  config.append([rid,rid,dst.name,'','preparation_'+str(row['source_replica']),0,guide,s[:20],'',1,s,edited['sequence']])
 with (out/'amplican.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['ID','Barcode','Forward_Reads','Reverse_Reads','Group','Control','guideRNA','Forward_Primer','Reverse_Primer','Direction','Amplicon','Donor']);w.writerows(config)
 (out/'evidence/INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2));(out/'evidence/CONSTRUCT.json').write_text(json.dumps(construct,indent=2));(out/'evidence/REFERENCE.json').write_text(json.dumps(ref,indent=2));(out/'evidence/REFERENCE_PROGRAMMED.json').write_text(json.dumps(edited,indent=2))
 (out/'evidence/FREEZE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_sha256':hashlib.sha256((Path(__file__).parent/'ANALYSIS_PLAN.md').read_bytes()).hexdigest(),'N':sum(x['raw_reads'] for x in manifest),'Q30':sum(x['Q30'] for x in manifest)},indent=2))
 print('Prepared',len(runs),'libraries; raw',sum(x['raw_reads'] for x in manifest),'Q30',sum(x['Q30'] for x in manifest),flush=True)
if __name__=='__main__':prepare(sys.argv[1])

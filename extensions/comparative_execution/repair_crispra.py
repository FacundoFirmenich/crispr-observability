"""Reexecute archived variant stages with a typed multi-row return contract."""
from pathlib import Path
import concurrent.futures, hashlib, json, subprocess, sys, tarfile, time, urllib.request

base='https://github.com/FacundoFirmenich/crispr-observability/releases/download/comparison-evidence-43493655364/'
upstream='https://api.bitbucket.org/2.0/repositories/synbiolab/crispr-a_nextflow/src/e5c6c8df1711df6de1f9c3c832ee297d5c48d18d/bin/VC_parser-cigar.R'
image='msanvicente/crisprgareq@sha256:1ff4d7fdcff07a4e4da6eb88ed1676284abe96d4444ee292822e08e33f6acf0d'
w=Path(sys.argv[1]).resolve();w.mkdir();archive=w/'source.tar.gz'
urllib.request.urlretrieve(base+'crispra-export-repair-evidence.tar.gz',archive)
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='34ea9b0669526602ed8dbb5c428f8cbf00a0b02e881092885c3b3ef64bf3cf13'
inputs=w/'inputs';inputs.mkdir();out=w/'result';out.mkdir()
with tarfile.open(archive) as tf:
    for member in tf:
        name=Path(member.name).name
        if member.isfile() and (name.endswith(('.sorted.bam','.sorted.bam.bai','_summary.csv','_reference-correctOrient.fasta','_template.fasta','_edits.csv')) or name=='CONSTRUCT.json'):
            dest=inputs/name
            b=tf.extractfile(member).read()
            if dest.exists(): assert dest.read_bytes()==b
            else: dest.write_bytes(b)
original=urllib.request.urlopen(upstream).read();text=original.decode()
for kind in ['delin','ins-out','ins-in','dels-out','dels-in']:
    old='c(rowReads$count, "'+kind+'", strsplit(rowReads$ids, ","))'
    new='list(sum(rowReads$count), "'+kind+'", unlist(strsplit(rowReads$ids, ","), use.names=FALSE))'
    assert text.count(old)==1; text=text.replace(old,new)
line='write.csv(t_info[3:length(t_info)],file=paste0(results_path, "_template-reads.csv"))'
assert text.count(line)==1
text=text.replace(line,'saveRDS(t_info, file=paste0(results_path, "_template-object.rds"))\n    '+line)
script=out/'VC_parser-cigar_typed.R';script.write_text(text)
(out/'REPAIR.json').write_text(json.dumps({'upstream':upstream,'original_sha256':hashlib.sha256(original).hexdigest(),'patched_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'image':image,'scope':'Typed aggregation of multiple matching alignment rows; no alignment or template-matching criterion changed'},indent=2))
subprocess.run(['docker','pull',image],check=True)
guide=json.loads((inputs/'CONSTRUCT.json').read_bytes())['spacer_RNA_5to3'].replace('U','T')
ids=sorted(p.name.removesuffix('_summary.csv') for p in inputs.glob('*_summary.csv'));assert len(ids)==9
def execute(rid):
    d=out/rid;d.mkdir()
    cmd=['docker','run','--rm','-v',str(w)+':'+str(w),'-w',str(d),image,'Rscript',str(script),str(inputs/(rid+'.sorted.bam')),rid,str(inputs/(rid+'_reference-correctOrient.fasta')),guide,rid,str(inputs/(rid+'_template.fasta')),'no',str(inputs/(rid+'_summary.csv')),'-3','false']
    t=time.monotonic()
    with (d/'run.log').open('w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    return {'run':rid,'exit':p.returncode,'seconds':time.monotonic()-t,'command':cmd,'control_table_unchanged':(d/(rid+'_edits.csv')).read_bytes()==(inputs/(rid+'_edits.csv')).read_bytes() if p.returncode==0 else None}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool: rows=list(pool.map(execute,ids))
(out/'EXECUTION.json').write_text(json.dumps(rows,indent=2))
assert all(x['exit']==0 for x in rows),rows
assert all(x['control_table_unchanged'] for x in rows if x['run'] not in ['SRR29623270','SRR29623271','SRR29623272'])
print(json.dumps(rows,indent=2))

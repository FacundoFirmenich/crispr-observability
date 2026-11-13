from pathlib import Path
import json,subprocess,sys,time,os,shutil
ROOT=Path(__file__).resolve().parents[2]
IMAGE={'crispresso':'pinellolab/crispresso2:v2.3.4','amplican':'quay.io/biocontainers/bioconductor-amplican@sha256:d32eebb4c7ebe251d81b24ad142f23f65699ccfc74b98efd24486e6021b2e4c1','crispra':'msanvicente/crisprgareq@sha256:1ff4d7fdcff07a4e4da6eb88ed1676284abe96d4444ee292822e08e33f6acf0d'}
def main(engine,work):
 w=Path(work).resolve();e=w/'evidence';image=IMAGE[engine];start=time.monotonic();steps=[]
 def run(cmd,cwd=None):
  begin=time.monotonic();p=subprocess.run(cmd,cwd=cwd);steps.append({'command':cmd,'seconds':time.monotonic()-begin,'exit':p.returncode});(e/'EXECUTION.json').write_text(json.dumps({'engine':engine,'image':image,'elapsed_seconds':time.monotonic()-start,'steps':steps},indent=2));p.check_returncode()
 run(['docker','pull',image]);(e/'IMAGE.json').write_bytes(subprocess.check_output(['docker','inspect',image]))
 container=['docker','run','--rm','-v',str(w)+':'+str(w),'-w',str(w),image]
 if engine=='crispresso':
  con=json.loads((e/'CONSTRUCT.json').read_bytes());ref=json.loads((e/'REFERENCE.json').read_bytes());seq=ref['sequence'].translate(str.maketrans('ACGT','TGCA'))[::-1]
  for row in json.loads((e/'INPUT_MANIFEST.json').read_bytes()):
   rid=row['run_accession'];run(container+['CRISPResso','-r1',str(w/'inputs'/(rid+'_R1.fastq.gz')),'-a',seq,'-an','observed_mutant','-n',rid,'-o',str(e/'runs'),'--n_processes','1','--write_detailed_allele_table','--fastq_output','--suppress_report','--suppress_plots','--prime_editing_pegRNA_spacer_seq',con['spacer_RNA_5to3'],'--prime_editing_pegRNA_extension_seq',con['RTT_RNA_5to3']+con['PBS_RNA_5to3'],'--prime_editing_pegRNA_scaffold_seq',con['scaffold_RNA_5to3'],'--prime_editing_nicking_guide_seq',con['nsgRNA_spacer_RNA_5to3']])
 elif engine=='amplican':
  shutil.copy2(Path(__file__).parent/'run_amplican.R',w/'run_amplican.R');run(container+['Rscript',str(w/'run_amplican.R'),str(w)])
 else:
  source=w/'upstream';run(['git','clone','--quiet','https://bitbucket.org/synbiolab/crispr-a_nextflow.git',str(source)]);run(['git','checkout','--quiet','e5c6c8df1711df6de1f9c3c832ee297d5c48d18d'],source)
  for f in (source/'bin').rglob('*'):
   if f.is_file():f.chmod(f.stat().st_mode|0o111)
  rows=json.loads((e/'INPUT_MANIFEST.json').read_bytes());ids=[r['run_accession'] for r in rows];con=json.loads((e/'CONSTRUCT.json').read_bytes());guide=con['spacer_RNA_5to3'].replace('U','T')
  values={'outDir':str(e/'native'),'input':str(w/'inputs/*_{R1,R2}.fastq.gz'),'r2file':False,'inSize':-1,'umiClustering':False,'umiType':'TV','minlen':55,'maxlen':57,'id':.99,'ubs':1,'error':.2,'simGE':False,'simGEfilesPath':str(source/'bin/SimGE'),'simGEid':'unused','simGEsequence':'','cutSite':'','referenceFasta':str(w/'refs/*'),'indexHuman':'unused','indexMouse':'unused','refOrganism':[[r,'other'] for r in ids],'gRNAseq':[[r,guide] for r in ids],'selfRef':[[r,False] for r in ids],'aligner':'minimap2','template_seq':str(w/'templates/*'),'mock':False,'template':True,'spike':'no','samplesNames_file':str(source/'samples.php'),'protCutSite':[[r,-3] for r in ids]}
  # Native DSL1 stages both index paths even with selfRef=false. Distinct empty
  # directories prevent a staging-name collision; no genome index is consumed.
  for key in ['indexHuman','indexMouse']:
   unused=w/('disabled_'+key);unused.mkdir(exist_ok=True);values[key]=str(unused)
  cfg='params {\n'+'\n'.join(k+' = '+json.dumps(v) for k,v in values.items())+'\n}\n'+"process.container = '"+image+"'\ndocker.enabled = true\nsingularity.enabled = false\nprocess.maxForks = 2\nexecutor.queueSize = 2\n"
  (e/'native.config').write_text(cfg)
  try:run([str(w/'nextflow'),'run','crispr-a.nf','-c',str(e/'native.config'),'-work-dir',str(w/'nf-work'),'-with-trace',str(e/'process-trace.tsv')],source)
  finally:
   if (source/'.nextflow.log').exists():shutil.copy2(source/'.nextflow.log',e/'nextflow.log')
 (e/'SUCCESS.json').write_text(json.dumps({'engine':engine,'status':'COMPLETE','elapsed_seconds':time.monotonic()-start},indent=2))
if __name__=='__main__':main(*sys.argv[1:])

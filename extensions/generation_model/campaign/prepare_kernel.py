"""Build a self-contained Kaggle source from the exact repository files."""
import base64,hashlib,io,json,zipfile
from pathlib import Path
R=Path(__file__).resolve().parents[1];OUT=R/'campaign/kaggle';OUT.mkdir(exist_ok=True)
paths=[R/'pyproject.toml',R/'protocols/MP27.json',R/'campaign/input_binding.json',R/'campaign/heldout_generator.py',R/'campaign/run_campaign.py']+sorted((R/'src/crispr_extended').glob('*.py'))
buf=io.BytesIO()
with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as z:
    for p in paths:
        zi=zipfile.ZipInfo(p.relative_to(R).as_posix(),date_time=(2026,11,10,0,0,0));zi.compress_type=zipfile.ZIP_DEFLATED;z.writestr(zi,p.read_bytes())
payload=buf.getvalue();encoded=base64.b64encode(payload).decode();digest=hashlib.sha256(payload).hexdigest()
source=f'''# MIT. Frozen source bundle, no credentials or manuscript.
import base64,hashlib,io,pathlib,subprocess,sys,zipfile
blob=base64.b64decode({encoded!r})
assert hashlib.sha256(blob).hexdigest()=={digest!r}
root=pathlib.Path('/kaggle/working/extended_source');root.mkdir(exist_ok=True)
with zipfile.ZipFile(io.BytesIO(blob)) as z:
 for name in z.namelist():
  dest=root/name
  assert dest.resolve().is_relative_to(root.resolve())
  dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(name))
try:import pysam
except ImportError:subprocess.run([sys.executable,'-m','pip','install','--no-cache-dir','pysam==0.23.3'],check=True)
subprocess.run([sys.executable,str(root/'campaign/run_campaign.py')],check=True)
'''
(OUT/'kernel.py').write_text(source,encoding='utf8',newline='\n')
metadata={'id':'fjfmad/crispr-extended-mp27-generation-v1','title':'CRISPR Extended MP27 Generation V1','code_file':'kernel.py','language':'python','kernel_type':'script','is_private':True,'enable_gpu':False,'enable_internet':True,'dataset_sources':[],'competition_sources':[],'kernel_sources':['fjfmad/crispr-p1-mp24-context-20261110','fjfmad/crispr-p1-mp25-heldout-20261110']}
(OUT/'kernel-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
receipt={'bundle_sha256':digest,'kernel_source_sha256':hashlib.sha256(source.encode()).hexdigest(),'files':[{'path':p.relative_to(R).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]}
(OUT/'SOURCE_FREEZE.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'bundle_sha256':digest,'files':len(paths)}))

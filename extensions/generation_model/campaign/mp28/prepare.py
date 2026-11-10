"""MIT. Deterministic, self-contained MP28 Kaggle source bundle."""
import base64,hashlib,io,json,zipfile
from pathlib import Path
R=Path(__file__).resolve().parents[2];OUT=R/'campaign/mp28/kaggle';OUT.mkdir(exist_ok=True)
binding=json.loads((R/'campaign/input_binding.json').read_text())
binding={k:v for k,v in binding.items() if k!='MP24'}
binding['MP27']={'manifest_sha256':'e9943f12e41ac0d9233dc41fbbb013febc251671554c7110d497f7df8a914c81','selection_sha256':'89649c750b9aaa49c59ac741ffa161f894eddcb97a68a586e7a3df9be61724ce'}
(OUT.parent/'input_binding.json').write_text(json.dumps(binding,indent=2)+'\n',encoding='utf8',newline='\n')
paths=[R/'protocols/MP28.json',OUT.parent/'input_binding.json',OUT.parent/'generator.py',OUT.parent/'run.py']+sorted((R/'src/crispr_extended').glob('*.py'))
b=io.BytesIO()
with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
    for p in paths:
        zi=zipfile.ZipInfo(p.relative_to(R).as_posix(),date_time=(2026,11,10,0,0,0));zi.compress_type=zipfile.ZIP_DEFLATED;z.writestr(zi,p.read_text(encoding='utf8').replace('\r\n','\n').encode())
raw=b.getvalue();digest=hashlib.sha256(raw).hexdigest()
source=f'''# MIT. Frozen MP28 bundle; no manuscript or credentials.
import base64,hashlib,io,pathlib,subprocess,sys,zipfile
blob=base64.b64decode({base64.b64encode(raw).decode()!r})
assert hashlib.sha256(blob).hexdigest()=={digest!r}
root=pathlib.Path('/kaggle/working/mp28_source');root.mkdir(exist_ok=True)
with zipfile.ZipFile(io.BytesIO(blob)) as z:
 for name in z.namelist():
  dest=root/name
  assert dest.resolve().is_relative_to(root.resolve())
  dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(name))
try:import pysam
except ImportError:subprocess.run([sys.executable,'-m','pip','install','--no-cache-dir','pysam==0.23.3'],check=True)
subprocess.run([sys.executable,str(root/'campaign/mp28/run.py')],check=True)
'''
(OUT/'kernel.py').write_text(source,encoding='utf8',newline='\n')
metadata={'id':'fjfmad/crispr-mp28-locus-calibration-v1','title':'CRISPR MP28 Locus Calibration V1','code_file':'kernel.py','language':'python','kernel_type':'script','is_private':True,'enable_gpu':False,'enable_internet':True,'dataset_sources':[],'competition_sources':[],'kernel_sources':['fjfmad/crispr-p1-mp25-heldout-20261110','fjfmad/crispr-extended-mp27-generation-v1']}
(OUT/'kernel-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
receipt={'bundle_sha256':digest,'kernel_source_sha256':hashlib.sha256(source.encode()).hexdigest(),'normalization':'UTF8 LF','files':[{'path':p.relative_to(R).as_posix(),'sha256':hashlib.sha256(p.read_text(encoding='utf8').replace('\r\n','\n').encode()).hexdigest()} for p in paths]}
(OUT/'SOURCE_FREEZE.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'bundle_sha256':digest,'files':len(paths)}))

"""Run the HG003/HG004 replica with its own frozen inputs and checkpoints."""

import zipfile
from pathlib import Path


bundle=Path('/content/hg003_hg004_inputs.zip')
with zipfile.ZipFile(bundle) as archive:
    runner=archive.read('run_panel.py')
namespace={'__name__':'frozen_runner'}
exec(compile(runner,'run_panel.py','exec'),namespace)
namespace['ROOT']=Path('/content/hg003_hg004_homopolymer')
namespace['main'](bundle)

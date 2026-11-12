"""Embed the frozen CFTR coordinate audit and existing modeled atomic coordinates."""
from collections import defaultdict
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
prior = ROOT.parent / 'CRISPR_ATOMIC_DNA_REVIEW_20261112'
audit_path = ROOT.parent / 'CRISPR_VISUAL_REVIEW_20261112/CFTR_VISUAL_SOURCE_RECEIPT.json'
model = json.loads((prior / 'MODEL_INPUTS.json').read_bytes())
ccd = model['CCD']
audit = json.loads(audit_path.read_bytes())

def molecule(path):
    points, residue_atoms, strands = [], {}, defaultdict(list)
    for line in path.read_text(encoding='ascii').splitlines():
        if not line.startswith('ATOM  ') or line[16] not in ' A':
            continue
        name, residue, strand, res_no = line[12:16].strip(), line[17:20].strip(), line[21], line[22:27]
        element = line[76:78].strip() or name[0]
        if residue not in ccd or element not in {'C', 'N', 'O', 'P'}:
            continue
        if name in {'O1P','O2P','O3P'}:
            name = 'OP' + name[1]
        key = (strand, res_no)
        if key not in residue_atoms:
            residue_atoms[key] = (residue, {})
            strands[strand].append(key)
        idx = len(strands[strand]) - 1
        residue_atoms[key][1][name] = len(points)
        points.append([*[float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))],
                       list(strands).index(strand), idx, element, name])
    assert len(strands) == 2
    a,b = strands.values()
    seq = [''.join(residue_atoms[key][0][-1] for key in chain) for chain in (a,b)]
    assert seq[0].translate(str.maketrans('ACGT','TGCA'))[::-1] == seq[1]
    bonds = []
    for chain in (a,b):
        for i,key in enumerate(chain):
            residue, atoms = residue_atoms[key]
            for x,y in ccd[residue]['bonds']:
                if x in atoms and y in atoms:
                    bonds.append([atoms[x], atoms[y]])
            if i+1 < len(chain):
                bonds.append([atoms["O3'"],residue_atoms[chain[i+1]][1]['P']])
    assert (len(points),len(bonds)) == (1716,1924)
    assert seq[0] == audit['cases'][0]['segment'][:42] == model['CFTR']['sequence']
    return {'points':points,'bonds':bonds,'seq':seq[0],
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

payload={'audit':audit,'molecule':molecule(prior/'CFTR_42.pdb'),
         'audit_sha256':hashlib.sha256(audit_path.read_bytes()).hexdigest()}
(ROOT/'data.json').write_text(json.dumps(payload,separators=(',',':'))+'\n',encoding='utf-8')

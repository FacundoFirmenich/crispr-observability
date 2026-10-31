"""Match three named contextual compounds by frozen structure keys, not name rank.

Registry identity only: no docking, affinity, efficacy or dose inference.
"""
from pathlib import Path
import argparse,csv,json

ROOT=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/chemical');args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    rows=[];excluded=[]
    for name in ['elexacaftor','tezacaftor','ivacaftor']:
        p=json.loads((ROOT/'sources'/f'pubchem_{name}.raw.json').read_text())['PropertyTable']['Properties']
        assert len(p)==1
        p=p[0]
        candidates=json.loads((ROOT/'sources'/f'chembl_{name}.raw.json').read_text())['molecules']
        matches=[c for c in candidates if (c.get('molecule_structures')or{}).get('standard_inchi_key')==p['InChIKey']]
        assert len(matches)==1,(name,len(matches))
        c=matches[0]
        exact=json.loads((ROOT/'sources'/f'chembl_exact_{name}.raw.json').read_text())
        assert exact['molecule_chembl_id']==c['molecule_chembl_id']
        assert exact['molecule_structures']['standard_inchi_key']==p['InChIKey']
        assert exact['molecule_properties']['full_molformula']==p['MolecularFormula']
        assert exact['pref_name'].lower()==name
        rows.append({'compound':name,'PubChem_CID':p['CID'],'ChEMBL_ID':exact['molecule_chembl_id'],'molecular_formula':p['MolecularFormula'],'InChIKey':p['InChIKey'],'PubChem_MW_g_mol_reported':p['MolecularWeight'],'ChEMBL_MW_g_mol_reported':exact['molecule_properties']['full_mwt'],'identity_keys_agree':True,'PubChem_source':f"https://pubchem.ncbi.nlm.nih.gov/compound/{p['CID']}",'ChEMBL_source':f"https://www.ebi.ac.uk/chembl/api/data/molecule/{exact['molecule_chembl_id']}.json"})
        for other in candidates:
            if other['molecule_chembl_id']==exact['molecule_chembl_id']:continue
            excluded.append({'search_name':name,'excluded_ChEMBL_ID':other['molecule_chembl_id'],'pref_name':other.get('pref_name'),'same_reported_formula':(other.get('molecule_properties')or{}).get('full_molformula')==p['MolecularFormula'],'same_standard_InChIKey':False,'standard_InChIKey':(other.get('molecule_structures')or{}).get('standard_inchi_key'),'decision':'Different returned structure key; not substituted for the named compound. Search results are not an exhaustive chemical universe.'})
    with (args.output/'CHEMICAL_COMPARATOR_IDENTITIES.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report={'scope':'Post hoc frozen chemical registry identity context for three named compounds; not experimental pharmacology or molecular validation.','matched_compounds':len(rows),'matched_by':'Exact full standard InChIKey plus named record and molecular formula; ChEMBL parent record retrieved separately. Database agreement is not independent experimental validation.','reported_mass_precision':'Retained per provider; displayed rounding differences are not a chemical identity discrepancy.','excluded_search_records':excluded,'new_docking_or_simulation':False,'efficacy_affinity_dose_or_genetic_correction_inferred':False}
    (args.output/'CHEMICAL_IDENTITY_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','matched_compounds':len(rows),'excluded_nonmatching_search_records':len(excluded)}))

if __name__=='__main__':main()

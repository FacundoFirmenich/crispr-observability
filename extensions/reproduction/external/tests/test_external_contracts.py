"""Non-simulated contract tests using the acquired complete external catalogs."""
import json,subprocess,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'software/src'))
from observed_contracts.tabular import audit,require_claim

def rows(case):return json.loads((ROOT/'results'/f'{case}_CATALOG.json').read_text())['rows']

class ExternalContracts(unittest.TestCase):
    def test_unmodified_nhej_does_not_imply_complete_template(self):
        c=audit(rows('nhej'),'source_status_and_reference','exact_declared_template')
        self.assertFalse(c['complete_catalog_authorized'])
        w=c['counterexamples'][0]
        self.assertEqual(w['representation_value'],['UNMODIFIED','Reference'])
        self.assertEqual(w['target_values'],[False,True])

    def test_unresolved_allele_cannot_be_certified_by_homogeneous_resolved_rows(self):
        c=audit(rows('allele_specific'),'observed_sequence','exact_declared_template')
        self.assertTrue(c['factors_on_resolved_subcatalog'])
        self.assertEqual(c['target_unresolved_weight'],572)
        with self.assertRaises(ValueError):require_claim(c,'finite_external_assignment_catalog')

    def test_constant_target_pass_exposes_its_lack_of_discrimination(self):
        c=audit(rows('base_editor'),'source_status_and_reference','exact_declared_template')
        self.assertTrue(c['complete_catalog_authorized'])
        self.assertTrue(c['target_constant_on_resolved_domain'])
        self.assertEqual(c['distinct_resolved_target_values'],1)
        with self.assertRaises(ValueError):require_claim(c,'biological_accuracy')

    def test_cli_accepts_finite_task_and_rejects_truth_promotion(self):
        base=[sys.executable,'-m','observed_contracts.tabular',str(ROOT/'results/nhej_CATALOG.json'),'--representation','observed_sequence','--target','source_status','--require-scope']
        import os
        env=dict(os.environ,PYTHONPATH=str(ROOT/'software/src'))
        a=subprocess.run(base+['finite_external_assignment_catalog'],env=env,capture_output=True,text=True)
        b=subprocess.run(base+['biological_accuracy'],env=env,capture_output=True,text=True)
        self.assertEqual(a.returncode,0,a.stderr)
        self.assertEqual(b.returncode,2,b.stderr)
        self.assertEqual(json.loads(a.stdout)['claim_authorization'],'ACCEPTED_WITHIN_ASSIGNMENT_CATALOG')
        self.assertEqual(json.loads(b.stdout)['claim_authorization'],'REJECTED')

if __name__=='__main__':unittest.main()

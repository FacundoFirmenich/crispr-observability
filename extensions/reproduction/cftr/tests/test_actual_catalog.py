"""Contract tests use only actual coding rows and their real counterexamples."""
from pathlib import Path
import json, unittest
from observed_contracts.core import factorization, compile_partition, require_scope, ScopeError

R = Path(__file__).resolve().parents[1]
DATA = json.loads((R/'results/CODING_CATALOGS.json').read_text())
ROWS = DATA['catalogs']['witness']['ALL']
EXACT = DATA['exact_programmed_coding']

class ActualCatalogContracts(unittest.TestCase):
    def test_coarsest_joint_partition_uses_three_real_classes(self):
        p=compile_partition(ROWS,['native_event','exact_programmed_event'],EXACT)
        self.assertEqual(p['minimum_observed_classes'],3)
        sizes={tuple(g['target_signature']):g['records'] for g in p['classes']}
        self.assertEqual(sizes,{(False,False):157745,(True,False):128,(True,True):12852})

    def test_peptide_preserves_native_but_not_exact_dna(self):
        good = factorization(ROWS,'peptide_or_frame_status','native_event',EXACT)
        bad = factorization(ROWS,'peptide_or_frame_status','exact_programmed_event',EXACT)
        self.assertTrue(require_scope('finite_observed_catalog',good))
        self.assertFalse(bad['factors_on_observed_catalog'])
        self.assertTrue(bad['counterexamples'])
        with self.assertRaises(ScopeError): require_scope('finite_observed_catalog',bad)

    def test_valid_finite_certificate_cannot_authorize_biological_claim(self):
        good = factorization(ROWS,'observed_coding','native_event',EXACT)
        for scope in ['full_allele','biological_function','future_sufficiency']:
            with self.assertRaises(ScopeError): require_scope(scope,good)

    def test_exact_binary_label_does_not_preserve_native_binary_label(self):
        bad = factorization(ROWS,'exact_programmed_event','native_event',EXACT)
        self.assertFalse(bad['factors_on_observed_catalog'])
        a,b = bad['counterexamples'][0]['actual_records']
        self.assertNotEqual(a['coding'],b['coding'])

if __name__ == '__main__': unittest.main()

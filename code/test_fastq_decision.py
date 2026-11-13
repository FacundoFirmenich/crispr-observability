"""Full authentic-data regression; no invented biological observations."""
from pathlib import Path
import json,unittest
from fastq_decision import count,ROOT
class AuthenticFASTQ(unittest.TestCase):
    def test_all_nine_libraries(self):
        c=json.loads((ROOT/'examples/observation-report/cftr_fastq_config.json').read_bytes())
        rows=json.loads((ROOT/'raw_replay/expected/M1_TARGET_COUNTS.json').read_bytes())
        selected=[r for r in rows if r['policy']=='common_raw_meanQ30' and r['event']=='CTT_plus_full_SE2' and r['trim']==0]
        self.assertEqual(len(selected),9)
        for r in selected:
            got=count(ROOT/'raw_replay/inputs/fastq'/(r['run_accession']+'.fastq.gz'),c)
            self.assertEqual({k:got[k] for k in ['N','E','NE','U']},{k:r[k] for k in ['N','E','NE','U']},r['run_accession'])
    def test_invalid_configuration(self):
        from fastq_decision import validate
        with self.assertRaises(ValueError):validate({'left_anchor':'','right_anchor':'A','exact_insert':'A','max_span':3,'min_mean_phred':30,'reference':'test'})
if __name__=='__main__':unittest.main()

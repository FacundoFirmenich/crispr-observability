import copy, json, unittest
from fractions import Fraction
from observability_report import analyse, demo_input, utility

class ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.data=demo_input()
    def test_real_frozen_contrast(self):
        r=analyse(self.data)
        self.assertEqual(Fraction(r['contrast']['lower']['exact']),Fraction(3563,20573)-Fraction(702,17177))
    def test_denominator_mutation_rejected(self):
        d=copy.deepcopy(self.data);d['treated']['N']+=1
        with self.assertRaises(ValueError):analyse(d)
    def test_boolean_not_count(self):
        d=copy.deepcopy(self.data);d['control']['E']=True
        with self.assertRaises(ValueError):analyse(d)
    def test_same_universe_rejected(self):
        d=copy.deepcopy(self.data);d['control']['universe_id']=d['treated']['universe_id']
        with self.assertRaises(ValueError):analyse(d)
    def test_event_mismatch_rejected(self):
        d=copy.deepcopy(self.data);d['comparability']='different_event'
        with self.assertRaises(ValueError):analyse(d)
    def test_all_unresolved_preserves_bounds(self):
        d=copy.deepcopy(self.data)
        for arm in ['treated','control']:d[arm].update(E=0,NE=0,U=d[arm]['N'])
        r=analyse(d);self.assertEqual(r['contrast']['lower']['exact'],'-1');self.assertEqual(r['contrast']['upper']['exact'],'1')
        self.assertIsNone(r['arms'][0]['resolved_only_event_fraction'])
    def test_local_partition_and_abstention_loss(self):
        u=utility();self.assertEqual(sum(r['N'] for r in u['local_rows']),241)
        self.assertEqual(u['totals'],{'initial':{'C':207,'D':26,'U':8},'final':{'C':94,'D':1,'U':146}})
        self.assertEqual(u['decision_loss']['break_even_c_D_over_c_U'],'138/25')
        empty=[r for r in u['local_rows'] if r['final']['C']+r['final']['D']==0]
        self.assertTrue(empty);self.assertTrue(all(r['discordance_among_retained_after'] is None for r in empty))

if __name__=='__main__':unittest.main()

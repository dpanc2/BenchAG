"""Regression tests for annotation exclusions versus genuine prediction errors."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import importlib.util
import json
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HAS_CLIENT = importlib.util.find_spec('alphagenome') is not None

@unittest.skipUnless(HAS_CLIENT, 'Requires the benchmark AlphaGenome environment')
class AnnotationSkipTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = ROOT / 'benchmarks/expression/scripts'
        sys.path.insert(0, str(folder))
        spec = importlib.util.spec_from_file_location('benchmark_under_test', folder / 'Alphagenome_prediction_benchmarking_json.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
    def fake_result(self, gene=None):
        return SimpleNamespace(n_obs=0 if gene is None else 1,
            X=np.empty((0,1)) if gene is None else np.array([[4.]]),
            obs=pd.DataFrame() if gene is None else pd.DataFrame({'gene_id':[gene+'.1']}),
            var=pd.DataFrame({'strand':['+'],'name':['track']}))
    def test_missing_annotation_saves_available_and_resumes_without_requery(self):
        m=self.module
        calls=[]
        def query(interval, **kwargs):
            calls.append(interval.name)
            return [self.fake_result('available' if interval.name=='available' else 'other')]
        model=SimpleNamespace(score_interval=query)
        intervals={g:SimpleNamespace(name=g,strand='+') for g in ['available','absent']}
        with tempfile.TemporaryDirectory() as tmp, patch.object(m.dna_client,'create',return_value=model), patch.object(m,'ONTOLOGY_LIST',['track']), patch.object(m,'EXPERIMENT_ONTOLOGY_DICT',{'track':'sample'}):
            cp=str(Path(tmp)/'checkpoint.csv')
            result=m.get_fold_intervals_for_intervals('FOLD_1','unused',intervals,None,checkpoint_path=cp,checkpoint_every=1)
            self.assertEqual(result.columns.tolist(),['available'])
            self.assertEqual(float(result.loc['sample','available']),4.)
            self.assertEqual(json.loads(Path(cp+'.outcomes.json').read_text())['absent']['status'],'annotation_missing')
            self.assertTrue(Path(cp).exists())
            calls.clear()
            again=m.get_fold_intervals_for_intervals('FOLD_1','unused',intervals,None,checkpoint_path=cp)
            self.assertEqual(calls,[])
            self.assertEqual(again.columns.tolist(),['available'])
    def test_empty_annotation_without_gene_id_is_skipped(self):
        m=self.module;outcomes={}
        model=SimpleNamespace(score_interval=lambda *a,**k:[self.fake_result()])
        self.assertIsNone(m.get_gene_scores_intervals('absent',SimpleNamespace(strand='+'),model,['track'],outcomes))
        self.assertEqual(outcomes['absent']['status'],'annotation_missing')
    def test_real_api_error_is_not_annotation_exclusion(self):
        m=self.module;outcomes={}
        def fail(*a,**k):raise RuntimeError('simulated service failure')
        model=SimpleNamespace(score_interval=fail)
        self.assertIsNone(m.get_gene_scores_intervals('gene',SimpleNamespace(strand='+'),model,['track'],outcomes))
        self.assertEqual(outcomes['gene']['status'],'error')

if __name__=='__main__': unittest.main()

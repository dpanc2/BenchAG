import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import json
import pandas as pd
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
launcher = load('launcher', ROOT / 'scripts/run_expression.py')
scorer = load('scorer', ROOT / 'benchmarks/expression/scripts/score_ct_specificity.py')

class ExpressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.cfg = yaml.safe_load((ROOT / 'configs/expression.human.example.yaml').read_text())
        self.cfg.update(annotation_release='synthetic', truth_normalization='synthetic', split_provenance='synthetic')
        for split, suffix, offset in [('test', 'test', 0), ('valid', 'val', 10)]:
            for direction, strand, num in [('forward', '+', 1), ('reverse', '-', 2)]:
                gid = f'ENSG{num+offset:05}'
                path = self.base / f'{split}.{direction}.tsv'
                pd.DataFrame(dict(gene_id=[gid+'.1'],gene_id_unversioned=[gid],chrom=['chr1'],TSS=[1000000],gene_strand=[strand])).to_csv(path, sep='\t', index=False)
                self.cfg[f'genes_{direction}_{suffix}'] = str(path)
            path = self.base / f'truth.{split}.csv'
            pd.DataFrame({'gene_id':[f'ENSG{1+offset:05}.1',f'ENSG{2+offset:05}.1'],'ENCFFa':[1,2],'ENCFFb':[3,1]}).to_csv(path,index=False)
            self.cfg[f'true_{split}_file'] = str(path)
        mapping = self.base / 'mapping.json'
        mapping.write_text(json.dumps({'track1':'ENCFFa','track2':'ENCFFb'}))
        self.cfg['ontology_json'] = str(mapping)
        selected = self.base / 'selected.csv'
        pd.DataFrame({'gene_id':['ENSG00001.1,ENSG00002.1'], 'assembly':['human'],'on_targets':['ENCFFa'],'off_targets':['ENCFFb']}).to_csv(selected,index=False)
        self.cfg['path_to_selected'] = str(selected)
        self.config = self.base / 'config.yaml'
    def check(self):
        self.config.write_text(yaml.safe_dump(self.cfg))
        return launcher.prepare(self.config, check_only=True)
    def test_valid_synthetic_inputs_and_configured_species(self):
        cfg = self.check()
        self.assertEqual(cfg['organism'],'HOMO_SAPIENS')
        self.assertIn('results/expression',cfg['output_dir'])
    def test_wrong_strand_rejected(self):
        path = Path(self.cfg['genes_reverse_test'])
        df = pd.read_csv(path,sep='\t');df['gene_strand']='+';df.to_csv(path,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'strand'):
            self.check()
    def test_missing_file_rejected(self):
        self.cfg['true_test_file']=str(self.base/'missing.csv')
        with self.assertRaises(FileNotFoundError): self.check()
    def test_duplicate_mapping_rejected(self):
        Path(self.cfg['ontology_json']).write_text(json.dumps({'a':'ENCFFa','b':'ENCFFa'}))
        with self.assertRaisesRegex(ValueError,'duplicate'): self.check()
    def test_scorer_aligns_shuffled_ids(self):
        truth = pd.DataFrame({'gene_id':['ENSG00001.1','ENSG00002.1'],'ENCFFa':[1.,4.],'ENCFFb':[3.,1.]})
        pred=truth.iloc[::-1].copy()
        score=scorer.score_predictions(truth,pred,self.cfg['path_to_selected'],need_log=False)
        self.assertAlmostEqual(score['deviation_r'],1.)
    def test_scorer_uses_common_finite_rows(self):
        path = Path(self.cfg['path_to_selected'])
        df = pd.read_csv(path);df.loc[0,'gene_id'] += ',ENSG00003.1';df.to_csv(path,index=False)
        truth = pd.DataFrame({'gene_id':['ENSG00001.1','ENSG00002.1','ENSG00003.1'],'ENCFFa':[1.,4.,np.nan],'ENCFFb':[3.,1.,2.]})
        pred = truth.fillna(3).iloc[::-1]
        score=scorer.score_predictions(truth,pred,path,need_log=False)
        self.assertEqual(score['genes_number'],2)
        self.assertAlmostEqual(score['deviation_r'],1.)
    def test_pending_task_fails_honestly(self):
        import subprocess
        proc=subprocess.run(['bash',str(ROOT/'run.sh'),'mpra'],capture_output=True,text=True)
        self.assertEqual(proc.returncode,2)
        self.assertIn('no contributed implementation',proc.stderr)

if __name__=='__main__': unittest.main()

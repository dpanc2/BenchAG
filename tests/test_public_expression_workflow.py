import importlib.util
from pathlib import Path
import json
import tempfile
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
mapping=load('public_mapping',ROOT/'scripts/prepare_expression_mapping.py')
prediction=load('public_prediction',ROOT/'scripts/predict_expression.py')
resolver=load('encode_resolver',ROOT/'benchmarks/expression/scripts/build_alphagenome_ontology_map.py')
comparison=load('public_comparison',ROOT/'scripts/compare_expression.py')

class PublicWorkflowTests(unittest.TestCase):
 def test_plain_sample_list(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'ids.txt';p.write_text('ENCFF000AAA\nENCFF000AAA\nENCSR000BBB\n')
   self.assertEqual(resolver.read_sample_rows(p),[{'id':'ENCFF000AAA'},{'id':'ENCSR000BBB'}])
 def test_parent_experiment_metadata(self):
  payload={'dataset':{'biosample_ontology':{'term_id':'CL:1','term_name':'cell'},'assay_title':'total RNA-seq'}}
  self.assertEqual(resolver.resolve_ontology(payload,{},1,0)[0::3],('CL:1','total RNA-seq'))
 def test_matching_keeps_assays_distinct_and_strands_in_same_group(self):
  samples=[{'id':'ENCFF1','ontology_curie':'CL:1','assay':'total RNA-seq'},
           {'id':'ENCFF2','ontology_curie':'CL:1','assay':'polyA plus RNA-seq'},
           {'id':'ENCFF3','ontology_curie':'CL:9','assay':'total RNA-seq'}]
  metadata=[{'name':'CL:1 '+assay,'ontology_curie':'CL:1','Assay title':assay,'strand':strand}
            for assay in ['total RNA-seq','polyA plus RNA-seq'] for strand in ['+','-']]
  groups,report=mapping.match_samples(samples,metadata)
  self.assertEqual(len(groups),2);self.assertEqual(groups['CL:1 total RNA-seq'],['ENCFF1'])
  self.assertEqual(report[-1]['status'],'no_matching_ontology_assay')
 def test_gene_input_does_not_need_truth(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'genes.tsv'
   pd.DataFrame({'gene_id':['ENSG00000001.3'],'chrom':['chr1'],'TSS':[1000000],'gene_strand':['+']}).to_csv(p,sep='\t',index=False)
   df=prediction.read_genes(p,'human')
   self.assertEqual(df.gene_id_unversioned.iloc[0],'ENSG00000001')
 def test_canonical_name_takes_priority_over_broad_assay_title(self):
  samples=[{'id':'ENCFF1','ontology_curie':'CL:1','assay':'polyA plus RNA-seq'}]
  metadata=[{'name':'CL:1 polyA plus RNA-seq','ontology_curie':'CL:1','Assay title':'RNA-seq'}]
  groups,report=mapping.match_samples(samples,metadata)
  self.assertEqual(groups,{'CL:1 polyA plus RNA-seq':['ENCFF1']})
 def test_wrong_species_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'genes.csv';pd.DataFrame({'gene_id':['ENSG00000001'],'chrom':['chr1'],'TSS':[0],'gene_strand':['+']}).to_csv(p,index=False)
   with self.assertRaisesRegex(ValueError,'organism'):prediction.read_genes(p,'mouse')
 def test_named_prediction_columns_and_truth_group_averaging(self):
  m={'CL:1 total RNA-seq':['sample1','sample2']}
  gt=pd.DataFrame({'sample1':[1.,3.],'sample2':[3.,5.]},index=['g1','g2'])
  pred=pd.DataFrame({'CL:1 total RNA-seq':[2.,4.]},index=gt.index)
  pd.testing.assert_frame_equal(comparison.group_values(gt,m),comparison.group_values(pred,m,True))
 def test_finite_pair_pearson_matches_numpy(self):
  x=np.array([[1.,2.,3.,4.,5.],[2.,3.,4.,5.,6.],[3.,4.,6.,8.,9.],[5.,6.,8.,10.,11.],[7.,9.,10.,12.,14.]])
  y=x*2+3;y[0,0]=np.nan
  for axis in [0,1]:
   actual=comparison.pairwise_corr(x,y,axis)
   np.testing.assert_allclose(actual,np.ones(5))
 @unittest.skipUnless(importlib.util.find_spec('alphagenome'),'Needs installed SDK')
 def test_public_predict_saves_ontology_names_and_skips_absent(self):
  from alphagenome.models import dna_client
  def score(interval,**kwargs):
   available=interval.name=='ENSG00000001'
   return [SimpleNamespace(n_obs=int(available),X=np.array([[4.]]) if available else np.empty((0,1)),
        obs=pd.DataFrame({'gene_id':['ENSG00000001.3']}) if available else pd.DataFrame(),
        var=pd.DataFrame({'strand':['+'],'name':['CL:1 total RNA-seq']}))]
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);g=d/'genes.csv';m=d/'mapping.json';out=d/'predictions'
   pd.DataFrame({'gene_id':['ENSG00000001.3','ENSG00000002.1'],'chrom':['chr1','chr1'],'TSS':[1000000,2000000],'gene_strand':['+','+']}).to_csv(g,index=False)
   m.write_text(json.dumps({'CL:1 total RNA-seq':'ENCFF1'}))
   argv=['predict_expression.py','--genes',str(g),'--mapping',str(m),'--out',str(out)]
   with patch.object(sys,'argv',argv),patch.dict('os.environ',{'ALPHAGENOME_API_KEY':'unused'}),patch.object(dna_client,'create',return_value=SimpleNamespace(score_interval=score)):
    prediction.main()
   df=pd.read_csv(out/'predictions.csv')
   self.assertEqual(df.columns.tolist(),['gene_id','CL:1 total RNA-seq']);self.assertEqual(len(df),1)
   self.assertEqual(json.loads((out/'checkpoint.csv.outcomes.json').read_text())['ENSG00000002']['status'],'annotation_missing')

if __name__=='__main__':unittest.main()

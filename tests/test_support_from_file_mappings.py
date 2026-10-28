import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import time
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('support_builder',ROOT/'benchmarks/expression/scripts/build_alphagenome_ontology_map.py')
builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)

class SupportPreparationTests(unittest.TestCase):
 def test_select_species_from_original_mapping_genome_column(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'file_mappings.csv';p.write_text('id,genome\nENCFF000AAA,GRCh38\nENCFF000BBB,mm10\n')
   self.assertEqual(builder.read_file_mappings(p,'human')[0]['id'],'ENCFF000AAA')
   self.assertEqual(len(builder.read_file_mappings(p,'mouse')),1)
 def test_cache_contains_raw_encode_records_not_resolved_table(self):
  with tempfile.TemporaryDirectory() as tmp:
   cache=Path(tmp)/'cache.json';cache.write_text(json.dumps({'ENCFF000AAA':{'dataset':{'biosample_ontology':{'term_id':'CL:1','term_name':'test cell'},'assay_title':'total RNA-seq'}}}))
   resolved=builder.resolve_rows([{'id':'ENCFF000AAA'}],cache,sleep=0)
   self.assertEqual(resolved[0]['ontology_curie'],'CL:1');self.assertEqual(resolved[0]['query_accession'],'ENCFF000AAA')
 def test_ontology_support_is_distinct_from_exact_assay_support(self):
  rows=[{'id':'a','ontology_curie':'CL:1','assay':'total RNA-seq','errors':''},
        {'id':'b','ontology_curie':'CL:1','assay':'polyA plus RNA-seq','errors':''}]
  metadata=[{'name':'CL:1 total RNA-seq','ontology_curie':'CL:1','strand':strand} for strand in ['+','-']]
  result=builder.support_rows(rows,metadata)
  self.assertEqual(result[0]['alphagenome_rna_track_count'],2)
  self.assertEqual(result[0]['alphagenome_ontology_assay_supported'],'yes')
  self.assertEqual(result[1]['alphagenome_supported'],'yes')
  self.assertEqual(result[1]['alphagenome_ontology_assay_supported'],'no')
 def test_permanent_encode_error_is_not_retried(self):
  session=SimpleNamespace(get=lambda *a,**k:SimpleNamespace(status_code=404))
  with patch.object(builder,'encode_session',return_value=session):
   with self.assertRaisesRegex(ValueError,'HTTP 404'):
    builder.encode_get('ENCFF000AAA',{},3,0)
 def test_complete_pipeline_starts_with_original_file_mappings(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);raw=d/'file_mappings.csv';raw.write_text('id,original_id,genome,csv\nENCFF000AAA,ENCSR000AAA,GRCh38,expression.csv\n')
   cache=d/'encode.json';cache.write_text(json.dumps({'ENCFF000AAA':{'biosample_ontology':{'term_id':'CL:1','term_name':'test cell'},'assay_title':'total RNA-seq'}}))
   meta=d/'ag.csv';meta.write_text('name,ontology_curie,strand\nCL:1 total RNA-seq,CL:1,+\nCL:1 total RNA-seq,CL:1,-\n')
   out=d/'qnorm_alphagenome_support_summary.csv'
   result=builder.run(raw,'human',out,cache,meta)
   self.assertEqual(result[0]['alphagenome_ontology_name'],'CL:1 total RNA-seq')
   self.assertTrue(out.exists())
   self.assertTrue((d/'qnorm_alphagenome_ontology_map.csv').exists())

if __name__=='__main__':unittest.main()

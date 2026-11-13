from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_official_expression as runner


class OfficialExpressionTests(unittest.TestCase):
    def setUp(self):
        self.metadata = pd.DataFrame({'name': ['A', 'A', 'B', 'Padding'],
                                      'strand': ['+', '-', '.', '.'],
                                      'ontology_curie': ['CL:1', 'CL:1', 'CL:2', '']})
        self.record = {'dna_sequence': np.broadcast_to(np.uint8(0), (1048576, 4)),
                       'start': 100, 'end': 1048676, 'chromosome': 'chr1',
                       'rna_seq': np.tile([2., 5., 7., 0.], (16, 1)),
                       'rna_seq_mask': np.array([[True, True, True, False]]),
                       'rna_seq_strand': np.array([[1, 2, 3, 0]])}
        self.start = 100 + (1048576 - 16) // 2

    def test_crop_strand_mask_and_same_base_prediction_aggregation(self):
        a = self.start
        assignments = {'g1': {'strand': '+', 'regions': [(a, a+4)], 'exon_bases': 4},
                       'g2': {'strand': '-', 'regions': [(a+5, a+10)], 'exon_bases': 5}}
        prediction = SimpleNamespace(resolution=1, interval=SimpleNamespace(chromosome='chr1', start=100, end=1048676),
                                     metadata=self.metadata.iloc[:3], values=np.broadcast_to([4., 10., 14.], (1048576, 3)))
        obs, pred = runner.aggregate_record(self.record, self.metadata, assignments, prediction)
        self.assertEqual(obs, {'g1': {'A': 2., 'B': 7.}, 'g2': {'A': 5., 'B': 7.}})
        self.assertEqual(pred, {'g1': {'A': 4., 'B': 14.}, 'g2': {'A': 10., 'B': 14.}})

    def test_missing_target_is_not_zero_and_strand_mapping_is_verified(self):
        self.record['rna_seq_mask'][0, 2] = False
        self.record['rna_seq_strand'][0, 2] = 0
        obs, _ = runner.aggregate_record(self.record, self.metadata,
                                         {'g': {'strand': '+', 'regions': [(self.start, self.start+2)], 'exon_bases': 2}})
        self.assertTrue(np.isnan(obs['g']['B']))
        self.record['rna_seq_strand'][0, 0] = 2
        with self.assertRaisesRegex(ValueError, 'strand IDs'):
            runner.validate_record(self.record, self.metadata)

    def test_wrong_coordinates_and_channels_rejected(self):
        self.record['end'] += 1
        with self.assertRaisesRegex(ValueError, 'input window'):
            runner.validate_record(self.record, self.metadata)
        self.record['end'] -= 1
        with self.assertRaisesRegex(ValueError, 'channel count'):
            runner.validate_record(self.record, self.metadata.iloc[:3])

    def test_onehot_exact_dna_and_unknown_bases(self):
        dna = np.r_[np.eye(4), np.zeros((1, 4))]
        self.assertEqual(runner.sequence_from_onehot(dna), 'ACGTN')
        dna[0] = 1
        with self.assertRaisesRegex(ValueError, 'one-hot'):
            runner.sequence_from_onehot(dna)

    def test_same_named_stranded_and_unstranded_measurements_stay_separate(self):
        self.metadata.loc[2, 'name'] = 'A'
        obs, _ = runner.aggregate_record(self.record, self.metadata,
                                         {'g': {'strand': '+', 'regions': [(self.start, self.start+2)], 'exon_bases': 2}})
        self.assertEqual(obs['g'], {'A | stranded': 2., 'A | unstranded': 7.})

    def test_only_transient_api_errors_are_retried(self):
        from unittest.mock import Mock, patch
        class TransientError(Exception):
            def code(self): return SimpleNamespace(name='UNAVAILABLE')
        model = SimpleNamespace(predict_sequence=Mock(side_effect=[TransientError(), 'prediction']))
        with patch.object(runner.time, 'sleep'):
            self.assertEqual(runner.predict_with_retries(model, 'ACGT'), 'prediction')
        self.assertEqual(model.predict_sequence.call_count, 2)
        model.predict_sequence = Mock(side_effect=ValueError('Invalid input'))
        with self.assertRaisesRegex(ValueError, 'Invalid input'):
            runner.predict_with_retries(model, 'ACGT')
        self.assertEqual(model.predict_sequence.call_count, 1)


if __name__ == '__main__':
    unittest.main()

import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


evaluation = load('evaluate_gene_expression')
aggregation = load('aggregate_gene_exons')


class GeneExpressionTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec('pyBigWig'), 'Requires optional pyBigWig')
    def test_bigwig_to_figure_end_to_end(self):
        import pyBigWig
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            gtf = []
            for i in range(12):
                strand = '+' if i % 2 == 0 else '-'
                for a, b in [(i*50+1, i*50+10), (i*50+6, i*50+15), (i*50+21, i*50+25)]:
                    gtf.append(f'chr1\tx\texon\t{a}\t{b}\t.\t{strand}\t.\tgene_id "g{i}";\n')
            (d/'genes.gtf').write_text(''.join(gtf))
            (d/'windows.bed').write_text('chr1\t0\t600\n')
            rows = []
            expected = {}
            for j in range(4):
                for strand in ['+', '-']:
                    path = d/f'track{j}_{strand}.bw'
                    values = np.array([1 + ((i//50 + 2)*(j+3)) % 19 + (100 if strand == '-' else 0)
                                       for i in range(600)], dtype=float)
                    with pyBigWig.open(str(path), 'w') as bw:
                        bw.addHeader([('chr1', 600)])
                        bw.addEntries('chr1', 0, values=values.tolist(), span=1, step=1)
                    rows.append({'track_id': f't{j}', 'strand': strand,
                                 'observed_bigwig': str(path), 'predicted_bigwig': str(path)})
                    if strand == '-': expected[f't{j}'] = values[50]
            pd.DataFrame(rows).to_csv(d/'tracks.csv', index=False)
            subprocess.run([sys.executable, str(ROOT/'scripts/aggregate_gene_exons.py'), '--gtf', str(d/'genes.gtf'),
                            '--windows', str(d/'windows.bed'), '--tracks', str(d/'tracks.csv'), '--assembly', 'synthetic',
                            '--annotation-release', 'fixture', '--split-provenance', 'fixture', '--out', str(d/'aggregated')], check=True, capture_output=True)
            obs = pd.read_csv(d/'aggregated/observed.csv', index_col=0)
            for track, value in expected.items(): self.assertEqual(obs.loc['g1', track], value)
            genes = pd.read_csv(d/'aggregated/genes.csv')
            self.assertTrue(genes.exon_bases.eq(20).all())
            subprocess.run([sys.executable, str(ROOT/'scripts/evaluate_gene_expression.py'),
                            '--predictions', str(d/'aggregated/predictions.csv'), '--observed', str(d/'aggregated/observed.csv'),
                            '--prediction-transform', 'log1p', '--observed-transform', 'log1p', '--input-kind', 'exon_coverage',
                            '--out', str(d/'evaluation')], check=True, capture_output=True)
            summary = pd.read_csv(d/'evaluation/summary.csv')
            np.testing.assert_allclose(summary['mean'], 1.)
            self.assertGreater((d/'evaluation/gene_expression.png').stat().st_size, 1000)

    def test_quantiles_known_answer_with_ties(self):
        x = np.array([[5, 4, 3], [2, 1, 4], [3, 4, 6], [4, 2, 8.]])
        expected = [[17/3, 31/6, 2], [2, 2, 3], [3, 31/6, 14/3], [14/3, 3, 17/3]]
        np.testing.assert_allclose(evaluation.quantile_normalize(x), expected)
        np.testing.assert_allclose(evaluation.normalize_expression(x).mean(axis=1), 0, atol=1e-14)

    def test_correlations_match_scipy_with_missing_and_large_offsets(self):
        rng = np.random.default_rng(17)
        x = rng.normal(size=(12, 5)) + 1e8
        y = x * 1.3 + rng.normal(size=x.shape)
        y[1, 2] = np.nan
        for axis in [0, 1]:
            actual, _ = evaluation.correlations(x, y, axis)
            for i, r in enumerate(actual):
                a, b = (x[:, i], y[:, i]) if axis == 0 else (x[i], y[i])
                mask = np.isfinite(a) & np.isfinite(b)
                self.assertAlmostEqual(r, pearsonr(a[mask], b[mask]).statistic, places=12)

    def test_perfect_predictions_and_alignment(self):
        rng = np.random.default_rng(9)
        x = pd.DataFrame(rng.normal(size=(20, 6)), index=[f'g{i}' for i in range(20)])
        summary, _, coverage = evaluation.evaluate(x.iloc[::-1, ::-1], x)
        np.testing.assert_allclose(summary['mean'], 1.)
        self.assertEqual(len(coverage['used_genes']), 20)

    def test_missing_tracks_and_genes_reported(self):
        rng = np.random.default_rng(9)
        x = pd.DataFrame(rng.normal(size=(20, 6)))
        y = x.copy(); y[0] = np.nan; y.iloc[0, 1] = np.nan
        summary, _, coverage = evaluation.evaluate(y, x)
        self.assertEqual(coverage['dropped_all_missing_tracks'], [0])
        self.assertEqual(coverage['dropped_incomplete_genes'], [0])
        self.assertEqual(summary.n_genes.iloc[0], 19)

    def test_gene_baseline_cannot_create_specificity_performance(self):
        rng = np.random.default_rng(29)
        baseline = np.arange(100)[:, None]
        observed = pd.DataFrame(baseline + rng.normal(size=(100, 12)))
        prediction = pd.DataFrame(baseline + rng.normal(size=(100, 12)))
        summary, _, _ = evaluation.evaluate(prediction, observed)
        values = summary.set_index('metric')['mean']
        self.assertGreater(values.cross_gene_raw, .995)
        self.assertLess(abs(values.cross_gene_normalized), .15)

    def test_gtf_coordinate_conversion_union_and_half_exon_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'genes.gtf'
            path.write_text('chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "g1";\n'
                            'chr1\tx\texon\t16\t25\t.\t+\t.\tgene_id "g1";\n'
                            'chr1\tx\texon\t101\t110\t.\t+\t.\tgene_id "g1";\n'
                            'chr1\tx\texon\t201\t210\t.\t+\t.\tgene_id "g1";\n')
            genes = aggregation.read_exons(path)
            assignments = aggregation.assign_windows(genes, [('chr1', 0, 30), ('chr1', 0, 120)])
            self.assertEqual(assignments['g1']['window_index'], 1)
            self.assertEqual(assignments['g1']['regions'], [(10, 25), (100, 110)])
            self.assertEqual(assignments['g1']['exon_bases'], 25)
            self.assertEqual(assignments['g1']['exon_fraction'], .75)
            self.assertEqual(aggregation.assign_windows(genes, [('chr1', 0, 15)]), {})

    def test_base_weighted_mean_missing_and_negative_coverage(self):
        signal = np.array([1, 1, 9, 9, 9, 9.])
        reader = SimpleNamespace(values=lambda chrom, a, b, numpy: signal[a:b])
        self.assertAlmostEqual(aggregation.mean_coverage(reader, 'chr1', [(0, 2), (4, 6)]), 5)
        signal[0] = np.nan
        with self.assertRaisesRegex(ValueError, 'Missing'): aggregation.mean_coverage(reader, 'chr1', [(0, 2)])
        self.assertEqual(aggregation.mean_coverage(reader, 'chr1', [(0, 2)], 'zero'), .5)
        signal[0] = -1
        with self.assertRaisesRegex(ValueError, 'nonnegative'): aggregation.mean_coverage(reader, 'chr1', [(0, 2)])

    def test_api_strand_selection_and_union_exons(self):
        # SDK mock is local; this test makes no network requests.
        import sys
        from unittest.mock import patch
        values = np.column_stack([np.arange(16), np.arange(16) + 100])
        interval = SimpleNamespace(chromosome='chr1', start=0, end=16)
        output = SimpleNamespace(rna_seq=SimpleNamespace(values=values, resolution=1, interval=interval,
                   metadata=pd.DataFrame({'name': ['track', 'track'], 'strand': ['+', '-']})))
        model = SimpleNamespace(predict_interval=lambda *a, **kw: output)
        modules = {'alphagenome': SimpleNamespace(), 'alphagenome.data': SimpleNamespace(genome=SimpleNamespace(Interval=lambda **kw: interval)),
                   'alphagenome.models': SimpleNamespace(dna_client=SimpleNamespace(OutputType=SimpleNamespace(RNA_SEQ='RNA_SEQ')))}
        genes = {'g': {'strand': '-', 'regions': [(0, 2), (10, 12)], 'exon_bases': 4}}
        tracks = pd.DataFrame({'track_id': ['track', 'track'], 'strand': ['+', '-']})
        with patch.dict(sys.modules, modules):
            result = aggregation.api_gene_means(model, ('chr1', 0, 16), genes, tracks, None)
        self.assertEqual(result['g']['track'], 105.5)


if __name__ == '__main__':
    unittest.main()

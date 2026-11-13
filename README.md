# BenchAG

BenchAG runs independent AlphaGenome benchmarks using public experimental data.

```bash
git clone https://github.com/dpanc2/BenchAG.git
cd BenchAG
```

| Benchmark | Start here |
| --- | --- |
| Gene expression | [Prepare ENCODE TPMs, predict expression, and compare](benchmarks/expression/README.md) |
| Gene expression against official targets | [Evaluate released RNA-seq coverage](docs/GENE_EXPRESSION_PAPER.md) |
| CAGI5 MPRA | [Run zero-shot DNase variant scoring](benchmarks/mpra/README.md) |
| ATAC-seq | [Prepare held-out regions and compare with ChromBPNet](benchmarks/atac/README.md) |
| eQTL | [Rebuild GTEx datasets, predict, and evaluate](benchmarks/eqtl/README.md) |
| Enhancer–promoter | [Score enhancer shuffles and calculate auPRC](benchmarks/pep/README.md) |

Run commands from the repository root. API benchmarks use
`ALPHAGENOME_API_KEY`. MPRA and enhancer–promoter inference use a CUDA GPU
and the AlphaGenome ALL_FOLDS checkpoint. See [data setup](data/README.md)
for shared inputs. Downloaded data and outputs go under `data/` and `results/`.

[License](LICENSE) · [Source attribution](docs/UPSTREAM.md)

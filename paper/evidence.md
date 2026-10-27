# Evidence needed before submission

| Claim/item | Current evidence | Next action |
|---|---|---|
| Expression experiments completed | Workbook checkpoints rows 18/20 and main method-labeled blocks, transcribed in paper/results | Recover original predictions/configs and exact filtering; no independent rerun yet |
| MPRA and ATAC experiments completed | Contributor report via lead author | Obtain code, data recipes, run manifests, Methods |
| Without-qnorm versus qnorm sensitivity | checkpoints AC:AR groups for AlphaGenome rows 18/20 | Verify qnorm group interpretation, preprocessing recipe, identical predictions and populations |
| Split mismatch with Borzoi | Lead-author observation | Quantify genomic interval overlap for exact model folds |
| Two expression definitions explored | Both code pathways plus main rows 102/103/107/108 | Recover historical configs and verify API semantics and matched populations |
| Strand/version problems | Lead-author report; hardcoded mouse call observed in current source | Recover historical artifacts; do not infer historical errors from current code alone |
| No ready-made protocol for our target | Team experience | Describe narrowly; do not claim official benchmarks do not exist |
| Quantitative superiority or robustness | Descriptive expression summaries now transcribed, without uncertainty or matched-run audit | Do not infer causal effects or cross-model ranking |
| One-command execution | Entry point implemented | Validate with synthetic inputs; real run requires data and key |
| Reference-paper style | Links received; OpenReview blocked | Read supplied PDFs when available |

Minimum next inputs: authors/affiliations, ground-truth preprocessing scripts and
accessions, historical expression results/configs, MPRA and ATAC contributions,
model fold region files, and GitHub owner/name/visibility for a hosted remote.

This revision adds expression only, following the lead author's explicit scope.
MPRA, ATAC-seq and eQTL results remain excluded even though other workbook tabs
contain such experiments. Historical specificity scores have multiple transform
variants and must not be merged. See paper/results/README.md for exact cell ranges.

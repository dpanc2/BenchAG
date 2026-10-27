# Evidence needed before submission

| Claim/item | Current evidence | Next action |
|---|---|---|
| Expression experiments completed | Lead-author report | Import result tables and exact run settings |
| MPRA and ATAC experiments completed | Contributor report via lead author | Obtain code, data recipes, run manifests, Methods |
| Raw versus qnorm sensitivity | Lead-author report of both analyses | Verify normalization recipe and paired results |
| Split mismatch with Borzoi | Lead-author observation | Quantify genomic interval overlap for exact model folds |
| Two expression definitions explored | Both code pathways exist in imported source | Recover historical configs and verify API semantics |
| Strand/version problems | Lead-author report; hardcoded mouse call observed in current source | Recover historical artifacts; do not infer historical errors from current code alone |
| No ready-made protocol for our target | Team experience | Describe narrowly; do not claim official benchmarks do not exist |
| Quantitative superiority or robustness | No verified consolidated results | Do not assert |
| One-command execution | Entry point implemented | Validate with synthetic inputs; real run requires data and key |
| Reference-paper style | Links received; OpenReview blocked | Read supplied PDFs when available |

Minimum next inputs: authors/affiliations, ground-truth preprocessing scripts and
accessions, historical expression results/configs, MPRA and ATAC contributions,
model fold region files, and GitHub owner/name/visibility for a hosted remote.

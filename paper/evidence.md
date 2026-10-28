# Evidence needed before submission

| Claim/item | Current evidence | Next action |
|---|---|---|
| Expression experiments completed | Workbook checkpoints rows 18/20 and main method-labeled blocks, transcribed in paper/results | Recover original predictions/configs and exact filtering; no independent rerun yet |
| MPRA and ATAC experiments completed | Contributor report via lead author | Obtain code, data recipes, run manifests, Methods |
| Without-qnorm versus qnorm sensitivity | checkpoints AC:AR groups for AlphaGenome rows 18/20 | Verify qnorm group interpretation, preprocessing recipe, identical predictions and populations |
| Borzoi / AlphaGenome FOLD_1 split correspondence | Borzoi README: test fold3, valid fold4; AG_versions.ipynb and installed AlphaGenome 0.7.0: valid fold3, test fold4; SDK loads Borzoi BED files | Keep Borzoi table labels explicit; audit resized prediction windows separately |
| Two expression definitions explored | Both code pathways plus main rows 102/103/107/108 | Recover historical configs and verify API semantics and matched populations |
| Strand/version problems | Lead-author report; hardcoded mouse call observed in current source | Recover historical artifacts; do not infer historical errors from current code alone |
| No ready-made protocol for our target | Team experience | Describe narrowly; do not claim official benchmarks do not exist |
| Quantitative superiority or robustness | Descriptive expression summaries now transcribed, without uncertainty or matched-run audit | Do not infer causal effects or cross-model ranking |
| One-command execution | Entry point implemented | Validate with synthetic inputs; real run requires data and key |
| Reference-paper style | Supplied ModernGENA PDF read and first page visually inspected on 2026-10-28 | Use concise research framing, concrete examples, contributions, and experimental sections; final venue template still required |

Minimum next inputs: authors/affiliations, ground-truth preprocessing scripts and
accessions, historical expression results/configs, MPRA and ATAC contributions,
model fold region files, and GitHub owner/name/visibility for a hosted remote.

This revision adds expression only, following the lead author's explicit scope.
MPRA, ATAC-seq and eQTL results remain excluded even though other workbook tabs
contain such experiments. Historical specificity scores have multiple transform
variants and must not be merged. See paper/results/README.md for exact cell ranges.

## Revision of 28 October 2026

- Removed the former sections 2.1 (shared provenance/software details) and 2.2
  (extended ground-truth discussion) from the paper body.
- Reserved one explicit placeholder for ground-truth construction.
- Replaced the introduction with scientific motivation, concrete evaluation
  problems, and a concise contribution statement, following the supplied example's
  rhetorical structure rather than its wording.
- Established the exact fold-label swap from two independent sources. The SDK
  points to Borzoi BED files at commit 5c9358222b5026abb733ed5fb84f3f6c77239b37.
- The paper keeps Borzoi labels for existing input files and result tables.
  FOLD_0's scores on these folds must not be framed as held-out improvements.
- This revision does not change numerical result cells or add MPRA/ATAC results.

TSS membership check on the current benchmark inputs (all genes matched):
human test 2,781/2,781 fold3, human validation 3,038/3,038 fold4;
mouse test 2,724/2,724 fold3, mouse validation 2,998/2,998 fold4.
No test TSS overlapped fold4 and no validation TSS overlapped fold3 in this check.
This point-membership check is not a full resized-window overlap audit.
The machine-readable check is saved as paper/results/fold_tss_membership_20261028.json.

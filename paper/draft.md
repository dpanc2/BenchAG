# BenchAG: Towards Reproducible Evaluation of AlphaGenome Across Gene Expression, Reporter Activity, and Chromatin Accessibility

**Working draft — revised 28 October 2026.**

**Authors:** [lead author: name and affiliation required], [MPRA contributor], [ATAC-seq contributor], [additional authors]. Author order and corresponding author remain to be confirmed.

## Abstract

Sequence-to-function models predict molecular phenotypes from DNA, but evaluating these predictions on an external dataset requires more than choosing an accuracy metric. We study AlphaGenome expression prediction in human and mouse, focusing on ground-truth normalization, gene-level scoring, input sequence length, and genomic partitions. Our expression experiments show that interval-based and variant-based scoring can produce substantially different across-gene correlations, while normalization changes the agreement across cellular contexts. We also establish the correspondence between the released Borzoi partitions and AlphaGenome FOLD_1: the same two folds are held out, but their test and validation labels are exchanged. BenchAG brings these comparisons into a common evaluation framework, with MPRA and ATAC-seq contributions to be added. The goal is to distinguish differences in predictive performance from differences in the quantity and population being evaluated.

## 1. Introduction

Predicting gene regulation from DNA is a central problem in computational genomics. Sequence-to-function models such as Borzoi [2] and AlphaGenome [1] address this problem by predicting experimental measurements across tissues and cell types. AlphaGenome extends this approach across multiple molecular modalities, providing a common model for studying gene expression, chromatin accessibility, and the consequences of sequence variation. These capabilities make independent evaluation important for understanding when a prediction can support a biological application.

Adapting a model to a new benchmark, however, requires choices that are not resolved by a published aggregate score. RNA-seq coverage must be converted into a gene-level expression estimate, experimental samples must be matched to predicted tracks, and genomic regions must be assigned to a suitable model version. A mismatch at any of these steps can change the interpretation of an otherwise valid correlation.

Our expression experiments illustrate several concrete difficulties. First, two available scoring pathways provide different definitions of a gene-level prediction: interval-based aggregation and a statistic returned by the variant-scoring interface. In the historical human results on the Borzoi test partition, their across-gene correlations are 0.713 and 0.548, respectively, whereas the across-context correlations are much closer (Section 3.3). Second, the representation of the experimental target matters: applying quantile normalization can affect across-context agreement differently from across-gene agreement. Third, dataset labels do not have the same meaning across models. Borzoi's test regions are AlphaGenome FOLD_1's validation regions, and Borzoi's validation regions are AlphaGenome FOLD_1's test regions [5, 6]. Selecting a model or a result file by its name alone can therefore produce an incorrect comparison. Gene identifiers, strand-specific tracks, and changes in the client or model version introduce further compatibility questions.

These issues motivate an evaluation that states what is being predicted, how it is measured, and which genomic regions are held out. We introduce BenchAG as a framework for consolidating expression, MPRA, and ATAC-seq benchmarks. This draft reports the expression component; the other task sections remain reserved for the corresponding contributions.

**Contributions.** We (i) compare alternative gene-level scoring pathways and ground-truth representations for AlphaGenome expression prediction; (ii) evaluate human and mouse expression across input sequence lengths using both across-gene and across-context correlations; and (iii) resolve the Borzoi–AlphaGenome split correspondence and provide an implementation for repeating the expression evaluation.

## 2. Experimental Setup

**Ground-truth construction.** **[To be added: describe the source expression datasets, how experimental measurements were converted to gene-level ground-truth values, how qnorm and non-qnorm versions were obtained, and how samples were matched to AlphaGenome outputs.]**

### Gene-level expression prediction

We obtain gene-level predictions by aggregating AlphaGenome RNA-seq output over annotated exonic regions around each gene's transcription start site. The retained pathway uses interval-based scoring with mean aggregation and selects tracks matching the gene's strand. We examine 16 kb and 1 Mb input configurations; the completed 128 kb experiment provides an intermediate context length. Input context and the region used for aggregation are separate settings and are kept explicit in the experimental configurations.

We also compare the interval pathway with the historical variant-scoring pathway, which uses identical reference and alternative alleles to obtain an activity statistic. These pathways summarize model output differently and are treated as distinct definitions of a gene-level prediction. Genes not returned in the AlphaGenome annotation output are recorded as exclusions; predictions for the remaining genes are retained.

### Genomic partitions and model versions

The genomic partition definitions come from the released Borzoi human and mouse BED files [5]. Borzoi's released model replicates share the same held-out partitions: **fold3 is test and fold4 is validation**. AlphaGenome's fold-interval API uses these Borzoi interval files and specifies which folds belong to each model version. We obtain this model-to-fold mapping using `get_fold_names`; the code and recorded output are provided in `AG_versions.ipynb` [6].

For **AlphaGenome FOLD_1**, the held-out folds are the same as for Borzoi, but their roles are reversed:

| Genomic fold | Borzoi subset | AlphaGenome FOLD_1 subset |
|---|---|---|
| fold3 | Test | Validation |
| fold4 | Validation | Test |

Both partitions are therefore held out from FOLD_1 training. We retain the Borzoi naming convention for the benchmark input files and result tables: **“test” denotes fold3 (AlphaGenome validation), and “validation” denotes fold4 (AlphaGenome test)**. We checked the input gene TSS coordinates against the released BED files: all test TSSs fall in fold3 and all validation TSSs fall in fold4, for both organisms. The correspondence concerns the genomic partitions; prediction context lengths and eligible gene populations remain separate experimental choices.

Other AlphaGenome model versions use different held-out folds. In particular, fold3 and fold4 are training folds for FOLD_0, so scores obtained with that version on the Borzoi evaluation regions do not constitute the same held-out comparison. The ALL_FOLDS model is trained across all eight folds and has no independent validation or test partition. We use FOLD_1 for the expression benchmark and treat results from other versions as separate diagnostic experiments.

### Evaluation metrics

We report two complementary Pearson correlations. **Gene correlation** is calculated across genes within each matched output group and then averaged over groups. **Cell correlation** is calculated across matched groups within each gene and then averaged over genes. These metrics measure agreement in gene abundance within a context and agreement in expression patterns across contexts, respectively. Only paired finite observations are used, and undefined correlations are excluded.

The historical results additionally contain a selected-target specificity score. We report that score separately where available, because applying a log transformation inside its scoring function is distinct from quantile-normalizing the experimental targets. Sample counts and eligible output groups accompany the comparisons; uncertainty estimates remain to be added for the final analysis.

### Additional benchmark tasks

**[MPRA: contributor to add assay description, prediction statistic, evaluation protocol, and results.]**

**[ATAC-seq: contributor to add data sources, prediction target, genomic partitions, and results.]**

The inclusion of an eQTL benchmark remains undecided.

## 3. Expression Results

### 3.1 Source and evaluation panels

We consolidate the expression results reported in the team's *Expression* workbook, retrieved on 27 October 2026 [4]. Values below are transcribed from the workbook rather than recomputed by a new benchmark run. The linked `checkpoints` tab separates a 14-cell-type panel from a broader AlphaGenome-matched panel. For the latter, the headers report 770 human source tracks, 240 unique ontologies and 285 ontology–assay heads, and 135 mouse source tracks, 55 unique ontologies and 68 ontology–assay heads. These are distinct levels of track matching, not interchangeable numbers of independent cell types. The exact aggregation and eligibility rules require the original preprocessing records.

Throughout this section, r_gene and r_cell preserve the workbook's `gene_corr` and `cell_corr` labels, corresponding to the across-gene and across-context correlation axes defined in Section 2. This correspondence must be checked against the historical evaluation scripts before final submission. No confidence intervals or significance tests are supplied by these summary cells. Split labels follow the Borzoi convention defined in Section 2; they must be exchanged when interpreted as AlphaGenome FOLD_1 subset names.

### 3.2 Input length and ground-truth normalization

**Table 1. AlphaGenome expression correlations in the broader matched panel, with and without qnorm ground truth.** The qnorm column group is headed “AlphaGenome” in the source; we identify it as qnorm from its contrast with “AlphaGenome (gt wo qnorm)” and the lead author's description. Values are rounded to four decimals. Source: `checkpoints!AC18:AR18` and `AC20:AR20`; input-length labels in `A18` and `A20`.

| Organism | Split | Input | r_gene, qnorm | r_cell, qnorm | r_gene, without qnorm | r_cell, without qnorm |
|---|---|---|---|---|---|---|
| Human | Validation | 16 kb | 0.5479 | 0.2211 | 0.5450 | 0.2700 |
| Human | Test | 16 kb | 0.5593 | 0.2051 | 0.5580 | 0.2590 |
| Mouse | Validation | 16 kb | 0.5930 | 0.2800 | 0.6035 | 0.4306 |
| Mouse | Test | 16 kb | 0.6030 | 0.2640 | 0.6135 | 0.4169 |
| Human | Validation | 1 Mb | 0.7040 | 0.3160 | 0.6990 | 0.4260 |
| Human | Test | 1 Mb | 0.7140 | 0.3080 | 0.7100 | 0.4130 |
| Mouse | Validation | 1 Mb | 0.7100 | 0.3900 | 0.7254 | 0.5768 |
| Mouse | Test | 1 Mb | 0.7200 | 0.3800 | 0.7349 | 0.5608 |

The reported 1 Mb configuration has higher r_gene than the 16 kb configuration in each organism and split under both target representations. With qnorm truth, human test r_gene increases from 0.5593 to 0.7140 and mouse test r_gene from 0.6030 to 0.7200. The corresponding r_cell values increase from 0.2051 to 0.3080 and from 0.2640 to 0.3800. These are descriptive comparisons of reported configurations. Context length and gene-mask width differ in the supported implementations, so the comparison does not isolate the causal contribution of additional sequence alone.

Normalization affects the correlation axes differently. For the reported 1 Mb human test configuration, r_gene is similar with qnorm and without qnorm (0.7140 and 0.7100), whereas r_cell is 0.3080 and 0.4130. For mouse test, the corresponding values are 0.7200 versus 0.7349 for r_gene and 0.3800 versus 0.5608 for r_cell. This pattern illustrates why a single correlation cannot summarize both questions. It does not establish that one normalization is universally preferable: the truth recipe, eligible population and reuse of the same predictions across the two column groups must first be confirmed.

The separate 14-cell-type panel reports 1 Mb validation/test r_gene of 0.7044/0.7163 and r_cell of 0.2380/0.2102. At 16 kb, the corresponding values are 0.5457/0.5718 and 0.1655/0.1570 (`checkpoints!O18:T18`, `O20:T20`). These numbers should remain separate from Table 1 because the target panel differs.

### 3.3 Alternative gene-level prediction definitions

The historical `main` tab explicitly distinguishes interval and variant scoring for FOLD_1. Table 2 transcribes the later block with method labels on every row, rather than combining it with earlier rounded entries.

**Table 2. Historical FOLD_1 human expression results by scoring pathway.** Source: `main!A102:J103` and `A107:J108`. `score_predictions` is preserved as a source label because its transformation setting is not specified in these rows. These runs are not assumed to share the broader panel or exact configuration in Table 1.

| Split | Scoring pathway | r_gene | r_cell | score_predictions |
|---|---|---|---|---|
| Validation | Interval | 0.7044 | 0.2383 | 0.5213 |
| Validation | Variant | 0.5816 | 0.2317 | 0.5203 |
| Test | Interval | 0.7129 | 0.4004 | 0.5449 |
| Test | Variant | 0.5482 | 0.4048 | 0.5705 |

Interval scoring yields higher reported r_gene than variant scoring on both validation (0.7044 versus 0.5816) and test (0.7129 versus 0.5482). The difference is much smaller for r_cell: 0.2383 versus 0.2317 on validation, and 0.4004 versus 0.4048 on test. Thus the choice of gene-level statistic is associated with a substantial change in the across-gene metric in these historical runs, without a uniform improvement across every metric. Confirming identical populations and model/client settings remains necessary before treating this as a controlled ablation.

The workbook also distinguishes applying a log transformation inside the specificity-scoring function from quantile normalization of the ground truth. In the FOLD_1 interval rows explicitly labeled “without log transformation in score function,” `score_predictions` is 0.6343 on test and 0.4731 on validation (`main!F111:F112`); the corresponding r_gene and r_cell entries match the interval correlations in Table 2. These values differ from 0.5449 and 0.5213 in the earlier explicitly labeled interval rows (`F107`, `F102`). We report the blocks separately because the earlier rows do not state the transform setting. The current repository's `need_log=False` setting must not be retroactively assigned to every historical score.

The same method-labeled historical block also records fold-dependent results. For interval scoring, FOLD_0 validation/test r_gene is 0.8030/0.8104 (`main!J100`, `J105`), compared with 0.7044/0.7129 for FOLD_1 (`J102`, `J107`). As established in Section 2, fold3 and fold4 are training folds for FOLD_0 but held-out folds for FOLD_1. The higher FOLD_0 scores are therefore not evidence of better held-out prediction and are excluded from the primary model comparison.

### 3.4 Remaining checks

These results establish a concrete descriptive basis for the expression section, while preserving several unresolved protocol questions. The 14-cell panel, broader matched panel, and historical `main` blocks have different reported values and must not be collapsed into a single run. The FOLD_1 correspondence is established at the level of released genomic folds. Before submission, attach the exact ground-truth construction procedure, gene/track counts after filtering, and original run settings; a fold-label correspondence does not by itself audit every resized prediction window. Paired uncertainty estimates require per-gene or per-region predictions and cannot be recovered from the summary cells alone.

MPRA, ATAC-seq and eQTL results are intentionally not included in this revision.

## 4. Discussion and Limitations

The main issue is defining what a benchmark measures. A model-to-truth correlation combines predictive ability with decisions about experimental targets, normalization, track identity, gene annotation, and evaluation population. Making these decisions executable can help distinguish a biological limitation from an evaluation mismatch. The reported expression summaries show different sensitivities across correlation axes: the historical scoring pathways differ substantially in across-gene agreement, while changing the truth representation is associated with larger differences in across-context agreement for the 1 Mb runs. Recovering matched run artifacts is needed to attribute these differences to individual protocol choices.

The initial repository has substantial limitations. It does not bundle experimental data or reproduce ground-truth construction. Supplied ontology mappings require review against exact sample metadata. The model-to-fold correspondence is established, while a full overlap audit of resized prediction contexts remains separate. Expression tables are transcribed from a shared workbook and have not been independently rerun. The MPRA and ATAC-seq contributions and quantitative tables are pending. Historical API/model changes may prevent exact reconstruction unless the original configurations and predictions are retained. These limitations motivate a provenance-centered release, but do not establish that other published benchmarks are incorrect.

## 5. Conclusion

BenchAG organizes an independent evaluation of AlphaGenome around explicit target definitions, strand-aware expression scoring, partition provenance, and recorded software/model choices. The initial draft and software provide a foundation for consolidating expression, MPRA, and ATAC-seq results. The expression summaries motivate reporting both correlation axes and recording the target representation and scoring pathway. Stronger causal claims and model rankings will depend on matched-run verification and the remaining experiments.

## References

[1] Avsec, Ž., Latysheva, N., Cheng, J., et al. 2026. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649:1206–1218. https://doi.org/10.1038/s41586-025-10014-0

[2] Linder, J., Srivastava, D., Yuan, H., et al. 2025. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. *Nature Genetics* 57:949–961. https://doi.org/10.1038/s41588-024-02053-6

[3] AIRI-Institute. GENA_LM, expression benchmark source, commit 6476f509c056ee2252fb86eb0f3e65b0ea76d8fe. https://github.com/AIRI-Institute/GENA_LM/tree/6476f509c056ee2252fb86eb0f3e65b0ea76d8fe/downstream_tasks/expression_prediction/benchmarks/alphagenome_benchmark

[4] BenchAG team. 2026. *Expression*, shared experiment workbook, `checkpoints` and `main` tabs. Retrieved 27 October 2026. https://docs.google.com/spreadsheets/d/1qdmZc7KqsdZN7jOuDYwjVRg_kPs0-17DbkmWZ-qpuec/edit?gid=1307246850#gid=1307246850

[5] Calico. *Borzoi: model availability and released genomic partitions*. The repository specifies test = fold3 and validation = fold4. https://github.com/calico/borzoi#model-availability. Human and mouse interval files used by AlphaGenome client 0.7.0: https://github.com/calico/borzoi/tree/5c9358222b5026abb733ed5fb84f3f6c77239b37/data.

[6] AlphaGenome client 0.7.0, `alphagenome.data.fold_intervals`, and BenchAG fold-inspection notebook, `benchmarks/expression/notebooks/AG_versions.ipynb`. Model/subset assignments verified on 28 October 2026.

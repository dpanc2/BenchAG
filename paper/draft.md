# BenchAG: Towards Reproducible Evaluation of AlphaGenome Across Gene Expression, Reporter Activity, and Chromatin Accessibility

**First working draft — 27 October 2026. Not submission ready.**

**Authors:** [lead author: name and affiliation required], [MPRA contributor], [ATAC-seq contributor], [additional authors]. Author order and corresponding author remain to be confirmed.

## Abstract

Sequence-to-function models offer a common interface for predicting molecular phenotypes from DNA, but independent evaluation requires decisions that are not uniquely specified by model outputs. We present BenchAG, a benchmark framework being assembled for AlphaGenome across gene expression, massively parallel reporter assays (MPRA), and ATAC-seq. Our expression workflow makes explicit the relationship between RNA-seq tracks and ground-truth samples, normalization, genomic partitions, gene-level aggregation, strand selection, and model and client versions. We describe a protocol for separating these evaluation choices from model performance and for documenting their consequences through controlled comparisons. The initial software release incorporates an expression-prediction workflow and contribution interfaces for MPRA and ATAC-seq. Reported expression results show that the 1 Mb configuration has higher across-gene correlations than the 16 kb configuration on both human and mouse evaluation splits. Ground-truth normalization and gene-level scoring definitions affect the correlation axes differently. These descriptive results are transcribed from the team's experiment workbook; matched-run verification and cross-task analyses remain pending. Our central methodological position is that a reported score should be accompanied by an executable protocol defining the target, population, and prediction statistic being evaluated.

## 1. Introduction

Predicting molecular measurements from DNA can connect regulatory sequence to cellular function. AlphaGenome predicts multiple functional genomic modalities, including RNA-seq and chromatin accessibility [1]. Borzoi provides a related sequence-based approach to RNA-seq coverage prediction [2]. For an independent evaluator, the model prediction is only one component of a benchmark: experimental measurements must also be selected, transformed, matched to output tracks, and assigned to evaluation partitions.

Our work focuses on this evaluation interface. During expression benchmarking, we encountered uncertainty about the appropriate ground-truth representation, differences between quantile-normalized and alternative expression targets, the correspondence of validation and test partitions across models, alternative definitions of a gene-level prediction, strand handling, and software/model versions. These are observations from our benchmarking process; their quantitative effects are described below from the experiment workbook, with controlled verification still pending. They should not be interpreted as evidence that official evaluations are absent or invalid. A publication and software interface do not necessarily provide a directly executable protocol for our particular target dataset.

BenchAG is intended to make these choices inspectable across three complementary tasks: endogenous gene expression, reporter activity measured by MPRA, and chromatin accessibility measured by ATAC-seq. These tasks interrogate different experimental quantities, so they require individual protocols rather than a single pooled accuracy score. The current contribution is an imported and configurable expression workflow, a shared provenance convention, and a structured plan for completing the multi-task evaluation. We ask: (i) how sensitive are expression conclusions to truth preprocessing and prediction definition, and (ii) which protocol details must be fixed to make independent evaluations reproducible?

**[Before submission: confirm the final multi-task scope and link the reported expression summaries to their original run artifacts.]**

## 2. Evaluation Framework

### 2.1 Task-specific targets and shared provenance

Each task specifies its experimental target, input sequence population, prediction statistic, partition, and metric. Every reported run should identify dataset accessions and checksums, genome assembly and annotation, preprocessing, output-track mapping, model selection, client version, code revision, and configuration. BenchAG records configuration and input hashes and isolates expression outputs and checkpoints by run identity. This supports auditing a local run; exact API replay additionally depends on the stability and identification of the remotely served model.

The expression code originates from the GENA_LM expression branch at commit `6476f509c056ee2252fb86eb0f3e65b0ea76d8fe` [3]. It includes human and mouse track mappings, prediction scripts, and comparison notebooks. Historical configurations are preserved for provenance. The supported entry point makes organism and sequence settings explicit and checks required data before inference. MPRA and ATAC-seq implementations will be supplied by the corresponding contributors; their presence in the repository structure does not imply that these benchmarks are already reproduced.

### 2.2 Expression ground truth and track correspondence

Let Y[g,t] denote the observed expression of gene g in experimental sample or context t. To compare this matrix with AlphaGenome predictions, each selected RNA-seq output track must be matched to the appropriate experimental target. Existing mappings connect output tracks to representative sample identifiers; auxiliary mappings retain multiple candidate identifiers. The choice between a representative replicate and aggregation of multiple measurements changes the target definition and must be documented.

We consider ground truth with and without quantile normalization, following the analyses undertaken by the lead author. The final Methods must specify original quantification units, any earlier normalization or log transformation, the axis of quantile normalization, tie handling, the reference population, and replicate aggregation. Until these details are established, we label alternatives by their exact preprocessing recipes rather than assuming that “raw” and “qnorm” identify universally comparable expression quantities. Transformations whose parameters are estimated from data must be defined without using held-out outcomes to tune the evaluation.

**[Required: dataset accessions, preprocessing code, annotation release, gene/sample counts, and provenance of each truth matrix.]**

### 2.3 Gene-level predictions and strand handling

The retained expression pathway uses interval scoring with RNA-seq output and a GeneMaskScorer mean aggregation. Intervals are centered on an annotated transcription start site; the selected output track strand is matched to the gene strand. The default configuration requests a 1,048,576 bp input interval and a 200,001 bp scorer width. A separate supported configuration uses 16,384 bp input and a 10,001 bp width. Input context and scorer width are distinct parameters and must be reported separately.

The historical code also contains a variant-scoring pathway using a reference-equals-alternate construction and GeneMaskActiveScorer. We retain it as an archival implementation of the second prediction definition explored during development. We do not assume that this statistic is equivalent to the interval-based expression score. To interpret its historical results as a controlled comparison, we must verify the exact client semantics, selected score field, gene matching, allele/coordinate convention, and historical configuration. The initial public runner exposes the retained interval pathway.

Strand selection is a correctness condition, not a hyperparameter to tune against test outcomes. A deliberately mismatched or merged-strand analysis may quantify the consequences of an implementation error, but must not be presented as an equally valid benchmark definition.

### 2.4 Partitions and model versions

Labels such as “test,” “valid,” or “fold 1” are insufficient to establish the same evaluation population across models. We therefore require an explicit mapping from evaluation windows to each model's training and held-out regions. The software checks disjointness of supplied gene identifiers, but this does not establish genomic separation: long sequence contexts may overlap even when gene lists do not.

Comparisons involving Borzoi should use a justified common held-out population, with context overlap and output alignment audited separately. If a common held-out population cannot be established, results should be described as performance on the specified external dataset rather than a controlled held-out model ranking. Client package version and model identifier must also be recorded separately. The initial dependency pin follows the imported workflow (AlphaGenome client 0.7.0); it is not evidence that all historical runs used the same served model.

### 2.5 Expression metrics

For prediction matrix P and aligned truth Y, we report Pearson correlation across genes for each track and average across eligible tracks. Separately, we compute correlation across tracks for each gene and average across eligible genes. These quantities answer different questions: agreement in gene abundance within a context and agreement in relative expression across contexts. Undefined correlations, including constant vectors or insufficient finite observations, are excluded with explicit denominators.

The imported workflow additionally computes a selected-target contrast correlation. With its current `need_log=False` setting, the contrast for contexts i and j is d(x_i,x_j) = (x_i − x_j)/2, pooled over selected genes and context pairs before correlating prediction and truth contrasts. This is distinct from centering each gene across all tissues. Its selected population and missing-data handling must be reported. We retain the statistic for continuity and avoid treating it as interchangeable with the other correlations.

For final comparisons we propose paired uncertainty estimates on a fixed eligible population, with resampling units chosen to respect genomic dependence. **[Specify resampling blocks, seeds, and replicate count after implementing this analysis.]** We will report mapping coverage, failed predictions, excluded genes/tracks, and the number of observations contributing to each metric alongside aggregate scores.

### 2.6 MPRA and ATAC-seq protocols

The MPRA contribution will define the reporter assay, sequence construction, activity target and normalization, replicate aggregation, biological context, prediction statistic, and split unit. Reporter activity is an experimentally distinct quantity from endogenous RNA-seq expression; mapping a genomic prediction to reporter activity must therefore be justified explicitly. **[Contributor to add assay accessions, Methods, metrics, and verified results.]**

The ATAC-seq contribution will define track selection, genome assembly, evaluation regions, resolution, aggregation, normalization, and genomic holdouts. Coverage regression and peak classification require different targets and metrics and must be distinguished. **[Contributor to identify the actual evaluation target and add Methods and results.]**

An eQTL evaluation is outside the initial scope. It may be added only after its phenotype, variant–gene mapping, allele orientation, tissue selection, and evaluation population are documented.

## 3. Expression Results

### 3.1 Source and evaluation panels

We consolidate the expression results reported in the team's *Expression* workbook, retrieved on 27 October 2026 [4]. Values below are transcribed from the workbook rather than recomputed by a new benchmark run. The linked `checkpoints` tab separates a 14-cell-type panel from a broader AlphaGenome-matched panel. For the latter, the headers report 770 human source tracks, 240 unique ontologies and 285 ontology–assay heads, and 135 mouse source tracks, 55 unique ontologies and 68 ontology–assay heads. These are distinct levels of track matching, not interchangeable numbers of independent cell types. The exact aggregation and eligibility rules require the original preprocessing records.

Throughout this section, r_gene and r_cell preserve the workbook's `gene_corr` and `cell_corr` labels, corresponding to the across-gene and across-context correlation axes described in Section 2.5. This correspondence must be checked against the historical evaluation scripts before final submission. No confidence intervals or significance tests are supplied by these summary cells.

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

The same method-labeled historical block also records fold-dependent results. For interval scoring, FOLD_0 validation/test r_gene is 0.8030/0.8104 (`main!J100`, `J105`), compared with 0.7044/0.7129 for FOLD_1 (`J102`, `J107`). These differences reinforce the need to establish how each fold relates to the evaluation regions. They do not independently demonstrate leakage or justify selecting a fold by its test score.

### 3.4 Remaining checks

These results establish a concrete descriptive basis for the expression section, while preserving several unresolved protocol questions. The 14-cell panel, broader matched panel, and historical `main` blocks have different reported values and must not be collapsed into a single run. Validation and test labels also do not establish exclusion from each model's training regions. Before submission, attach exact input matrices, gene/track counts after filtering, fold-region audits, normalization recipes, model/client versions, and run manifests. Paired uncertainty estimates require per-gene or per-region predictions and cannot be recovered from the summary cells alone.

MPRA, ATAC-seq and eQTL results are intentionally not included in this revision.

## 4. Discussion and Limitations

The main issue is defining what a benchmark measures. A model-to-truth correlation combines predictive ability with decisions about experimental targets, normalization, track identity, gene annotation, and evaluation population. Making these decisions executable can help distinguish a biological limitation from an evaluation mismatch. The reported expression summaries show different sensitivities across correlation axes: the historical scoring pathways differ substantially in across-gene agreement, while changing the truth representation is associated with larger differences in across-context agreement for the 1 Mb runs. Recovering matched run artifacts is needed to attribute these differences to individual protocol choices.

The initial repository has substantial limitations. It does not bundle experimental data or reproduce ground-truth construction. Supplied ontology mappings require review against exact sample metadata. Gene-list checks do not certify training-region independence. Expression tables are transcribed from a shared workbook and have not been independently rerun. The MPRA and ATAC-seq contributions and quantitative tables are pending. Historical API/model changes may prevent exact reconstruction unless the original configurations and predictions are retained. These limitations motivate a provenance-centered release, but do not establish that other published benchmarks are incorrect.

## 5. Conclusion

BenchAG organizes an independent evaluation of AlphaGenome around explicit target definitions, strand-aware expression scoring, partition provenance, and recorded software/model choices. The initial draft and software provide a foundation for consolidating expression, MPRA, and ATAC-seq results. The expression summaries motivate reporting both correlation axes and recording the target representation and scoring pathway. Stronger causal claims and model rankings will depend on matched-run verification and the remaining experiments.

## References

[1] Avsec, Ž., Latysheva, N., Cheng, J., et al. 2026. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649:1206–1218. https://doi.org/10.1038/s41586-025-10014-0

[2] Linder, J., Srivastava, D., Yuan, H., et al. 2025. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. *Nature Genetics* 57:949–961. https://doi.org/10.1038/s41588-024-02053-6

[3] AIRI-Institute. GENA_LM, expression benchmark source, commit 6476f509c056ee2252fb86eb0f3e65b0ea76d8fe. https://github.com/AIRI-Institute/GENA_LM/tree/6476f509c056ee2252fb86eb0f3e65b0ea76d8fe/downstream_tasks/expression_prediction/benchmarks/alphagenome_benchmark

[4] BenchAG team. 2026. *Expression*, shared experiment workbook, `checkpoints` and `main` tabs. Retrieved 27 October 2026. https://docs.google.com/spreadsheets/d/1qdmZc7KqsdZN7jOuDYwjVRg_kPs0-17DbkmWZ-qpuec/edit?gid=1307246850#gid=1307246850

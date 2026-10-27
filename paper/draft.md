# BenchAG: Towards Reproducible Evaluation of AlphaGenome Across Gene Expression, Reporter Activity, and Chromatin Accessibility

**First working draft — 27 October 2026. Not submission ready.**

**Authors:** [lead author: name and affiliation required], [MPRA contributor], [ATAC-seq contributor], [additional authors]. Author order and corresponding author remain to be confirmed.

## Abstract

Sequence-to-function models offer a common interface for predicting molecular phenotypes from DNA, but independent evaluation requires decisions that are not uniquely specified by model outputs. We present BenchAG, a benchmark framework being assembled for AlphaGenome across gene expression, massively parallel reporter assays (MPRA), and ATAC-seq. Our expression workflow makes explicit the relationship between RNA-seq tracks and ground-truth samples, normalization, genomic partitions, gene-level aggregation, strand selection, and model and client versions. We describe a protocol for separating these evaluation choices from model performance and for documenting their consequences through controlled comparisons. The initial software release incorporates an expression-prediction workflow and contribution interfaces for MPRA and ATAC-seq. Quantitative results and cross-task analyses remain to be inserted after verification of the underlying experiments. Our central methodological position is that a reported score should be accompanied by an executable protocol defining the target, population, and prediction statistic being evaluated.

## 1. Introduction

Predicting molecular measurements from DNA can connect regulatory sequence to cellular function. AlphaGenome predicts multiple functional genomic modalities, including RNA-seq and chromatin accessibility [1]. Borzoi provides a related sequence-based approach to RNA-seq coverage prediction [2]. For an independent evaluator, the model prediction is only one component of a benchmark: experimental measurements must also be selected, transformed, matched to output tracks, and assigned to evaluation partitions.

Our work focuses on this evaluation interface. During expression benchmarking, we encountered uncertainty about the appropriate ground-truth representation, differences between quantile-normalized and alternative expression targets, the correspondence of validation and test partitions across models, alternative definitions of a gene-level prediction, strand handling, and software/model versions. These are observations from our benchmarking process; their quantitative effects require controlled analysis. They should not be interpreted as evidence that official evaluations are absent or invalid. A publication and software interface do not necessarily provide a directly executable protocol for our particular target dataset.

BenchAG is intended to make these choices inspectable across three complementary tasks: endogenous gene expression, reporter activity measured by MPRA, and chromatin accessibility measured by ATAC-seq. These tasks interrogate different experimental quantities, so they require individual protocols rather than a single pooled accuracy score. The current contribution is an imported and configurable expression workflow, a shared provenance convention, and a structured plan for completing the multi-task evaluation. We ask: (i) how sensitive are expression conclusions to truth preprocessing and prediction definition, and (ii) which protocol details must be fixed to make independent evaluations reproducible?

**[Before submission: replace provisional contribution wording with the scope actually completed; insert numerical findings only from verified runs.]**

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

The historical code also contains a variant-scoring pathway using a reference-equals-alternate construction and GeneMaskActiveScorer. We retain it as an archival implementation of the second prediction definition explored during development. We do not assume that this statistic is equivalent to the interval-based expression score. Before including its results, we must verify the exact client semantics, selected score field, gene matching, allele/coordinate convention, and historical configuration. The initial public runner exposes the retained interval pathway.

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

## 3. Results and Planned Sensitivity Analyses

No numerical performance results have yet been verified within this repository. The lead author reports completed expression experiments and colleagues report MPRA and ATAC-seq experiments, but the corresponding final tables, manifests, and executable contributions have not yet been consolidated. We therefore reserve quantitative conclusions until their provenance is checked.

**Table 1. Results to be populated from verified runs; “pending” is not a measured outcome.**

| Task | Dataset and eligible population | Prediction/target definition | Main result |
|---|---|---|---|
| Expression | [accessions; genes; tracks; audited holdout] | strand-matched interval mean; [truth recipe] | pending |
| MPRA | [assay; sequences; split] | [activity target and prediction statistic] | pending |
| ATAC-seq | [accessions; regions; split] | [coverage or classification target] | pending |

Expression sensitivity analyses should change one justified choice at a time on a common population: normalization recipe; representative versus aggregated target mapping; interval versus historical alternative score; and model/client version where matching runs exist. The 16 kb versus 1 Mb comparison changes both context and supported scorer width, so it cannot alone isolate a context-length effect. Partition discrepancies should be evaluated by reporting overlap and recomputing on an audited population, not by selecting whichever split yields a favorable score.

**[Insert Table 2: paired metric differences, uncertainty, population sizes, and run identifiers. Insert one figure only after deciding which sensitivity has the clearest supported scientific consequence.]**

## 4. Discussion and Limitations

The main issue is defining what a benchmark measures. A model-to-truth correlation combines predictive ability with decisions about experimental targets, normalization, track identity, gene annotation, and evaluation population. Making these decisions executable can help distinguish a biological limitation from an evaluation mismatch. Our framework is designed to support that distinction; the size of any resulting score changes remains an empirical question.

The initial repository has substantial limitations. It does not bundle experimental data or reproduce ground-truth construction. Supplied ontology mappings require review against exact sample metadata. Gene-list checks do not certify training-region independence. The MPRA and ATAC-seq contributions and quantitative tables are pending. Historical API/model changes may prevent exact reconstruction unless the original configurations and predictions are retained. These limitations motivate a provenance-centered release, but do not establish that other published benchmarks are incorrect.

## 5. Conclusion

BenchAG organizes an independent evaluation of AlphaGenome around explicit target definitions, strand-aware expression scoring, partition provenance, and recorded software/model choices. The initial draft and software provide a foundation for consolidating expression, MPRA, and ATAC-seq results. Claims about robustness or comparative model performance will depend on the verified experiments added to this framework.

## References

[1] Avsec, Ž., Latysheva, N., Cheng, J., et al. 2026. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649:1206–1218. https://doi.org/10.1038/s41586-025-10014-0

[2] Linder, J., Srivastava, D., Yuan, H., et al. 2025. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. *Nature Genetics* 57:949–961. https://doi.org/10.1038/s41588-024-02053-6

[3] AIRI-Institute. GENA_LM, expression benchmark source, commit 6476f509c056ee2252fb86eb0f3e65b0ea76d8fe. https://github.com/AIRI-Institute/GENA_LM/tree/6476f509c056ee2252fb86eb0f3e65b0ea76d8fe/downstream_tasks/expression_prediction/benchmarks/alphagenome_benchmark

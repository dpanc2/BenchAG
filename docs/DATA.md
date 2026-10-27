# Expression data contract

External files are required; no biological ground truth is bundled.

| Input | Format and required fields |
|---|---|
| Forward/reverse intervals for both valid and test | Tab-separated despite `.csv` suffix: gene_id, gene_id_unversioned, chrom, TSS, gene_strand |
| Ground truth for each requested split | CSV: unique gene_id, one numeric expression column per mapped sample ID |
| Track mapping | JSON object: AlphaGenome RNA-seq track name -> one truth sample ID |
| Selected targets | CSV: gene_id (comma-separated allowed), assembly, on_targets, off_targets (comma-separated sample IDs) |
| Optional tissue mapping | CSV: original_id,id; normally null when JSON values already match truth |

Intervals must use the selected genome assembly and nonnegative integer TSS
coordinates consistent with AlphaGenome's interval convention. Forward tables
contain `+`, reverse tables `-`. Record exact annotation release and coordinate
conversion; neither is known from filenames. Gene IDs in truth must agree with
versioned IDs restored by interval tables. Review chromosome naming and boundary
handling against the actual API annotation. The launcher cannot infer these.

The current upstream workflow reads all four interval tables even for a single
split, to restore gene-ID versions. Reference FASTA is not needed for interval
scoring. The archival neutral-variant pathway uses FASTA and is not exposed by
the supported launcher.

## Ground truth decisions still to document

Record source accessions, assay, quantification units, library normalization,
replicate handling, downloaded file checksums, gene filtering, and an executable
preprocessing recipe. Distinguish raw/unnormalized from already normalized units.
For qnorm record the axis, implementation, ties, and fit/reference population;
avoid fitting transforms using held-out outcomes. Existing filenames alone do
not establish these choices. This repository does not yet reproduce truth creation.

Mappings include representative (`first`) and many-to-one (`all`) alternatives.
The supported runner uses the representative mapping. It does not automatically
aggregate truth using the all mapping. Match truth preprocessing to the selected
mapping and report coverage and excluded biosamples.

Gene-list disjointness is weaker than genomic-window disjointness. Audit actual
prediction windows against each model fold's train/validation/test intervals,
including overlapping context; do not equate equally named folds across models.

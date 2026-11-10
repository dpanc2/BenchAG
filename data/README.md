Sample-mapping CSVs listed below are included in Git. Large experimental data
belongs here or at DATA_ROOT and is excluded from Git. See ../docs/DATA.md.

## Included sample mappings

- `file_mappings_borzoi_and_v1_human.csv`: human sample mappings.
- `file_mappings_borzoi_all.csv`: human and mouse sample mappings; select with `--organism`.

These CSVs are inputs to `benchmarks/expression/scripts/build_alphagenome_ontology_map.py`.
Paths inside them refer to the original dataset; the ontology script only uses
accessions and species labels, so those external files are not needed for mapping.

## Shared human reference

Use `data/genomes/hg38.fa` for human benchmarks requiring a local genome.
ATAC and MPRA use this path. On aglab it points to the full GRCh38 reference
in `CAGE_index_build/reference/`. The ATAC runner downloads full UCSC hg38
only if the shared FASTA is missing. Expression and eQTL API predictions do not require a local FASTA.
Mouse benchmarks use their own mm10 assembly.

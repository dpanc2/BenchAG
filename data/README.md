Sample-mapping CSVs listed below are included in Git. Large experimental data
belongs here or at DATA_ROOT and is excluded from Git. See ../docs/DATA.md.

## Included sample mappings

- `file_mappings_borzoi_and_v1_human.csv`: human sample mappings.
- `file_mappings_borzoi_all.csv`: human and mouse sample mappings; select with `--organism`.

These CSVs are inputs to `benchmarks/expression/scripts/build_alphagenome_ontology_map.py`.
Paths inside them refer to the original dataset; the ontology script only uses
accessions and species labels, so those external files are not needed for mapping.

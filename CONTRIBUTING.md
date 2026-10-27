# Adding a benchmark

Work in a branch and submit a pull request once a shared remote is configured.
Each contributor owns one directory under benchmarks and adds an example config
under configs. Keep reusable code separate from exploratory notebooks.

For MPRA and ATAC, supply an executable entry point, a task-specific dependency
file with versions, and one documented command. Add dispatch to run.sh only after
a real implementation is present. Do not return placeholder scores.

Every contribution must include:

1. Input accessions, download/preprocessing steps, schemas and checksums.
2. Genome assembly, coordinate convention, split definition and overlap audit.
3. AlphaGenome client/model version, prediction/scoring definition and settings.
4. A tiny synthetic offline check for alignment, strand/allele orientation,
   aggregation, and metric behavior where relevant.
5. Machine-readable metrics with sample counts and missing/failure counts;
   prediction IDs must align explicitly to truth IDs.
6. A results manifest identifying code revision, input/config hashes, environment,
   and a brief Methods paragraph for the paper.

Keep API keys in environment variables. Do not commit large data, generated
predictions, or credentials. Add paper evidence to paper/evidence.md and edit the
appropriate Methods/Results subsection without inventing cross-task comparability.

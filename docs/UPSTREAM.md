# Source provenance

Imported from https://github.com/AIRI-Institute/GENA_LM, branch `expression`,
commit `6476f509c056ee2252fb86eb0f3e65b0ea76d8fe`, directory
`downstream_tasks/expression_prediction/benchmarks/alphagenome_benchmark`.
Retrieved 2026-10-27. The upstream MIT license is retained in LICENSE.

The imported notebooks, data mappings, legacy runner, configs, and Slurm scripts
are retained. Upstream README and Slurm paths describe the GENA_LM checkout;
use the root README and run.sh here. Legacy scripts are archival, not the public entry point.

Changes to the JSON runner: repository root adjusted; organism, sequence length
and gene mask width configured; annotation-output exclusions are recorded and skipped, while genuine
prediction errors stop scoring after saving a checkpoint; empty checkpoint concatenation guarded; duplicate gene IDs rejected;
loaded prediction paths resolved against output_dir; missing truth fails;
correlations use finite pairs and expose unambiguous axis names.

The new launcher validates inputs and isolates runs by configuration/input hash.
No ground-truth data or experimental performance claims were imported as verified results.

The selected-target scorer now aligns identical finite gene rows before computing
contrasts and returns undefined correlation for constant contrasts. Synthetic
tests cover shuffled IDs and asymmetric missingness. Run identity also includes
launcher/scorer source hashes and the installed AlphaGenome client version.

Annotation exclusions are persisted beside each checkpoint in an outcomes JSON.
Empty annotation responses without gene_id are treated as exclusions. Successful
predictions are preserved, and validation continues after test exclusions.

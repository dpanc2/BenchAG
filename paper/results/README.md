# Expression results from the shared workbook

Source: https://docs.google.com/spreadsheets/d/1qdmZc7KqsdZN7jOuDYwjVRg_kPs0-17DbkmWZ-qpuec/edit?gid=1307246850#gid=1307246850
Retrieved 2026-10-27. The spreadsheet was read without modifying it.

`expression_sheet_snapshot.json` preserves source headers, method labels, exact
unrounded expression values and cell locations. `expression_source_cells.csv`
contains the same selected cells in a flat representation. The snapshot includes
SHA-256 hashes of the downloaded XLSX and independently retrieved linked-tab CSV.
All selected checkpoints cells were cross-checked against the CSV export at its
displayed precision (the CSV uses decimal commas and rounded display values).
The XLSX numeric cells provide the unrounded values retained in these extracts.
Only expression results are included; no MPRA, ATAC, PEP or eQTL result columns
were retained. Values are reported historical summaries, not new benchmark runs.

Paper table 1: checkpoints AC18:AR18 and AC20:AR20.
14-cell panel: checkpoints O18:T18 and O20:T20.
Paper table 2: main A102:J103 and A107:J108 (A/F/H/J only).
Scoring transform discussion: main A110:J112 (A/F/H/J only).
Broader-panel headers: checkpoints AC2:AR3.

The unqualified “AlphaGenome” group is interpreted as qnorm in contrast to the
explicit “gt wo qnorm” group, consistent with the lead author's description.
Exact normalization and historical run settings remain to be verified.
The files are curated source extracts, not an independently verified results release.

# Benchmark visualizations

LOSAT keeps benchmark visualizations separate from ordinary CI and from the v0.1.0 certification authorities.

The workflow `.github/workflows/benchmark-visualizations.yml` always generates fresh outputs from the checked-out LOSAT commit and the official NCBI BLAST+ 2.17.0 Linux archive. It never reads the historical PNG files under `LOSAT/tests/plots/` as input.

## Automatic main snapshot

A relevant push to `main` automatically runs the bounded `representative` visual scope. The representative set is declared in `LOSAT/tests/generate_benchmark_visuals_v010.py` and currently covers:

- BLASTN task-blastn and megablast;
- BLASTP serial pairwise output;
- ordinary TBLASTX, non-default query genetic code, and the approved non-default local-subject `db_gencode` behavior.

The automatic job produces fresh LOSAT and NCBI raw outputs, `hit_distribution_data.tsv`, `visual_case_results.tsv`, `benchmark_metadata.json`, `summary.md`, and `hit_distribution.png`. The entire directory is uploaded as a GitHub Actions artifact.

This workflow is intentionally not part of ordinary pull-request CI. Expensive NCBI comparison runs should not turn every documentation or maintenance PR into a long-running parity campaign.

## Manual full visual snapshot

Use `workflow_dispatch` with:

- `scope=full`
- `mode=visual`

This executes every case currently declared by the three v0.1.0 manifests and generates the same raw tables and hit-distribution plot from fresh outputs.

The manifests remain authoritative for case selection. The visualization script does not maintain a second full case list.

## Manual performance snapshot

Use `workflow_dispatch` with `mode=performance` or `mode=both`. Timing uses:

- the same GitHub runner for NCBI and LOSAT;
- the same release profile and output sink style;
- one warmup by default;
- at least three timed samples by default;
- the median wall time as the plotted statistic.

The artifact contains `performance_data.tsv` with every raw timing sample and `execution_time.png` with medians.

GitHub-hosted timing is **diagnostic only**. Runner hardware and load are not controlled tightly enough to make small timing changes a release claim. A release-facing performance claim still requires the repository performance contract: same machine/toolchain/profile, fixed fixtures, warmup, repeated samples, medians, and output validation.

## Output and provenance

`benchmark_metadata.json` records the candidate SHA, selected scope/cases, NCBI executable versions and SHA-256 values, runner identity, generation timestamp, and performance sampling policy. Raw NCBI and LOSAT outputs are retained below the artifact `raw/` directory.

No line-ending normalization, row sorting, or semantic repair is applied before plotting. Known approved differences remain visible in the raw data and in `visual_case_results.tsv`.

These plots are not `CROSS_PLATFORM_NATIVE_CERTIFIED`, `RC_HANDOFF_READY`, or any other release authority. They are current visual diagnostics tied to an exact commit.

## README publication policy

Do not embed an old tracked PNG merely because the file still exists. A README benchmark figure should be promoted only from a successful fresh workflow artifact whose metadata identifies the exact LOSAT SHA and NCBI version used to generate it.

Generated artifacts are not automatically committed back to `main`. This prevents CI from mutating release history and avoids an automated commit loop. Promoting a figure into the README remains an explicit reviewed change.

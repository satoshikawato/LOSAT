# LOSAT benchmark snapshot

This directory is a frozen, reproducible compute-program benchmark comparing
LOSAT with NCBI BLAST+ 2.17.0+. It characterizes program output and wall-clock
execution; it is not a biological interpretation of the input sequences.

## Snapshot summary

- LOSAT revision: `428c83cbc7ec5eb37e9afe2f1c4782d5e4d68d86`
- Hit distributions: 45 cases and 1,565,630 normalized rows (782,815 rows per
  implementation). Every retained LOSAT/NCBI output pair is byte-identical.
- Execution time: 108 retained samples: 6 cases × 6 modes × 3 repetitions.
- Statistic: median of three wall-clock samples; whiskers show the complete
  three-sample minimum–maximum range.

The hit-distribution and timing datasets intentionally use different NCBI target
protocols:

| Dataset | NCBI target | LOSAT target |
| --- | --- | --- |
| TBLASTX hit distribution | Prebuilt `-db`, with the matching `-db_gencode` | Local `-subject` |
| BLASTN/BLASTP hit distribution | Local `-subject` | Local `-subject` |
| Execution time, every program | Prebuilt `-db` | Local `-subject` |

NCBI database construction is completed before the timed search and is excluded
from the plotted wall time. This lets NCBI's requested thread count take effect;
NCBI reduces local-subject searches to its minimum thread count. The distribution
protocol retains local-subject NCBI output for BLASTN and BLASTP, while TBLASTX
uses a database so non-default subject genetic codes, including code 4, are
applied.

The BLASTP timing target is intentionally timing-only: the NCBI database search
returns a different result set from LOSAT's local-subject search. It is not used
as hit-distribution evidence. This follows the benchmark's compute-timing scope;
the exact BLASTP distribution comparison uses the separate NCBI `-subject` run.

## Timing cases and modes

| Case ID | Program/task |
| --- | --- |
| `PesePMNV.MjPMNV.task_blastn` | BLASTN (`blastn`) |
| `Sakai.MG1655.megablast` | BLASTN (`megablast`) |
| `WSSV.PajaWSV.blastp` | BLASTP |
| `p03_mela_pemojnva` | TBLASTX |
| `d04_ap027131_ap027133_code4` | TBLASTX, query/db genetic code 4 |
| `p11_avclpv_psclpv` | TBLASTX |

Each case has NCBI n1/n8, LOSAT native n1/n8, LOSAT serial Wasm n1, and LOSAT
threaded Wasm requested n8 measurements. The protocol is one untimed warmup plus
exactly three retained timed repetitions for each case and mode. Bars do not
select the fastest sample. Additional repetitions are collected only when
explicitly requested or when the three retained samples are demonstrably
inconclusive.

## Environment

- Intel Core i9-14900HX, 32 logical CPUs
- Ubuntu 24.04.3 LTS under WSL2
- Rust/Cargo 1.92.0
- Node.js 26.8.2, V8 14.6.202.34-node.28
- NCBI BLAST+ 2.17.0+

Process startup, Node startup, and Wasm compilation are included in each search
measurement. Output is written to regular files. Thread labels state the
requested configuration; they do not claim that every worker is continuously
busy for every input.

## Files

- `metadata.json`: source revision, exact commands, target semantics, environment,
  checksums, and collection protocol.
- `alignment_results.tsv.gz`: normalized rows used by the distribution plot.
- `execution_times.tsv`: all 108 retained wall-clock samples.
- `plot_data.json`: deterministic histogram data and timing summaries.
- `hit_distribution.png`: distributions using the target policy above.
- `execution_time.png`: median timing bars and full three-sample ranges.
- `render_manifest.json`: input/output checksums from deterministic rendering.

## Reproduction

Build a snapshot from completed comparison runs:

```bash
python scripts/build_benchmark_snapshot.py \
  --distribution-run /path/to/distribution-run \
  --timing-root /path/to/timing-root \
  --snapshot benchmarks/v0.1.0
```

Re-render the figures without running a search:

```bash
python scripts/render_benchmark_plots.py \
  --snapshot benchmarks/v0.1.0 \
  --output benchmarks/v0.1.0
```

The renderer reads only the frozen local snapshot and verifies dataset checksums
before drawing the figures.

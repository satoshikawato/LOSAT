# LOSAT v0.1.0 benchmark snapshot

This directory is a durable, data-only snapshot. Its alignment dataset was
recovered from retained LOSAT evidence, and its current timing dataset was
collected once under the protocol below. Rendering it does not build or execute
LOSAT, execute an NCBI BLAST+ program, collect a benchmark, or access the
network.

## Data completeness

| Figure data | Status | What is available | Final-v0.1.0 interpretation |
| --- | --- | --- | --- |
| Alignment output | `AVAILABLE_EXACT` | 1,371,516 normalized rows from all 43 retained PR 5 native contracts; exact source-output and command hashes are recorded in `metadata.json` | Representative of the certified Linux x86_64 runtime. LOSAT SHA `5845d22ed9842449628a647f29b8c6762511ca59`; NCBI BLAST+ `2.17.0+` |
| Historical wall time | `AVAILABLE_HISTORICAL` | 219 retained log values that reproduce the committed historical execution-time plot | Not attributable to the final v0.1.0 SHA; producing LOSAT SHA and complete NCBI identity are unknown |
| PR 93 wall time | `PARTIAL` | Four reported scalar values for base `3ba0024e88444adc2554eebae26aae7703e38ae0` and candidate `d92905569317ff9c6bde4e9dabcc66f6e4c15f81` | Candidate algorithm is in the certified lineage, but raw samples, statistic definition, thread count, and commands are unavailable |
| Final-v0.1.0 timing samples | `AVAILABLE_CURRENT` | 60 retained wall-clock samples: six representative cases, two tools, five timed repetitions per tool/case | Collected directly from final main SHA `1ef13ed09e9cf6e3d6c3b1a1ed9527a389c3992d` on one machine/session; medians are used in the current figure |

The TBLASTX alignment plot includes six certified, approved local-subject
non-default `db-gencode` deviation contracts. The BLASTN inputs include the one
unchanged source-undetermined accepted contract described by the release
certification. The figures visualize recorded outputs; they are not themselves
a parity score or a new certification.

## Current timing protocol

The current series uses NCBI BLAST+ 2.17.0 and the final-main LOSAT binary on
the same WSL2 Ubuntu 24.04.3 session with an Intel Core i9-14900HX. Both tools
use one thread and regular-file output on the same local filesystem. Each
tool/case group has one warmup followed by five retained timed repetitions;
points show the median and whiskers show the full five-sample range.

The six cases are `PesePMNV.MjPMNV.task_blastn`,
`Sakai.MG1655.megablast`, `pairwise_default_serial`,
`p03_mela_pemojnva`, `d06_ap027131_ap027133_db4`, and
`p11_avclpv_psclpv`. The p11 case is the heavier TBLASTX representative because
it is exact-certified and has both retained historical timing context and the
additional PR 93 scalar lineage. Sakai retains its
`SOURCE_UNDETERMINED_ACCEPTED` classification, and d06 retains the approved
local-subject non-default `db-gencode` deviation.

## Files

- `metadata.json`: snapshot identity, checksums, completeness classifications,
  per-contract source paths, exact command records, versions, SHAs, dates, and
  known provenance gaps.
- `alignment_results.tsv.gz`: normalized output rows. Columns are `program`,
  `case_id`, `implementation`, `classification`,
  `primary_for_distribution`, `source_row`, `qseqid`, `sseqid`, `pident`, and
  `length`. Duplicate BLASTP thread-4 contracts remain recoverable but are not
  counted twice in the distributions.
- `execution_times.tsv`: the 60 individual current samples and output hashes,
  plus the separately identified 219 recovered historical timing values and
  four PR 93 scalar-lineage values. Historical values are not used as current.
- `plot_data.json`: deterministic numerical plot product used by CI. It records
  the renderer-source hash, bins, weighted histograms, sampled-scatter hashes,
  all current timing samples, and current median/mean/min/max summaries.
- `hit_distribution.png` and `execution_time.png`: rendered figures.
- `render_manifest.json`: input and output checksums for the committed render.

## Render

With Python, NumPy 2.2.5, and Matplotlib 3.10.3 installed:

```bash
python scripts/render_benchmark_plots.py \
  --snapshot benchmarks/v0.1.0 \
  --output benchmarks/v0.1.0
```

The renderer accepts only the snapshot and output paths. It verifies the two
source-data checksums before reading them. It contains no process-launching,
binary-discovery, build, download, or network path.

## CI and collection boundary

The `Benchmark plot rendering` workflow installs only the plotting dependency,
runs renderer tests, renders into runner-temporary storage, and compares the
resulting `plot_data.json` with the committed product. PNG byte hashes remain in
the render manifest for local traceability, but CI avoids treating platform font
rasterization as benchmark-data drift.

Benchmark collection is intentionally absent from the renderer and its CI job.
The final-v0.1.0 samples were collected once by a separate, temporary manual
collector; no collector was added to the repository. Any future collection must
remain separate and must replace snapshot data only through an explicit,
provenance-reviewed change.

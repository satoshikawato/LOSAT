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
| Prior final-v0.1.0 timing samples | `SUPERSEDED_HISTORICAL` | 60 retained wall-clock samples: six representative cases, two one-thread tools, five timed repetitions per tool/case | Preserved unchanged as historical evidence and excluded from the current figure |
| Full final-v0.1.0 timing samples | `AVAILABLE_CURRENT` | 180 retained wall-clock samples: six representative cases, six modes, five timed repetitions per mode/case | Collected from exact SHA `af3e2ea837afdb8a00cf19920f68be4f0bf3bfb5`; medians and min–max ranges use all five samples |

The TBLASTX alignment plot includes six certified, approved local-subject
non-default `db-gencode` deviation contracts. The BLASTN inputs include the one
unchanged source-undetermined accepted contract described by the release
certification. The figures visualize recorded outputs; they are not themselves
a parity score or a new certification.

## Current timing protocol

The current series uses NCBI BLAST+ 2.17.0 and LOSAT binaries built from exact
SHA `af3e2ea837afdb8a00cf19920f68be4f0bf3bfb5` on the same physical WSL2
Ubuntu 24.04.3 machine with an Intel Core i9-14900HX. The six modes are NCBI
BLAST+ n1/n8, LOSAT native n1/n8, LOSAT serial `wasm32-wasip1`, and LOSAT
threaded `wasm32-wasip1-threads` requested n8. All use regular-file output on
the persistent evidence filesystem. Each case/mode group has one protocol
warmup followed by five retained timed repetitions; points show the median and
whiskers show the full five-sample range.

This is a same-machine controlled benchmark collected in two execution
segments, with binary/toolchain identities revalidated after restart. Segment 1
uses boot ID `d95508ef-366f-4b1c-862c-b613b0caa163` and contains 207 completed
protocol invocations; segment 2 uses boot ID
`942079cc-0a9f-43b4-9b08-7c377d67a6d6` and contains only the nine previously
missing timed samples. Before the first resumed sample of each affected mode,
one additional restart-conditioning `resume_warmup` was run. These six runs are
fully recorded in the external evidence and excluded from the 36 protocol
warmups and all timed statistics.

Only the six p11 case/mode groups span both segments. Independent review found
acceptable continuity: the segment-median shifts were evenly split between
faster and slower, their median shift was -0.225%, and the combined medians
differed from segment-1 medians by -1.720% to +1.270%. No p11 replacement was
required, no samples were discarded or adjusted, and the full 216-invocation
campaign was not restarted. Threaded-Wasm n8 is the requested configuration;
the per-case effective-thread classification and probe evidence remain in the
snapshot data. This performance characterization does not expand the certified
feature or target scope.

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
- `execution_times.tsv`: the 180 individual current samples with output hashes,
  boot-segment provenance, and effective-thread evidence; the 60 superseded
  one-thread samples, 219 recovered historical timing values, and four PR 93
  scalar-lineage values remain separately identified and are not used as
  current data.
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
The final-v0.1.0 samples were collected by a separate manual collector and a
narrow external resume wrapper; neither collector was added to the repository.
Any future collection must remain separate and must replace snapshot data only
through an explicit, provenance-reviewed change.

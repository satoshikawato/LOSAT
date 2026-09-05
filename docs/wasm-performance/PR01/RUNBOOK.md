# PR01 baseline harness

The 2026-09-05 follow-up uses separate
[600-second](followup-20260905/manifest.json) and
[bounded concurrent](followup-20260905/manifest-concurrent.json) manifests, with the original
120-second manifests and raw evidence retained. See the follow-up
[dependency boundaries](followup-20260905/REQUIREMENTS.md).

This package records local baseline work only. It does not certify a release or
authorize later PRs. The production Rust sources and existing canonical
authority are unchanged by PR01.

The current frozen measurement contract is [contract/manifest.json](contract/manifest.json).
`frozen/` belongs to the initial 43/41 audit. `baseline/` identifies the earlier warm smoke validation, not an alternative
expected-output authority. Unused development snapshots were removed. The TBLASTX no-hit trial on the original LC pair produced
hits even at `1e-180`; the final additional no-hit contract uses the compact
280-byte nucleotide fixture. Neither trial updated canonical output.

Run commands from the repository root. Output directories must be new; the
harness rejects existing directories and conflicting controlled input files.
It never chooses the newest artifact or rewrites a frozen expected hash.

```bash
python3 LOSAT/tests/wasm_performance.py freeze \
  --native /tmp/losat-pr01-20260905/native/release/LOSAT \
  --serial /tmp/losat-pr01-20260905/serial/wasm32-wasip1/release/LOSAT.wasm \
  --threaded /tmp/losat-pr01-20260905/threaded/wasm32-wasip1-threads/release/LOSAT.wasm \
  --oracle-dir /home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin \
  --source-archive /home/kawato/tools/ncbi-blast-source/ncbi-blast-2.17.0+-src.tar.gz \
  --timeout 120 --output /tmp/pr01-new-contract

python3 LOSAT/tests/wasm_performance.py run \
  --manifest /tmp/pr01-new-contract/manifest.json \
  --phase audit --output /tmp/pr01-new-audit

python3 LOSAT/tests/wasm_performance.py run \
  --manifest /tmp/pr01-new-contract/manifest.json \
  --phase benchmark --output /tmp/pr01-new-cold

python3 LOSAT/tests/wasm_performance.py run \
  --manifest /tmp/pr01-new-contract/manifest.json \
  --phase warm --output /tmp/pr01-new-warm

python3 LOSAT/tests/wasm_performance.py run \
  --manifest /tmp/pr01-new-contract/manifest.json \
  --phase profile --targets native --threads 1 \
  --case p01_ap027280_self --output /tmp/pr01-new-profile

python3 -m unittest discover -s LOSAT/tests -p test_wasm_performance.py -v
node --check LOSAT/tests/run_wasm_performance_warm.js
```

`--targets`, `--threads`, and repeated `--case` select explicitly reported
subsets. Their success is not full-matrix success. Audit preserves original
thread arguments, including the two native BLASTP thread-4 contracts. Serial
Wasm has 41 directly applicable rows. Benchmark changes only the requested
thread count and output sink. Empty manifest fields continue to omit options.

All commands use shell-free ordered argv and the PR5 lexical input root
`/tmp/losat-pr5-runtime-cert-5845d22/LOSAT`. Inputs are copied from the current
worktree with exact bytes and hashes, without overwriting conflicting files.
The canonical hash and command selectors are independently reconstructed from
the existing certification catalog. Official NCBI executable hashes are fixed
to the existing 2.17.0 Linux authority. Platform-native Gate B remains NOT_RUN
on this Linux host; its Windows/macOS fingerprints never become LOSAT expected
output.

Cold samples use one untimed warmup, five fresh processes, then a separate
diagnostic run. Wall time includes process startup, Node compilation,
instantiation, search, output, and process exit. An independent watchdog bounds
the blocking process wait. `/usr/bin/time` supplies CPU seconds and peak RSS;
missing measurements remain null. Output hashes are checked on every repeat.
Failed oracle gates prevent benchmark execution. A failed repeat stops that
series and invalidates its median. Timeouts are not byte mismatches or zero-time
samples.

Warm serial samples compile one Module and execute six distinct instances,
using Node WASI `returnOnExit=true`. Each instance gets fresh linear memory and
WASI state. The first run is warmup. The measured boundary includes creating
WASI, instantiating, and returning from the full command; compile time is
separate. This is not a resident search engine and does not call `_start` twice
on one instance. Raw hash failure prevents the next job. RSS is labeled either
post-job RSS or cumulative process-lifetime peak RSS, not a per-job peak.
Threaded warm-instance lifecycle and direct API remain NOT_RUN.

Existing `LOSAT_TIMING` output is collected only in diagnostic runs. Stage times
can overlap and can sum work across threads; they must not be subtracted from
wall time to infer a serial fraction. BLASTP disables these timers on Wasm.
Threaded diagnostics record configured effective computation threads per
stage and actual helper spawns. Active worker utilization and native effective
N remain unknown where no measurement exists. Requested N is not proof of
effective N.

The frozen 4 GiB RSS budget is a rejection threshold. The threaded runner
enforces a 1 GiB linear-memory maximum; the serial runner has no corresponding
host limit. Its limit is module-defined. No browser-memory budget or browser
performance is certified.

Each execution retains raw output, stdout, stderr, ordered command JSON, and
resource usage. `samples.jsonl` and `summary.json` accompany each run. On a
byte difference, `ordered.diff` is decisive; `sorted-diagnostic.diff` is only a
diagnostic. Build logs, source identities, input identities, oracle archive
checks, and the final evidence record are adjacent to this runbook.

To revalidate the retained evidence and regenerate the review tables on this
same checkout (explicit artifact and oracle paths must still exist):

```bash
python3 LOSAT/tests/summarize_wasm_performance.py docs/wasm-performance/PR01
```

The exporter rehashes retained outputs against their recorded byte lengths and
frozen canonical SHA. It checks successful process measurements against saved
`result.json`, requires every warmup/timed/diagnostic row for a cold series,
recomputes cold medians/ranges/RSS, and rejects inconsistent saved summaries.
`RUN_MANIFEST.jsonl` is filled from the supplied per-run template; its thread
observations reference the separate diagnostic output. Its whole-series status
cannot be replaced by the first successful warmup. Missing metadata remains
null with an explicit reason. Warm results retain their distinct boundary.

The final TBLASTX V8 CPU profile is a separate diagnostic command after all
baseline searches finish. `cpu-profile/search/command.json` records the exact
Node invocation (`--cpu-prof`, 1 ms requested sample interval), and
`cpu-profile/gate.json` records the frozen raw-output gate. This profiling run
is excluded from all latency medians. The raw `.cpuprofile` is inspectable in
DevTools. Samples may include runtime startup/compilation and profiling can
perturb execution; percentages are sampled CPU attribution, not benchmark
wall time or an additive decomposition of existing engine timers.

## Follow-up execution and reconciliation

`followup-20260905/retry_audit.py` imports the unchanged test harness. It runs
only the incomplete initial-audit oracle/process executions and references
successful results after SHA/byte-count verification. In particular, p10
native is reused despite the old row's ORACLE_FAIL, and d04's completed oracle
is reused. Each missing process gets one 600-second attempt; retry timings
are diagnostics only and never enter cold/warm medians. The driver rejects
an existing `retry-audit` destination. Re-execution needs a new evidence
directory and separately declared manifest, not removal of old outputs.

The sequential driver completed the p11 oracle timeout and was interrupted
during p11 native to switch audit scheduling. The owned processes and reason
are recorded in `followup-20260905/sequential-interruption.json`. Partial
files are retained and are neither TIMEOUT nor successful evidence.
`retry_concurrent.py` uses `manifest-concurrent.json`: 21 remaining processes,
each with 600 seconds, at most three searches concurrently. It reuses the
completed p11 oracle TIMEOUT without another attempt, the successful p10
native result and d04 oracle. Contended audit durations do not support any
performance claim or timeout-speed comparison.

`followup-20260905/reconcile.py` produces `combined-audit.json` after the
`concurrent-audit/summary.json` exists. It rechecks manifest identity, exact selector sets, saved
ordered argv/cwd, raw SHA/byte counts, and the existing exception validators.
Every row names its oracle and LOSAT process evidence. The original 84 rows
and their timeout statuses remain unchanged on disk.

The old exporter was run against `followup-20260905/retained-revalidation`,
whose read-only input directories link to the original runs. Generated
tables go only into that new directory. Thus this retained-data check does
not regenerate the original top-level tables or convert a timeout to PASS.
All baseline cold and warm measurements are reused; an oracle timeout still
prevents a new p11 timing series. Warm's watchdog is a whole-series limit of
the per-search allowance times the number of jobs; it is not independently
enforced on each warm job.

# Wasm threading remediation evidence, 2026-09-13

See [the implementation report](../../wasm_threading_remediation_results_20260913.md)
for scope, results, artifact identities, limitations, and reproduction commands.

- `cold-performance.tsv`: 63 recorded cold-process conditions with correct raw output;
  median/min/max seconds, peak RSS bytes, and raw output SHA-256.
- `reuse-performance.tsv`: all 35 candidate reuse conditions, each with five
  timed calls; includes invocation, worker startup, and host exit wait intervals.
- `index.json` / `SHA256SUMS`: final archive identity and scope.
- `evidence.tar.gz`: raw outputs, inputs, ordered argv, environments, stderr,
  result hashes, host events, tool/build logs, before/after source records,
  independent review records, CRT evidence, and the local validation bundle.

Verify the archive before extracting it into a new directory:

```sh
sha256sum --check SHA256SUMS
mkdir extracted
tar -xzf evidence.tar.gz -C extracted
```

`SHA256-MANIFEST.json` inside the archive binds every other member by size and
SHA-256. The final roots are `matrix-repaired-formats/`, `frozen-deferred.json`, `frozen-final/`, `frozen-timeout-retries/`, `frozen-continuation/`,
`lc-thresholds/`, `performance-final/`, `format-artifacts/`,
`format-artifacts-reverse/`, and `validation-package-smoke-final/`.
`frozen-inputs/` and `frozen-authority/` retain the original Git blob bytes.
The inner validation bundle has its own SHA256SUMS and usage README.

Of the 63 cold conditions, 35 candidate conditions pass the current output and
execution contracts. The other 28 are output-correct baseline measurements;
the baseline can silently shrink requested pools. Those measurements describe
prior behavior and repair costs, not acceptance of the current thread contract
or a valid like-for-like parallel speedup.

`frozen-final/` preserves the original 98 successful cases, five deadline
failures, and three interrupted searches. One serial p11 retry passed.
At the user's request, the seven remaining active retry/continuation searches
were stopped and further long regressions deferred. `frozen-deferred.json`
verifies 99/145 distinct successful cases and retains all 114 attempts
(99 successes, five timeouts, ten interruptions). The full gate is INCOMPLETE.
Its 46 pending records include exact argv, environment and expected hashes for
later execution; none of the six retained-Linux oracle or twelve repeatability
cases has completed in this frozen run. Original failures are never relabeled.
`frozen-deferred-stop.json` records the user-directed stop. The full acceptance
script remains strict and has not produced a `frozen-accepted.json`.

Historical diagnostic roots and failed intermediate logs are retained for
traceability; they are not additional passing tests. `matrix-all-formats/`
records the formatting discrepancies before repair. `frozen-regression/`
contains metadata for the superseded, incomplete run. Its p09 artifact binding
was ambiguous and is excluded from acceptance. `old-audit/` preserves the prior
402-file audit by its original hashes. `cost-profile/` is a separately
instrumented build predating the final BLASTP formatting changes. Its counters
are bounded diagnostics, not final production speed or complete heap profiles.

Absolute paths in raw records are the paths actually used. Reproduction on a
different host requires the repository and exact registered input blobs at the
recorded lexical fixture path where output formats include those paths. Do not
rewrite headers or normalize outputs to make a comparison pass.

This is local validation of the recorded dirty source tree. It is not a release
candidate or a substitute for the existing hosted platform Gate B and a new
release certification lineage.

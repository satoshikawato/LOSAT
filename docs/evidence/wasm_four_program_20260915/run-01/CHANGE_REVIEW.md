# Change review boundaries

The initial dirty working tree is frozen in `baseline-inputs.tar.gz` and
`initial-user.diff`. Existing user changes are not the implementation diff for
this task. Measurement has ended at the user's instruction. The integrated
implementation is adopted by the user's subsequent explicit instruction;
see [ADOPTION.md](ADOPTION.md). Six final time-guard failures remain recorded.

## Production Rust

| Owner | Task changes |
|---|---|
| `algorithm/blastn/alignment/greedy.rs` | M0 preserves nonaffine scratch writes on a fence return; M1 addresses existing base/pool rows without round-trip copies. The affine path retains its behavior. |
| `algorithm/blastn/blast_engine/run.rs` | N0 uses the existing global query-context search space for subject-set cutoffs. N1 carries prepared offsets with speculative DP results. N2 reuses batch storage. N3 supplies context-correct exact four-base lookup slices. |
| `algorithm/blastn/extension.rs` | N3 precomputes the existing BLASTNA expression and uses it in the approximate word>=11 extension path. Exact short-word extension retains its behavior. |
| `algorithm/blastp/blast_engine.rs` | P1 shares lazily initialized immutable initial-matrix values. P2 shares single-query preparation and reuses privately owned mutable scratch. Ordered result/error consumption is preserved. |
| `algorithm/tblastx/sum_stats_linking/linking.rs` | X2 defers coordinate reads after the sum predicate for non-traced rejected helpers, preserving helper visits, jumps, comparisons, predecessor selection and floating-point evaluation order. |

All paths are under `LOSAT/src/`. Each modification carries its NCBI source
owner/snippet in the code. M0 and N0 are correctness corrections and are included
in C0; performance comparisons do not attribute their effect to an optimization.
Individual candidate patches and the independent review are in this directory.

## Tests and measurement harnesses

- Embedded Rust regressions cover fence/reused scratch, speculative batch16 and
  shifted-subject boundaries, the exact four-base expression, and subject-set
  cutoffs. The new word7 golden file is
  `LOSAT/tests/unit/helpers/ncbi_reference_data/blastn_word7_subject_set.out`.
- `LOSAT/tests/benchmark_wasi_reuse.js` corrects its child count to `N−1`, matching
  the existing caller-included runtime contract. A same-artifact before/after
  reproduction retains the original assertion failure and corrected PASS.
- The run's Python/Rust diagnostic harnesses, body timer builds and browser
  staging tools are evidence tools. They do not add runtime counters, host
  flags, external BLAST execution or fallback behavior to production LOSAT.
- Final Cargo test, clippy and format logs apply to the current integrated
  source. Existing formal release fingerprints/expected outputs are unchanged.

## Documentation

The requested plan's status links to this run. `REPORT.md`, measurement
declarations, `INDEPENDENT_REVIEW.md` and `REPRODUCE.md` distinguish candidate
decisions, completed checks, failed probes and unmeasured scope.

## Generated evidence

Build bindings, exact process argv, raw outputs, state transcripts, timings and
browser records are generated artifacts. Diagnostic counts are not adoption
timings. Native requested allocation bytes are not Wasm linear memory or RSS.
The evidence archive is produced after all measurement jobs are stopped.

## Commit representation

`commit-source-bindings.json` binds the committed source to the measured working
files. The commit retains HEAD's LF line endings for `linking.rs`; the measured
CRLF working file and archive remain intact. This removes pre-existing line-end
churn from the commit. The other four owned production files match byte for
byte. No source logic changes or new benchmark runs accompany this conversion.

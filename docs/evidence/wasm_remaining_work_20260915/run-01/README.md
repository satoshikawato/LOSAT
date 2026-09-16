# Wasm remaining-work evaluation, 2026-09-15

This is a new evaluation against commit `8f23f774` **plus the working-tree inputs
frozen at the start of this run**. It does not reopen the previous adoption.

The [implementation plan](../../../wasm_remaining_work_implementation_plan_20260915.md)
is the task authority. See [STATUS.md](STATUS.md) for the final decision and
[DIAGNOSIS.md](DIAGNOSIS.md) for eligibility decisions and diagnostic limits.

## Reproduction and evidence layout

The experimental working directory is `/tmp/losat-remaining-20260915`. All build,
run and measurement commands are recorded with input and artifact hashes. Native,
serial command and threaded command have separate Cargo target directories.

- `baseline-manifest.json` and `baseline-inputs.tar.gz` freeze the initial Rust,
  host, configuration, lockfile and fixture inputs.
- `supplemental-test-manifest.json` and `supplemental-tests.tar.gz` add the nested
  Rust tests and already-existing boundary fixtures needed for focused checks.
  `supplemental-integration-fixtures.json` and its archive add the unchanged
  megablast integration expected file omitted from the initial isolated snapshot.
- `archives.json` lists all archive sizes/hashes. `bundle-manifest.json` lists
  every member hash in the four result/source/artifact bundles;
  `archive-verification.json` records their successful integrity checks.
  Run `python3 verify_bundle.py` in this directory to repeat those checks.
- Experimental sources are stored as overlays relative to the frozen baseline.
  `H1-source/LOSAT/tests` additionally contains all host dependencies. Build
  manifests enumerate the exact source inputs used for each artifact.
- Source preparation and comparison scripts are retained, including unsuccessful
  preparation/trace attempts. Their later corrections are recorded explicitly;
  failed attempts are never counted as successful comparisons.

For replay, use a fresh copy of `/tmp/losat-remaining-20260915`:

1. Extract `baseline-inputs.tar.gz` and `supplemental-integration-fixtures.tar.gz`
   at the experiment root. Both contain a `baseline-source/` prefix.
2. In `supplemental-tests.tar.gz`, extract `LOSAT/...` members under
   `baseline-source/`, and `boundary-fixtures/...` members at the experiment root.
3. For each named experimental `*-source` overlay, first copy the complete
   `baseline-source` tree to that name, then extract `sources-and-diagnostics.tar.gz`
   at the experiment root. Its changed files replace only the corresponding
   baseline copy. H1 additionally carries all host JavaScript dependencies.
4. Extract the verification, measurement and experimental-artifact archives at
   the experiment root. The bundle manifest preserves every original file path.

Commands contain absolute fixture, oracle, checkout and output paths. Adapt
machine-specific checkout/oracle locations when replaying elsewhere; do not
pretend relocated runs are the original observations. Use a fresh output
location: comparison scripts intentionally refuse to overwrite results.
Install the recorded Rust toolchain/targets, Node, Cargo dependencies and NCBI
oracle binaries separately. NCBI binaries are comparison tools only.

## Adoption criterion clarification

The Wasm/Native ≤1.20 ratio is an aspirational target, not a prerequisite for
adopting a useful speedup. The user explicitly confirmed this during execution.
H1 is evaluated on its measured improvement and correctness/regression guards;
X1 is rejected because it is slower than its own baseline.

## Interpretation limits

Timing runs are exclusive. Counter, code-print and correctness runs are separate
from performance samples. The main condition is ordinary Node without compiler
flags; TurboFan-fixed results are supplementary diagnostics. All warmups and
three alternating A/B pairs are retained. A failed or inconclusive candidate is
not extended with an automatic search for better samples.

These comparisons use fresh NCBI BLAST+ 2.17.0 output for the listed valid-query
fixtures. They are not PR 6 native certification: frozen PR 5 Gate A and the
registered official-platform Gate B remain separate contracts. The long code4
fixture uses a fresh NCBI database oracle for the approved local-subject genetic
code behavior. No expected file or platform fingerprint is changed.

Node/WASI results do not certify browser, Wasmtime, or wasm32-unknown-unknown
performance. No browser package is updated by this evidence directory. Bounded
reactor checks do not establish unlimited-repeat memory stability. The existing
invalid-TBLASTX-query status and reactor-memory issues remain outside this change.

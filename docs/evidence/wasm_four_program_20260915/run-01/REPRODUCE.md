# Reproduction and evidence boundaries

Status: adopted by user; measurements remain stopped with exactly three complete pairs. See [ADOPTION.md](ADOPTION.md).

## Fixed identities

- Original working-tree baseline: `baseline-inputs.tar.gz` and the file hashes in
  `baseline-inputs.json`. HEAD alone is insufficient because the starting tree
  contained user changes.
- C0: original baseline plus the separately tested scratch-fence and global
  BLASTN subject-set-cutoff corrections. See `C0.patch`.
- Integrated candidate: C0 plus N1/N2/N3/M1/P1/P2/X2. The candidate remains
  in the working tree and is adopted by user instruction, accepting the six
  final time-guard failures as documented performance tradeoffs.
- Final production and coarse-body artifacts: `final-build-bindings.json`.
  Body artifacts have two outer timer reads and are diagnostic only.
- Runtime: Node 26.8.2, native/serial WASI/threaded WASI. Command-Wasm threaded
  builds use `wasm32-wasip1-threads` with the `wasm-threads` feature. Plain
  `wasm32-wasip1` remains serial. `wasm32-unknown-unknown` is not measured.
- NCBI source and executable identities are recorded in `REPORT.md`, build and
  gate metadata. NCBI executables are comparison oracles only.

## Scripts and original layout

Exact executed argv and working directories are retained per process. Scripts
use the original repository path `/mnt/c/Users/genom/GitHub/LOSAT` and evidence
workspace `/tmp/losat-four-program-20260915`. Relocation requires updating those
paths while preserving source, fixture, runner and artifact hashes. Use fresh
output directories; the harness rejects existing evidence/output paths.

`build.py VERSION` runs locked offline release builds in the corresponding
`VERSION-source/LOSAT` directory with independent target directories. The final version records
source and artifact hashes. Initial individual builds recorded artifact hashes only. Run from the snapshot crate so its Cargo target
configuration applies. Required Rust targets and already cached Cargo
dependencies must be available.

`gate.py`, `long_gate.py`, the edge drivers and `final_validation.py` reproduce
the raw and runtime checks. `review_validation.py` verifies explicit record
counts/statuses and preserves the known failed probes; a driver exit code alone
is not proof of a complete matrix. `validation-review.json` is the completed
review for this run.

`final_measurements.py` declares the cold/body/reuse commands. Do not run other
benchmarks, compilers, browser checks or archives concurrently with adoption
samples. Affinity is CPUs0–7. The primary runtime uses default Node flags and
the same runner paths for both versions. Each cold/body condition has one
warmup plus five alternating A/B pairs. No sample is selected by speed.

`reuse_benchmark.py` collects two independent AB/BA sessions per case, each with
one warmup and five measured invocations per condition. Reactor n1/n8 share a
per-case instance. An invocation includes API input/result copies and host
worker completion; source FASTA loading and module compilation are separately
reported setup costs. Failed fixed sessions are retained, and other declared
independent processes continue. The timeout remains3,600 seconds per process.

`final_turbofan_sequence.py` separately declares the TBLASTX cold condition with
equal `--no-liftoff --no-wasm-tier-up` flags. It starts after the default sequence.
This secondary condition does not change the primary runtime or imply measured
TurboFan body/reuse/browser behavior.

`summarize_final.py` recomputes medians, all samples, regression guards and
same-version Native/Wasm ratios. Body values come from one finite positive
`[BODY_SCOPE_SECONDS]` record in each eligible stderr, not the process-time
summary. Reuse reports session-wide memory and RSS, complete-window status,
individual sessions and the combined two-session medians separately.

## Retained failures and limits

See `REPORT.md` for the original code4 oracle timeout, invalid scratch comparison,
word7 baseline difference, old child-count assertion, unsupported `-strand`
probe and statistically invalid TBLASTX query exit-status defect. Replacement
evidence has distinct directories; original failures are not overwritten or
relabelled. Existing reactor-memory and TBLASTX-reuse failures retain their
original criteria and remain open unless independently resolved.

Browser checks serve the actual gbdraw source locally with staged artifact
overlays. They verify the consumer lifecycle and offline boundary; they do not
certify a rebuilt/distributed gbdraw package or browser speed.

This run does not modify frozen PR5 outputs, PR6 platform fingerprints or formal
release certification.

## User stop and durable archive

`measurement-stop.json` supersedes the planned five-repeat count for the final
cold summary. The five-repeat/body/reuse/TurboFan commands above are the original
protocol, not instructions to resume the canceled queue. Run only
`python3 summarize_three_pairs.py /tmp/losat-four-program-20260915` to regenerate
the three-pair summary from preserved records; this launches no benchmarks.

`execution-evidence.tar.gz` stores the execution workspace without Cargo object
caches. `execution-evidence-manifest.json` records every retained path, SHA256,
size and mode, including duplicate-content hardlinks. Original RUNNING checkpoint
files remain unmodified evidence of interruption; `measurement-stop.json` is
the final status. C0.patch is derived from the frozen baseline archive; final
source snapshots and `final-build-bindings.json` are the build identities.
Historical initial source snapshots are retained as observed and are not claimed
to be contemporaneous build-time bindings.


Archive verification initially compared full `stat` mode (including file type)
with tar's permission-only mode and stopped. The original archive was retained
unchanged; the completed verification compares permission bits (`07777`) plus
every content hash and size. `archive-verification-initial.log` preserves the
failed verifier run; the manifest records the correction.

## Archive transport in Git

The execution archive is committed as four ordered `.part01`–`.part04` files.
`archive-parts.json` records each part's SHA256 and the unchanged full archive
SHA256. Reconstruct from this directory before reading or extracting:

```bash
cat execution-evidence.tar.gz.part0[1-4] > execution-evidence.tar.gz
sha256sum execution-evidence.tar.gz
```

Expected archive SHA256: `5fe6dfd0b9ef33beded1be5b95978641a83c0de2c1889a90404d969793022c53`.
The reconstructed archive is locally ignored; the parts and full content
manifest are tracked. Splitting changes no archive contents or evidence.

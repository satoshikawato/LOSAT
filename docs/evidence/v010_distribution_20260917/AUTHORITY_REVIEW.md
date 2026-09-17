# v0.1.0 release authority conflict

Status: **authority revision approved; recertification pending**. The maintainer
selected “認証仕様を改訂して再認証する” in this implementation session. Runtime
is unchanged. Gate A and the current BLASTN classification are revised as below;
publication remains blocked until all fresh certification gates pass.

At `89a85d6664a44c2dd41fe738158222c88c15b1e8`, a fresh Rust 1.92.0 locked
release build reproduces CI run 35089176534. With the original fixture bytes,
lexical paths, outfmt 7 and one thread, current LOSAT and the official Linux
NCBI BLAST+ 2.17.0 oracle both emit SHA-256
`ac8177764f35b2fb92777f8b1d01ea646cc55ae21b05f1f486ee231348a34a45`.
The oracle executable is
`33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5`.
The retained PR5 official oracle has this same raw output hash.

Frozen PR5 LOSAT emits
`3f27c1f1396b59ecb7e78827b5da390e9c97a276c05cbc7a6c82f56bd0a03262`.
Both files have 536,922 bytes and 6,476 data rows. Five rows differ in identity,
length, mismatch and/or gap-open fields; headers, coordinates, E-values,
bit scores and ordering agree. See `frozen-current.diff`; `runs.json` preserves
both exact commands and executable hashes. `fixture-staging.json` binds inputs.

The source-backed correction in `f5955c5952998e50c1262186212d8cfc1527eb1d`
sets traditional megablast's initial X-drop to zero, then uses the subject word
cutoff. NCBI source owners:

- `c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54`: zero-initialized
  initial-word options.
- `c++/src/algo/blast/api/blast_nucl_options.cpp:115-133,163-174`: megablast
  sets the window but does not set blastn's 20-bit X-drop.
- `c++/src/algo/blast/core/blast_parameters.c:218-221,380-383`: zero initial
  X-drop selects the subject cutoff.

The prior diagnostic owner is
`docs/evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/megablast-first-divergence.md`.
The five changes have an upstream cause; the equal-HSP allowance cannot explain
or authorize them. Reverting the correction solely to satisfy frozen Gate A
would restore a demonstrated NCBI divergence.

## Approved disposition

Keep the corrected Rust algorithm. Retain PR5 evidence as historical. After an
explicit decision, revise the BLASTN classification and Gate A authority through
the existing Product Decision and certification owners, characterize the full
current matrix on every required platform, and freeze the resulting reviewed
LOSAT outputs. Linux evidence supports the hash above for this case only; it
cannot certify other platforms. Keep Gate B's registered platform-specific
NCBI fingerprints independent and unchanged unless new evidence separately
requires their established review procedure. Do not normalize output, add a
runtime compatibility switch, or claim old evidence measured the new commit.

The approved decision supersedes the plan’s original-hash preservation condition
for this one source-corrected fixture. Other expectations and Gate B remain
unchanged. Phases 3–5 require new measured evidence.

## Additional observed CI failure

The CI artifact contains 104 frozen runs: 99 passed, four Sakai executions
(native, threaded, two repeats) mismatched, and threaded TBLASTX
`p11_avclpv_psclpv` timed out after 3,600 seconds. The exception printed only
the first failure. `ci-failures.json` records all five. Threshold regressions
after that failing command did not run. The timeout requires separate diagnosis;
it must not be treated as a successful parity check or silently extended.

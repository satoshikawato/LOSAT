# Product Decision: Pure-Rust runtime authority

- Decision ID: `PD-PURE-RUST-RUNTIME-AUTHORITY`
- Version: 1.0
- Date: 2026-08-30
- Status: Accepted

## Scope

This decision defines who owns the implementation of LOSAT's declared local
BLAST-compatible behavior. It governs project-authored production source,
project build logic, native and serial Wasm execution, and the sorting migration
recorded in [`docs/architecture/pure_rust_runtime_authority.md`](../architecture/pure_rust_runtime_authority.md).

It does not claim that Rust's standard library, allocator, compiler runtime,
filesystem, or thread implementation avoids operating-system ABIs.

## Authorities

Official NCBI BLAST+ source is the sole behavioral authority for supported
source-defined behavior. Official NCBI executables are external validation
oracles only. They may be used for tests, diagnostics, comparisons, and release
certification, but never as a LOSAT runtime, build, feature, fallback, FFI, or
subprocess implementation dependency.

LOSAT production behavior is owned by LOSAT Rust source. Project-authored
production and build code must not directly delegate supported BLAST behavior
to a host C/system function, native algorithm library, dynamically loaded
implementation, NCBI executable, or bundled C/C++ build.

Cargo `links` metadata is a review signal, not proof of delegation. Every
resolved normal/build dependency with `links` metadata must be traced and
classified. A new or changed finding blocks automatic acceptance until review;
only a demonstrated implementation path for supported BLAST behavior violates
the zero-delegation product boundary.

## Sorting policy

The NCBI source comparators, their fields and directions, and their invocation
timing remain authoritative. The migration changes only the project-owned sort
mechanism:

```text
existing NCBI-compatible comparator
+ existing invocation timing
+ Rust stable sorting
+ incoming-order preservation when the comparator returns Ordering::Equal
+ no additional comparison key
```

After the staged migration, native and serial Wasm paths use the same Rust
stable-sort semantics. Stable incoming-order preservation is a LOSAT
target-neutral implementation policy; it is not a claim that NCBI C `qsort`
defines a canonical order for comparator-equal elements.

No edit script, identity, alignment length, mismatch count, gap count, pointer,
address, hash iteration, accession, fixture, insertion ordinal, target, or
platform value may be added as an extra comparator key. Emulation of glibc,
MSVCRT, Apple libc, or another platform's `qsort` partition behavior is not
authorized.

## Existing compatibility decisions

[`PD-BLASTN-HSP-CANONICALIZATION`](PD-BLASTN-HSP-CANONICALIZATION.md) remains
controlling for BLASTN common-endpoint comparator-equal HSPs. Stable sorting
does not convert that narrow source-underdetermined class into an NCBI-defined
canonical representation, and a migration-induced representation change still
requires explicit review.

TBLASTX explicit non-default query genetic-code behavior remains exact against
NCBI. The only accepted TBLASTX product deviation is the existing local
`-s/--subject` non-default `--db-gencode` behavior, where LOSAT applies the
subject genetic code during translation, search, and reporting. This decision
does not broaden that deviation or permit any timing, ordering, scoring,
filtering, statistics, pruning, or formatting difference.

## Cross-platform output

LOSAT output across supported native platforms is a raw-byte contract. LOSAT
output must not be normalized or reordered after generation to make a
comparison pass.

During representative platform-native NCBI checks, a source-defined oracle
difference may be recorded as `LINE_ENDING_ONLY` only after proving that the
sole difference is `CRLF` versus `LF`, with unchanged row order and every
non-newline byte identical. This is a diagnostic classification, not a LOSAT
product exception and not permission to normalize LOSAT output or hashes.

## Enforcement and migration

The production-boundary checker in
`LOSAT/tests/check_pure_rust_runtime_boundary.py` freezes the live native-sort
adapters exactly. PR 2 removes both BLASTN-owned adapter groups. PR 3 removes
the shared `Hit`/`SubjectGroup` output adapter. PR 4 removes the TBLASTX
`UngappedHit` search/linking adapter. The temporary production-delegation
allowlist must reach zero after PR 4, and ordinary CI must continue to reject
new delegation paths and stale entries.

Reviewed non-algorithm Cargo `links` entries are dependency-review baselines,
not production-delegation exceptions and are recorded in a separate reviewed
dependency inventory. A version, `links` value, dependency kind, or
dependency-path change requires review.

Changing this ownership boundary, equal-element policy, or exception scope
requires a new explicit Product Decision. A production fallback or a new
platform-specific sort path may not be introduced as an incidental bug fix.

## Non-goals

- No production algorithm, comparator, pruning, linking, formatting, or output
  behavior changes in this decision.
- No generic BLAST, database-search, remote-search, or threaded-Wasm claim.
- No claim that LOSAT or Rust uses no operating-system ABI internally.
- No change to manifests, biological fixtures, expected outputs, or parity
  classifiers.

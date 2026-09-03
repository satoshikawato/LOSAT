---
name: verify-ncbi-parity-and-speed
description: Implement, diagnose, or verify LOSAT behavior against the NCBI BLAST source and executable oracle while preserving performance. Use for LOSAT TBLASTX, BLASTN, BLASTP, native or Wasm parity, output differences, hit-count discrepancies, coordinate or score mismatches, NCBI source ports, benchmarks, SIMD, threading, or release parity evidence.
---

# Verify NCBI parity and speed

Read `AGENTS.md` before any LOSAT work. It is authoritative. Then read:

- [references/exceptions.md](references/exceptions.md) for accepted parity exceptions;
- [references/evidence.md](references/evidence.md) for the comparison record;
- [references/commands.md](references/commands.md) for repository entry points.

NCBI source defines behavior. NCBI executables are comparison oracles only. Never add an NCBI runtime, build, fallback, FFI, or subprocess dependency to LOSAT.

## Freeze one reproducible difference

1. Name the program, task, input pair, options, target, thread count, output format, LOSAT commit, and NCBI BLAST+ version.
2. Build the current release binary before diagnosing old output.
3. Run both implementations into a fresh temporary directory.
4. Compare bytes first. When they differ, produce a structured field diff for coordinates, score, bit score, E-value, identity, gaps, frames, ordering, and formatting.
5. Reduce the fixture only if the reduced case retains the same first divergence.

Do not diagnose from hit-count percentages alone. Identify the first stage or record where outputs diverge.

## Trace the NCBI owner

Find the exact NCBI call path for the behavior. Record:

- source file, function, and current line range;
- caller and call timing;
- input encoding and coordinate system;
- sort, pruning, or tie-break order;
- floating-point type and rounding;
- output conversion step.

Trace the corresponding LOSAT path from CLI options to formatted output. A matching local function is insufficient when it runs at a different time or receives different state.

## Port behavior faithfully

Copy the relevant NCBI snippet into a Rust comment with the source path and lines, as required by `AGENTS.md`. Preserve timing, state transitions, sentinel handling, precision, comparison order, and edge cases.

Make Rust-only changes for ownership, bounds safety, data layout, or parallel scheduling only when they preserve the same observable behavior and asymptotic complexity. Unsupported behavior must fail explicitly.

Do not mix a parity correction with an unrelated optimization. Establish parity first unless the defect itself is caused by an optimized path.

## Prove the parity gate

Run the focused fixture, then the program/task comparison suite. Require exact output for every field covered by the declared gate. Apply only the named exception in `references/exceptions.md`; do not broaden it by analogy.

For native and Wasm paths, compare each path with the same input and options. Preserve NCBI result order after parallel reduction. Test serial, requested thread count, and repeated execution when scheduling changes.

Keep old resolved discrepancies as regression fixtures. Do not reopen a resolved hypothesis without a new current diff.

## Measure performance without changing the oracle

Use a release build, fixed fixture, fixed thread count, and the same output sink for baseline and candidate. Record warmup policy, repetitions, wall time, CPU time when available, peak memory when relevant, and output checksum.

Compare at least:

- current baseline versus candidate LOSAT;
- candidate LOSAT versus NCBI for output, not necessarily speed;
- native versus Wasm when the optimized code is shared.

Reject a faster change when output differs. Treat noisy timing as inconclusive and rerun; do not report the best single sample.

## Finish with an evidence record

Complete the template in `references/evidence.md`. Report the first divergence, NCBI owner, LOSAT owner, code change, exact parity result, benchmark result, targets tested, and remaining unsupported cases.

Ask the `ncbi_parity_auditor` custom agent for an independent read-only pass when the change is release-facing, alters pruning or ordering, changes floating-point behavior, or claims a resolved parity discrepancy.

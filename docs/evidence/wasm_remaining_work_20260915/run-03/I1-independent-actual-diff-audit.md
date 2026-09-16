# I1 actual-diff audit before hypothesis screening

**Status: no static pre-screen blocker found.** Complete the running I1 raw/thread gate before exclusive timing. This review supports only the bounded hypothesis screen described in README.md; it is not performance, full parity, tier-mechanism or adoption acceptance.

Audit mode: read-only source/evidence inspection; no builds, tests or timing run by the auditor. Only this note was written. Existing user changes and production source were preserved.

## Bound source and exact comparison

- Candidate: `work/I1-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`.
- SHA256 `88aac9cb7d35b419537e3a8fa3f16a48c7ffe28332756ef45740c7756b77a9d7`, freshly recomputed and matched to `I1-source-manifest.json`.
- Baseline: `work/baseline-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`, the frozen current working-tree source described in `I1-independent-preflight.md`; HEAD `8f23f774b44d6812d0943149877d93835eab1d52` alone is not its identity.
- Reviewed `I1-vs-baseline.diff`, actual helper at lines 910–990, caller around 1758, and both tests at 3082 and 3138.

Independent lexical comparison preserved string literals while excluding comments/whitespace. The complete original scan and extracted scan each contain 277 tokens; they match allowing only optional trailing macro-argument commas from rustfmt. Replacing the candidate's tuple-assignment call with the original scan recreates **all 5,432 tokens of link_hsp_group_ncbi exactly**. This proves the extraction did not silently alter another group-state statement, condition or arithmetic expression. It does not assert identical machine code or source-location panic text.

NCBI authority remains `c++/src/algo/blast/core/link_hsps.c`, `s_BlastEvenGapLinkHSPs`: scan 827–861, prior-best setup 812–823, downstream updates 863–894. The exact current source and wider call path were inspected in the preceding preflight.

## Caller and helper correctness

- The helper runs only within the unchanged `!can_skip_ncbi` path and strict `i_score > cutoff_big` branch. Previous-best `sum - 1` initialization, cutoff counters and header diagnostics happen first. Cached and filtered paths still bypass it. There is one call per scan-eligible HSP per existing DP pass, including an empty predecessor range; no synthetic warmup or extra search was added.
- Inputs are the existing helper/link/hit slices, cursor bound, frame-relative trim ends, trace flag, and four scalar selection values. The returned tuple has the same `(i32, i16, f64, usize)` ordering as the caller assignment. Immutable borrows expire before pool writes; no allocation, payload clone, mutable pool alias or persistent owner was introduced.
- Save current_idx/helper, read sum/next_larger, compute b0, decrement, jump and optional early continue occur in the original order. Non-traced b0 rejects avoid coordinate reads. Traced b0 rejects still compute and print the same predicates using the original helper/index, not the destination cursor.
- Strict sum and trim comparisons, helper-to-HSP index mapping, selected num/sum/xsum/link copies and trace argument order are unchanged. The helper adds no f64 arithmetic.
- The caller retains new_sum, the left-associative `(h_xsum + score * lambda) - logK`, all helper/link writes, next_larger construction, best update, xsum assignment and linked_to increments. Outer chain selection, preliminary/final linking timing, group partitioning and indexed parallel reduction are unchanged.
- Wasm requests no-inline and Native requests inline(always). These are compiler directives to inspect in emitted output, not proof of optimized callee entry or preserved native performance. The revised source comment correctly describes a repeated function-entry boundary without promising a measured speedup.

## Added boundary tests: actual coverage and limits

`I1-unit.log` shows **2 Native release library tests passed, 0 failed, 432 filtered out**. This is not a full suite or Wasm execution claim.

1. `large_gap_scan_preserves_initial_state_and_zero_one_jumps` verifies an empty range; exact return of nonzero initial sum/num/link and negative-zero f64 bits; real b0 jumps to sentinel 0 and sentinel 1; trace off/on. Using empty link/hit arrays and an invalid helper HSP index also ensures those rejected predecessors are not selected or dereferenced through that index. It closes the production jump-to-0 gap identified in the old C-X3 audit.
2. `large_gap_scan_preserves_strict_ties_mapping_and_xsum_bits` verifies that the first visited equal-score predecessor remains selected, helper index maps to a different HSP index, query/subject trim equality independently rejects a predecessor, and distinct selected xsum bits survive. Each case runs trace off/on.

These are meaningful focused checks, but not the complete planned differential proof. In particular, helper 3's `next_larger=2` is not a demonstrated nonadjacent interior jump: in the tested branch it reaches helper 2 by the normal decrement. Negative initial sums, a jump that skips multiple live helpers, broad mixed helper/link layouts, complete captured trace equality and whole-group multi-round state remain for the full proof if the screen is promising. The unit log interleaves test-harness output and trace text; it demonstrates branch execution, not a byte-perfect standalone trace transcript. Preserve per-target f64 bits and any existing cross-target delta in subsequent state checks.

## Gate and screening boundary

The parent reports the frozen baseline six-fixture Native/threaded n1/n8 gate completed 24 PASS; I1's corresponding gate was still running at inspection. This note does not count unfinished I1 cases as passed. `run_gate.py` reads exact oracle output bytes, records artifact/fixture hashes, and validates thread evidence in diagnostic mode; a complete candidate run is required before screening.

The gate fixtures are the two primary TBLASTX gencode-1 inputs, LC738874/LC738875 thresholds 10/100/10000 and valid-query no-hit. No non-default local-subject genetic-code exception applies to those conditions. Existing NCBI executable/oracle identities must remain recorded in this resumed run; no lost C8d/C9 raw evidence supplies acceptance.

README.md prospectively separates normal-runtime n8/n1 hypothesis screening from adoption. That is appropriate: a promising timing result must still be followed by the emitted-function/actual-executed-tier diagnosis, complete state/FP/order proof, Native/serial controls, module/reactor reuse, broader parity and actual browser validation. Fixed warmup and AB/BA/AB samples must remain intact, with no concurrent audit/build/test/profile work during measurement. Do not claim the n8 reversal resolved unless candidate n8 is faster than its own n1 for the same input; report each primary input separately.

Prior large-gap C-X3 extraction is acknowledged in both README.md and the preflight. This revision differs in Native inlining, retained current b0 fast path and explicit complete initial state. No old rejected experiment or unavailable round-02 result is relabelled as new evidence.

Audit completed: 2026-09-15T23:08:50.637373+00:00 UTC. All auditor tools have exited; no background work remains.

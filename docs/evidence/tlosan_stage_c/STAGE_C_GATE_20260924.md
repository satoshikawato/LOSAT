# TLOSAN Stage C gate: fixed-input search and traceback parity

Branch `feature/tlosan-tblastn-v0.2.0`, resumed from `f9832b6ce5fc8cb8dbde16286542186525819786`. Pinned NCBI C/C++ source is `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison-only BLAST+ `tblastn` SHA-256 is `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`. Read the [latest masking continuation](CONTINUATION_20260924_MASKS.md), [integrated chunk continuation](CONTINUATION_20260924_INTEGRATED.md), and their linked predecessors for exact commands, input checksums, full traces and first corrected differences.

**Stage C's fixed-input search and traceback gate passed on 2026-09-24.** This is the Stage C contract in `docs/tlosan_tblastn_v0.2.0_plan.md:102`, with the NCBI-recorded context cutoffs, x-drop, matrix and HSP limits supplied at their function boundaries. It is not a claim that default composition adjustment, sum statistics, linking, formatting, all 27 genetic codes, every CLI option, threading, or TBLASTN certification passes. Those belong to later gates or require their own input fixtures. The public TBLASTN CLI continues to fail explicitly with `TBLASTN local search is unimplemented (Stages C-E)` until outfmt 0/6/7 is complete.

| C function order / source | Required positive/boundary evidence | Result |
| --- | --- | --- |
| `blast_engine.c:747-850`, `aa_ungapped.c:492-505`: six-frame translation and candidate scanning by chunk | 15,000,962-nt subject, >5,000,000 translated residues, 12 scan calls, three query contexts, all 394 ordered candidate pairs | Full saved order, call states and context cutoffs match. |
| `aa_ungapped.c:516-614`, `blast_gapalign.c:3827-3927`: WordFinder and gapped score before the next chunk scan | Four long-case initial/gapped HSPs, source-indexed per-context cutoff and distinct x-drop intervention | Full initial/gapped HSP fields and event order match. |
| `blast_engine.c:539-586,840-850`, `blast_hits.c:2455-2537,2809-3050`: chunk endpoint purge, sort, offset adjustment, merge, frame append | Real-call 3→2 positive endpoint deletion and six append snapshots with cap three; chunk overlap and NO_RANGE | Removed HSP, merged HSP fields, cap and frame order match. |
| `blast_traceback.c:380-693`: full-translation retry, start-offset failure, interval-tree containment, re-evaluation, `Blast_HSPTest` deletion | Positive real-path start failure, positive test deletion, 12 natural positive containment events, query-indexed retry | Saved inputs, predicates, score, identity, survivor order and deletion timing match. |
| `blast_filter.c:1241-1255`, `blast_setup.c:614-638`: hard/soft SEG and query/subject lowercase mask state | Hard/soft SEG crossing, hard/soft query lowercase, SEG/lowercase overlap and overlap+soft real calls | Working and identity query bytes, merged `[0,45)` mask, candidate order, HSPs and traceback match. |

The first newly reproduced defects in this continuation were the per-context cutoff read, hard-SEG gapped/traceback working query, soft-mask WordFinder/gapped/traceback query, and query lowercase mask location. Each was fixed at its NCBI source-defined function boundary and has a retained fixture. No divergence remains in the Stage C fixtures listed above.

## Gate checks

- `cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn`: **61 passed, 0 failed**.
- `cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings`: passed.
- `cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check` and `git diff --check`: passed.
- All seven retained hard/soft SEG and lowercase mask SHA-256 manifests passed. The five newly added mask traces/output streams reproduced byte for byte in fresh pinned NCBI runs, excluding manifests whose command paths differ by design. Earlier long candidate, purge/cap, x-drop and start/deletion manifests/replays are recorded in their linked continuations.
- An independent read-only NCBI parity auditor reviewed the final function/input/order mapping and mask code. It found no concrete Stage C mismatch in the fixed-input scope and confirmed that default composition, sum statistics and linking calculations are Stage D work.

Stage D must generate NCBI-equivalent *computed* cutoffs and traceback input lists for default mode. Compare those actual generated values back across the C/D boundary; reopen C if that comparison reveals an upstream search or traceback mismatch. Do not start Stage E until D's statistics, linking, composition/Kappa, ranking and deletion-order checks pass. Code 32 stays a comparison-only NCBI C++ API oracle; NCBI `-db` statistics or headers are never local `-subject` expected bytes. NCBI is not linked or invoked by LOSAT runtime, build or fallback code.

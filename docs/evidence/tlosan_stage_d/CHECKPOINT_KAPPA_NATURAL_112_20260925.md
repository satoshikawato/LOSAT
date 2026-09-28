# TBLASTN Stage D: natural 112-subject preliminary stream through Kappa results

Branch: feature/tlosan-tblastn-v0.2.0; starting LOSAT commit:
2419e2db489d6e9c33afd7a370276f17b13744c7. Pinned NCBI source:
598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4. This checkpoint extends
[the early termination comparison](CHECKPOINT_KAPPA_EARLY_20260925.md) and
[the result hitlist comparison](CHECKPOINT_KAPPA_RESULT_20260925.md).
**Stage D remains partial. Stage E has not started. The public TBLASTN CLI
still returns an explicit unimplemented error.**

The focused internal Rust test
algorithm::tblastn::stage_d_linking::tests::natural_heap_replacement_fixture_matches_ncbi_preliminary_stream
now carries its own Stage C preliminary HSPs for every subject through
preliminary link/reap, the 27-list Kappa stream, early termination, 11 redo
calls, postredo containment, link, E-value calculation and reap, score
normalization, translated-subject identity, composition heap insert/pop,
Blast_HitListUpdate, and result reversal. Its input is the retained
[120-aa query and 112 nucleotide subjects](kappa_heap_rejection_20260925);
the local subject total is 40,320 nt, code 1, composition mode 2, sum
statistics enabled, max_target_seqs=2, one thread. The Rust initial
cutoffs and preliminary lists come from the fixture and Rust parameter and
Stage C code. The Kappa gapping/x-drop/cutoff and matrix setup in this
focused test still uses values recorded at the NCBI call boundary. Thus
the test does **not** establish a general runtime calculation of every
Kappa call parameter.

Pinned blast_kappa.c:3383-3427,3525-3736,2494-2515,
composition_adjustment/redo_alignment.c:1560-1582,
link_hsps.c:1765-1810, and blast_hits.c:3243-3297,3420-3437
define the call sequence. The test checks the saved NCBI
[preliminary/early call trace](kappa_heap_rejection_20260925/natural_c_d_early_20260925/manifest.txt),
[mode-2 redo and composition trace](kappa_heap_rejection_20260925/natural_mode2_20260925/manifest.txt),
and [Kappa heap/result trace](kappa_heap_rejection_20260925/result_order_20260925/manifest.txt).
It checks all 12 redone alignments across 11 calls for scores, rules, query indices, frames and endpoints;
all 11 postredo link, E-value and reap event triples, input HSP fields, link
multiplicities and full double E-value bits; 11 heap decisions, one rejected
candidate, ten insertions; every inserted HSP's normalized raw score,
bit-score bits, E-value bits, identities, frame, coordinates and within-list
order; and all ten hitlist updates with each retained HSP's score, E-value,
coordinates, frame and link multiplicity. The two-HSP subject OID 5 is
included. Bit score and identity are checked at heap input, while the
hitlist type does not yet retain those fields or edit scripts. The test
reconstructs postredo HSPs from alignments with placeholder gapped starts,
so it does not establish a complete report-ready result payload. The final
reversed result OID order is 10,11. This is a function/input/order comparison on one local fixture;
the saved outfmt 6 is a diagnostic trace source, **not Stage E byte parity**.

Three comparison-only NCBI runners were replayed into fresh /tmp directories.
All three directories matched their retained evidence byte for byte with
diff -rq exit 0. Each runner pins the executable SHA-256 and confirms
tracing preserves NCBI stdout and ordinary stderr. Exact commands and
focused Rust/clippy/fmt outcomes are in
[VALIDATION_KAPPA_NATURAL_112_20260925.txt](VALIDATION_KAPPA_NATURAL_112_20260925.txt).
No NCBI binary or library entered LOSAT runtime, build, or fallback paths.

Next Stage D boundaries are a general local orchestration path using
Rust-derived Kappa parameters and query/subject state; complete natural
multi-query, SEG and option cases; all 27 genetic codes (code 32 through
the CLI-calibrated comparison-only FindGeneticCode(32) API oracle);
natural positive postredo containment and composition-heap replacement
fixtures where reached; result post-pipes and edit-script-bearing payloads.
Identify the first mismatch from new comparison traces and pinned source.
Do not use NCBI CLI code-32 rejection or -db statistics/headers as a local
-subject oracle. Only after C and D pass all required cases may Stage E
compare local outfmt 0/6/7 bytes. Keep the public CLI gate until all three
formats pass; do not claim completion or certification from this checkpoint.

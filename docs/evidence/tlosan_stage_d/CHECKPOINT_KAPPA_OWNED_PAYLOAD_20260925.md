# TBLASTN Stage D: owned Kappa HSP scripts and report fields

Branch `feature/tlosan-tblastn-v0.2.0`; starting LOSAT commit
`4e229304f3387cd29027df8de3214631b0997ff6`; pinned NCBI source
`598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`.
This is an **internal Stage D checkpoint, not a Stage D pass**. Stage E has
not started. The public TBLASTN CLI still returns its explicit unimplemented
error.

The natural 112-subject internal comparison now moves each Kappa alignment's
edit script into a converted HSP, following pinned `blast_kappa.c:305-358`.
NCBI passes the two `unknown_value=0` gapped-start arguments to
`Blast_HSPInit`; the Rust conversion uses those source-defined zeros.
The converted HSPs are score-sorted, containment-survivor positions carry
their scripts into the link input, and link `source_index` carries each
script through reap, normalization, heap insertion, hitlist update/replacement,
HSP sort, and reversal. The former score/frame/endpoint/context lookup for
recovering alignment ownership is gone from this natural path. A focused
duplicate-numeric-key unit test confirms that two alignments with identical
score, frame, context and endpoints retain their distinct scripts.

All **12** redo scripts are checked against the saved
[NCBI Kappa return/edit trace](kappa_heap_rejection_20260925/result_order_20260925/ncbi.trace)
as sorted (raw score, query endpoints, subject endpoints, script) records at
conversion. This includes both OID 5 HSPs. The final two hitlist scripts are
still compared after result reversal. The other ten are not yet individually
compared after hitlist updates or eviction.

After NCBI's postredo link/reap and score normalization timing
(`blast_kappa.c:3687-3713,515-526`), Rust computes identities, positives,
alignment length, mismatches, gap letters, and gap opens from the transferred
script and translated subject. The source rule is
`blast_hits.c:745-835,966-989`. These values travel in each
`KappaHspPayload`. The retained [comparison-only NCBI report-field
fixture](kappa_heap_rejection_20260925/report_payload_20260925/manifest.txt)
uses local `-subject`, code 1, one thread. Its 2-target output checks the
two final hitlist HSPs. Its 112-target output checks all **11 HSPs in ten
accepted composition-heap insertions**, including gapped alignments and
the two-HSP OID 5. The comparison checks raw score, identity, positives,
length, mismatch, gaps, gap opens, frame, and query/subject output
coordinates. The coordinate conversion follows pinned
`blast_seqalign.cpp:1348-1383`.

The 112-target output has a different `-max_target_seqs` option from the
2-target call-state trace. Each compared HSP is confirmed by subject, within
subject order, score, frame, and both coordinate pairs; this establishes the
reported fields for those shared HSPs. It does not establish identical
candidate selection or heap behavior across those two option values. The
new runner pins the NCBI source and executable, records the exact commands
and input checksums, and its saved files replay byte for byte. Details and
checks are in [VALIDATION_KAPPA_OWNED_PAYLOAD_20260925.txt](VALIDATION_KAPPA_OWNED_PAYLOAD_20260925.txt).

The independent read-only NCBI parity audit found no concrete source
mismatch in the owned conversion, source-defined zero starts, containment
index path, or report-field arithmetic. Its scope is this internal
checkpoint. The natural 112-subject path is still test-only. A general
local-subject orchestration path, natural positive postredo containment,
composition-heap replacement, and multi-HSP E-value sorting with payload
reordering remain. So do remaining local fixtures/options, all 27 genetic
codes (ID 32 through the CLI-calibrated `FindGeneticCode(32)` comparison-only
API oracle), result post-pipes, and Stage E outfmt 0/6/7 bytes. Do not use
NCBI CLI's code-32 rejection or `-db` statistics/headers as local
`-subject` truth. No NCBI code or executable entered LOSAT runtime, build,
or fallback paths.

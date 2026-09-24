# TBLASTN Stage D Kappa continuation — 2026-09-25

Branch: `feature/tlosan-tblastn-v0.2.0`. Started from `490ea5b8f117f952dee164551da9582e0c9e7f5d`; pinned NCBI source: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. The [Stage C fixed-input gate](../tlosan_stage_c/STAGE_C_GATE_20260924.md) remains the accepted C boundary. **Stage D is partial and has not passed. Stage E has not started.** The public TBLASTN CLI still returns `TBLASTN local search is unimplemented (Stages C-E)`.

## New function/input/order evidence

The local `-subject` fixtures have positive `BlastSeqSrcGetTotLen`; `BLAST_OneSubjectUpdateParameters` does not run. The initial setup computes the cutoffs. [The parameter/linking checkpoint](CHECKPOINT_PARAMETERS_LINKING_20260924.md) records this upstream call state. This continuation uses those same NCBI initial inputs and does not treat earlier composition matrices or SEG bytes as end-to-end Kappa parity.

| NCBI function and call order | Rust comparison and observed boundary |
| --- | --- |
| `blast_kappa.c:1677-1704,1475-1552`; `redo_alignment.c:1172-1208` | `translated_subject_get_range` selects the query range, copies U as C, translates the selected subject window, then checks `subject_maybe_biased` and near identity before SEG. The saved 19 get-range calls, multiple query indices, near-identity decisions and seven changed SEG byte sequences match. The code-32 API trace additionally confirms a 120-aa translated subject from the selected table. |
| `redo_alignment.c:1219-1254`; `blast_kappa.c:1896-1957` | The TBLASTN callback redoes alignments with NCBI local gapped starts, scaled gap costs/x-drop and mode-2 matrix. All 21 retained code-1 traceback calls match query/subject bytes, adjusted 28×28 matrices, score, coordinates and edit script. The integrated Rust redo on the hard-SEG and multi-query fixtures matches NCBI alignment order and fields. The code-32 API case matches all 784 matrix entries, preliminary 656, redo 20493, rule 4, full-length coordinates and edit script. |
| `blast_kappa.c:224-273,3658-3689` | Converted alignments are score sorted, tested for containment, then relinked and reaped. On the multi-query fixture, NCBI has 10/9 post-redo HSPs and 7/6 survivors. Rust matches the saved before/after order, scaled raw scores and full-precision E-value bits. A source-based positive containment unit test covers score, frame and endpoint conditions; the retained natural Kappa fixtures do not show a positive post-redo containment deletion. |
| `blast_kappa.c:102-114,459-532,3689-3736` | Rust normalizes the 13 surviving multi-query HSPs with `BLAST_Nint` and the scaled lambda. Every raw score, bit-score `double` bit pattern, E-value bit pattern and translated-subject identity count matches the NCBI heap input trace. Code 32 matches raw 640, bit-score `251.13596086715009`, E-value `3.3904733561272406e-93`, identity 120. |
| `composition_adjustment/compo_heap.c:89-101,252-275,330-391,439-466` | The Rust heap matches the saved 11-subject local fixture's 11 decisions, insertions and pop order. A separate 112-subject local fixture with hitlist size 2 records 11 `WouldInsert` calls, one rejected candidate, ten insertions and ten pops; Rust matches state, capacity, E-value bits, rejection and pop order. A source-based unit test covers the comparator tie and replacement branch. No natural positive replacement event has yet been retained. |

The comparison-only probe is [ncbi_kappa_traceback_trace.c](ncbi_kappa_traceback_trace.c); [its runner](run_ncbi_kappa_traceback_trace.py) pins the NCBI executable SHA-256 and verifies that traced stdout and ordinary stderr equal unprobed NCBI bytes. The saved [traceback and heap inputs](kappa_traceback_20260924) include source/command/checksum manifests. The new [heap rejection fixture](kappa_heap_rejection_20260925) has a reproducible runner, FASTA inputs, raw output, trace and checksums.

## Code 32 local-subject oracle

The [comparison-only C++ API oracle](tblastn_code32_local_oracle.cpp) calls `FindGeneticCode(32)` and supplies that code on subject `BlastSeqSrc` retrieval, the input consumed by `blast_engine.c:1460-1466` throughout search. This is necessary because the stock local `MultiSeqBlastSeqSrcInit` reconstructs `SSeqLoc` at `seqsrc_multiseq.cpp:140-164`, losing the subject code metadata. A trial that set only the options and `SSeqLoc.genetic_code_id` still searched with code 1: preliminary 581 and final raw 590, while changing only the displayed identity. That trial is **not** the code-32 oracle.

The corrected [API runner](run_code32_local_api_oracle.sh) archives and verifies the pinned source, compiles the oracle outside LOSAT, compares plain and traced API stdout/stderr for traceback, mode-2 and composition probes, and requires code-1 API outfmt-6 bytes and the complete Kappa call trace to equal the same local `-subject` NCBI CLI run. The saved [code-32 local API evidence](code32_local_api_20260925_full) shows preliminary 656, mode-2 784-element matrix equality, redo 20493, final raw 640, identity 120, and full-precision bit/E-value equality with Rust. NCBI CLI code-32 rejection and database statistics/headers were not used as the local oracle. This oracle is comparison-only and is not linked, called or built by LOSAT.

## Gate status and next boundary

The current Rust implementation is an internal Stage D path exercised at function boundaries with retained NCBI call inputs. The Stage C generated preliminary HSPs, per-query Kappa redo, containment, relink/reap, normalization, identities, heap, result post-pipes and formatting are not yet wired together as one public TBLASTN execution path. The fixed-input C gate does not establish that integration. Full D acceptance needs that source-order path compared on all retained local cases, including code 32 and a positive natural containment or heap replacement fixture if that branch is reached. Resolve every new discrepancy from a fresh fixture and pinned source, without a guessed fallback.

Only after C and D pass as a whole may Stage E compare local `-subject` outfmt 0/6/7 bytes. The public CLI's explicit unimplemented error remains until all three formats are complete. No D pass, E pass, TBLASTN completion, certification or performance claim is made.

Focused verification and saved checksum results: [VALIDATION_KAPPA_REDO_20260925.txt](VALIDATION_KAPPA_REDO_20260925.txt). Reproduce the comparison runs from repository root with fresh output directories:

```bash
python3 docs/evidence/tlosan_stage_d/run_ncbi_kappa_traceback_trace.py /tmp/tlosan-kappa-trace-replay
python3 docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/run_ncbi_heap_rejection_trace.py /tmp/tlosan-heap-replay
bash docs/evidence/tlosan_stage_d/run_code32_local_api_oracle.sh /tmp/tlosan-code32-local-replay
```

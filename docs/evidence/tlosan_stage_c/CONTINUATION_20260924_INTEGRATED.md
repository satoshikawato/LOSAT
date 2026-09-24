# TLOSAN Stage C continuation: integrated chunk order and SEG query state

Branch `feature/tlosan-tblastn-v0.2.0`, after `bd4770e9`. Read the [chunk and multi-query checkpoint](CONTINUATION_20260924_CHUNK_MULTI.md) and its linked predecessors. Authority remains pinned NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; the comparison-only `tblastn` SHA-256 remains `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

**Stage C acceptance is still under review; D and E have not started.** These are local `-subject`, one-thread, comparison-only function-state checks. The public TBLASTN CLI still returns its explicit unimplemented error. Neither the NCBI executable nor its source is used by LOSAT runtime, build, or fallback code.

## Function, input, and execution order

| Pinned NCBI source and call order | Retained function input/output | Rust comparison |
| --- | --- | --- |
| `blast_engine.c:478-493,804-850`; `aa_ungapped.c:492-505`: translate each frame, scan every available chunk, then WordFinder | [Long three-query fixture](long_multi_query_20260924/) has a 15,000,962-nt subject, >5,000,000 translated residues, 12 WordFinder calls, three query contexts, and all 394 ordered candidate pairs. | `scan_protein_words_by_chunk` and `find_protein_init_hsps_by_chunk` emit each chunk before its next scan. The integrated test compares every ordered candidate, all four initial HSPs, and context cutoffs. |
| `blast_engine.c:489-586`: within each chunk, GetGappedScore, endpoint purge, score sort, offset adjustment, and HSPListsMerge follow WordFinder | The long fixture retains four gapped HSPs and two frame merges. [Positive purge/cap probe](chunk_purge_cap_real_path_20260924/) changes an actual purge input from three HSPs to two, then records all six real HSPListAppend outputs with cap three. | `preliminary_protein_hsps_in_ncbi_order` now executes one combined event stream. The integrated tests compare candidate → initial → gapped → purge → merge → append order and full HSP fields, including the positive deletion and cap case. |
| `blast_engine.c:840-850`; `blast_traceback.c:401-409,675-693`: append each frame before the next scan, then query-indexed full translation traceback and interval-tree containment | Long and [natural BLOSUM45 containment](alternate_matrix_word2_20260924/) traces retain append snapshots and positive containment events. Earlier [start failure](CONTINUATION_20260924_START_FAILURE.md) and [Blast_HSPTest deletion](CONTINUATION_20260924_CAND_HSP.md) traces remain valid. | The integrated append output is fed to traceback in NCBI query order. Tests compare six frame snapshots, start offsets, containment inputs/predicates, final HSPs, and the retained real-call deletion cases. |
| `blast_setup.c:614-625`; `blast_engine.c:484-525`; `blast_gapalign.c:2410-2442`: hard SEG masks the working query before WordFinder and GetGappedScore | [Hard SEG prefix](seg_hard_query_20260924/) records the full 160-byte query at six WordFinder and one GetGappedScore call. | The integrated path passes the same SEG-masked query frame to seed extension and gapped scoring. |
| `blast_traceback.c:380-391,583-596`: masked working query drives traceback score; `sequence_nomask` drives identity counting before Blast_HSPTest | [SEG crossing](seg_cross_traceback_20260924/) inserts 18 K residues at query offset 60 and matching subject codons. The HSP spans query 0..138. The probe records six WordFinder query buffers, one gapped buffer, and the unmasked identity input/result. NCBI raw traceback score is 638 and identity 138/138. | The integrated HSP and masked traceback score equal NCBI; the shared protein edit-script walker counts 138/138 on the unmasked query before Blast_HSPTest. |

## First reproduced difference and correction

The initial integrated hard-SEG path masked lookup and WordFinder input but passed the unmasked query to GetGappedScore. The hard-prefix NCBI probe showed that both calls receive the same masked working sequence. The integrated path now passes the masked frame to both calls. The SEG-crossing fixture then exposed a second concrete input difference at query offset 60: NCBI traceback uses masked `X` (NCBISTDAA 21), while the old Rust path used unmasked `K` (10). On the saved full-span HSP, that produced raw score 746 instead of NCBI's 638. Rust now traces the masked sequence and computes identity on the separately retained unmasked sequence, matching NCBI 638 and 138/138. These are function input and timing corrections, with no output exception.

## Reproduction and verification

From the repository root:

    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-seg-hard-replay --seg-hard
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-seg-cross-replay --seg-cross
    (cd docs/evidence/tlosan_stage_c/seg_hard_query_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/seg_cross_traceback_20260924 && sha256sum -c outputs.sha256)
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn
    cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings
    cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check

Both fresh NCBI replays completed; every trace/output file besides the manifest and its checksum was byte-identical to its retained counterpart. The retained manifests preserve their original command paths and now correctly state that the unprobed bytes equal **five** probe outputs. The Stage C Rust suite passed 56/56 before the temporary diagnostic print was replaced with a fixed 746 regression assertion; final post-edit verification is recorded by the checkpoint commit. The independent bounded auditor found no concrete mismatch in the integrated, purge/cap, or SEG-crossing cases before this final gate review.

The fixture coverage remains bounded by its saved matrices, code IDs, query/subject inputs, and one native thread. `-soft_masking true` and query-side `-lcase_masking` are followed in the [masking continuation](CONTINUATION_20260924_MASKS.md); the full v0.2.0 search matrix still requires source-backed scope resolution before complete Stage C parity. No D statistics/linking/composition/Kappa mapping or E outfmt comparison may use this bounded checkpoint as a C pass.

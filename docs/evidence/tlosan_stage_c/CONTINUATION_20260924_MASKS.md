# TLOSAN Stage C continuation: query masking function states

Branch `feature/tlosan-tblastn-v0.2.0`, after the [integrated chunk/SEG checkpoint](CONTINUATION_20260924_INTEGRATED.md). Authority remains pinned NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`, comparison-only `tblastn` SHA-256 `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

**Subsequent gate:** [Stage C fixed-input gate passed](STAGE_C_GATE_20260924.md) after the final 61-test sweep and independent review. D and E had not started at this masking checkpoint. The public TBLASTN CLI still returns its explicit unimplemented error. All NCBI binaries and probes here are comparison oracles, not LOSAT runtime/build/fallback dependencies.

## Fixed source, input, and order

| NCBI source/function state | Saved real-call comparison | Rust boundary |
| --- | --- | --- |
| `blast_fasta_input.cpp:486-502`, `blast_setup_cxx.cpp:655-657`, `blast_filter.c:1241-1255`: `-lcase_masking` adds query lowercase spans to filter locations, alongside SEG intervals | [Hard lowercase query](lcase_query_20260924/) has query positions 40..59 lowercase, SEG off. NCBI passes `X` (NCBISTDAA 21) at those positions to six WordFinder calls, one GetGappedScore call, and traceback; identity reads the original query. | `encode_tblastn_lookup_query` merges query lowercase and SEG intervals into the lookup mask, retains `aa_seq_nomask`, and hard-masks the working query. The integrated path compares all 154 ordered candidate pairs, two initial HSPs, gapped HSP raw score 528, six append calls, and traceback identity 120/120. |
| `blast_args.cpp:389-390`, `blast_setup.c:614-638`: `-soft_masking true` sets `mask_at_hash`, skips `BlastSetUp_MaskQuery`, but still constructs the complement lookup segments | [Soft SEG crossing](seg_soft_traceback_20260924/) uses the same 138-aa sequence as the [hard SEG crossing](seg_cross_traceback_20260924/). Both have the same 183 ordered candidate pairs. NCBI passes unmasked `K` (10) at positions 60..77 through WordFinder, GetGappedScore and traceback; soft raw score is 746 and hard raw score is 638. | Lookup remains restricted by the SEG intervals. Under soft masking, WordFinder, gapped scoring and traceback use the unmasked query. The integrated soft test compares all candidate pairs, initial HSP, gapped HSP, six frame append events, raw score 746, and identity 138/138. |
| `blast_filter.c:1241-1255`, `blast_setup.c:614-638`: lowercase filter locations plus `mask_at_hash` | [Soft lowercase query](lcase_soft_query_20260924/) and the hard lowercase fixture have identical 154-candidate streams. The soft working query retains the original residues; NCBI's initial/gapped/traceback raw score is 656 versus hard 528. | The combined soft/lowercase integrated fixture matches the NCBI gapped HSP, candidate trace, six append calls, traceback score 656, and identity 120/120. |
| `blast_filter.c:1241-1255`: `BlastSeqLocAppend(filter_out, lcase_mask_slp); BlastSeqLocCombine(filter_out, 0);` | [Overlapping SEG/lowercase query](seg_lcase_overlap_20260924/) masks a 40-K prefix with SEG and lowercase positions 35..44. NCBI's six WordFinder and one GetGappedScore working-query buffers contain one merged `X` interval `[0,45)`; 168 candidates and the gapped/traceback raw score 623 are retained. | The shared helper merges both intervals to `[0,45)`. The integrated test checks the full buffer bytes, all 168 ordered candidates, initial HSP, gapped HSP, six append calls and traceback identity 115/115. |
| `blast_setup.c:614-638`: `mask_at_hash` keeps the merged lookup restriction but skips hard query replacement | [Soft overlapping SEG/lowercase query](seg_lcase_overlap_soft_20260924/) has the same 168 ordered candidates as the hard overlap case. All six WordFinder buffers, GetGappedScore and identity input retain the original residues; NCBI's initial/gapped/traceback score is 656 with 120/120 identities. | The combined soft test compares all saved input buffers, candidate pairs, initial HSP, gapped HSP, six append snapshots and traceback score/identity. |

## First reproducible differences

The auditor found that the initial integrated path had no `mask_at_hash` input, even though TBLASTN argument parsing accepts `-soft_masking true`. The saved soft SEG trace made the first function-state difference concrete: at query offset 60, NCBI WordFinder sees `K` (`0x0a`) while the old Rust hard-mask path saw `X` (`0x15`). Candidate equality concealed that difference; subsequent NCBI gapped scores are 746 soft and 638 hard. The new input flag preserves the masked lookup but passes an unmasked working query through WordFinder, GetGappedScore and traceback.

The auditor also identified unported query-side `-lcase_masking`. The hard lowercase trace records `X` at offsets 40..59 in the six actual WordFinder inputs. Earlier lowercase fixtures masked only the subject. The shared TBLASTN query helper now joins SEG and query lowercase masks for lookup and preserves the unmasked identity sequence. Hard and soft lowercase fixtures both match their saved NCBI states.

## Reproduce and verify

From the repository root, use fresh output directories:

    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-soft-seg-replay --seg-soft
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-lcase-replay --lcase-query
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-lcase-soft-replay --lcase-soft
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-seg-lcase-overlap-replay --seg-lcase-overlap
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-seg-lcase-overlap-soft-replay --seg-lcase-overlap-soft
    (cd docs/evidence/tlosan_stage_c/seg_soft_traceback_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/lcase_query_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/lcase_soft_query_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/seg_lcase_overlap_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/seg_lcase_overlap_soft_20260924 && sha256sum -c outputs.sha256)
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn
    cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings
    cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check

All seven retained masking SHA-256 manifests passed. Five new comparison fixtures (soft SEG, hard lowercase, soft lowercase, hard overlap and soft overlap) reproduced their 12 trace/output files byte for byte on fresh pinned NCBI runs, excluding only command-path-specific manifests. Their focused Rust comparisons passed. The pre-final-overlap Stage C suite passed 60/60, and its lint, format and diff checks passed. The final Stage C sweep passed 61/61; lint, format and diff checks passed. Independent review found no concrete mismatch in the fixed-input Stage C scope. The precise acceptance statement and D boundary are in the subsequent gate record. The C acceptance boundary still requires review of the declared option/data scope and exact input cutoffs under default composition settings; a saved bounded match does not certify the entire TBLASTN program.

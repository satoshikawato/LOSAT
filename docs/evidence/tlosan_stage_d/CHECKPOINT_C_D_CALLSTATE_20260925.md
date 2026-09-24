# TBLASTN Stage C→D call-state correction — 2026-09-25

Branch: `feature/tlosan-tblastn-v0.2.0`; starting commit `490ea5b8f117f952dee164551da9582e0c9e7f5d`; pinned NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. The [Stage C fixed-input gate](../tlosan_stage_c/STAGE_C_GATE_20260924.md) remains the accepted C boundary. **Stage D remains partial and has not passed. Stage E has not started.** The public TBLASTN CLI still returns its explicit unimplemented error.

## Observed oracle discrepancy and resolution

The earlier [code-32 API diagnostic](code32_local_api_20260925_full) used the default `CLocalDbAdapter(subjects, options)` argument `dbscan_mode=false`. The source at `blast_app_util.cpp:203-210` passes `true` for the CLI local `-subject` path unless `BL2SEQ_LEGACY` is set; `seqsrc_multiseq.cpp:175-180,290-297` gives the adapter positive total length only in that mode. The old API source returned `BlastSeqSrcGetTotLen=0`, so `blast_engine.c:1407,1434-1443` called `BLAST_OneSubjectUpdateParameters`. Comparison probes measured initial hit cutoff 9, updated cutoff 19, and one preliminary 656-point HSP. Code-1 API outfmt-6 and Kappa trace happened to equal the CLI, but the D call state did not. The old bundle is retained only to diagnose that mistake.

The [C++ API oracle](tblastn_code32_local_oracle.cpp) now passes `dbscan_mode=true` and supplies the registered `FindGeneticCode(32)` table on `BlastSeqSrcGetSequence`. This keeps NCBI as a **comparison-only** oracle and supplies the selected code to the search and formatting path. The [runner](run_code32_local_api_oracle.sh) rebuilds against the pinned source outside LOSAT, pins the CLI SHA-256, checks plain versus probed stdout and ordinary stderr for Kappa, D-call and sequence-source probes, and compares code-1 API versus the local `-subject` CLI bytes for outfmt 6, Kappa, D-call and sequence-source traces. The resulting [CLI-calibrated evidence](code32_local_api_20260925_cli_calibrated) records `BlastSeqSrcGetTotLen=360` for code 1 and 32, no one-subject update, initial word/gapped cutoff 9, and two code-32 preliminary HSPs in NCBI order: score 656/frame +1 and score 16/frame -3. No NCBI CLI code-32 result or database statistical/header output is treated as this oracle.

## Rust function, input and order comparison

- NCBI `blast_engine.c:870-899` links the two code-32 preliminary HSPs before E-value reap. The Rust Stage C generated both from FASTA input with computed initial cutoffs, in the same order and with the same frame, score, query/subject endpoints and gapped starts. Rust preliminary link/reap retained both and matched NCBI `num` and full `double` E-value bit patterns (`1.4156323507376273e-95`, `8.7911823288453839`).
- NCBI `blast_kappa.c:3577-3595` passes both retained HSPs into `Blast_RedoOneMatch`. The fixed API trace now shows two translated get-range/composition/traceback calls in that order. The weak frame -3 HSP redoes to score 48, then the full frame +1 HSP to 20493; only the latter appears in the redone alignment list. The Rust internal C→D code-32 test feeds its **own** Stage C HSPs into Kappa and matches that final list, the 120-residue identity count, normalized raw score 640, bit score `251.13596086715009` and E-value `3.3904733561272406e-93`. The traceback test compares both calls' query/subject input bytes, all 784 adjusted-matrix entries each, score, coordinates and edit script.
- For the saved multi-query code-1 fixture, Rust Stage C creates 21 HSPs. Preliminary link/reap retains 20 in two groups of ten; all incoming Kappa fields match the NCBI trace before the Rust-produced groups enter per-query redo. Existing comparisons cover the redone alignment order, containment, relink/reap, 7+6 survivors, normalized raw/bit/E-value and identities. The saved local six-case Stage C gate and this focused integration have different scopes; neither is an all-case public Stage D result.

## Gate status and next work

The corrected code-32 call state and the two focused natural C→D tests pass. The other retained local fixtures have not all been carried through the same Rust source-order D pipeline. A subsequent [13-subject continuation](CHECKPOINT_MULTI_SUBJECT_STAGE_D_20260925.md) now connects one natural all-subject Stage C→Kappa→heap path and hard-SEG preliminary search to redo. The general runtime path, natural positive postredo containment and heap replacement fixtures if reached, result post-pipes and all outfmt 0/6/7 byte comparisons remain open. Resolve new differences against fresh comparison-only NCBI traces and the pinned source. Do not open the public CLI until C, D and all three Stage E formats pass; do not claim completion or certification from this checkpoint.

Reproduce the corrected code-32 oracle with a new output directory:

```bash
bash docs/evidence/tlosan_stage_d/run_code32_local_api_oracle.sh /tmp/tlosan-code32-local-replay
```

Focused test and checksum commands/results are recorded in [VALIDATION_C_D_CALLSTATE_20260925.txt](VALIDATION_C_D_CALLSTATE_20260925.txt).

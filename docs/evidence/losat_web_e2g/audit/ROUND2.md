# E2g independent audit, round 2 (after the round-1 findings were fixed)

Follow COMMON.md (same rules, read-only). FINAL_BINARY for this round: given in your prompt.
Work dir: /home/kawato/.cache/losat-web-gui-target/e2g-audit/r2/ ; append results to CASES.tsv (case_id, item, command, env, ncbi_exit, losat_exit, stdout_same, stderr_same, first_difference) and NOTES.md as you go, and give the findings in your final message (the harness may refuse a FINDINGS.md file; then report in text).

Round 1 found, and these commits changed (read each diff with `git show`):
- 5dac71f72: CHUNK_SIZE=1000 without BATCH_SIZE now reproduces NCBI (outfmt 0 prolog, then "BLAST engine error: Empty CBlastQueryVector", exit 3); a negative CHUNK_SIZE with a negative OVERLAP_CHUNK_SIZE is rejected only when the chunk size is above the overlap; R2 rejection messages reworded; V1 comment.
- 993310891: rejection of the CHUNK_SIZE/OVERLAP_CHUNK_SIZE values where NCBI's setup of a chunk would split it again (NCBI: CCoreException eNullPtr) — `chunk_would_be_split` in LOSAT/src/algorithm/blastn/query_split.rs; check its condition against NCBI's source (split_query_aux_priv.cpp:99-146,190-201; blast_aux_priv.cpp:206-208; split_query_cxx.cpp:134-180) and against the oracle on both sides of the boundary for several chunk sizes, query lengths (one long query, several queries whose boundaries fall inside chunks, a batch of many short queries with BATCH_SIZE large), both tasks.
- 19a2f8397 / earlier docs: INVENTORY.tsv additions (blast_options.c validators, GetDisplayIds, ~CBlastFormat) and reclassifications.
- AUDIT_B_FIXES (filled in by the orchestrator if audit (b) leads to fixes).

Tasks:
1. Verify each fix above (source, then oracle; at least 150 compared cases for the T7 branches together).
2. Re-run a sample of the round-1 cases (at least 300, spread over all items T1–T14, R1–R3, V1, including multi-batch inputs with BATCH_SIZE and title/invalid-query warnings, outfmt 0/6/7, other organisms from LOSAT/tests/fasta, IUPAC-heavy queries) with FINAL_BINARY, capturing stdout and stderr separately (never 2>&1), to show nothing regressed.
3. Verdict for the round: supported / unsupported / inconclusive, with every remaining difference classified (defect / explicit rejection whose reason holds / approved exception).

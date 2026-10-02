# E2g independent audit (ncbi_parity_auditor role) — common brief

You are an independent, read-only auditor of LOSAT's BLASTN (Rust reimplementation of NCBI BLAST+ blastn) against NCBI BLAST+ 2.17.0. Do not edit, commit, or create worktrees/clones in any git repository. Write everything you produce under your own work dir (given in your angle brief), appending results to files as you go.

- Repository (read-only): /mnt/c/Users/genom/GitHub/LOSAT-web-gui (branch feature/losat-web-gui). LOSAT crate in LOSAT/; BLASTN engine LOSAT/src/algorithm/blastn/ (run.rs is the main engine), reports in LOSAT/src/report/.
- NCBI source (fixed commit 598d8ae6): /mnt/c/Users/genom/GitHub/ncbi-blast/c++ — the only source of truth. Cite file:line for every claim.
- NCBI binary (comparison oracle only): /home/kawato/micromamba/bin/blastn
- LOSAT binary under audit: FINAL_BINARY (given in your angle brief), run as `LOSAT blastn <args>`.
- The inventory: docs/evidence/losat_web_e2g/INVENTORY.tsv (1006 rows: NCBI function/branch, LOSAT location, status faithful/n-a/divergent/unported/rejected, e2g_action, e2g_result) built by build_inventory.py; the review of the inventory: docs/evidence/losat_web_e2g/stage2/REVIEW.md; the session's transpile items T1–T14, R1–R3, V1 are summarized in build_inventory.py RESULTS and the git log (`git -C <repo> log --oneline 9c810a3d7..HEAD`).
- Approved exceptions (not defects): AGENTS.md "Approved CLI exceptions" (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES: argument-parser syntax errors and -help text/exit 2; -num_threads with -subject without NCBI's thread warnings; outfmt 6/7 write failures; allocation failure aborts). Explicit rejections ("not supported by LOSAT", "which LOSAT does not reproduce") are not parity defects when the reason holds.
- Comparison discipline: run NCBI and LOSAT with identical arguments from the same cwd, with an environment without BATCH_SIZE, CHUNK_SIZE, OVERLAP_CHUNK_SIZE, BL2SEQ_LEGACY, CTOOLKIT_COMPATIBLE, PRE_FETCH_SEQS_LIMIT, DIAG_*, NCBI_CONFIG_*, LOSAT_*, RAYON_* unless a case sets one on purpose; compare stdout bytes, stderr bytes and exit status. Hit counts alone are never evidence of parity.
- Verdict format: for each finding, severity (high/medium/low), the NCBI file:line, the LOSAT file:line, a reproducing command if any, and whether it is a defect, an accepted exception, or a justified rejection. End with an overall verdict for your angle: supported / unsupported / inconclusive.
- Use at most 8 concurrent processes. Record every command you run and its outcome in a TSV in your work dir.

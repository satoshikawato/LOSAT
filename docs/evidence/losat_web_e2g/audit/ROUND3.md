# E2g independent audit, round 3 (final binary, after the maintainer decisions and the audit (b) fixes)

Follow COMMON.md (same rules, read-only). FINAL_BINARY for this round: /home/kawato/.cache/losat-web-gui-target/e2g-bins/AC7/LOSAT (engine at the commit given in your prompt). Inventory: 1012 rows now. Use at most 8 concurrent processes (another auditor runs at the same time).

Write results to your work dir as you go (CASES.tsv with case_id, item, command, env, ncbi_exit, losat_exit, stdout_same, stderr_same, first_difference, classification; NOTES.md). The harness may refuse a FINDINGS file: then give the findings in your final message.

## What changed since rounds 1 and 2 (read each diff with `git show <sha>`)

Maintainer decisions of 2026-10-02 (docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md, committed in 56454a208; read it):
- d846e9bbe: APPROVED EXCEPTION 1. A query chunk that NCBI would split again (NCBI: CCoreException null pointer, exit 3): LOSAT searches each chunk once (exit 0). Not a defect. Fixtures env.resplit_* compare with NCBI at the largest overlap that does not split a chunk again.
- c72452236: APPROVED EXCEPTION 2. outfmt 0 subject deflines made only of punctuation that end in a run of spaces and separators (e.g. ", ,"): NCBI's x_CleanAndCompress reads past the string and crashes when such a subject has hits; LOSAT stops the cleanup at the end of the string (title ", "). Without hits on such subjects NCBI runs and LOSAT must equal it.
- 30713884f: PORT. -evalue +inf, -nan, +nan(1), 1e999 (NCBI accepts them and searches as with the largest e-value): LOSAT now searches too; a bare "inf"/"nan" stays an argument error in both.
- e37099f44: FIX of audit (b) D3 and D4. Warnings are written between the query reports as NCBI posts them (cerr tied to cout): a batch's title warnings before the report of its first query, an invalid query's warning before its preamble; outfmt 0 write failure (-out /dev/full, stdout to /dev/full) prints only "BLAST failed to write output" (exit 6), no warning. Merged streams (2>&1) are now in scope: compare them byte for byte.
- ad5fa9c85: a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE is rejected only when it splits the batch (calculate_num_chunks > 1); otherwise LOSAT searches as NCBI (e.g. -1 / -2147483648). Kept as explicit rejections by maintainer decision when the batch splits (NCBI either throws or searches chunk ranges with gaps).
- Kept by maintainer decision (valid, not findings): CHUNK_SIZE=1000 without BATCH_SIZE reproduced (empty batch, exit 3); -max_target_seqs 2^30..2^31-51 reproduced (prelim size 10); explicit rejections: K-A table miss with an invalid query first in its batch; -max_target_seqs > 2^31-51; reward 32767 / penalty -32768 and rewards that wrap to <= 0; subjects >= 2^31 letters; megablast gap costs > 32767.
- Closed pipe: outfmt 0 exception 5 (LOSAT exit 6, NCBI SIGPIPE); outfmt 6/7 exception 3 (LOSAT reports the write error; when its writes complete before the reader closes, LOSAT exits 0 — record such rows as "exception 3/5 timing", not defects).

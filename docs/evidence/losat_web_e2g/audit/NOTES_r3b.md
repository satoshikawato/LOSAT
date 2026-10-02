# E2g round 3, part B (auditor notes). FINAL_BINARY AC7 (engine ad5fa9c85), NCBI 2.17.0 oracle.
Layout: CASES.tsv (all rows: T7 1823, EV 520, PUNCT 432, R1 sample 484, R1 merged-stream 224; classification column);
T7_r3*.tsv raw T7 rows (r2 harness re-run, 198 rows replayed from r2/CASES.tsv with 90 s timeout because the r2 generator no longer reproduces their ids),
T7_classified.tsv, EV_r3.tsv, R1_r3.tsv, R1merged_r3.tsv, sets_CASES.tsv (54 punctuation sets), exhaustive_all.tsv (6196 deflines x hit/nohit),
exc1_results.jsonl / exc1b.json (approved exception 1 comparison against NCBI at the largest non-throwing overlap), empty_chunk_probe.tsv,
model.py (independent port of x_CleanAndCompress with the surrounding trims), h/ (harness copy, h/diff_*/ holds stdout/stderr of every differing case).

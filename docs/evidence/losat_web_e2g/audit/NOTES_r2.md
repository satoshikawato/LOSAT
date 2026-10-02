# E2g audit round 2 (binary AC2 = engine 993310891), notes
Files: CASES.tsv (T7 branches: 1823 rows), CASES_R1.tsv (round-1 re-run, 1354 rows), CASES_R1b.tsv (241 T9/T8 rows), CASES_ALL_classified.tsv (all 3418, last column = classification),
CASES_R1_excluded_merged_streams.tsv (224 round-1 rows captured with 2>&1; excluded: the brief forbids merged capture).
Generators: cases_r2_t7.py (grid + CHUNK_SIZE=1000 + negative pairs + independent model of NCBI's chunk setup), cases_r2_t7b.py (explicit BATCH_SIZE, model follows CBlastInput batching),
cases_r2_neg.py, cases_r2_big.py (default chunk sizes, 1.5M-11M queries), run_r1.py (round-1 sample), final_classify.py.
Model (cases_r2_t7.py model_ncbi_fails, written from split_query_aux_priv.cpp:99-146, split_query_cxx.cpp:145-171, :190-201 with size_t wrap): agrees with the oracle on 1152/1152 comparable cases
(782 no failure, 370 CCoreException); the 9 other rows are -ungapped (NCBI never splits: SB-1082, LOSAT rejects the option) and batches made only of an all-N query (no valid context, NCBI sets up no lookup table).
Round-1 sample: 1578 selected (all 963 that differed in round 1 + a stratified 615 of the identical ones), 1354 kept after dropping merged-stream cases; outcome changed vs round 1 in 35 only:
34 T7 (chunk re-split cases that LOSAT used to search and NCBI refuses: now rejected, 993310891) and T14_pipe_f7_head100 (closed pipe race).

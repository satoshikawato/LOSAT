# Range C result notes (core setup and parameters)

Final code read: `LOSAT/src` as copied to the scratch directory (identical to the working tree of `/mnt/c/Users/genom/GitHub/LOSAT-web-gui` at the time; the git repository returned I/O errors during the run, so I could not re-check that the tree is commit 90c5f0181; the binary used for every run is `LOSAT-90c5f0181`). Line numbers in `result_C.tsv` are from that copy.

## GAP rows
None. No input I tried (about 600 NCBI/LOSAT pairs, listed below) showed a byte difference in stdout, stderr or exit status, apart from approved exceptions (clap text and exit 2 for `-template_type 1`, and unsupported `-outfmt 5` / extra `-outfmt 6` columns, which are other ranges).

## UNSURE
None blocking. Weak evidence, stated in the rows:
- Row 21: no positive control for the array-size rule. The sizing difference (2^k-40, 2^k] can only show through diagonal aliasing with a diagonal difference of exactly 2^k, which needs 2L+1+gap >= 4096 (L >= 2039 at 4096, with the second hit within 40 of the first). I built such inputs (exact 24-nt segments, both orders, L = 2040, 2043, 2046, 2047); NCBI and LOSAT both give 0 rows. I could not reproduce the old (wrong) size to prove that the inputs would have distinguished it.
- Row 5: the tandem-repeat tests (42 cases) agree, but NCBI has no option to switch the 6 to 50, so sensitivity to min_diag_separation is not demonstrated; the value is covered by the unit test `task_defaults_follow_ncbi_handles`.
- Row 29: query batching and splitting tests agree, but a difference would show only in rare chunk-border HSPs.

## Window size 40 consequences (diag array/hash), what I checked
- Array size: `run.rs:6692-6698` uses `query_concat_length + config.window_size`; hit_level/hit_len arrays are allocated only when window_size > 0 (`run.rs:7772-7790`); initial offset = window (`run.rs:3807`); hash created with window = offset = 40 (`run.rs:3806`); hash insert stale window = window + min(scan_range, window - word_length) + 1 = 41 for dc (`run.rs:893`; word_length is the template length).
- Per-subject (per chunk) advance by length + window (`run.rs:652`, calls at 7818 and 10739); reset at offset >= INT4_MAX/4 sets offset = window, last_hit = -window, flag 0, hit_len 0 (array) or occupancy 1, offset = window, backbone 0 (hash).
- Deviation without output effect (X4): a whole subject shorter than the word size (11/12, 4-9 for short) returns at `run.rs:7394` without the offset advance; NCBI advances by length + window. Stale entries are always at least 40 behind the next subject's offset in either case, so only the wrap point moves.
- For dc `word_length == lut_word_length == template_length`, so `s_TypeOfWord` returns 1 and `extended` = 0; hit_len values are at most 21 (fits NCBI's Uint1); off-diagonal code is dead because `-off_diagonal_range` is rejected (Delta = 0).
- Container choice: concat length of K queries is 2*sum(L)+2K-1 (always odd), array when <= 8000. Boundary pair 3999/4000 nt (concat 7999/8001) and several multi-query batches agree.

## Tests run (all with NCBI 2.17.0 vs LOSAT-90c5f0181, outfmt 6 unless noted; scripts and outputs in the scratch directory `t/`)
- dc, single query L = 3997-4002 and 7998-8003 (mutated unit copies, both strands, 60 kb subject): 12+5 SAME.
- NCBI `-window_size 0` gives 526 rows against 157 rows with window 40 on the same input, so these inputs are sensitive to the window.
- Crafted aliasing inputs for L = 2040-2047 (see UNSURE): SAME, 0 rows.
- 3000 subjects of 1-600 nt (many shorter than the template) with query 800 / 3999 / 4000 / 4600 nt: 4 SAME.
- Long subjects: EDL933 (5.53 Mb, 2 chunks) with queries across 4999900; 10.3 Mb single record (3 chunks, overlaps at 4999900 and 9999800) with 60 mutated pieces around the overlaps (one batch = hash, six single = array), 3 templates x 3 lengths: all SAME.
- Offset wrap: 124 records x 4.6 Mb (570 Mb), `-num_threads 1`, planted diagonals repeated in every record, query 1200 nt (array) and 4500 nt (hash): SAME (484 and 489 rows; hits in all 124 records, 118 passes the wrap point). The 570 Mb file was deleted afterwards.
- Random stress: 180 cases (1-8 queries of 30-8100 nt, 1-400 subjects of 1-6000 nt, dc with random word size / template length / type / evalue / dust, and blastn-short with word size 4,7,8,11): all SAME.
- Tandem repeats (period 9-48, 12% divergence, HSPs on diagonals 6-49 apart): 42 SAME.
- dc with 7 scoring variants (-reward/-penalty/-gapopen/-gapextend), outfmt 0 and 7, lowercase masking (5 variants): SAME.
- blastn-short: default, -dust yes/no/"20 64 1"/"10 32 2", -lcase_masking, -evalue 10 / 0.01, word 4 and 11, outfmt 0 and 7 on repeat-rich queries (t/shq x t/shs): SAME. Word 4..9 with 700-1500 queries of 20-49 nt around the 32767 table switch: SAME for 700 (4-9), 1000 (4,5,7,8,9), 1200 and 1500 (7,8,9). Word 4-6 at 1200+ queries and 6 at 1000 not run (up to 1.4 M rows per run; machine load).
- dc with a 2.6 Mb single query and 100 x 12 kb queries (query chunk size 5 M vs 1 M): SAME; the same inputs with `-task blastn` SAME.
- Rejected options (`-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`, `-window_size`, `-off_diagonal_range`, `-no_greedy`, `-ungapped`, `-soft_masking`, `-use_index`, also `-strand`, `-min_raw_gapped_score`, `-culling_limit`, `-index_name`): each gives "error: the NCBI BLAST+ option -X is not supported by LOSAT's BLASTN", exit 2, with both dc-megablast and blastn-short.

## Inventory errors / outdated text
- Row 25 note "LOSAT advances by subject_length + 0" and row 22 "window 0" described the baseline; no longer true.
- Row 3 is listed under "short" only; for dc the table-based `s_NuclUngappedExtend` is used (word_length = template length >= 11) and LOSAT takes the same branch (`run.rs:8854`).
- Row 21 says the array is sized from `query->length`; precisely it is `buflen - 2` = last context offset + length (X2), same value as `query_concat_length`.

## Extra rows
X1 container_type choice (blast_parameters.c:226-231), X2 query->length for BlastExtendWordNew, X3 subject chunking with Blast_ExtendWordExit per chunk, X4 Blast_ExtendWordExit on an empty scan range (deviation without output effect), X5 discontiguous scan_step 1 / stride, X6 CDiscNucleotideOptionsHandle defaults (scan_range left 0, hit saving not overridden).

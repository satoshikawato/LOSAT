# Auditor B (BLASTP) log
Started. Counts below are cumulative at the end of each batch (comparison = one NCBI-vs-LOSAT run of one command line).

## Progress log
- Batch 1 (c1: basic outfmt 0/6/7, tasks, threads): 18 cases; 0 diffs (blastp-short = pre-existing rejection; NCBI thread warning = approved exception).
- Batch 2 (c2: every accepted custom field, outfmt 6/7, query/subject/both loc): 110 cases; 0 diffs.
- Batch 3 (c3: 6 input pairs incl. ambiguity letters/low-complexity/pipes ids/multi-record, 13 range combos, outfmt 0/6/7): 234 cases; 183 same, 51 explicit rejections (R2 start==len+1 and empty defline), 0 diffs.
- Batch 4 (c4: 60 range spellings x query_loc/subject_loc/both): 180 cases; 112 same, 68 rejections (66 R1 verified NCBI rc 255 CStringException, 2 R2), 0 diffs.
- Batch 5 (c5: 7 multi-record sets, skipped records at batch boundaries/ends, BATCH_SIZE 1/300/800/1300/2500/10000, 5 range sets, outfmt 0 and 6): 490 cases; all same incl. 82 'Empty CBlastQueryVector' and 84 'Invalid from coordinate'.
Cumulative: 1032 comparisons, 0 findings so far.
- Batch 6 (c6/c7: query splitting, long (22000-32000 aa) queries, intervals 19799/19800/19801/20001/9999/10000/10001, multi-record with long+short, long subject_loc): 60 cases; all same (3 R2 rejections).
- Batch 7 (c8: SEG variants, comp_based_stats, matrices, max_target_seqs, max_hsps, evalue, word_size/threshold/window_size, tasks, 4 input pairs x 5 range sets): 1056 cases; 672 same, 384 pre-existing LOSAT rejections (comp_based_stats 0/1/3, use_sw_tback, word_size 2/4/6, non-default matrix gap costs: verified also rejected without ranges), 0 diffs.
- Batch 8 (c9/c10: O/U/B/Z/J/X/* letters at range edges, skipped records before them, BATCH_SIZE variants): 506 cases; all same.
Cumulative: 2662 comparisons, 0 findings so far.
- Batch 9 (c11: 53 option-spelling/ordering/error-priority cases: -query_loc=..., duplicates, missing values, both ranges invalid, missing files, bad -evalue/-task/-matrix/-seg with bad ranges): 53; 35 same, 2 R1/pre-existing rejections, 16 parser-syntax differences (NCBI USAGE exit 1 vs LOSAT exit 2: duplicate option, missing value, --query_loc, -query_loc10-100, abbreviation, bad -evalue/-task/-max_target_seqs/-num_threads/-word_size values) = approved exception.
- Batch 10: stdin query (-query -), -out FILE: 24; empty subject intervals (start == length+1) in first/middle/last/all subjects x 9 ranges x 5 option sets x 3 formats: 810; CRLF/blank lines/dup ids/trailing '*' inputs: 96; title warnings (>=50 aa tail) with batches and skipped records: 240; CHUNK_SIZE/OVERLAP_CHUNK_SIZE/BATCH_SIZE env: 560; tiny intervals 2..40 residues: 468; all NCBI tabular field names: 208 (80 unsupported-field rejections); -num_threads 2/3/4: 292 (54 = duplicated -num_threads parser case, approved); outfmt 7 batches: 245; 46-record default-batch input: 19; large LvMJNV vs MeenMJNV: 8.
Cumulative: 5631 comparisons, 0 findings so far.
- Batch 11: random fuzz (3400 cases: random multi-record queries/subjects with low-complexity runs, ambiguity letters, random query/subject ranges, outfmt 0/6/7 + fields, tasks, evalue/max_target_seqs/max_hsps/window/threshold/ungapped/cbs, BATCH_SIZE/CHUNK_SIZE envs, threads): 3400 comparisons, all same (incl. ~300 'Empty CBlastQueryVector' and ~300 'Invalid from coordinate' outcomes); SEG-focused fuzz (ranges cutting masked runs, SEG parameters, outfmt 0 lowercase): 1500, 1 R2 rejection, rest same; threshold/window/word_size sweep: 320 (128 pre-existing word_size 6/7 rejections); matrices/gap costs: 53 (all pre-existing rejections except BLOSUM62 11/1); output formats 1-18: 20 (15 pre-existing rejections); empty-defline/Query_N-titled queries: 270 (180 pre-existing empty-defline rejections); pairwise error-priority matrix (32 error kinds x 32): 522 (46 = duplicated-option parser syntax, approved).
Cumulative: ~11700 comparisons, 0 findings so far.
- Batch 12: 138 exotic range spellings (en dash, full-width/Arabic digits, whitespace/control chars, leading zeros, +/- signs, 2^31 limits; query/subject/both): 0 diffs, 66 R1 rejections (NCBI rc 255 CStringException verified for each); X/U/O/J/*/B/Z-only and low-complexity intervals: 192 (12 R2); title warnings x range errors x BATCH_SIZE: 352 (64 R2); chunk-boundary splitting with hits at 9700-10300/19600-20300/29600-30300 relative to the interval start and CHUNK_SIZE/OVERLAP/BATCH_SIZE variants, 15-20k-residue intervals vs 100-record subject set (NCBI output verified to change with chunking; LOSAT matched): 83; high thread counts: 10 (2 = NCBI 'reduced to 32 CPUs' + 'ignored with subject' warnings = approved thread-warning exception, same without ranges).

## Totals
About 12,550 NCBI-vs-LOSAT comparisons (each = stdout + stderr + exit status for one command line).
Breakdown of non-identical outcomes: ~1,600 explicit LOSAT rejections (R1 range parts that NCBI cannot convert, NCBI rc 255 verified; R2 start == length+1 for queries; pre-existing rejections: comp_based_stats 0/1/3, use_sw_tback, word_size 2/4/6/7, non-default matrices/gap costs, outfmt 1-5/8-12/15-18, unsupported tabular fields, empty defline, headerless FASTA); ~160 argument-parser syntax differences (NCBI USAGE exit 1 vs LOSAT exit 2: duplicated options, missing option value, `--query_loc`, `-query_loc10-100`, abbreviations, bad -evalue/-task/-max_target_seqs/-num_threads/-word_size values) and NCBI thread warnings with -subject = approved exceptions.
All other comparisons (~10,800) were byte-identical on stdout, stderr and exit status.

## Findings
None. No difference outside R1, R2, the approved exceptions and the pre-existing rejection classes was found for BLASTP.

## Conclusion (perspective B, BLASTP)
supported

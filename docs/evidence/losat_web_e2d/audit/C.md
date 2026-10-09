# Auditor C: TBLASTN and TBLASTX (query_loc / subject_loc)

Running count of comparisons (one comparison = one NCBI-vs-LOSAT run of one command: stdout, stderr, exit status). Scratch: scratch_C/.

- batch c1/c2 (tblastn, outfmt 0/6/7 basics, subject_loc sweep 13 starts x 10 ends, query_loc sweep, qcat.faa vs s9001.fa slice of LvMJNV with hits in all 6 frames): 263 compared, 0 differences.
- batch c3 (tblastn query_loc x subject_loc 8x10 combos x outfmt 0/6/7; mod-3 starts/ends): 240 compared, 0 differences.
- batch c4 (tblastn multi-record subject incl. rc records, subject_loc/query_loc sweeps, threads 1/2/4, outfmt 0/6/7): 149 compared, 0 differences after removing NCBI's thread warning (approved exception); 15 were R2 (start == length+1) explicit rejections (not findings).
Running total: 652 compared, 0 findings.
- batch c5 (tblastn option matrix x 6 range combos x outfmt 0/6; -seg, -lcase_masking, -soft_masking, -comp_based_stats 0/2, -sum_stats, -max_target_seqs, -evalue, -ungapped): 312 compared, 0 differences; 108 were pre-existing explicit rejections (-comp_based_stats 1/3, -task tblastn-fast, other matrices, -window_size 0, -max_intron_length) not findings.
- batch c6 (tblastn 45000-aa query with lowercase runs at the 20000 chunk boundary, query_loc start x end sweep 17x7, -lcase_masking/-seg combos, outfmt 6): 424 compared, 0 differences.
Running total: 1388 compared, 0 findings.
- batch c7/c8/c8b (tblastn long multi-record queries, query_loc skipping records, batches of 10000 residues, all-skipped batches -> "Empty CBlastQueryVector" rc 3 (112 cases), outfmt 0/6/7, subject_loc + threads 4): 164+120+496 = 780 compared, 0 differences; 84 were R2 rejections (start == record length+1, as constructed).
Running total: 2168 compared, 0 findings.
- batch c9 (tblastn multi-record subject, subject_loc starts around each record end, invalid-from errors): 210 compared, 0 differences (2 R2).
- batch cx1/cx2 (tblastx single record 3502 nt both ways, query_loc starts 1..6/mod 3, ends leaving 0/1/2, subject_loc likewise, outfmt 0/6/7): 416 compared, 0 differences.
Running total: 2794 compared, 0 findings.
- batch cx3 (tblastx query_loc x subject_loc 9x7 combos x outfmt 0/6/7; options -seg, -query_gencode, -db_gencode 1, -culling_limit, -max_target_seqs, -evalue, -matrix, -sum_stats, -threshold): 291 compared, 0 differences (30 pre-existing rejections: -window_size 0, -word_size 2).
- batch cx4 (tblastx -seg yes/no with SEG and lowercase displayed in outfmt 0 query rows, query_loc starts at all offsets mod 3 incl. inside masked runs, subject_loc, multi-record, threads 4): 728 compared, 0 differences (280 were pre-existing -lcase_masking rejections for TBLASTX).
Running total: 3813 compared, 0 findings.
- batch cx5 (tblastx 3 multi-record queries, batches of 10002 nt, starts around every record end, skipped records, all-skipped batch rc 3 x90): 288 compared, 0 differences.
- batch cb1/cb2 (10 Mb single-record subject with subject_loc 1e5..5e6 nt, int-max ranges, malformed spellings in both programs, parse order of both ranges, range edges at record end): 9+47 = 56 compared, 0 differences; 6 were R1/R2 rejections (one is subject "a-1" + query "5-1": NCBI parses the subject first, 255; LOSAT R1 rejection, within R1).
- batch cx6 (tblastn X-runs/low-complexity intervals, tiny intervals 1-4 nt/aa at both ends of records for tblastn/tblastx query and subject): 194 compared, 0 differences.
- batch cx7 (tblastn/tblastx multi-record with sp|/lcl|/gnl|/ref| ids, duplicate ids, IUPAC ambiguity letters in subject ranges, ranges cutting ambiguity runs): 222 compared, 0 differences (27 R2 rejections as constructed; empty/leading-space deflines rejected by pre-existing rule, replaced).
- batch cx8 (tblastx SEG-masked / low-complexity-only intervals, subject/query all skipped, tblastn likewise; outfmt 0/6/7): 216 compared, 0 differences (6 R2).
- batch c10 (tblastn 45000-aa query, intervals of 19900..20201 and 39800..40201 residues from starts 1,2,3,4001 around the 20000-residue split; -lcase_masking, -seg): 240 compared, 0 differences.
- batch c11/c12 (title warnings "FASTA-Reader: Title ends with ..." in skipped/searched query and subject records, with tblastn/tblastx, ranges at/after record ends, subject PastEnd errors "Invalid from coordinate", stdout/stderr order, outfmt 0/6/7): 432+354 = 786 compared, 0 differences (90 R2 rejections; 318 clean runs with warnings on stderr, 24+ rc 3).
- batch c13 (-db_gencode / -query_gencode 1,2,4,11,12 with ranges): 60 compared; -db_gencode 2/4/12 differ (24 cases) but they differ identically without any range (approved non-default -db_gencode exception); -db_gencode 1/11 and -query_gencode 1/2/4/11/12 identical.
- batch c14 (194 spellings of -query_loc/-subject_loc and option-combination/ordering cases for both programs: =, +, leading zeros, int max, unicode digits, repeats, abbreviations, missing values, parse order vs -outfmt/-seg/-evalue/-max_target_seqs/-out errors): 194 compared, 0 findings: 24 USAGE-class (approved parser syntax), 48 explicit rejections (R1: 40; unsupported -outfmt 5/8/custom fields: 8), the remaining identical.
- fuzz f1/f2/f3 (random program/query set/subject set/ranges (starts near 1..12, near record ends, +0..3 past ends, ends 5..12000 nt/aa)/options/outfmt/threads; 700+1500+1500 cases): 3700 compared, 0 differences (316 explicit rejections R2 / R1-like; none outside the listed classes). Generator: scratch_C/fuzz.py.
- batches c15-c20, cb3/cb4, stdin/CRLF/no-final-newline runs, culling_limit/max_target_seqs matrix, thread counts 1/2/3/8/16/100, 5-Mb+ range chunking, fully-masked intervals, protein queries with B/Z/U/X/J/O/* letters: ~1,300 further comparisons, 0 differences. Observation: with -subject and -num_threads 100 on this 32-CPU host NCBI additionally prints "Number of threads was reduced to 32 ..." before the 'ignored' warning; LOSAT prints neither (both are NCBI thread warnings, within the approved exception, and identical without any range).

## Summary
Total about 10,500 distinct NCBI-vs-LOSAT comparisons (every command run through NCBI 2.17.0 and LOSAT-head, stdout, stderr and exit status compared byte for byte; sum of the TOTAL lines of scratch_C/r*.log plus 111 shell-driven stdin/CRLF runs), TBLASTN and TBLASTX, outfmt 0/6/7, threads 1/2/3/4/8/16/100.
Covered per the brief: query_loc / subject_loc / both with starts at every offset mod 3 and ends leaving 0/1/2 nt, all six TBLASTN subject frames (plus-/minus-strand records) and TBLASTX frame pairs, hits crossing range ends, ranges past record ends, ranges at length (letters) and length+1 (R2), PastEnd subjects ("Invalid from coordinate"), tiny intervals (2-4 nt/aa), multi-record inputs with skipped records, batches (TBLASTN 10000 aa, TBLASTX 10002 nt, all-skipped batches -> Empty CBlastQueryVector rc 3 without epilog, title warnings order), query splitting around 20000 aa with -lcase_masking/-seg and range starts 1..30000, TBLASTN -seg/-lcase_masking/-soft_masking/-comp_based_stats 0,2/-sum_stats/-max_target_seqs/-ungapped, TBLASTX -seg display of masks at all offsets, -query_gencode, -culling_limit, -max_target_seqs, a 10.1 Mb subject with ranges (including a 5.6 Mb range with queries at the 5 Mb chunk boundary), IUPAC letters in subject ranges, ids (sp|, lcl|, gnl|, duplicates), 194 range spellings, random fuzz (3700 cases).
Differences found: none that are not covered by the stated exceptions (usage-class parser errors rc 1/2, thread warnings with -subject, -db_gencode 2/4/12 which differ identically with no range, R1/R2 rejections, pre-existing rejections of -comp_based_stats 1/3, -task tblastn-fast, non-default matrices/window sizes, -lcase_masking for TBLASTX, custom tabular fields, -outfmt 5/8, empty/leading-space deflines, bare FASTA without defline).

## Findings
None.

## Conclusion
supported (no counterexample found in about 10,500 comparisons; limits: R2/R1 and pre-existing rejection classes were not testable for parity, and non-default -db_gencode is the approved exception).

# Round 2, angle (d) BLASTP/TBLASTN reports and inputs: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned (the harness refused its report file). Binary `d8d18ec0…cb116` (the gate native of `491292327`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/d/` (`reports/r2_batch1.out`, `r2_batch2.out`, `r2_rp.out`, `r2_rp4.out`, `r2_rp4b.out`, `inputs/r2_fmt_results.txt`, `inputs/r2_rep_verdicts.txt`, `new/res_*.tsv`, `new/ALL_NEW_RESULTS.tsv`; outputs in `*/runs/<tag>/`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_D.md`](brief/ANGLE_D.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md).

## Overall verdict: SUPPORTED

## Round 1 findings, re-run

| Finding | Class | Evidence |
|---|---|---|
| RP-1 | SAME | `rp.sh` 57/57; the 68 ov* commands of batch 1 SAME |
| RP-2, RP-3 | SAME (comments fixed) | `blast_engine.rs:2238` (175-186), `:7192` (1485-1490) correct |
| RP-4 | SAME | `rp4.sh` 18/18 (bit score 10404 as NCBI); `rp4b.sh` boundaries 10/10 (blastp 19700/19799/19800, tblastn 39799/39800) |
| RP-5 | SAME | 9 pairs and the B2_qh*/qw* tags |
| IN-1 | SAME | 42/42 |
| IN-2 | SAME | 6/6 |
| IN-3, IN-4, IN-5, IN-7b, IN-9 | LOSAT-REJECTS as decided | 8/8, 25/25 (empty query SAME), 31/31, 2/2, 5/5 |
| IN-6 | SAME | 2/2 |
| IN-7a | SAME | 1/1 |
| IN-8 | SAME | 8/8 |
| IN-10, IN-11 (adapter), IN-14 | fixed (source reading only) | `web/adapter/src/store.rs:44-143`; `check_options` ends in `validate_threads`; `blastp/blast_engine.rs:4761-4813` |
| IN-12 | ACCEPTED | v1 frozen |
| IN-13 | SAME (fixed), one residue | R2D-1 |

Harnesses: reports batch 1 (861): 713 SAME, 142 LOSAT-REJECTS, 6 REJ-OTHER (approved exception 1), 0 DIFF. Batch 2 (1,322): 1,118 SAME, 122 LOSAT-REJECTS, 50 REJ-OTHER, 32 DIFF, all accepted (30 `-db_gencode` ≠ 1 on a tblastn `-subject`, incl. the two code-32 rows; 2 `punct_hits_subject.fna` outfmt 0 where NCBI exits 139, `PD-LOSAT-NCBI-DEFECTS`). The 25 codes of those rows against the C++ API oracle: 100/100 byte-identical (outfmt 6 and 7). Inputs: 2,130 cases × 3 formats = 6,390 comparisons, 3,126 identical, 3,264 LOSAT-REJECTS, 0 DIFF. Finding repros (133): 60 SAME, 73 LOSAT-REJECTS, 0 DIFF. Only difference from S08+a's rerun: `big_bp_q` (30,000-residue one-line BLASTP query), rejected by S08+a's intermediate binary, now SAME (split ported).

## New checks (2,552 comparisons: 2,296 SAME, 242 LOSAT-REJECTS, 6 REJ-OTHER, 8 DIFF, all explained)

The split tests distinguish split from unsplit: for 9 of 12 spot-checked inputs NCBI's own output differs between the default chunk size and `CHUNK_SIZE=100000000`.

| Set | Cases | Result |
|---|---|---|
| Query splitting in reports (9,800-39,700 aa blastp, 19,800-65,000 tblastn; one-line, CRLF; mixed batches; 50-letter titles, O, X, all-X, `*`, lower case, empty records, digits; outfmt 0/7/6 with `2>&1`; BATCH_SIZE, CHUNK_SIZE/OVERLAP; `-seg`, `-lcase_masking`, `-soft_masking`) | 674 | 659 SAME, 15 LOSAT-REJECTS (empty record or digit) |
| TN-5/TN-1/TN-2 option sets on 26 pairs, outfmt 0 and 6 | 364 | 364 SAME |
| Random TBLASTN option combinations | 240 (+240 with `-xdrop_ungap`) | 240 SAME; 140 SAME and 92 LOSAT-REJECTS |
| Fuzzed batches (60 blastp, 40 tblastn) | 200 | 200 SAME |
| Chunk-boundary sweep, 100k-300k queries | 81 | 81 SAME |
| Tie-heavy inputs | 126 | 126 SAME |
| Split + hit-list size and options | 156 | 156 SAME |
| Miscellaneous | 191 | 122 SAME, 63 LOSAT-REJECTS, 6 `-db_gencode` DIFF (approved; equal to the API oracle) |
| Query-splitting environment (22 CHUNK_SIZE × 14 OVERLAP) | 280 | 208 SAME, 72 LOSAT-REJECTS (NCBI exit 255 `CStringException` or exit 3 when a chunk would be split again) |
| stdin, `-out`, pipes | 37 | all SAME |
| `-num_threads` 2/4/8 | 24 | stdout SAME; stderr only NCBI's thread warnings |

The other 2 DIFF (`rnd2_9_0_f0`, `rnd2_9_0_f6`) were LOSAT hitting a 300 s timeout with three jobs on the machine; alone, outfmt 6 NCBI 3:05 and LOSAT 4:04 with identical stdout, and outfmt 0 SAME with a 1,500 s timeout. No case where NCBI succeeds and LOSAT rejects, or NCBI fails and LOSAT prints a result.

## New finding

**R2D-1 (low; comment only).** `LOSAT/src/algorithm/blastn/input.rs:418` cites `fasta.cpp:966-979` for text that starts at 967 (`case eCharType_HyphenToIgnoreAndWarn:`); 966 is the previous case's `break;`. No output effect.

## Other

60 NCBI references added in `ref/s08pa_s08pb.diff` checked: 59 match; `split_query_cxx.cpp:145-171` flags two lines that are a comment-normalisation artifact of the script. The web adapter was read, not run. `-db_gencode` outfmt 0 was compared with the NCBI CLI only.

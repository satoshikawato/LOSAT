# Re-run of the inputs audit (round 1), LOSAT-head vs NCBI BLAST+ 2.17.0

Binaries: NCBI `/home/kawato/micromamba/bin/{blastp,tblastn,tblastx}`; LOSAT `/home/kawato/.cache/losat-web-gui-target/s08pa/bin/LOSAT-head`.
Work dir: `/home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/inputs/` (inputs copied from `$D/work/inputs`, `$D/LOSAT` replaced by LOSAT-head; nothing under `$D` or /mnt/c touched; no source edited).
vperf.lock was absent at every check (the wrapper scripts wait on it before every run).
Compared: stdout, stderr, exit status and `2>&1` (outfmt 0, 6, 7 for the harness cases; the outfmt given in the report for the finding repros, plus 0/6/7 for IN-1/IN-2).

## Verdict per finding (CLI repros, files `rep/<tag>/`, `repro_out.txt`, `rep_verdicts.txt`)

| Finding | ROUND1 decision | Verdict | One-line evidence |
|---|---|---|---|
| IN-8 | fixed | SAME (8/8) | `blastp -query in/min_qO.faa -subject in/min_sX.faa -outfmt 6`: both `qO sX 99.760 417 1 0 ...`; outfmt 0 both `Identities = 416/417 (99%), Positives = 416/417 (99%)`; qX/sO, qO/sO, TBLASTN O-vs-NNN also SAME |
| IN-2 | fixed | SAME (6/6) | `s_empty_only.faa`, `s_empty_two.faa`, outfmt 0/6/7: LOSAT prints `Effective search space used: 0` (8 blocks), byte-identical to NCBI |
| IN-4 | explicit rejection | LOSAT-REJECTS (25/25 subject cases); empty-query case SAME | `blastp -query $F/e2e_protein_query.faa -subject in/gt_only.faa -outfmt 6` NCBI rc 0 (warning), LOSAT rc 1 `Error: subject record 1 has a defline that is empty; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTP (use ASCII deflines without control characters)`. Same text with TBLASTN/TBLASTX; le/p_gtonly_first.fa, n_gtonly_first/file/file_nonl/mid/last/nonl all rejected. `-query le/p_ws.fa -subject le/p_gtonly_first.fa` (empty query) is SAME (`Query is Empty!`) |
| IN-9 | fixed (explicit rejection) | LOSAT-REJECTS (8/8) | `tblastx -query - -subject - -outfmt 6 < in/n_base.fna`: NCBI rc 0 no output; LOSAT rc 1 `Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX`. Same for outfmt 0/7, `-subject -` with query omitted, and BLASTP/TBLASTN (`...LOSAT's BLASTP` / `...TBLASTN`) |
| IN-1 | fixed | SAME (42/42) | q_O3, q_xO, q_allo, q_allOmulti, q_Oshort, q_ox, q_o25 x BLASTP/TBLASTN x outfmt 0/6/7: LOSAT stderr `Warning: [blastp] Query_1 O3: Could not calculate ungapped Karlin-Altschul ... filtering options One or more O characters replaced by X ...` identical to NCBI |
| IN-3 | explicit rejection | LOSAT-REJECTS (8/8) | `blastp -query in/q_tt_sp1.faa -subject in/s_base.faa -outfmt 6`: NCBI rc 0 silent; LOSAT rc 1 `Error: query record 1 has a defline that ends with white space after 50 amino-acid letters (...); NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTP`. Also sp2, sptab, TBLASTN query, BLASTP subject (s_ttl_trailing_space, s_ttl_trailing_space_only) |
| IN-5 | explicit rejection (at read time) | LOSAT-REJECTS (31/31) | `blastp -query ws/p_seqlast_nbsp.fa ...`: NCBI warns `On line 3: 358-359`, rc 0; LOSAT rc 1 `Error: query record 1 has byte 0xc2 in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by LOSAT's BLASTP (use letters and '*')`. p_seqtrail_* and p_lineonly_* for emsp, ideo, lsep, nbsp, nel x BLASTP query/subject, TBLASTN query: all rejected |
| IN-6 | fixed | SAME (2/2) | `blastp -query in/q_ttl_multi.faa -subject in/s_base.faa -outfmt 0 2>&1` (and TBLASTN): title warnings now after the prolog, byte-identical |
| IN-7 | (a) fixed, (b) explicit rejection | (a) SAME, (b) LOSAT-REJECTS | (a) `blastp -query le/p_ws.fa -subject in/s_empty_mid2.faa -outfmt 6`: both only `Query is Empty!`. (b) `-subject pc/p_mid_31.fa`: NCBI `FASTA-Reader: Ignoring invalid residues ... On line 3: 41` + `Query is Empty!`; LOSAT rc 1 `subject record 1 has '1' in a sequence line; ... not supported by LOSAT's BLASTP (use letters and '*')`. NBSP subject (`ws/p_seqlast_nbsp.fa`) likewise rejected |
| IN-10 | fixed (adapter) | adapter, not re-run | no CLI repro in the report; CLI-side rejections themselves verified in the IN-3/4/5/7b rows above |
| IN-11 | fixed (adapter validate) | adapter, not re-run; CLI half unchanged: LOSAT-REJECTS | `blastp\|tblastx -num_threads 100000 ...`: NCBI rc 0 (warning), LOSAT rc 1 `Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT` |
| IN-12 | accepted (web ABI v1 frozen) | adapter/web ABI, not re-run (ACCEPTED) | no CLI repro |
| IN-13 | fixed (comments) | not a CLI repro; source read | `report/query_warnings.rs` no longer says O messages come first; it cites blast_aux.cpp:1043-1054 (RemoveDuplicates sort). Minor: `blastn/input.rs:418` still cites `fasta.cpp:966-979` (report said it starts at 967); 376-385 and 1651-1673 are now correct |
| IN-14 | fixed (adapter run_local) | adapter, not re-run | no CLI repro |

## Harness re-run (the "Checked OK" corpus)

Driver: copy of `c.sh` with the binary changed to LOSAT-head (`run1.sh`), cases rebuilt from the old `r/` names with `job.sh`/`job2.sh` logic for the bulk sets (`le`, `ws`, `pc`, `nt`) and by hand for the other 396 (`gen.py`/`gen2.py`, `cases.tsv`); run 8 in parallel, total about 5 minutes (limit of 40 minutes per script not approached). Per-case results: `final_results.txt`, outputs `r/<case>/{n,l}.<fmt>.{out,err,rc,both}`.

Reconstruction check: the NCBI fmt 6/7 stdout+stderr+rc of the re-run equal the old NCBI results (path-normalised) for 2129 of 2130 cases (`base_p` is not reproduced: it is not a finding repro, an e2e fixture pair guess; LOSAT-vs-NCBI result for it is SAME). The fuzz sets fz/fz2 needed `-evalue 1000` and oo_b45/b80/pam30 `-evalue 1000` (b45 with `-comp_based_stats 0`) to match the old NCBI output; the inputs of a few small groups (stdin, tn_q) were re-derived from the old outfmt 0 headers.

Counts over 2130 cases x 3 outfmt = 6390 comparisons (stdout+stderr+rc+`2>&1` all compared):

| Class | outfmt comparisons | cases (all 3 fmts alike) |
|---|---|---|
| SAME | 3123 | 1041 |
| LOSAT-REJECTS ("not supported by LOSAT", rc != 0) | 3267 | 1089 |
| ACCEPTED | 0 | 0 |
| DIFF | 0 | 0 |

By set (cases): le 234 (SAME 108, REJECT 126), ws 552 (SAME 54, REJECT 498), pc 837 (SAME 477, REJECT 360), nt 111 (SAME 66, REJECT 45), other 396 (SAME 336, REJECT 60). The 3 outfmts agree for every case (no case is SAME in one format and REJECT in another).

Transition against the old binary's results (`old/r`, $D/LOSAT): DIFF -> SAME 111 (55 cases: the IN-1, IN-2, IN-6, IN-7a, IN-8 sets, 50-letter title/ttl and batch-split cases `tl_bpq_*`, some fz fuzz sets), DIFF -> REJECT 227, REJECT -> REJECT 3024, SAME -> SAME 3009, SAME -> REJECT 7 (3 cases, below), no REJECT -> DIFF, no SAME -> DIFF, no DIFF left.

Cases that were identical before and are now rejected (not listed in ROUND1.md as a table row, listed here for review):
- `big_bp_q` (BLASTP, 30000-residue single-line query in/q_bigline.faa): LOSAT-head `Error: a query batch of 30000 residues, which NCBI BLAST+ searches in 3 query chunks (the query chunk size 10000, overlap 100), is not supported by LOSAT's BLASTP`; the base binary ran it and matched NCBI. This is the `check_protein_query_split` rejection documented in head-src `docs/evidence/losat_web_e2e/s08pa/NOTES.md:74`.
- `bpq_q_ttl_cr_crlf`, `bpq_q_ttl_tabend` (BLASTP query defline ending in CR / TAB after 50 letters, fmt 6 and 7): now rejected by the IN-3 rule (`defline that ends with white space after 50 amino-acid letters`); NCBI warns for these, the base binary matched.

## DIFFs

None. No repro and no harness case classifies as DIFF.

## Notes
- IN-10, IN-11 (adapter half), IN-12, IN-14 concern web/adapter or web ABI v1 and were not re-run.
- `ws2_*` cases have no outfmt 7 in the old results, so 12 comparisons (3 SAME, 9 REJECT) have no old-binary reference; they are included in the counts above.

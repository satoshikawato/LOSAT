# Rerun of the BLASTP / TBLASTN reports audit (round 1) against LOSAT-head

Work dir: audit-rerun/reports/ (cmp2.sh = harness, runs/<tag>/{n,l}.{out,err,rc,both}). NCBI 2.17.0 vs LOSAT-head, compare stdout, stderr, exit status and 2>&1. Classes: SAME, LOSAT-REJECTS ("not supported by LOSAT"), ACCEPTED, DIFF.
vperf.lock was absent for the whole run (checked before every single run inside cmp2.sh).

## Batch 1: every command line in the original logs (861 pairs, replayed from work/reports/*.log; quoting reconstructed)
- 713 SAME, 142 LOSAT-REJECTS, 6 REJ-OTHER (both fail, exit 1 vs 2, different syntax-error text for -window_size -1/2147483648/abc: approved exception 1 -> ACCEPTED), 0 DIFF.
- All 68 runs that were DIFF in the original (RP-1 O-queries, ov*) are now SAME.
- 47 runs that were identical before are now LOSAT-REJECTS, as documented in ROUND1.md: blastp/tblastn -evalue +inf / 1e999 (D12, 46) and blastp -window_size 2147483647 (D11/TX-1, 1).

## RP-1 (O query warning order): SAME
RP-1 repro (`printf '>oo\nOO\n' > q.faa`, blastp and tblastn, outfmt 0/6/7, stdout/stderr/2>&1/exit) and the inputs o1 o2 o3 o4 o5 o6 o7 o8 o9 mixE mixF x outfmt 0/6/7 (blastp): 57 pairs, all SAME (runs RP1_*). The 27 original ov* runs (blastp and tblastn, ov1..ov27t) that were DIFF are SAME as well (batch 1). NCBI order is now reproduced: `Warning: [blastp] Query_1 oo: Could not calculate ungapped Karlin-Altschul parameters ... filtering options One or more O characters replaced by X ... positions 0, 1 `.

## RP-2, RP-3 (comments only): not testable at run time (no output effect); ROUND1.md lists them as fixed.

## RP-5 (query-title warning placement): SAME
qh11, qw1, qw2 x blastp/tblastn x outfmt 0/6/7 with 2>&1 (runs RP5_*, plus B2_qh0..16, B2_qh_all, B2_qw1/2 in batch 2): all SAME. In RP5_qw2_blastp_0 the two FASTA-Reader warnings now come after the `Database:` block (line ~23) as in NCBI.

## RP-4 (60000-aa TBLASTN query): explicit rejection, and SAME with CHUNK_SIZE=100000
Inputs: round1/rp4/*.gz (decompressed to rp4/, byte-identical to the audit inputs).
- Default env, tblastn bigq vs bigs.fna (outfmt 6/0/7), vs bigs4.fna, with -comp_based_stats 0, with -seg no -evalue 1e-5: all LOSAT-REJECTS (NCBI exit 0). Message: `Error: a query batch of 60000 residues, which NCBI BLAST+ searches in 3 query chunks (the query chunk size 20000, overlap 100), is not supported by LOSAT's TBLASTN` (exit 1).
- Default env, blastp bigq vs bigs.faa (protein version, outfmt 6/0/7; = Checked OK item 10): LOSAT-REJECTS: `Error: a query batch of 60000 residues, which NCBI BLAST+ searches in 6 query chunks (the query chunk size 10000, overlap 100), is not supported by LOSAT's BLASTP`.
- CHUNK_SIZE=100000 for both NCBI and LOSAT-head: tblastn bigs.fna outfmt 6/0/7, bigs4.fna outfmt 6, -comp_based_stats 0, -seg no -evalue 1e-5, and blastp bigs.faa outfmt 6/0/7: 9 of 9 SAME. sbig4 first HSP bit score is 10412 in both (NCBI default-env value 10404 was the chunked result; LOSAT-head with CHUNK_SIZE gives 10412 = NCBI).
- Boundary (query = first N aa of bigq): blastp N=19700, 19799 SAME (outfmt 6 and 0); N=19800 LOSAT-REJECTS. tblastn N=39799 SAME; N=39800 LOSAT-REJECTS. Matches the stated thresholds.

## Other numbered findings: none in the report beyond RP-1..RP-5.

## Checked OK items 1-10: replay
No driver scripts exist under work/reports (only cmp.sh, stat.sh, fz.py, fz2.py, verify_refs.py); the commands were taken from the 16 logs (861 command lines; shell quoting reconstructed, e.g. `-seg "8 1.5 2.0"` and `-outfmt "6 qseqid sseqid"` as one argument) and the unlogged tag families of runs/ were regenerated from the tag names and the report text (1322 further pairs, gen2.py). All pairs ran in minutes; no script came near the 40-minute limit.

Batch 1 (logged, 861 pairs): 713 SAME, 142 LOSAT-REJECTS, 6 REJ-OTHER, 0 DIFF.
 - 68 of the original DIFF (RP-1) now SAME; 4 `-seg 8 1.5 2.0` rows that were REJ in the original (the log had lost the quotes) are SAME with the argument quoted.
 - 47 rows that were identical before are now LOSAT-REJECTS as ROUND1.md says: -evalue +inf / 1e999 (D12: blastp "an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP"; tblastn likewise) and blastp -window_size 2147483647 (D11/TX-1: "a -window_size of N with a query of N letters ... is not supported by LOSAT's BLASTP").
 - LOSAT-REJECTS otherwise (unchanged from the original REJ rows): TBLASTN matrix/gap/word/threshold/window/composition combinations other than BLOSUM62 defaults, blastp BLOSUM45/80/PAM30 and non-default gaps, -comp_based_stats 0 for blastp, TBLASTN custom outfmt fields, outfmt 10 / delim=, empty-record queries ("has no residues"), '-' in sequence lines.
 - REJ-OTHER (6, ep3_43..48): both fail (NCBI exit 1, LOSAT exit 2 clap error), different syntax-error text for -window_size -1 / 2147483648 / abc = approved exception 1 (ACCEPTED).
Batch 2 (regenerated, 1322 pairs): 1118 SAME, 122 LOSAT-REJECTS, 50 REJ-OTHER, 32 DIFF that are all ACCEPTED (below), 0 unexplained DIFF.
 - Families: qh0..16/qw1-2/nh/nh2/rq1-12 (items 6, 11) x blastp/tblastn x outfmt 0/6/7; max_target_seqs 18 values x many (300) and big (720) subjects x blastp/tblastn x outfmt 0/6 + defaults (item 2); max_hsps; duplicate subjects 1..8; seg (10 value sets x 2 queries x 2 subjects x outfmt 0/6, bad forms) and masking combos (item 3); u1..u7 / t1..t6 / fz f1..f70,g1..g40 title subjects (item 1; 110 files); sh, tw/pwidth, huge (3200 subjects), long-path-with-spaces (items 9, 11); threshold / window / matrix / evalue sweeps (items 4, 5); fixture matrix of queries x subjects (e2e_*, punct_hits, pwidth/twidth) x outfmt 0/6/7; db_gencode 1..33.
 - REJ-OTHER (50): max_target_seqs -1/0/abc/1.5 (32), threshold inf/nan/-1 (6): both fail, NCBI exit 1 vs LOSAT exit 2, different text = approved exception 1 (item 2 and 4 of the report say so); db_gencode 7, 8, 17-20 (12): both reject the unassigned code with different text, same class.
 - LOSAT-REJECTS (122): max_hsps for TBLASTN, blastp -lcase_masking / -soft_masking, -evalue inf (D12), TBLASTN non-default matrix/threshold/window rows, query/subject deflines with control characters / non-ASCII / empty / leading white space, invalid UTF-8 FASTA. All explicit "not supported by LOSAT" messages with exit 1.
 - ACCEPTED (32): db_gencode != 1 on TBLASTN -subject (30 rows: gencode 2-6, 9, 10, 12-16, 21-31, 33 score differences; 32 accepted by LOSAT, rejected by NCBI = PD-TLOSAN-LOCAL-GENCODE-32); report item 7 and AGENTS.md list this as an approved exception.
   The other 2 DIFF rows are B2_fx_punct_hits_query_punct_hits_subject_0 and B2_fx_e2e_protein_query_punct_hits_subject_0 (tblastn on punct_hits_subject.fna, outfmt 0): NCBI hangs/crashes (killed by timeout 120, exit 139; x_CleanAndCompress on `, ,` / `;~ ;` titles), LOSAT prints a result. This is approved exception 2 of PD-LOSAT-NCBI-DEFECTS (see make_inputs.py); with the stand-in subject punct_hits_standin.fna outfmt 0/6/7 are SAME, and outfmt 6/7 on the real subject are SAME. ACCEPTED.
Item 9 stdin query (`-query -`) and `-out file`, blastp and tblastn x outfmt 0/6/7: 6/6 SAME (io_out.txt).
Item 10 (60000-aa query): LOSAT-REJECTS for blastp and tblastn as expected (see RP-4).

## Summary table
| Finding | Verdict | Evidence |
|---|---|---|
| RP-1 O-query warning order | SAME (fixed) | 57 repro pairs + 27 original ov* runs (blastp/tblastn, outfmt 0/6/7): byte-identical incl. 2>&1 |
| RP-2 comment line range | n/a (comment only) | listed as fixed in ROUND1.md |
| RP-3 comment line range | n/a (comment only) | listed as fixed in ROUND1.md |
| RP-4 60000-aa TBLASTN query | LOSAT-REJECTS (default) and SAME (CHUNK_SIZE=100000) | "a query batch of 60000 residues, which NCBI BLAST+ searches in 3 query chunks ... is not supported by LOSAT's TBLASTN"; CHUNK_SIZE=100000: 9/9 SAME (sbig4 bit score 10412 in NCBI and LOSAT); thresholds 19800/39800 confirmed |
| RP-5 query title warning placement | SAME (fixed) | qh11, qw1, qw2 x 2 programs x outfmt 0/6/7 with 2>&1: SAME |
| Checked OK 1-9, 11 | SAME / LOSAT-REJECTS / ACCEPTED | 2183 pairs replayed: 1831 SAME, 264 LOSAT-REJECTS, 56 REJ-OTHER (approved exception 1), 32 ACCEPTED (30 db_gencode, 2 punct-title NCBI crash), 0 unexplained DIFF (see below) |
| Checked OK 10 (60000-aa blastp) | LOSAT-REJECTS | "a query batch of 60000 residues, which NCBI BLAST+ searches in 6 query chunks ... is not supported by LOSAT's BLASTP"; SAME with CHUNK_SIZE=100000 |

## DIFF list
None unexplained. The only non-SAME, non-rejection rows are the 32 ACCEPTED ones above (30 db_gencode rows, 2 punct_hits_subject.fna outfmt 0 rows where NCBI hangs/crashes).
Caveats: (1) batch 2 commands were regenerated from tag names, not the original command lines, so coverage of the original ~1,000 pairs is approximate; (2) NCBI timeout was 120 s per run (1200 s for RP-4); (3) vperf.lock never existed during the runs.

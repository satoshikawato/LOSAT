# Angle (d): BLASTP and TBLASTN reports, and inputs of all three programs

Read `COMMON.md` first. Work dir: `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/d/`.

1. Re-run round 1: every reproducing command and harness in `src/docs/evidence/losat_web_e2e/audit/round1/reports.md` (RP-1..RP-5, ~1000 commands) and `audit/round1/inputs.md` (IN-1..IN-14, ~2100 inputs × outfmt 0/6/7 for blastp, tblastn, tblastx). S08+a's copies: `/home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/reports/` (`cmp2.sh`, `drive.py`, `cmds*.json`) and `.../audit-rerun/inputs/` (`c.sh`, `job*.sh`, `cases*.tsv`, `gen*.py`, input dirs `in`, `le`, `nt`, `pc`, `ws`, `fz*`). Copy them into `d/`, point LOSAT at `$FINAL`, and write only under `d/`. Give a verdict per round 1 finding and totals per class.
2. New, around S08+a changes that reach reports and inputs:
   - Query splitting in reports (BLASTP chunk 10000, TBLASTN 20000; `s08pa/NOTES.md` "RP-4"): outfmt 0 and 7 of batches with long queries (above ~9,800 residues for blastp, ~19,800 for tblastn) mixed with short ones; queries with title warnings (50-letter protein titles), O residues, invalid residues, lower case, empty records in the same batch as a long query; `BATCH_SIZE` that moves the long query to another batch; `2>&1` ordering of warnings and reports. Long inputs: `audit/round1/rp4/*.gz` (gunzip into `d/`), `LOSAT/tests/fasta/outfmt0/e2e_split_query.faa`, `e2e_split_two_query.faa`.
   - Long single-line queries (30000+ residues on one line), CRLF line ends in long queries, long queries with `-lcase_masking` (tblastn) and `-seg yes`.
   - TBLASTN outfmt 0 of the TN-5/TN-1/TN-2 fixtures' option sets on other inputs (`-matrix BLOSUM45 -word_size 2 -comp_based_stats 0 -evalue 1000`; `-comp_based_stats 0 -xdrop_gap_final 1e10 -evalue 1000`; `-comp_based_stats 0 -lcase_masking -seg no -max_target_seqs 2`): Sbjct rows, descriptions, footers (Lambda/K/H lines, effective search space).
3. Report per COMMON.md.

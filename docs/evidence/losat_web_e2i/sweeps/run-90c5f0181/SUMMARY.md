# E2i sweeps of `-task dc-megablast` and `-task blastn-short` (LOSAT-90c5f0181 vs NCBI BLAST+ 2.17.0)

LOSAT: /home/kawato/.cache/losat-web-gui-target/sd/bin/LOSAT-90c5f0181. NCBI: /home/kawato/micromamba/bin.
Scripts: ../../sweeps.sh (runs all), ../../check_inputs.py (new), and the `--tasks` option of
docs/evidence/losat_web_e2c/scoring_sweep.py and word_size_sweep.py. Nothing under LOSAT/ was touched.

| sweep (file) | cases | not same / same-error / expected |
|---|---|---|
| scoring-sweep-fmt0.tsv (dc-megablast, blastn-short) | 880 | 0 (260 same, 620 same-error) |
| scoring-sweep-fmt6.tsv | 880 | 0 (260 same, 620 same-error) |
| scoring-sweep-fmt7.tsv | 880 | 0 (260 same, 620 same-error) |
| word-size-sweep.tsv | 256 | 0 differing |
| batch-sweep-task-dc-megablast.tsv | 120 | 0 differing |
| batch-sweep-task-blastn-short.tsv | 120 | 0 differing (e.g. 8321 output lines on a case at e-value 1000, equal) |
| batch-sweep-task-blastn-short--evalue-1e-5.tsv | 120 | 0 differing |
| slice-sweep-task-dc-megablast.tsv (pool viral) | 120 | 0 differing |
| slice-sweep-task-blastn-short.tsv | 120 | 0 differing |
| slice-sweep-task-blastn-short--evalue-1e-5.tsv | 120 | 0 differing |
| check-inputs.tsv | 774 (237 `.dc`, 237 `.short`) | 23 differ from the E2g expectation, none is a LOSAT defect except possibly 2 timeouts (see below) |
| default-scoring-sweep-fmt6.tsv (unchanged invocation) | 880 | 0 (300 same, 580 same-error), exit 0 |
| default-word-size-sweep.tsv (unchanged invocation) | 256 | 0 differing, exit 0 |

The word-size sweep with dc-megablast covers NCBI's error for word sizes other than 11 and 12
("Invalid discontiguous template parameters: word size must be either 11 or 12", exit 1): stdout, stderr and
exit status equal.

## check-inputs.tsv: the 23 rows whose result is not the expectation

Result counts: same 370, losat-rejects 146, same-error 212, arg-error 37, exception-2 7, timeout 2 (the E2g
classes; the 23 below are those with result != expect).

1. 14 cases, expectation `same` (default task), result `same-error`: with the other task NCBI itself fails and
   LOSAT fails with identical stdout, stderr and exit status (no differing bytes):
   audit.reward_65538.{dc,short}, audit.large_divisible.{dc,short} ("Gap existence and extension values 5 and 2
   are not supported for substitution scores ..."), audit4.hex.word_size_16.dc and
   audit14.iupac_seed.word_size_24.dc ("word size must be either 11 or 12"), audit5.bare_0x.gaps.{dc,short}
   ("Greedy extension must be used if gap existence and extension options are zero"),
   audit7.bit_score_99.fmt{0,6,7}.{dc,short} (gap values 5 and 2 with scores 2 and -5).
2. 6 cases, expectation `losat-rejects`, result `same` (LOSAT no longer rejects; byte-identical to NCBI):
   audit.megablast_gap_max.{dc,short}, audit2.greedy_gap_limit.{dc,short}, and the two E2g cases
   audit6.task.dc_megablast and audit6.task.blastn_short (their E2g expectation predates the two tasks).
3. 1 case, `exception-2`: audit12.crash_title_no_hit.fmt0.short (query multi_query.fasta, subject with the defline
   ">, ,", `-task blastn-short`). NCBI dies with SIGSEGV (shell exit 139) after writing 17 bytes of stdout
   (`BLASTN 2.17.0+\n\n\n`, bytes 424c4153544e20322e31372e302b0a0a0a), no stderr; LOSAT exits 0 and its stdout starts with
   the same 17 bytes followed by `Reference: Stephen F. Altschul, Thomas L. M...` (11149 bytes). This is the
   known NCBI crash for a punctuation-only title with hits (hits exist with blastn-short at e-value 1000; the
   default-task case has none and is `same`), i.e. exception 2 of PD-LOSAT-NCBI-DEFECTS, the same class as
   audit10.title_10/11.
4. 2 cases, `timeout` (60 s limit of the case runner): audit15.small_evalue.1e-300.short and
   audit15.small_evalue.1e-297.short (AP027131 x AP027132, full genomes, `-evalue 1e-300`/`1e-297`, outfmt 6,
   blastn-short; 1e-295.short is `same`). Not classified: NCBI blastn-short on these needs about 14.6 GB RAM
   (killed by the orchestrator, which asked to drop the combination). A LOSAT-only run with an 8 GB address
   space limit aborted (SIGABRT, allocation failure) after 83 s and 4.7 GB RSS, so LOSAT is also memory heavy
   here; no output comparison exists for these two. `check_inputs.py` now leaves out the `.short` variants of
   audit15.small_evalue.* (docstring says so); check-inputs.tsv is the run made before that change (all 774
   cases, including those two).

## Notes
- All sweeps except check-inputs had 0 differences. The batch and slice sweeps ran six at a time with 2 jobs
  each (12 jobs); after the coordinator's request for at most 8 jobs `sweeps.sh` was changed to run three at a
  time (6 jobs max, same commands and results); that version was not rerun.
- A rerun of check_inputs.py with the audit15 `.short` skip failed repeatedly with "Input/output error" from the
  /mnt/c file system (the shared drive was failing; no result file of this directory contains such an error);
  the complete first run is kept as check-inputs.tsv.

> S08+a: report of the TBLASTN audit-rerun agent, run with `LOSAT-head` = `6fa06cb6a` (before the TN-1 residual fix `68154c73e`). Re-running the two-record repro prints 0.40 (not 0.39) with `LOSAT-head` and with the pre-change binary; the final binary prints NCBI's 0.36 (NOTES.md, "TN-1 の残り").

# TBLASTN audit re-run (LOSAT-head vs NCBI tblastn 2.17.0)

Work dir: /home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/tblastn/ (runner c.sh, per-run files r/<label>.{n,l}.{out,err,rc}, logs log*.txt)
$D=/home/kawato/.cache/losat-web-gui-target/s08p/audit ; F=$D/src/LOSAT/tests/fasta/outfmt0 ; generated inputs copied to work dir in/.

## TN-1
All 12 repro variants SAME (stdout+stderr+rc identical, non-empty output): -xdrop_gap 5 -xdrop_gap_final 5 cbs 0 (e2e_many_subject); -xdrop_gap_final 1e10/1e9/4e9/1e300/+inf cbs 0; gen_q/gen_s -xdrop_gap_final 10 and 25 cbs 0; -sum_stats false variant; sweep find `-xdrop_gap_final -5 -comp_based_stats F` on gen_q; cbs 2 -xdrop_gap_final 1e9.

## TN-2
LOSAT-head SAME for all: gen_q_lc/gen_s `-comp_based_stats 0 -lcase_masking -seg no -max_target_seqs {1,5,10,100,500} -outfmt 6`, and BLOSUM45 variant (`-matrix BLOSUM45 -word_size 2 ... -max_target_seqs 2` and without -max_target_seqs).
LOSAT-base (old behaviour, reproduced): DIFF, e.g. `-max_target_seqs 1` row 1 NCBI `9.33e-177` vs base `9.00e-177`; BLOSUM45 `4.31e-157` vs base `4.21e-157`.

## TN-3
LOSAT-head SAME (stdout outfmt 0) for `-xdrop_gap_final 100 -seg no`, `-xdrop_gap_final 100 -soft_masking true -sum_stats false`, and defaults; NCBI output does contain the two all-gap `Sbjct       ----...` rows (lines 2201, 2205), LOSAT prints them identically.
Note: LOSAT-base is byte-identical to LOSAT-head on this command, i.e. the base binary already had this fix (not reproduced as a difference with LOSAT-base).

## TN-4 (time, peak RSS via /usr/bin/time; timeout 120; output compared by cmp; all outputs byte-identical to NCBI, rc 0)
Subject e2e_tblastn_subject.fna / e2e_many_subject.fna, query e2e_protein_query.faa, `-outfmt 6`.
| args | NCBI tbs | LOSAT-head tbs | NCBI many | LOSAT-head many |
|---|---|---|---|---|
| cbs 0 -xdrop_gap_final 1e7 | 0.32s 70MB | 0.23s 8.7MB | 0.47s 75MB | 4.13s 9.1MB |
| cbs 0 -xdrop_gap_final 1e8 | 0.33s 71MB | 0.23s 8.4MB | 0.53s 75MB | 4.26s 9.3MB |
| cbs 0 -xdrop_gap_final 5e8 | 0.34s 71MB | 0.24s 8.5MB | 0.51s 88MB | 4.23s 9.1MB |
| -xdrop_gap 1e7 | 0.58s 70MB | 0.72s 9.2MB | 0.63s 74MB | 4.27s 9.3MB |
| -xdrop_gap 1e8 | 0.66s 221MB | 0.70s 9.1MB | 0.62s 218MB | 3.93s 9.6MB |
| -xdrop_gap 5e8 | 0.60s 70MB | 0.64s 9.1MB | 0.74s 76MB | 4.09s 9.4MB |
Reference, LOSAT-base (old): tbs cbs0 final 1e8 55.2s 4.07GB; tbs -xdrop_gap 1e8 4.34s 4.06GB; many cbs0 final 1e7 39.4s 424MB.
Verdict: SAME output; the huge-X-drop time/memory blow-up is gone (peak RSS ~9 MB for every case).
Side observation (not a TN-4 regression): LOSAT-head takes ~4 s on e2e_many_subject.fna with default x-drop (cbs 0: 3.64s; cbs 2: 4.10s) vs NCBI 0.18-0.39 s; LOSAT-base is the same (3.91s), so the ~4 s on that subject is a pre-existing constant cost unrelated to X-drop.

## TN-5
LOSAT-head SAME (stdout, stderr, rc) for BLOSUM45 `-word_size 2 -comp_based_stats 0 -evalue X -outfmt 6` on e2e_protein_query x e2e_tblastn_subject: X = 2000, 1e4 (29,900 rows), 5000, 1e5, 1e10 (69,454 rows); BLOSUM62 default profile `-evalue 1e10` SAME.
`-evalue +inf` (BLOSUM45): LOSAT-REJECTS (see TN-6; NCBI rc 0 on this subject).
LOSAT-base (old behaviour, reproduced): DIFF, -evalue 1e4: 4 of 29,900 rows differ (first: NCBI `BDT62620.1 LvMJNV_160001_165000 50.000 6 3 0 39 44 4712 4695 3879 10.3`, base `... 4710 4693 ...`); -evalue 1e5: 12 of 69,448.

## TN-6 (-evalue +inf / 1e999)
LOSAT-REJECTS in every case, rc 1 and stdout empty: `Error: an infinite -evalue (inf), with which NCBI BLAST+'s tblastn crashes on some inputs, is not supported by LOSAT's TBLASTN`.
- e2e_many_subject `-evalue +inf`, `1e999`, `+inf -sum_stats false`, `+inf -seg no`: NCBI rc 139 (SIGSEGV), LOSAT rc 1 (explicit rejection, ROUND1 decision D12).
- `+inf` with `-comp_based_stats 0`, with `-max_target_seqs 5` (e2e_many_subject) and on e2e_tblastn_subject: NCBI completes rc 0 (11878 / 149 / ... rows), LOSAT still rejects (by design: LOSAT cannot know beforehand whether NCBI would crash; D12).
- `-evalue 1e300` and `-evalue 1e308` (e2e_many_subject, 7221 rows): SAME.

## TN-7 (-dryrun)
LOSAT-REJECTS: rc 2, `error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's TBLASTN` (also with `-matrix FOO`). NCBI rc 0 with empty output. Explicit rejection (D9) in place of the old generic unknown-option error.

## TN-8 (ACCEPTED per ROUND1)
- `-num_threads 100000`: NCBI rc 0 (warning: threads reduced to 32 ...), LOSAT rc 1 `Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT` (wording without the program name, accepted).
- `NCBI_CONFIG__BLAST__BATCH_SIZE=3`: NCBI rc 0, LOSAT rc 1 `Error: the environment variable NCBI_CONFIG__BLAST__BATCH_SIZE, which sets an entry of NCBI BLAST+'s registry that LOSAT does not know to change no output (it accepts only the entries listed for registry files), is not supported by LOSAT`.

## TN-9 (ACCEPTED per ROUND1)
- `-evalue 0 -outfmt 5`: NCBI rc 1 `BLAST query/options error: expect value or cutoff score must be greater than zero`; LOSAT rc 1 `Error: output format 5 is not supported by LOSAT's TBLASTN` (both rc 1, explicit failure).
- subject `ACGT1ACGT...` with `-matrix FOO` / `-evalue 0`: NCBI rc 1 (FASTA warning + options error); LOSAT rc 1 `Error: subject record 1 (s) has '1' at residue 5, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not supported by LOSAT's TBLASTN ...`. (With no bad option NCBI rc 0 with a warning, LOSAT rc 1: same explicit rejection.)
- control `-evalue 0 -outfmt 6`: SAME.

## TN-10, TN-11 (comments, checked in head-src/LOSAT/src, NCBI tree $D/../ncbi)
TN-10: protein_options.rs:390 now cites blast_options.c:910-936 (line 910 is the quoted `if`), :492 cites 1518-1521 (matches), stage_e_report.rs:623 cites blast_format.cpp:2267 (line 2267 = `options.GetMatrixName()`), args.rs:373 now cites both 2874-2886 and 2960-2977 with the explanation and the snippet from the second range: FIXED.
TN-11: args.rs: the composition-mode comment now carries a snippet (blast_traceback.c:1486-1499) and sits above the `composition_mode2` match; the compressed-lookup comment no longer cites blast_aalookup.c as if it were the rejection (now says it is a LOSAT limit); longest_intron / num_threads / blast_format.cpp:1443-1451 comments now carry snippets: FIXED.
Note: running the audit's heuristic verify_refs.py over the head files reports 123 "MISS" of 201 references, but it is the same noisy heuristic (it also flags unfenced/paraphrase citations in files outside TN-10/11); no TN-10/TN-11 site is among the cases examined by hand above.

## Harness re-runs (copies in scripts/ with the binary repointed to LOSAT-head; outputs in work dir)
- pairs.py (650 ordered pairs of 26 NCBI-rejected options): 350 SAME, 194 PARSEERR (NCBI USAGE exit 1 vs LOSAT exit 2, approved), 106 REJECT (LOSAT "not supported by LOSAT's TBLASTN"), 0 DIFF. Per-pair classes identical to the old audit's pairs.out line by line.
- pairs2.py (18 LOSAT-limit options x 14 NCBI-rejected options, both orders, 512 runs): 392 SAME, 38 PARSEERR, 74 REJECT, 0 DIFF. Per-pair classes identical to old pairs2.out except `mtbig q_nothing` where NCBI now segfaulted (rc -11) instead of timing out (rc -9) in this run; LOSAT rejects both ways.
- grid.py (160 runs: lcase/seg/max_target_seqs/cbs/soft_masking; the TN-2 grid): 160 SAME, 0 DIFF (old audit: DIFF for lcase + seg no + soft_masking false + cbs 0 + mts <= 10).

## NEW DIFF found by sweep.py (seed 1) - NOT one of TN-1..TN-11 as reproduced; same family as TN-1 (cbs 0, e-value stat_length with a huge final X-drop) but only with a large -evalue
Sweep hit: gen_q_lc x gen_s `-comp_based_stats F -xdrop_gap_final -5 -xdrop_gap 1e10 -evalue 1e5` (outfmt 0): 51 lines differ, e.g. `gs37 len=2643 ... 19.2    9.0` (NCBI) vs `... 9.9` (LOSAT), `Expect = 21` vs `23`, `66` vs `70`.
Isolated (outfmt 6, `-comp_based_stats 0`): effective final X-drop >= 1e9 bits (`-xdrop_gap_final 1e9/1e10`, or `-xdrop_gap 1e9/1e10` with any final < prelim, or final -5/0) AND `-evalue` >= 1000 on the generated set (>= 1e5 on the fixtures). Output is identical for -xdrop_gap_final <= 1e8, for -evalue <= 100, for cbs 2, and not fixed by -sum_stats false. LOSAT-base gives byte-identical output to LOSAT-head on these, so the session's TN-1 fix does not cover it and it is not a regression.
Fixture-only repro (F=$D/src/LOSAT/tests/fasta/outfmt0):
  `tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_many_subject.fna -comp_based_stats 0 -xdrop_gap_final 1e10 -evalue 1e5 -outfmt 6` -> NCBI 12148 rows, LOSAT 12139: NCBI has 9 more rows (first `BDT62567.1 n159 57.143 7 3 0 363 369 423 403 438 15.8`, also E=2344, 6203 rows); remaining rows differ in E-value. With `-evalue 1e4` and below: SAME.
Two-record repro (work dir r/min_q.faa = gq4, r/min_s.fna = gs37): `-comp_based_stats 0 -xdrop_gap_final 1e10 -evalue 1000 -outfmt 6` -> NCBI row 1 `... 0.36 19.2`, LOSAT `... 0.39 19.2` (ratio ~1.08 like TN-1's constant per-subject scale).
Logs: r/m1.*, r/m_fixt.*, r/g_gen_q_*, log_d.txt.

- sweep.py / sweep_b..e.py random option sweeps (1-6 random options of 11 per run, 60 s timeout per run): 10 invocations, 1690 runs: 1636 SAME, 51 REJ (LOSAT "not supported by LOSAT's TBLASTN"), 0 PARSE, 3 DIFF.
  Seeds: sweep.py 1/2/3 n=150/160/160, sweep.py 4 n=120, sweep_b/c/d/e.py seed 5 n=200 each, plus two extra runs sweep.py seeds 11 and 12 n=150 each.
  The 3 DIFF: (a) sweep.py seed 1 and (b) seed 4 are the NEW DIFF above (huge final X-drop + large -evalue, cbs 0): `-outfmt 0 -comp_based_stats F -xdrop_gap_final -5 -xdrop_gap 1e10 -evalue 1e5` (gen_q_lc) and `-outfmt 0 -comp_based_stats x -evalue 1e10 -xdrop_gap_final +inf -sum_stats true` (gen_q_lc; reproduced as outfmt 6 `-comp_based_stats 0 -evalue 1e10 -xdrop_gap_final +inf`: DIFF). (c) sweep_c.py seed 5: `-evalue 1e5 -xdrop_gap 100 -soft_masking false -outfmt 7 -db_gencode 1 -seg '10 1.0 1.0...'` with LOSAT rc -9 (60 s timeout) -> not reproducible: re-run alone LOSAT-head 2.3 s, output byte-identical to NCBI (system load average 17-30 from other sessions at that moment): classed SAME (flake).
- cmp2.sh battery replay (battery1-6.out + bool.out of the old audit, 688 argv replayed through the repointed cmp2.sh as scripts/replay.py; args with embedded spaces had lost their quoting in the old logs, so a few -seg/-outfmt/empty-value cases are re-split heuristically and change class SAME<->PARSEERR for that reason only): new classes 485 SAME, 65 PARSEERR (NCBI USAGE exit 1 vs LOSAT exit 2), 133 REJECT-OK, 5 "DIFF!!!".
  The 5: `-db_gencode 4` (output differs: approved PD for local-subject db_gencode, AGENTS.md), `-db_gencode 32` (NCBI USAGE rc 1, LOSAT accepts: approved PD-TLOSAN-LOCAL-GENCODE-32), `-num_threads 2`, `-num_threads 4`, `-num_threads 0x2` (NCBI stderr warning "'num_threads' is currently ignored when 'subject' is specified." not printed by LOSAT: approved PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, AGENTS.md) -> all ACCEPTED, identical to old audit's battery4 results (the old report's "-num_threads 1..1000 identical" holds for stdout/rc only).
  Of the old non-SAME classes: old "-evalue +inf / 1e309" lines (SAME or NCBI-segfault DIFF) are now REJECT-OK (D12, intended).
  The huge-X-drop subset (-xdrop_gap/_final >= 1e7 spellings, 37 argv, run serially after the sweeps): 31 SAME, 4 PARSEERR (`inf`/`Inf`/`INF`/`infinity`), 2 REJECT-OK, 0 DIFF.
- xf_run.sh / t4.sh (cbs 0, -xdrop_gap_final 1e10..1e5 x e2e_tblastn_subject / e2e_many_subject, default -evalue): hit the 40 minute cap after the first 11 combos because of system stalls caused by other sessions (the 5e8/2e8/1e8 combos reported `DIFF(n=124 ...)` = NCBI itself hit the 60 s timeout after 300-490 s wall time; reproduced alone NCBI 0.26 s, so not a LOSAT result). Re-run cleanly afterwards (xf_rerun.out): all 12 combos (5e8, 2e8, 1e8, 1e7, 1e6, 1e5 x 2 subjects) SAME, LOSAT-head 0.3-0.6 s on e2e_tblastn_subject and 2.3-5.6 s on e2e_many_subject; together with the first run (1e10, 1e9 SAME): 16/16 SAME.
- t3.sh/cmp.sh/verify_refs.py: t3/cmp are drivers, used via c.sh equivalents above; verify_refs.py was run over the head sources (heuristic, see TN-10/11).

## Summary
| Finding | Verdict | Evidence |
|---|---|---|
| TN-1 | SAME (fixed) | all 12 repros byte-identical (xdrop 5/5, final 1e9..+inf, gen set final 10, sweep -5 F); caveat: NEW DIFF with -evalue >= 1e3..1e5 (below) |
| TN-2 | SAME (fixed) | gen_q_lc mts 1/5/10/100/500 + BLOSUM45 variant + 160-run grid: all SAME; base still DIFF 9.33e-177 vs 9.00e-177 |
| TN-3 | SAME | outfmt 0 all-gap Sbjct rows identical (base binary already identical) |
| TN-4 | SAME output, fixed time/memory | xdrop 1e7-5e8: LOSAT-head 0.2-0.7 s / 9 MB (tbs), 4 s / 9 MB (many, same 4 s as base with defaults); base 55 s 4.07 GB |
| TN-5 | SAME (fixed) | BLOSUM45 -evalue 2000..1e10 byte-identical (up to 69,454 rows); base 4/29,900 rows differ at 1e4 |
| TN-6 | LOSAT-REJECTS | `an infinite -evalue (inf), with which NCBI BLAST+'s tblastn crashes on some inputs, is not supported by LOSAT's TBLASTN` rc 1 (NCBI segfault rc 139 / rc 0 on some inputs); 1e300/1e308 SAME |
| TN-7 | LOSAT-REJECTS | `error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's TBLASTN` rc 2 (NCBI rc 0, empty) |
| TN-8 | ACCEPTED | thread cap / NCBI_CONFIG__* wording unchanged ("... which is not supported by LOSAT", rc 1; NCBI rc 0) |
| TN-9 | ACCEPTED | `-evalue 0 -outfmt 5`: both rc 1, different text; bad-residue subject: LOSAT rc 1 explicit rejection |
| TN-10 | fixed | cited lines now 910-936, 1518-1521, 2267, args.rs cites both ranges (checked against NCBI source) |
| TN-11 | fixed | the comments carry snippets and sit above the code they describe |
| harness | 0 unexplained DIFF except NEW DIFF | pairs 650: 350 SAME/194 PARSEERR/106 REJECT; pairs2 512: 392/38/74; grid 160 SAME; sweeps 1690: 1636 SAME/51 REJ/3 DIFF; battery replay 688: 485 SAME/65 PARSEERR/133 REJECT-OK/5 accepted exceptions |

## DIFF list
1. NEW (not in the first audit, pre-existing in LOSAT-base, not fixed by the TN-1 fix): `-comp_based_stats 0` with an effective final X-drop >= ~1e9 bits AND a large -evalue (>= 1e3 generated set, >= 1e5 on fixtures).
   `F=$D/src/LOSAT/tests/fasta/outfmt0; tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_many_subject.fna -comp_based_stats 0 -xdrop_gap_final 1e10 -evalue 1e5 -outfmt 6`: NCBI 12148 rows, LOSAT 12139. First differing line: diff `6795,6803d6794` = NCBI has 9 rows LOSAT lacks (`BDT62567.1 n159 57.143 7 3 0 363 369 423 403 438 15.8` ...). Two-record repro (gq4 x gs37, r/min_q.faa r/min_s.fna): `-comp_based_stats 0 -xdrop_gap_final 1e10 -evalue 1000 -outfmt 6` row 1 e-value NCBI 0.36 vs LOSAT 0.39. Other triggers: `-xdrop_gap 1e9|1e10` with final < prelim / -5 / 0; `-xdrop_gap_final +inf -evalue 1e10`.

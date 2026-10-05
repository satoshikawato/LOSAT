# Replay of the audit3 scripts with FIX3 (LOSAT-wip4) against stored FIX2 (LOSAT-wip3) results

Copy: audit3/replay/scratch/ (cp -a of scratch/; scratch/ untouched). h.py there has FIX=LOSAT-wip4.
Replay script: replay/scratch/replay.py (runs only FIX3; same argv, cwd = the copied scratch, same env filter as h.py;
3 worker threads, same as original order via ThreadPoolExecutor.map). Analysis: analyze2.py. Per-script FIX3 results: replay_*.pkl.
Compared (stdout, stderr, exit status) with the stored FIX2 result. Stored FIX2 results used: t1/t2/t6 w; t3/t4/t5 w; t7 eq and sep results;
t8, t9 direct (w) and U+FFFD twin (t); twin.pkl 4th element (twin of t1/t2/t6 non-path cases). a1.py / a6.py are analysis-only (no runs).
HEAD, FINAL, NCBI not re-run.

## Counts (each "case" = one command line run with FIX3)
| script | cases | identical | different | F-1 diffs | other diffs |
|---|---|---|---|---|---|
| t1 | 8736 | 8736 | 0 | 0 | 0 |
| t2 | 3125 | 3123 | 2 | 0 | 2 (harness file-state artefact, below) |
| t3 | 1998 | 1998 | 0 | 0 | 0 |
| t4 | 303 | 303 | 0 | 0 | 0 |
| t5 | 3500 | 3500 | 0 | 0 | 0 |
| t6 | 12400 | 12400 | 0 | 0 | 0 |
| t7 (2040 eq + 2040 sep) | 4080 | 4068 | 12 | 8 | 4 (harness file-state artefact, below) |
| t8 (11169 direct + 11169 twin) | 22338 | 21778 | 560 | 560 | 0 |
| t9 (6000 direct + 6000 twin) | 12000 | 11874 | 126 | 126 | 0 |
| twin.pkl | 22592 | 22592 | 0 | 0 | 0 |
| total | 91072 | 90372 | 700 | 694 | 6 |

(First pass, raw. The 6 "other" diffs are explained and removed below.)

## (a) F-1 cases: 694 diffs
All 694 stored FIX2 crashes (the only crashes in the stored data; every one is rc -6 = SIGABRT after the panic at
value_parsers.rs:273 "byte index 4 is not a char boundary") now give, in FIX3, exit 2 with
`error: invalid value '<value>' for '-<opt> <...>': expected a decimal number (other forms, which NCBI BLAST+ may read, are not supported by LOSAT's <PROG>)` -
clap parser error, text contains "not supported by LOSAT". 0 of the 694 crash in FIX3.
In every one of the 694, FIX3's result for the original argv is byte-identical (stdout, stderr, rc) to FIX3's result for the U+FFFD form of the same argv
(extra FIX3 run of the twin for each, 694 additional runs).
Breakdown: t7 8 (4 eq, 4 sep; non-UTF-8 argv); t8 260 direct non-UTF-8 + 20 direct valid-UTF-8 (value with real multibyte char, e.g. -xx€-type, pre-existing panic for UTF-8) + 280 twin (valid UTF-8 U+FFFD form);
t9 63 direct non-UTF-8 + 63 twin. The 260 non-UTF-8 t8 cases equal the "260 of 11169" in R.md.

## (b) other diffs: 6 raw, all caused by files in cwd whose content depends on earlier -out cases (not by the binary)
t2 idx 3039, 3045 (blastx; -query nq300<ff>.fa / -subject nq300<ff>.fa): the stored FIX2 run saw nq300<ff>.fa already overwritten by an earlier
run of the case `blastx ... -out nq300<ff>.fa` (idx 3051; the only case that changes that file; blastx output, 3964 bytes), so FIX2 stored a CFastaReader error (rc 1).
My copy started with the file clean (311 bytes), so FIX3 reads it as valid FASTA (rc 0). Evidence: (1) with the file in the overwritten state (the state the replay leaves) 4 further
t2 replays give 0 diffs vs stored; (2) control with a clean file in a separate dir: LOSAT-wip3 and wip4 give identical results on both cases (rc 0, same stdout).
t7 (blastn -query 1, -query=1 ; blastp -subject A, -subject=A): stored FIX2 saw "File is not accessible" because files "1" and "A" did not exist yet (t7's own -out 1 / -out A create them);
my copy already contains those files from the original run. Evidence: in a copy with the files created by t7/t8/t9 removed (replay/chk_t7), t7 replay gives 8 diffs, all F-1, 0 others.
No other differences.

## Crash check
FIX3 crashes (rc<0, 101/134/139, or "panicked" in stderr) over all 91072 + 694 runs: 0.
Exit >= 128: 800 cases with exit 255 `Error: Formatting choice is out of range` (-outfmt 99 etc.; NCBI-style error exit, no panic); all identical to the stored FIX2 results.

## Conclusion: supported
Only F-1 cases differ (694), none crashes in FIX3, and FIX3 equals its own U+FFFD-form result in all of them; the 6 other raw diffs are
cwd-file-state artefacts shown to vanish when the stored file state is reproduced (and wip3 == wip4 on a clean file). Total replayed command lines: 91072 (plus 694 twin checks, ~4200 extra runs for state diagnosis/controls).
Caveat: strict reading of "any other difference" would count the 6 raw diffs; they are listed above with evidence.

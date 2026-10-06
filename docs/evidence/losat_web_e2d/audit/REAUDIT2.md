# R3 audit of fix2 (fc9970539)
started Mon Oct  5 16:25:51 JST 2026
Interim (t1-t6 + twin check done, see below). FINDING so far: `blastn -query Q -subject S -evalue -$'\xc3\xc3'` -> FIX2 panics (rc 134, value_parsers.rs:273 "byte index 4 is not a char boundary"); LOSAT-final: explicit rejection. The panic is pre-existing for valid UTF-8 (HEAD/FINAL: `-evalue -xx€` panics), FIX2 makes it reachable from non-UTF-8 argv via U+FFFD substitution. Fuzzing for the full extent next (t8).

## Final (FIX2 = LOSAT-wip3, fc9970539)
Harness: scratch/h.py (HEAD, FINAL, FIX2, NCBI 2.17.0). t1 8736, t2 3125, t3 1998, t4 303, t5 3500, t6 12400 (focused matrix: help/parser/post-parse/typed/repeat/path), t7 2040 (eq vs sep fuzz), t8 11169 (numeric/odd-value fuzz, FIX2 vs lossy UTF-8 twin vs FINAL), t9 6000 (random argv, FIX2 vs twin), twin.pkl 22592 derived twin checks. 49271 command lines, ~187k process runs.
Results:
- Valid UTF-8 argv (t1 416, t2, t3 1998, t5 3500, t4 root/blastx): HEAD == FINAL == FIX2, 0 differences. blastx/root identical for non-UTF-8 too.
- Twin check (every non-UTF-8 argv of the four programs, non-path, 22592): FIX2 output == the output of the same command with U+FFFD in place of the bad bytes (parse-stage: help / parser error / -h,-version unsupported) in 17,000+ cases; otherwise (5,000 cases) the twin passes the parser and FIX2 gives the explicit "value of -X is not UTF-8 ... not supported by LOSAT's X" rejection. No third outcome.
- R-1 resolved (-help wins when first or when bad value's parse succeeds; same as UTF-8 `-evalue abc -help` behaviour when typed value fails parse first). R-2 resolved (inline == separate for 224 pairs, file-name rejection rc 1).
- FINDING F-1: panic (rc 134) for numeric options -evalue/-perc_identity/-threshold/-xdrop_gap/-xdrop_gap_final with a non-UTF-8 value that, after U+FFFD substitution, has a sign/"in"/"na" prefix and a multibyte char at byte 4 (e.g. -evalue $'-\xff\xff'): value_parsers.rs:273. Pre-existing for valid UTF-8 (HEAD/FINAL `-evalue -xx€` panics) but FINAL gave explicit rejection, HEAD a parser error. 260 of 11169 t8 cases.
Conclusion: unsupported (claim 1 has a counterexample).

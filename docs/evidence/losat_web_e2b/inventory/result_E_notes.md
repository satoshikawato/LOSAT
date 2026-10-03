# Range E result notes: what the TBLASTX report receives from the search (final code 0d533ba76)

Result TSV: `result_E.tsv` (60 inventory rows + extra rows `X1`, `X2`, ... for findings that no inventory row covers).
Scratch: `rscratch_E/` (inputs `rnd/`, `edge/`, `bs/`; outputs `out/`, `rnd/out/`; runner `diff.sh`, `diff2.sh`).
Binaries: LOSAT = `s08/native/release/LOSAT` (final), NCBI = `/home/kawato/micromamba/bin/tblastx` 2.17.0, pre-S08 LOSAT = `s08/bin/LOSAT-base`
(outfmt 6 only; used to tell what S08 changed).

## 0. How the rows were checked

* Read the NCBI function for each row, then the LOSAT final code, then ran both binaries. `diff.sh`/`diff2.sh` run the same
  argument list on NCBI and LOSAT and compare stdout, stderr and exit status byte for byte.
* Frozen fixtures re-run with the final binary: `tests/tblastx_regression_fixtures.py check` 40/40 same;
  `docs/evidence/losat_web_e2a/check_losat.py --programs tblastx` 26 fixtures: 24 same, 2 `exception` (code4 rows; approved), 0 differ.
* Differential runs of my own (all stdout + stderr + exit identical unless listed under X1, X2, row 29):
  * `rnd/q1..q160 x s1..s160`: random coding-like DNA with low-complexity inserts on both strands (CAG, poly-A, ...), N runs, IUPAC letters, 1-4 queries
    (some invalid: all-N, 1-2 nt, stop-only), 1-3 subjects (forward/reverse, lengths not multiples of 3), outfmt 0 and 7, random option mix
    (`-seg no`, `-query_gencode 2|4|5|11`, `-max_target_seqs 1|2|3|5`, `-evalue`, `-threshold`, `-window_size`): 320 runs.
  * `rnd/bq1..4 x bs1..4`: 26-42 queries, 15k-26k nt, default 10002-nt batches, invalid queries inside batches, outfmt 0, 6, 7.
  * `rnd/mq1..25 x ms1..25`: 8-16 subjects with exact duplicates (ties), `-max_target_seqs 1,2,3,4,6`, outfmt 0 and 7: 250 runs.
  * `rnd/bi1..40`: stop-free biased frames with `-seg no` (7 footers with non-ideal Karlin blocks), outfmt 0 and 7.
  * `edge/*`: 1, 2, 3-nt and all-N subjects; mixed valid/short/all-N subjects; 3-nt and 6-nt queries; poly-A, poly-CAG, poly-T queries; long query
    deflines in the invalid-query warning (35-byte title truncation); lowercase FASTA (rnd/lc_*); `-num_threads 4|8` (stdout same; the stderr thread warning is the approved exception).
  * Inventory inputs `qA,qT,qE,qS8,qS,qS6,qS10,qB,qC,qD,q1,qPA,qMulti,qBigplus2,qSmallplus2,q2,q4,qN,qStop` against their subjects, outfmt 0/6/7, with and
    without `-seg no`, incl. `qPA` (1440 comparator-equal HSPs) and `qBigplus2` (second batch fully invalid: two warnings).
* No input in these runs showed a difference except the ones in X1, X2 and row 29 below (the `-seg no` sweep over seeds 1-700 gave seeds 31, 113, 430, 668 = X2 and 158, 683 = row 29).

## 1. Findings (GAP rows)

### X1. `-culling_limit N` (N >= 1) changes or empties the report (outside E; belongs to G rows `CreateHspWriter / CreateHspPipe`)
`-culling_limit` is in the path scope of COMMON.md and is neither ported nor rejected. Same result in the pre-S08 binary.
```
cd rscratch_E
tblastx            -query qE.fa -subject sC.fa -seg no -culling_limit 1 -outfmt 6 | wc -l   -> 18   (NCBI)
LOSAT tblastx      -query qE.fa -subject sC.fa -seg no -culling_limit 1 -outfmt 6 | wc -l   -> 0    (LOSAT; also in LOSAT-base)
tblastx -query qA.fa -subject sA.fa -culling_limit 1 -outfmt 6 | wc -l -> 91 (NCBI) ; LOSAT -> 0
tblastx -query qA.fa -subject sA.fa -culling_limit 2 -outfmt 6 | wc -l -> 128 (NCBI) ; LOSAT -> 90
-culling_limit 100 -> equal (206 / 206); 2, 3, 5 on qE: equal; qA with 5: 148 vs 143.
```
`algorithm/tblastx/hsp_culling.rs:apply_culling` runs per subject inside the per-subject job (run_impl.rs ~3319), while NCBI's
culling writer (hspfilter_culling.c) culls across all subjects of a query per context and builds the hit lists in
`s_BlastHSPCullingFinal`; the E-range rows (30-33) describe the default collector only. outfmt 0/6/7 all inherit the wrong HSP set.

### X2. Chain (sum statistics) e-values differ from NCBI in the last bits, which changes the HSP order inside a subject
Reproduction (compact; outfmt 6 shows it too, and the pre-S08 binary behaves the same):
```
cd rscratch_E
tblastx       -query bs/gap_q.fa -subject bs/gap_s.fa -seg no -outfmt 6   # NCBI order of rows 1-4: 31-105, 195-263, 101-27, 264-214
LOSAT tblastx -query bs/gap_q.fa -subject bs/gap_s.fa -seg no -outfmt 6   # LOSAT:                   31-105, 101-27, 264-214, 195-263
```
NCBI (`-outfmt 8`, ASN) gives the two chains of 2 HSPs the sum_e values 3.2628699764908501e-18 (scores 114+65, plus frames) and
3.2628699764908802e-18 (scores 69+110, minus frames), both printed `3.26e-18`. Both chains have the same raw score sum (179) but NCBI's e-values
differ by about 9e-15 relative = 1 ULP of the normalized score. LOSAT (stage dump `LOSAT_DUMP_TBLASTX_STAGE`) gives all four HSPs the same e-value
3.26286997649087786e-18, so the comparator falls to the raw score: 114, 110, 69, 65 -> the order above. NCBI sorts by e-value first and keeps each chain together.
Across the 10 chains of that run the LOSAT e-value differs from NCBI's (ASN has 15 digits) by a median 7.1e-15 relative, whereas the 13 single HSPs agree to the 15-digit
resolution (max 3e-15): the chain path (BLAST_LargeGapSumE) is off by about 1 ULP of the normalized score.
Root cause (confirmed by experiment): `stats/sum_statistics.rs:793` (large_gap_sum_e) evaluates `xsum - num*ln(prod) + lnfact` left to right; NCBI evaluates
`xsum -= num*log(lcl_subject_length*lcl_query_length) - BLAST_LnFactorial((double) num)` (blast_stat.c:4560-4561), i.e. `xsum - (num*ln(prod) - lnfact)`.
Experiment: a copy of `LOSAT/src` in `rscratch_E/lcopy` (outside the repository, built into `rscratch_E/lcopy_target`) with that one line changed to
`xsum - (num*ln(prod) - lnfact)` makes seeds 31, 113, 430 and 668 byte-identical to NCBI in outfmt 6 (binary kept as `rscratch_E/lcopy_LOSAT_assoc_fix`); nothing else changes in 700 random `-seg no` inputs (seeds 1-700; seeds 158 and 683 stay different, see row 29).
Frequency: 4 of 700 random inputs with `-seg no` (seeds 31, 113, 430, 668); with default SEG 0 of 200 (outfmt 0). The small-gap path (`small_gap_sum_e`, sum_statistics.rs:687-690) already
subtracts in NCBI's order and is not used by tblastx (ordering_method 1). This belongs to the engine/linking inventory (not E), but E's order rows (18, 26-29) depend on exact double equality of e-values.

### Row 29 (GAP): comparator-equal HSPs of different query frames arrive in a different order
```
cd rscratch_E
diff <(tblastx -query bs/gap2_q.fa -subject bs/gap2_s.fa -outfmt 6 -seg no) <(LOSAT tblastx -query bs/gap2_q.fa -subject bs/gap2_s.fa -outfmt 6 -seg no)
< q158 s158_1 75.000 12 3 0 367 402 391 426 0.002 24.0    (NCBI line 89; LOSAT has it at line 90, after the 368-403 HSP)
```
outfmt 0: the two blocks `Frame = +1/+1` (Query 367 VGQQKKKKKKKK, Positives 11/12) and `Frame = +2/+1` (Query 368 LGSKKKKKKKKK, Positives 10/12) are swapped.
Both HSPs: raw 46, equal e-value, query aa offset 122, subject aa offset 130, equal ends: `s_EvalueCompareHSPs` returns 0, so the order is the incoming order of the list. The pre-S08
binary gives the same output as the final one. The single query q158 against the single subject s158_1 is enough (547 + 443 bytes, poly-K region); seed 158 of `rnd/`.


Second counterexample for row 29 (`bs/gap3_q.fa`, `bs/gap3_s.fa`, `-seg no`, outfmt 6; seed 683): the subject has a 9-N island. NCBI rows 4 and 6 are
`217 321 126 22` (2.84e-53) and `219 323 124 20` (1.37e-52), LOSAT has the two HSPs exchanged: the two HSPs have the same score (62.0 bits), and each is the partner of
a different 146/144-bit HSP in a 2-HSP chain; which partner a chain takes is a tie of the linking DP, i.e. again a matter of the incoming order. The pre-S08 binary equals NCBI here
(it did not resolve the N island with `CRandom`), the final one differs. Over the 514 runs (seeds 1-700, `-seg no|yes`) with an ambiguous subject, the final binary differs from the pre-S08 one in 236 and equals NCBI in 233 of those
(the pre-S08 binary equals NCBI in 1); the three exceptions are seeds 430 and 668 (X2) and 683 (this tie). The random resolution (`resolve_ncbi4na_to_ncbi2na`) is therefore right in all but possibly this one case; whether the 683 tie comes from the
resolved bases or from the DP order (UNSURE: cannot be told from the outside; with the N island replaced by A, C or ACGT letters all three binaries agree with NCBI).

## 2. Rows whose classification needs a remark

* Row 3 (`eSequenceComparison`): n/a (env var BL2SEQ_LEGACY, out of scope). Observation: BLASTN rejects BL2SEQ_LEGACY (`algorithm/blastn/blast_engine/run.rs:5624`);
  TBLASTX ignores it and prints the database-scan report while NCBI prints the legacy bl2seq report (`out/leg.n.out` vs `out/leg.l.out`).
* Row 6, 30, 32: `faithful`: the old behaviour (common.rs `ncbi_query_subject_groups`) is no longer on the TBLASTX path; S08 sorts all formats in
  `report.rs:final_hit_order` (shared `blastn::hsp::HitList`). Behaviour is NCBI's in all my runs.
* Row 28: `faithful` while NCBI's `qsort` is a stable merge sort (glibc 2.39 here). LOSAT sorts with a stable sort always; ties of comparator-equal HSPs
  (`qPA`) are identical. A libc with an unstable qsort could order such ties differently (not testable here).
* Row 31: was divergent (score of the first HSP after the e-value sort); fixed. Oracle input kept in `bs/keep_q1.fa` + `bs/keep_sXY1.fa`
  (`-seg no`): X and Y have the same best e-value (2.96e-21, a chain); X has the larger max raw score. NCBI: X, Y. LOSAT-base: Y, X. Final: X, Y.
* Row 20: Seq-align `num_ident`/`num_positives` are not consumed (tabular.cpp:1003-1021 overwrites num_ident from the displayed rows); the final code computes both from the displayed rows.
* Row 36/38/59: the footer is a new TBLASTX writer (`write_tblastx_query_footer`), not `write_blastx_query_footer` as the inventory proposed.
* Row 39: `write_tblastn_unsearched_query_footer_spacing(writer, false)`: the TBLASTN writer got a `trailing` flag in S08; the TBLASTN call keeps `true`.
* Row 41: LOSAT checks "average subject length > 0" (run_impl.rs:1901) before "any valid context" (1914); NCBI checks the score block first. Only a zero-length subject record can tell them apart (deferred: records without data).

## 3. Inventory corrections

* Row 17 note "Check LOSAT num after the second BLAST_LinkHsps": done; the engine links twice (prelink at run_impl.rs:3128-3139, final link ~3270), `num` comes from the final pass.
* Row 11 mentions "no culling here": true for the Seq-align builder, but `-culling_limit` itself is not implemented as NCBI does (X1).
* Row 28 says "differs only for comparator-equal HSPs": also true for chain e-values that differ only in the last bit (X2).
* Row 29 UNSURE: qPA is identical, but a second random input differs (GAP, see above); the UNSURE was justified.
* Row 33 says "LOSAT TBLASTX parses max_target_seqs but nothing applies it": that was the state before the port; now applied (`final_hit_order`, `process_options` warning).

## 4. Order of stderr and stdout for E (checked)

* `-max_target_seqs < 5` warning before the outfmt 0 prolog and before any query is read (run_impl.rs:759-767); unsearched-batch warnings come with
  their query's report, in query order, exit 0 (`QueryWarnings`); all identical with `2>&1` fixtures `warnings.*.merged`.

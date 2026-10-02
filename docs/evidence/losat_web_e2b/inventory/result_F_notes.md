# Result notes, range F: tblastx pairwise report residues vs the search translation, and the genetic codes

Final code: LOSAT-web-gui commit `0d533ba76` (`LOSAT/src`). LOSAT binary: `/home/kawato/.cache/losat-web-gui-target/s08/native/release/LOSAT`. NCBI oracle: `/home/kawato/micromamba/bin/tblastx` 2.17.0+. Nothing in the repositories was modified; all inputs and outputs are in `rscratch_F/` (comparison harness `cmp.sh TAG Q S [args]` runs NCBI and LOSAT for `-outfmt 0`, `6`, `7` and compares stdout, stderr and exit status).

Result TSV: `result_F.tsv` (36 rows). Counts: ported 16, faithful 8, n/a 7, reused 2, GAP 2, exception 1, rejected 0, deferred 0. No UNSURE rows.

## GAP rows

### Row 13 (and the cause shared with row 25): letters that NCBI's FASTA reader removes are kept by TBLASTX

NCBI never gives the translation a letter outside IUPAC-na (including `-`, `*`, `X`, `x`, digits, `.`, and protein letters such as `J`, `E`, `O`, `Z`): `CFastaReader` drops it and writes a warning on stderr. The final TBLASTX code reads the records with `bio::io::fasta` (`run_impl.rs:602 read_tblastx_fasta_records`) and keeps every byte; the letter then translates to `X` and shifts all later nucleotides by one. No warning is written. This is neither ported nor rejected, and is not in the "deferred" list of RESULT_COMMON.md (BLASTN had it as an S07p item, `docs/losat_web_gui_sessions/session_s07p_e2c_blastn_scoring_options.md` item 6; BLASTX has the NCBI warnings in `algorithm/blastx/native.rs:446-574`).

Reproduction (from `rscratch_F/`; `t4_qgap.fa`, `t4_qs.fa`, `t4_qx.fa` are a 544-nt query plus one hyphen, `*` or `X` at nt 303 (a hyphen file is 545 bytes of sequence); `lt_*.fa` are the same query with one extra letter inserted at nt 324; subject `t4_s.fa`):

```
tblastx -query t4_qgap.fa -subject t4_s.fa -seg no -outfmt 6        (NCBI)
LOSAT tblastx -query t4_qgap.fa -subject t4_s.fa -seg no -outfmt 6
```

- NCBI stderr: `CFastaReader: Hyphens are invalid and will be ignored around line 7`; LOSAT stderr: empty.
- NCBI `-outfmt 0`: `Length=544`; LOSAT: `Length=545`.
- NCBI outfmt 6 line 2: `qgap	s1	100.000	80	0	0	543	304	539	300	9.24e-126	203`; LOSAT: `qgap	s1	100.000	80	0	0	544	305	539	300	9.24e-126	203` (the lines for nt before the letter are equal; later lines differ in coordinates, and some HSPs, E-values and lengths).
- For `*` (`t4_qs.fa`), `X`, `x`, `J`, `E`, `O`, `Z`, `1`, `.` (`lt_star.fa`, `lt_X.fa`, `lt_x.fa`, `lt_J.fa`, `lt_E.fa`, `lt_O.fa`, `lt_Z.fa`, `lt_dig.fa`, `lt_dot.fa`): NCBI stderr `FASTA-Reader: Ignoring invalid residues at position(s): On line 7: 24` (position depends on the file), stdout differs from LOSAT's in the same way. `B` (valid IUPAC) is identical.
- The subject goes through the same reader in NCBI, so the same applies to it (not run separately). Letters tested for the query: `- * X x J E O Z 1 .` (dropped by NCBI) and `B` (kept, identical); the full set of letters NCBI drops was not enumerated.

The mapping `display_base` (`algorithm/tblastx/report.rs:56`) itself is a literal port of `sm_BaseToIdx`; the GAP is only that the branch is reachable in LOSAT.

### Row 27: `U` (RNA) on the minus strand of the search translation

NCBI reads `U` as `T` (the Bioseq holds T) so every frame of an RNA FASTA equals the DNA result. `generate_frames` (`algorithm/tblastx/translation.rs:100`) builds the minus frames with `bio::alphabets::dna::revcomp`, whose complement table (bio 1.6, `alphabets/dna.rs:43`, `AGCTYRWSKMDVHBN`) leaves `U` unchanged, and `GeneticCode::get`/`base_mask` then reads that `U` as `T` (mask 8) instead of `A`. The plus frames and the report rows (`display_base`: `U` -> 8, then bit-reversal complement) are right; the minus frames of the search are wrong, so the HSP list differs.

Reproduction: `e_core.fa` (600 random nt) and `e_core_rna.fa` (the same with every `T` replaced by `U`).

```
tblastx -query e_core_rna.fa -subject e_core.fa -outfmt 6      NCBI == output for e_core.fa vs e_core.fa
core	core	100.000	199	0	0	3	599	3	599	8.85e-141	481
core	core	100.000	199	0	0	599	3	599	3	5.53e-139	475
LOSAT:
core	core	100.000	199	0	0	3	599	3	599	8.85e-141	481
core	core	100.000	200	0	0	1	600	1	600	1.22e-124	427
```

An RNA subject differs the same way (`e_core.fa` vs `e_core_rna.fa`: NCBI 15 rows, LOSAT 12). `t4_qu.fa` (one `U` in a 545-nt query) differs likewise: NCBI first rows `qu s1 100.000 101 0 0 303 1 303 1 2.60e-126 251`, LOSAT `qu s1 100.000 100 0 0 300 1 300 1 9.24e-126 249`; stderr is empty in both. Impact low (RNA letters only). Fix sketch (not applied): map `U/u` to `T/t` before `dna::revcomp` in `generate_frames` and in `preliminary_subject_bases` consumers (the latter already maps U to 8).

## Rows resolved by oracle comparison (all byte identical in outfmt 0, 6, 7, stdout + stderr + exit status unless noted)

- `t5_amb.fa` vs `t5_conc.fa` and the reverse, `-seg no`: all 3375 IUPAC triplets, 27,007 codons, random lowercase, six frames (rows 8-12, 14, 23, 26, 34).
- `-query_gencode N` for all 26 accepted codes on `t5_amb.fa` vs `t5_conc.fa`, outfmt 0/6/7 (`codes.sh`, 78 runs, 0 differ) (rows 1, 3, 7, 12, 14, 22-25).
- `-db_gencode N` for all 26 codes, query `t5_conc.fa`, subject `t5_amb.fa`, outfmt 0: residue by residue comparison of the Query and Sbjct rows keyed by (frame, nucleotide position) over every HSP both programs print: 4,212,879 columns in common, 0 differ (`dbcodes.py`); the HSP sets differ slightly because of the approved exception (rows 3, 7, 31).
- `t6_q.fa` vs `t6_s.fa` (3375 queries, one IUPAC triplet each) `-seg no -max_target_seqs 5`, also with `-query_gencode 5`: identical (125,690 HSP rows) (rows 25, 34).
- `t7_q.fa` vs `t7_s.fa` (RAY/SAR/MTA/NNN/RAY vs SAR/RAY/HTA/NNN/NNN): identical, `97.500 120 3 0 ...` (rows 17, 19, 34); `t2_qa/qb/qc.fa` vs `t2_s.fa` identical.
- `pairs_q.fa` vs `pairs_s.fa`: all 625 ordered pairs of the letters `ARNDCQEGHILKMFPSTWYV B Z J X *` as aligned columns (the file builds each pair between identical scaffold codons; 8569 HSPs, 2.5 MB of outfmt 0): identical, so the midline `+`, Identities and Positives follow the display matrix for every pair (rows 17, 34).
- `t8_q.fa` vs `t8_sa/sb/sc/sd.fa` (N islands), `t9_s.fa`, `t9x_0..t9x_11.fa` (subject lengthened by 0..11 nt, the NCBI output takes 11 distinct values): identical, also with `-num_threads 3/4` (stdout; the stderr thread warning is the approved exception) (rows 28, 29). `tblastx_ambig_query/subject.fasta` identical.
- LC738874 vs LC738875 (default options, SEG lowercase present): identical in outfmt 0/6/7 (rows 4, 5, 35, 36). `EDL933` 6 kb piece vs `Sakai.fna` (5.5 Mb subject): identical.
- Edge inputs: 1- and 2-nt queries, all-N query and subject, lowercase queries/subjects, mixed-case, lengths mod 3 = 0/1/2, IUPAC mixes in both files, several queries and subjects (`e1`..`e10`, `lc1`, `lc3`): identical.
- `-db_gencode 5` on `t3_q.fa`/`t3_s.fa` (row 31, approved exception): NCBI first rows `q1 s1 88.889 180 20 0 540 1 540 1 5.57e-126 431`, `q1 s1 100.000 180 0 0 1 540 1 540 1.44e-125 430`; LOSAT `q1 s1 100.000 180 0 0 1 540 1 540 1.44e-129 443` (same extents, pident, mismatch; score, bit score, E-value and the HSP set differ). With `-db_gencode 1` identical.
- Parser: `-query_gencode`/`-db_gencode` 7, 8, 17, 32, 34, 0 rejected by both (NCBI exit 1, `Error: Argument "db_gencode". Illegal value, expected values between: 1-6, 9-16, 21-31, 33:  `7'`; LOSAT exit 2, clap text: approved exception).
- Defline modifiers `[gcode=4] [organism=...]` are ignored by NCBI tblastx (output equals the unmodified file), so row 22's "no Source descriptor" holds and LOSAT needs nothing.

## Inventory corrections and remarks

- F row 13 (`n/a`, "gap/X letters cannot reach the display from FASTA") is true for NCBI but not for the final LOSAT code: see GAP row 13.
- F row 27 / F_notes section 5 treat `generate_frames` as equal to NCBI's reverse strand translation; it is not for `U` (GAP row 27). Everything else about frames (order, sentinels, trailing 1-2 nt ignored) is equal.
- F row 31 (approved exception) and F rows 28/29/19/34-36 (divergent/missing before) are resolved; the three divergences in F_notes section 5 (pident/mismatch from engine `num_ident`, preliminary subject frames from 4na, rows not built) no longer reproduce (t7, t8, t9x are byte identical).
- F row 29 says `resolve_ncbi4na_to_ncbi2na` was "used by TBLASTN, not by TBLASTX"; now also used by TBLASTX through `preliminary_subject_bases` (`run_impl.rs:546`).
- F row 17 listed `is_positive_match` as `reusable`; the final code uses the shared `is_positive_match` for the midline but a new function `row_counts` (`algorithm/tblastx/report.rs:182`) for the counts, hence `ported`.
- F row 14: `display_codon` was made `pub(crate)` and is shared with BLASTX (`reused`); its pinned test is `display_codons_match_pinned_cpp_all_26_codes` (`algorithm/blastx/report.rs:1891`).
- F_notes open point about TBLASTN (`algorithm/tblastn/stage_e_report.rs:632-700` showing `X` where NCBI shows `B/Z/J`) is outside TBLASTX and was not examined.
- The LOSAT line numbers in `F.tsv` were taken before the port (the notes there say the tree also moved while the inventory ran); `result_F.tsv` gives the final ones.

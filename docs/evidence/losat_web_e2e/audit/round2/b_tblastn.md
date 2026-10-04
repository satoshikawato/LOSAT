# Round 2, angle (b) TBLASTN: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned. Binary `d8d18ec0…cb116` (the gate native of `491292327`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/b/` (`results.tsv` 4,221 new runs, `FINDINGS.txt`, `logs/`, saved DIFF outputs in `w/diff/`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_B.md`](brief/ANGLE_B.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md).

## Overall verdict: SUPPORTED

About 8,000 NCBI-vs-LOSAT runs, no output difference that is not an approved exception or an explicit rejection. Two new items, neither a wrong-output parity defect. Load average was 30-55 throughout, so timings are noisy. The round 1 `xf_run` harness runs NCBI under `ulimit -v 25000000`, where NCBI times out at X-drops 1e7-5e8; re-run without the limit, the outputs are identical.

## Round 1 findings, re-run

| ID | Verdict | Evidence |
|---|---|---|
| TN-1 | SAME | 9 original commands; X-drop sweeps XD (104) and S14 (151, 150 SAME) |
| TN-1 residual | SAME | `min_q.faa` × `min_s.fna` with `-comp_based_stats 0` (xf 1e10 ev 1000/1e5, xg 1e10 ev 1e5, xf +inf ev 1e10); window-end suite 1,208 runs SAME |
| TN-2 | SAME | `grid.py` 160/160; mts 1/5/10/100/500 and BLOSUM45; MK suite 396 runs SAME |
| TN-3 | SAME | 3 outfmt 0 commands |
| TN-4 | SAME output; time and RSS comparable | XD suite, 104 runs (below) |
| TN-5 | SAME; `+inf` LOSAT-REJECTS (D12) | BLOSUM45 w2 cbs 0 at ev 1e4/2000/5000/1e5/1e10/1e300; TIE45/TIE62 351 SAME, 27 rejections (BLOSUM62 `-word_size 2`) on 24 generated AT-rich and repetitive pairs |
| TN-6 | LOSAT-REJECTS as decided (D12) | NCBI 139 for `+inf`, `1e999` on `e2e_many_subject`; 1e308, 1e300 SAME; see R2b-1 |
| TN-7 | LOSAT-REJECTS as decided (D9) | `-dryrun` |
| TN-8, TN-9 | ACCEPTED | wording; both exit 1 |
| TN-10, TN-11 | FIXED | references checked by hand |
| RP-4 | SAME | 18/18 on bigq × bigs/bigs4 |

Harness totals (about 3,750 runs): `pairs.py` 650 (350 SAME, 194 parser, 106 reject, 0 DIFF); `pairs2.py` 504 (392, 38, 74, 0); `grid.py` 160 SAME; original random sweeps about 1,690 (1,637 SAME, 51 reject, 2 DIFF = LOSAT's 60 s timeout under load on about 10 MB outputs, SAME with a longer timeout); replay of 693 battery commands (485, 65, 133, 5 DIFF = approved `-db_gencode 4`/`32` and NCBI's thread warnings); 37 huge commands (31, 4, 2, 0); new random sweeps seeds 21-23, 450 (434 SAME, 16 reject, 0 DIFF).

## New checks around S08+a

- Query splitting: 98 SAME (19,800-20,100 not split, 39,700-40,100 two chunks, 59,699-59,701 three chunks; batches of 2-14 queries of 100-45,000 residues incl. all-X, 1- and 2-residue and lower-case queries; `BATCH_SIZE` 1-100000; outfmt 0/6/7; `-num_threads` 1/2/4/8).
- `CHUNK_SIZE`/`OVERLAP_CHUNK_SIZE` (about 40 values): every case NCBI fails is a justified rejection (null-pointer exception when a chunk would be split again, `StringToInt` exit 255). 150-value environment matrix: 56 SAME, 94 rejections.
- Masks at chunk boundaries: 560 SAME (lower-case runs of 1, 2, 3, 30 at and around chunk ends, `-soft_masking`, SEG windows across a boundary, X/B/Z/J/U/`*`).
- Long HSPs across boundaries: 162 SAME (300-6,000 residues; mts 2, cbs 0, `-sum_stats false`, ev 1000, X-drops).
- Option crosses: 80 SAME on split sets; about 190 random option and environment combinations SAME or both reject.
- Toolkit words in a value position: 177 runs, explicit rejections (pending the maintainer) or SAME/parser errors where the word is a plain value.
- TN-4 (`/usr/bin/time`): `-xdrop_gap` and `-xdrop_gap_final` 1e6…1e300 and `+inf`, cbs 0 and 2, all identical to NCBI. `e2e_tblastn_subject`: NCBI 0.1-0.9 s and 70 MB (427 MB at 1e6), LOSAT 0.1-1.7 s and about 9 MB. `e2e_many_subject`: NCBI 0.2-0.9 s, LOSAT 4-6 s (same gap with default options: the known cost). Long-HSP sweeps to 5e4 bits: NCBI up to 60 s, LOSAT up to 71 s; bigq × bigs `-xdrop_gap_final 5000`: NCBI 278 s, LOSAT 159 s.

## New findings

**R2b-1 (low-medium; D12 incomplete).** `-evalue` exactly DBL_MAX crashes NCBI (SIGSEGV), LOSAT prints a result. `tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_many_subject.fna -evalue 1.7976931348623157e308 -outfmt 6`: NCBI 139 every time (also `1.7976931348623158e308`, `-sum_stats false`, `-seg no`); no crash with `-comp_based_stats 0`, `-max_target_seqs 5`, on `e2e_tblastn_subject.fna`, or at `1.7976931348623156e308` (7,221 rows). LOSAT: 0 and 7,221 rows, identical to NCBI at `…156e308`. gdb: `blast_kappa.c:409` sets `*pbestEvalue = DBL_MAX` for a list without surviving HSPs; `blast_kappa.c:3687` `best_evalue <= expect_value` is then true; `blast_hits.c:3266` reads `hsp_array[0]->score` of the empty list (the TN-6 crash). LOSAT `tblastn/args.rs:497-503` rejected only infinite values. BLASTP shares `blast_kappa.c` (not tested).

**R2b-2 (low, informational; performance only).** Tiny `CHUNK_SIZE` makes LOSAT about 10× slower per chunk: `CHUNK_SIZE=1 OVERLAP_CHUNK_SIZE=0 tblastn -query b/in3/L40000.faa -subject b/in3/L40000.fna -outfmt 6` (40,000 one-residue chunks): NCBI 70-83 s with empty output, LOSAT still running at 500 s. `CHUNK_SIZE=7 OVERLAP_CHUNK_SIZE=0`: NCBI 8.7 s, LOSAT 89.6 s, same empty output. `CHUNK_SIZE=500 OVERLAP_CHUNK_SIZE=250`: NCBI 1.2 s, LOSAT 4.3 s, same output. Cost: the per-chunk `query_set_setup` in `tblastn/stage_d_pipeline.rs` `split_preliminary_hitlists`. Default settings make at most a few chunks.

## Limits

Two heavy random cases (`-evalue 1e10` with a 100k-residue query, an H16 BLOSUM45 case) timed out on both sides (inconclusive). Subject-side chunking (above about 15 Mb) not exercised. The composition-window sorts in `core/composition_adjustment/redo_alignment.rs` remain `sort_unstable_by`; 243 TIE62 runs including cbs 2 found no effect.

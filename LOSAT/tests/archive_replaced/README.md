# Archives and strategies with replaced bytes (line RG K4b)

`cases.tsv` lists the cases; `ncbi/` holds NCBI BLAST+ 2.17.0's output of each (the comments in `cases.tsv` say how).

`in/`: query and subject files of these cases; no generator or source of the letters is recorded. Added in commit
e8d3cee03 (2026-11-14, RG K4b): `n_ctrl.fa`, `n_ff.fa`, `n_lcl_del.fa`, `n_two_e.fa`, `p_e.fa`, `p_plain.fa`, `s_n_e.fa`,
`s_p_e.fa`; in commit 85446c15b (2026-11-14, RG K4m): `n_m3.fa`, `s_m8.fa`; in commit 6f9b46b5c (2026-11-14, RG K4b
review L7): `k_pallx.fa`, `s_pm24.fa`, `s_pmx.fa`. (`k_*.fa` are the queries of the strategy rows and `s_pm24.fa`, `s_pmx.fa`
the subjects whose titles are in the archive's request, as `cases.tsv` says.)

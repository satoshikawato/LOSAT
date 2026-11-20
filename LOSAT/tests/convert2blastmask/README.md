# convert2blastmask fixtures (line RG item K3)

`cases.tsv` lists the cases; `ncbi/` holds the official NCBI BLAST+ 2.17.0 `convert2blastmask` output of each.

`in/`: made by `remaining-gaps-20261112/K/k3/probe/gen_inputs.py` (a probe script kept with the session's files, not in
this repository), except `gapinv.fa` and `gapwrap*.fa` (by hand) and `id_del.fa` and `id_del_gnl.fa` (from a review's
probe inputs); see the comments in `cases.tsv`. The files came in with commits 8e9b162a8 (2026-11-13), 0df633f67 and
9fbaaf389 (2026-11-13), 499bce02c (2026-11-14) and d8bb13d2a (2026-11-15). No source of the letters is recorded.

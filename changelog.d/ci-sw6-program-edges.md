- CI: a ratchet (`LOSAT/tests/check_program_edges.py`, run by `rust lint`) fails when code under
  `LOSAT/src/algorithm/` gains a reference from one program to another (blastn, blastp, blastx,
  common, psiblast, pssm_search, tblastn, tblastx) that `LOSAT/tests/program_edges.tsv` does not
  list, when a reference goes up the layering `common < blastn < tblastx < blastp < tblastn <
  blastx`, listed or not, when a path in a program names the `algorithm` module itself instead of
  a program, and when a listed reference is gone without leaving the list. Rows are added only by
  review, for an edge down the layering; remove a row when its edge goes (fewer edges are welcome).
  (SWEBOK 2026-11-04 audit, KA3.)

# Inputs of the outfmt 11 (archive) rows: where the `arch_*` files come from

Inputs of the rows `sa.arch.*` and `db.sa.arch.*` of `LOSAT/tests/outfmt0_manifest.tsv` and `db_search_manifest` (outfmt 11
with its outfmt 8 and 5 siblings; the frozen outputs are NCBI BLAST+ 2.17.0's). No generator or source of the letters is
kept in the repository ("origin not recorded"); the commits that added the files, by program:

- `arch_alllc.fna`, `arch_crlf.fna`, `arch_lcase.fna`, `arch_nohit.fna`, `arch_qstar.fna`: BLASTN, commit 7ff889e94,
  2026-10-22; `arch_lcpart.fna`: commit c12b5c1a1, 2026-10-23.
- `arch_bp_*`: BLASTP, commit 3f475e72e, 2026-10-23.
- `arch_bx_*`: BLASTX, commit 7b4befdae, 2026-10-23.
- `arch_tbn_*`: TBLASTN, commit 8a53b9d5a, 2026-10-24.
- `arch_tbx_*`: TBLASTX, commit c6bd112a4, 2026-10-24.

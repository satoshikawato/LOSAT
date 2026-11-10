- Web ABI v1: `losat_web_run_pair` accepts `-searchsp` and `-dbsize` for `blastn`, `tblastx` and `blastp`
  (gbdraw request, 2026-11-10: a comparison's E-values take the scale of its genomes). The values are read
  by the CLI's parsers and give the CLI's report, which is NCBI BLAST+ 2.17.0's (shared vector
  `LOSAT/tests/fixtures/web_v1_search_space/`); an invalid value is refused as on the CLI. `blastx` and
  `tblastn` already took both. The other options and the version of the ABI are unchanged (`docs/web/abi_v1.md`).

- Web ABI v1: `losat_web_run_pair` and `losat_web_run_pair_handles` accept `tblastn` (a protein query
  against nucleotide subjects; Owner request 2026-11-04, for gbdraw's serial module). The request runs
  `losat tblastn` itself over the call's bytes, so the report and the error messages are the CLI's
  (shared vector `LOSAT/tests/fixtures/web_v1_tblastn/`); `-db` and the file and help words are
  refused. The other programs of ABI v1 are unchanged (`docs/web/abi_v1.md`).

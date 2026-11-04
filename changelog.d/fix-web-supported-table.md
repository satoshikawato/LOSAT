- LOSAT Web: a run with the default options is marked "Supported" in every mode and task, with one thread and
  with Auto (WA-RS-03). The default BLASTN task (megablast), default BLASTX and TBLASTN `-task tblastn-fast`
  gain outfmt 7 rows frozen from NCBI BLAST+ 2.17.0 (`def7.*` in `LOSAT/tests/outfmt0_manifest.tsv`; they had
  outfmt 0 and 6 only), V-BR checks the browser runtime at 1, 2, 4, 8, 16 and 32 threads in Chromium, Firefox
  and WebKit (was 1, 2 and 4), and a run's thread count is covered up to the largest count checked, since Auto
  picks half of the user's processors (up to 32).

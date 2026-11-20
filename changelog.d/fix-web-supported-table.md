- LOSAT Web: a run with the default options is marked "Supported" in every mode and task, with one thread and
  with Auto (WA-RS-03). The default BLASTN task (megablast), default BLASTX and TBLASTN `-task tblastn-fast`
  gain outfmt 7 rows frozen from NCBI BLAST+ 2.17.0 (`def7.*` in `LOSAT/tests/outfmt0_manifest.tsv`; they had
  outfmt 0 and 6 only). V-BR checks the browser runtime at every thread count that Auto (serial below 20,000
  bytes of input, then half of the processors: at most 8 below 1 MB, 32 from 1 MB) and the Threads list can
  choose: each of its 21 searches at 1, 2, 4, 8, 16 and 32 threads, one search of each program at every other
  count up to 32, in Chromium, Firefox and WebKit (was 15 searches at 1, 2 and 4). It also runs the defaults of
  megablast, dc-megablast, blastn-short, blastp-fast, BLASTX and tblastn-fast.

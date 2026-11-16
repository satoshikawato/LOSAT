# LOSAT Web: V-BR at every thread count and the tasks' defaults (line WA, W1 PR-b, WA-RS-03)

Auto picks a thread count from the user's processors (serial below 20,000 bytes of input, then half of
the processors: at most 8 below 1 MB, 32 from 1 MB; `web/app/src/infra/engine-worker/policy.ts`), and
the Threads list offers 1 to 32. V-BR had checked the browser runtime at 1, 2 and 4 threads only, and a
run's mark is Supported only at a thread count that was checked (`web/app/src/domain/verification.ts`),
so most runs on a machine with more than 8 processors were marked Not supported.

V-BR (`web/app/tests/e2e/engine.spec.ts`) now runs:

- its 15 earlier cases at 1, 2, 4, 8, 16 and 32 threads (`v-br-<browser>.json`);
- the defaults of the search form's tasks (multi.megablast, dc.default, short.default,
  e2e.blastp.fast, pyl.blastx and tfast.e2e.def, with the `def7.*` outfmt 7 rows of
  `LOSAT/tests/outfmt0_manifest.tsv`) at the same counts (`v-br-defaults-<browser>.json`);
- one case of each program whose work the thread count divides (multi.megablast, e2e.blastp.fast,
  blastx.batches, tblastn.batches, tblastx.threaded-direct) at every other count from 3 to 31
  (`v-br-counts-<case>-<browser>.json`).

Each group runs in a page of its own: in WebKit 26.5 a page that had run some 95 searches failed the
next dc-megablast search with "ran out of memory" (memory.grow refused) at every thread count.

- [`run-20261104T111500Z-vbr/`](run-20261104T111500Z-vbr/): the records of commit `98bac1391`
  (`environment.txt`; each `records/*.json` row is one search: argv, path, threads, engine build,
  output SHA-256 by format, the native CLI's hashes, the NCBI frozen hashes matched, hits, elapsed
  time; `e2e-engine-vbr-<browser>.log` the Playwright lines; `summary.md` the counts).
- `docs/web/verification_cells.tsv`: the browser rows (S09/S12, SXc, SXd3) name these cases and counts.

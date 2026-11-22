# makeprofiledb fixtures (row PS10a)

What they are: NCBI BLAST+ 2.17.0 `makeprofiledb` runs that `tests/makeprofiledb_fixtures.rs` repeats with
`losat makeprofiledb`, and NCBI `rpsblast` 2.17.0 output on two of the databases it built.

## Origin
- Program: NCBI BLAST+ 2.17.0 `makeprofiledb` and `rpsblast` (bioconda build, `makeprofiledb: 2.17.0+`, build Aug 11
  2025), run on the local development host (WSL2) by the PS10 survey oracle on 2026-11-21: `run_cases.sh` of
  `psiblast-20261025/PR/ps10/oracle/` (outside the repository), each case in a fresh directory holding `in/`, with
  `env -i PATH=/usr/bin:/bin TZ=UTC LC_ALL=C BLAST_USAGE_REPORT=false HOME=<run dir>`, under `prlimit --core=1`, the
  wall clock frozen at 2026-11-21 05:07:00 UTC by an `LD_PRELOAD` time shim.
- Inputs (`in/`): `pc.smp`, `pe.smp`, `pf.smp` are `psiblast 2.17.0 -num_iterations 2 -out_pssm` of 32-36 residue
  protein cuts (scores and frequency ratios, scale 1); `c1.smp`-`c3.smp` are CDD v3.22 cd00001-cd00003 cut to their
  first 22/18/20 columns (scale 100, word thresholds); `fc.smp`, `fe.smp` are frequency-ratio-only forms of the
  psiblast PSSMs (the IMPALA path); `t_source.smp`, `t_orgnotax.smp` are `pc.smp` with a `source` taxid and an
  organism without a taxid; `e_trunc.smp` and `e_unknownfield.smp` are `pc.smp` cut in its frequency ratios and with
  an unknown member, `e_notpssm.smp` a `Seq-entry`, `e_empty.smp` empty (the reader's errors). Lists (`*.lst`) name the PSSMs relative to the working directory.
- Then `freeze_fixtures.py` (same place) kept the 45 rows of `cases.tsv` and, per row, NCBI's standard output,
  standard error and exit status as written (`ncbi/<row>.{out,err,rc}`; the test maps NCBI's words to LOSAT's with
  `losat::wording`) and `ncbi/<row>.files`: for each file the run wrote (outside `in/`), its path, size and FNV-1a 64
  digest; in a `.pin` or `.pal` the build's date (its 4-byte length and padding in a `.pin`) is replaced by `<DATE>`
  first and the size is `*`.
- `rps/q3.faa` is the survey's three-protein query; `rps/R01.{6,0}` and `rps/R05.{6,0}` are NCBI rpsblast
  `-query q3.faa -db db/db -outfmt 6|0` on NCBI's databases of rows `O00` (score PSSMs) and `S_f2`
  (frequency-ratio-only PSSMs).

The whole oracle (206 cases) is compared by `compare.py` of the survey; these rows are the part kept in the
repository (about 545 KB).

# Fixtures of `losat blastdbcmd` (`LOSAT/src/blastdbcmd`)

What NCBI BLAST+ 2.17.0 blastdbcmd printed for `-info`, its argument errors, `-out` and
`-logfile`, for the tests in `LOSAT/tests/blastdbcmd_info.rs`. Frozen on 2026-10-29 by line BDC's
probe on the workstation (cases `I01`-`I38` and `B01`-`B61`; NCBI was not run on this tree).

- `db/`: databases of that probe. `p5` (version 5, `-parse_seqids`) and `p4` (version 4) are
  proteins with every Seq-id type, `n5n` a version 5 nucleotide database without parsed ids,
  `nL` a version 4 nucleotide database written by `losat makedb`, `pm5` a protein database
  with the seg masks of `makeblastdb -mask_data` (column files `.paa .pab .pac`), `palias` an
  alias file of `p5` and `pm5` (blastdb_aliastool). Kept: the files LOSAT reads and the LMDB
  file of each version 5 database, which must exist.
- `expected/ID.stdout`, `expected/ID.stderr`: NCBI's bytes for case `ID`, with `{R}` for the
  directory the probe ran in (here `tests/blastdbcmd`, the current directory of the runs) and
  `{T}` for the repository's `LOSAT/tests/` (the probe's copy of the seqdb fixtures).
  `B07.logfile` and `B08.logfile` are the files NCBI's `-logfile` wrote (`B08.logfile` after
  two runs, B08 and B38, into one file).

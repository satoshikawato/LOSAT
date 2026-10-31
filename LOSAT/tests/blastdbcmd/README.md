# Fixtures of `losat blastdbcmd` (`LOSAT/src/blastdbcmd`)

What NCBI BLAST+ 2.17.0 blastdbcmd printed for `-info`, its argument errors, `-out` and
`-logfile` (`LOSAT/tests/blastdbcmd_info.rs`), and for the retrieval of sequences
(`-entry`, `-entry_batch`, `-outfmt`, `-range`, `-strand`, `-line_length`, `-ctrl_a`,
`-target_only`, `-get_dups`, `-long_seqids`: `LOSAT/tests/blastdbcmd_retrieve.rs`). Frozen on
2026-10-29 and 2026-10-30 by line BDC's probe on the workstation (cases `I01`-`I38`,
`B01`-`B61`, `C01`-`C20`, `F01`-`F34`, `S01`-`S46`, `D01`-`D112`; NCBI was not run on this
tree).

- `db/`: databases of that probe. `p5` (version 5, `-parse_seqids`) and `p4` (version 4)
  are proteins with every Seq-id type, `n5` / `n4` the nucleotides, `p5n` / `p4n` / `n5n` /
  `n4n` the same without parsed ids, `pL` / `nL` version 4 databases written by LOSAT
  (`losat makedb`, now `losat makeblastdb`), `pm5` a protein database with the seg masks of
  `makeblastdb -mask_data` (column files `.paa .pab .pac`), `palias` an alias file of `p5`
  and `pm5` (blastdb_aliastool). Kept: the files LOSAT reads, the LMDB file of each version 5
  database (which must exist; the retrieval reads its list of named databases for
  `acc2oid`), and the ISAM files of `p4` and `n4` (whose presence makes a version 4 database
  one with parsed ids).
- `in/`: the `-entry_batch` files of the probe.
- `expected/ID.stdout`, `expected/ID.stderr`: NCBI's bytes for case `ID`, with `{R}` for the
  directory the probe ran in (here `tests/blastdbcmd`, the current directory of the runs) and
  `{T}` for the repository's `LOSAT/tests/` (the probe's copy of the seqdb fixtures).
  `B07.logfile`, `B08.logfile` and `D72.logfile` are the files NCBI's `-logfile` wrote
  (`B08.logfile` after two runs, B08 and B38, into one file), `D70.out` and `D84.out` those
  of `-out`. The stdout of `D82`, `D87`, `D97` and `D101` is not NCBI's but the correct
  output of an approved difference (line BDC DECISIONS 12: `-long_seqids` with `-range` on
  the minus strand, where NCBI writes the forward letters under its `:c` header): NCBI's
  header with the letters NCBI writes for the same range with short ids
  (`dropin-survey-20261029/BDC/probe/make_expected.py`).

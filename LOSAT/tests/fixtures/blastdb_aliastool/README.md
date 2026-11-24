# blastdb_aliastool fixtures (row RG-D2a)

What they are: NCBI BLAST+ 2.17.0 `blastdb_aliastool` runs that `tests/blastdb_aliastool_fixtures.rs` repeats with
`losat blastdb_aliastool`, each row in a fresh copy of `in/`.

## Origin
- Program: NCBI BLAST+ 2.17.0 `blastdb_aliastool` and `makeblastdb` (the official `ncbi-blast-2.17.0+-x64-linux`
  build), run on the local development host (WSL2) by the RG-D2a oracle on 2026-11-24: `run_cases.py` of
  `remaining-gaps-20261112/D/d2/oracle/` (outside the repository), under its `oracle.lock`, each row in a fresh
  directory holding `in/`, with `env -i PATH=/usr/bin:/bin TZ=UTC LC_ALL=C BLAST_USAGE_REPORT=false HOME=<run dir>`,
  standard input `/dev/null`, under `prlimit --as=8GiB --core=1` and a 60 s timeout, the wall clock frozen at
  2026-11-24 07:00:00 UTC by an `LD_PRELOAD` time shim.
- Inputs (`in/`), written by `make_src.py` (same place): `p.faa`, nine proteins whose ids are GI with `ref`
  (`NP_`, `XP_`), `gb`, `emb` and `sp`, `pir`, `prf` and `lcl`; `q.faa`, three proteins; `n.fna`, four nucleotides
  (`NM_`, `XM_`, `gb`, `emb`); `p.taxmap`, `n.taxmap`; text lists of GIs (`gi_*.txt`: valid, bad lines, CR LF line
  ends, duplicates, none found, one past 32 bits), of sequence ids (`ids*.txt`: accessions with and without
  versions, GI, FASTA-style ids, bad ids, CR LF, spaces, long ids, none found), of taxids (`tax*.txt`), of
  databases (`dbs*.txt`) and of VDB runs (`vdbs.txt`), and `empty.txt`. Then `build_in.sh` (same place) built with
  NCBI's `makeblastdb` (`-parse_seqids`): `p4`, `q4`, `n4` (version 4); `p5`, `q5`, `n5` (version 5; `p5` and `n5`
  with taxids); `pm` (`-oid_masks 1`), `pvm` (the same in two volumes, `-max_file_sz 400B`) and `pmv` (two volumes,
  `-max_file_sz 600B`); and with NCBI's `blastdb_aliastool` the binary lists `gi_ok.bgl`, `gi_big.bgl`, `ids.bsl`,
  `ids_db.bsl` (of `pvm`), `ids_n.bsl` (of `n5`) and `ids_zero.bsl`; `trunc.bsl` is the first 30 bytes of `ids.bsl`.
- Rows (`cases.tsv`, bash words with single quotes only): `H` help, `A` `-db` with `-gilist`, `S` `-seqidlist`,
  `T` `-taxidlist`, `O` `-oid_masks`, `D` `-dblist` and `-dblist_file`, `N` `-num_volumes`, `V` `-vdblist` and
  `-vdblist_file`, `G` `-gi_file_in` and `-gi_file_out`, `C` `-seqid_file_in` and its options, `I`
  `-seqid_file_info`, `X` argument errors and conflicts. `-version` is not a row: its `Package:` line is LOSAT's
  (`tests/ncbi_toolkit_version.rs` covers it).
- Rows `R01`-`R08` (review 1 of row RG-D2a, run the same way on 2026-11-24 20:5x, inputs from `add_review_inputs.py`,
  same place): a text `-gilist` converted onto its own name (`self.p.gil`), `-dblist_file` and `-vdblist_file` lines
  with `\v`, `\f` and a byte that is not UTF-8 (`dbs_vt.txt`, `vdbs_vt.txt`, `dbs_ff.txt`), and binary seqid lists
  whose count of ids is 2^31, 2^32 and 2^32 + 1 (`ids_n31.bsl`, `ids_n32.bsl`, `ids_n32p1.bsl`, from `ids.bsl`) or
  whose database part names no volume (`ids_dbnovol.bsl`, from `ids_db.bsl`).
- Then `freeze.py` (same place) kept, per row, NCBI's standard output, standard error and exit status as written
  (`ncbi/<row>.{out,err,rc}`, the run directory's absolute path replaced by `<RUN>`; the test maps NCBI's words to
  LOSAT's with `losat::wording`) and `ncbi/<row>.files`: for each file the run made or changed (in `in/` too), its
  path, size and FNV-1a 64 digest, and `REMOVED` for each file of `in/` it removed. The date of an alias file
  (`# Alias file created <date>`) and of a binary seqid list (its 8-byte file size, zeroed, and its date with the
  date's length byte) is replaced by `<DATE>` first and the size is `*`.

## Rows that wait for LOSAT's reader
The rows `S02`, `S03`, `S05`, `S07`, `S09`, `S10`, `T01`, `T02`, `T05` and `T06` make an alias file with a seqid or
taxid list over a version 5 database, which NCBI then opens to count its sequences. LOSAT's reader does not yet read
such lists (line RG-D, item D1c-2), so there the test expects LOSAT's refusal and no alias file left
(`AWAITING_V5_LISTS`); NCBI's frozen result stays for the day the reader reads them.

The row `R08` names a database by a `-dblist_file` line that is not UTF-8, which NCBI writes into the alias file as
it is (then the database does not open, and the alias file is removed); LOSAT reads the line by its bytes and gives
NCBI's result (line RG item K10d, 2026-11-25; before it, LOSAT refused the line).

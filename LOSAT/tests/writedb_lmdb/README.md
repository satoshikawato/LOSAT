# Fixtures of the LMDB writer (`LOSAT/src/writedb/lmdb.rs`, `LOSAT/src/writedb/mdb/`)

Files written by makeblastdb 2.17.0 (`-blastdb_version 5 -parse_seqids`, Linux x86-64, page size
4096) for the unit tests in `LOSAT/src/writedb/lmdb_tests.rs`, which make the same inputs in code and
compare every byte. `NAME.EXT` is `db.EXT` of the run `NAME` of the makedb-v5 survey A
(`gen_fixtures.py`, `run_fx.sh` in the line's task folder, outside the repository):

- `f1p`: `f1_tiny.fa` (3 proteins).
- `f4p200t7`: `f4_ids200.fa` with `MAX_LMDB_TXN_ENTRY=7`.
- `f5pdbt3`: `f5_pdb.fa` with `MAX_LMDB_TXN_ENTRY=3`.
- `f5d600t50`: `f5_dup600.fa` with `MAX_LMDB_TXN_ENTRY=50`.
- `t9`: `f9_multi.fa` with `-taxid_map f9_multi.map`.

The `.ptf` of a database without tax ids is one fixed file; the tests compare it, and the files of
the larger runs, by size and FNV-1a digest. The tests run no NCBI program.

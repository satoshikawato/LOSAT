- Library: the engine modules `common`, `stats`, `core`, `writedb` and `blastdbcmd` are private to the crate
  (`pub(crate)`) instead of public and hidden from the docs, so no other crate can call their functions (line LIB,
  P-V3b). The code that no build used once they were private is gone (about 2,100 lines): the old
  `common::write_output*` writers and the two `report` writers only they called, the packed-sequence types of
  `core::blast_encoding`, a second BLOSUM62 table in `core::blast_util` that only tests used (they use the shipped one
  now), unused `stats` functions, constants and fields, and the Seq-data choices LOSAT's FASTA reader never makes;
  with them go 25 tests that tested only that code. The command line's output is unchanged.

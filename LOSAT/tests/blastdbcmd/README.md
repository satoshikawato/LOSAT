# Fixtures of `losat blastdbcmd` (`LOSAT/src/blastdbcmd`)

What NCBI BLAST+ 2.17.0 blastdbcmd printed for `-info`, its argument errors, `-out` and
`-logfile` (`LOSAT/tests/blastdbcmd_info.rs`), and for the retrieval of sequences
(`-entry`, `-entry_batch`, `-outfmt`, `-range`, `-strand`, `-line_length`, `-ctrl_a`,
`-target_only`, `-get_dups`, `-long_seqids`: `LOSAT/tests/blastdbcmd_retrieve.rs`), and for
the taxonomy, the metadata, the lists of databases and the search path (`-taxids`,
`-taxidlist`, `-tax_info`, `-metadata`, `-list`, `-show_blastdb_search_path`, `-ipg`:
`LOSAT/tests/blastdbcmd_discovery.rs`), on databases with masks (`-info`,
`-mask_sequence_with`, `%m`: `LOSAT/tests/blastdbcmd_masks.rs`), and on several names in `-db`
(`LOSAT/tests/blastdbcmd_several_names.rs`). Frozen on
2026-10-29 and 2026-10-30 by line BDC's probe on the workstation (cases `I01`-`I38`,
`B01`-`B61`, `C01`-`C20`, `F01`-`F34`, `S01`-`S46`, `D01`-`D112`, on 2026-10-31 `H01`-`H119`,
`M01`-`M07` on 2026-10-29, on 2026-11-01 `I101`-`I1195`, `J101`-`J222`, `K01`-`K39` and
`M101`-`M1003`, on 2026-11-02 `N101`-`N2042`; NCBI was not run on this tree).

- `db/`: databases of that probe. `p5` (version 5, `-parse_seqids`) and `p4` (version 4)
  are proteins with every Seq-id type, `n5` / `n4` the nucleotides, `p5n` / `p4n` / `n5n` /
  `n4n` the same without parsed ids, `pL` / `nL` version 4 databases written by LOSAT
  (`losat makedb`, now `losat makeblastdb`), `pm5` a protein database with the seg masks of
  `makeblastdb -mask_data` (column files `.paa .pab .pac`), `palias` an alias file of `p5`
  and `pm5` (blastdb_aliastool), `d5` a version 5 database whose accessions `XP_800001` and
  `XP_800002` repeat on OIDs 1 and 256, 9 and 10 (the order of LMDB's duplicates), `e5` one
  whose `XP_900001` is on 600 sequences and `XP_900002` on 300, `d5b` the same FASTA as `d5`,
  `h_dd`, `h_p5p5n`, `h_p5np5` alias files over two LMDB files, `nm5` and `nm4` nucleotides
  with the dust masks 11 and 12, `m/nms4` the FASTA of `nm4` with masks 12 on one of its
  sequences, `m/pmg` proteins with GIs and seg masks whose `.ppi` and `.ppd` are copies of
  their `.pni` and `.pnd` (an index of PIGs). Version 4 with parsed ids
  besides `p4` and `n4`: `d4` (the FASTA of `d5`: the OIDs of a key in the text order of the
  ISAM, `10` before `9`), `i4` (pdb ids of four digits that read as GIs, an EMBL and a GenBank
  id in one definition line, CDD, GNOMON, PIG and trace ids as general ids, versions of one
  accession, locus names), `it4` (nucleotide trace, SRA and TSA ids), `dv4` the FASTA of `d5`
  in three volumes, with the alias files `i_dvr` (two volumes in reverse order) and `i_dv1`
  (the second volume). Kept: the files LOSAT reads, the LMDB file of each version 5 database,
  the ISAM files of GIs of `p5` and `n5`, and the string and GI ISAM files of the version 4
  databases (whose presence makes a version 4 database one with parsed ids).
- `in/`: the `-entry_batch` files of the probe (`i_p_all.txt`, `i_p_safe.txt`: every spelling
  of an id of `p4`, the second without the ids that stop the batch; `i_i_all.txt`,
  `i_i_range.txt`: those of `i4`; `n_p.txt`: accessions of `p5` and `pm5` with one that neither
  holds).
- `expected/ID.stdout`, `expected/ID.stderr`: NCBI's bytes for case `ID`, with `{R}` for the
  directory the probe ran in (here `tests/blastdbcmd`, the current directory of the runs) and
  `{T}` for the repository's `LOSAT/tests/` (the probe's copy of the seqdb fixtures).
  `B07.logfile`, `B08.logfile` and `D72.logfile` are the files NCBI's `-logfile` wrote
  (`B08.logfile` after two runs, B08 and B38, into one file), `D70.out` and `D84.out` those
  of `-out`. The stdout of `D82`, `D87`, `D97` and `D101` is not NCBI's but the correct
  output of an approved difference (line BDC DECISIONS 12: `-long_seqids` with `-range` on
  the minus strand, where NCBI writes the forward letters under its `:c` header): NCBI's
  header with the letters NCBI writes for the same range with short ids
  (`dropin-survey-20261029/BDC/probe/make_expected.py`). `H09` the same difference on a
  parsed database. `H107`, `H108` and `H109` are neither: NCBI stops on the accession of 600
  sequences (`PD-LOSAT-NCBI-DEFECTS` exception 6), and the expected bytes are the OIDs in the
  order NCBI gives a key of fewer sequences (`H110`). The stdout of the cases of the approved
  differences 3 to 6 in the table of `blastdbcmd_masks.rs` is not NCBI's either: NCBI's
  headers with the letters of NCBI's `-outfmt %s` run of the same words, or, for two `%m`, the
  runs with one `%m` each (`make_expected.py`, `PD-LOSAT-BLASTDBCMD-DEVIATIONS`).
- `j/`: whole copies (every file but the build logs; no line-end conversion, `j/.gitattributes`:
  the sizes are in `-list %U` and `-metadata`) of probe J's databases: `j/lst` holds `n4`,
  `n5`, `p4`, `p5`, `p5n`, `dv4` (three volumes and their alias file) and `pal`
  (blastdb_aliastool's alias of `p5` and `p5n`); `j/jd` copies of `p4` whose index file has
  another date (`am12`, `pm12`, `d01`, `d1`, `feb30`, `bad`) or a title with characters JSON
  escapes (`esc`); `j/jtax5/taxonomy4blast.sqlite3` a taxonomy over the taxids of `p5` (562,
  9606, 10090, 10116 under 9604, 10066 and 1, with 3000 more taxids so that both b-trees have
  interior pages; probe K, written by Python's `sqlite3`). The expected lines of `J160`,
  `J161`, `J163` and `J164` are NCBI's lines of these databases (the probe's directory held
  two more, `nal2` and `s18mv`). `J165`, `J166`, `J167` and `J176` are NCBI's lists of the
  probe's `db/j/lsx` (empty files, links, a directory named `d.pin`), which the test builds
  again on Unix. `K13` (`-taxids 1`, which never ends in NCBI: taxid 1 is its own parent) is
  the approved difference of DI-8 DECISIONS 10; its expected bytes are NCBI's `K25`, the same
  run with the explicit list of every taxid that taxid 1 reaches.

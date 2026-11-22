# makembindex fixtures (line RG item D5a)

`cases.tsv` lists the cases; `ncbi/` holds the official NCBI BLAST+ 2.17.0 `makembindex` output of each: its standard
output, standard error and exit status, and in `<name>.files` the path, size and FNV-1a 64 digest of each file the run
made. `tests/makembindex_fixtures.rs` runs every row in its own copy of `in/` and compares.

`in/`: `s1.fa`, `m.fa`, `empty.fa` and `crlf.fa` written by hand and by the probe of line RG item D5a; `ids.fa` and `mv.fa`
(the first 40 records of the generated `many.fa`) from the probe's generator. The databases were made by NCBI BLAST+
2.17.0: `s1db/s1` of `s1.fa`; `mdb/m` of `m.fa` with the masks of `dustmasker -outfmt maskinfo_asn1_bin` (algorithm id
11); `idsdb/ids` of `ids.fa` with `-parse_seqids`; `mvdb/mv` of `mv.fa` with `-max_file_sz 8KB` (three volumes).
`adb/zz` is a second database of `s1.fa`, and `ext/b.nin` and `ext/b-2` two more (`-blastdb_version 4`, as NCBI's
LMDB check refuses a version 5 volume named `b.nin`). `al.nal` is an alias of `mvdb/mv`, `s1db/s1` and `idsdb/ids`;
`al2.nal` of `s1db/s1` and `adb/zz`, whose volumes NCBI indexes in the order of their base names (`s1` before `zz`),
not of their paths; `al3.nal` of `ext/b.nin` and `ext/b-2`, where the base name drops the `.nin` (`b` before `b-2`).
`ol.nal` is an alias of `mvdb/mv` with `OIDLIST ol.msk` (OIDs 0, 3, 5, 6, 8, 10 and 12 of each volume), `olmissl.nal` and
`olmiss.nal` the same with a missing OID list file, with and without `NSEQ` and `LENGTH` (line RG item D1b-1, review B):
a new style index reads each volume whole, and an old style one every OID of the alias file, as NCBI does; NCBI's
open of the alias file scans the list without `NSEQ` and `LENGTH`, and the old style's first read builds it.

The `-idmap` titles (line RG item D5a-t): `titles.fa` (titles with HTML entities, end punctuation, leading spaces, a tab,
`TPA:`), `titlesdb/t` of it; `titlespdb/p` with `-parse_seqids` (`tpg`, `pat`, `lcl` and `gb` ids with and without
titles); `txdb/tx` with `-parse_seqids -taxid 9606`; and `ttdb/tt` of `tt.fa` with `-parse_seqids -taxid_map
tt_taxid.txt` (the taxids of `taxt5/`, a copy of `../seqdb/taxdb/t5`: a lower case name, a name with spaces around it, one
with entities, one starting `TPA: `, an empty name, and a taxid it does not hold). The rows with `BLASTDB=taxt5` find
that taxonomy database.

`gen/` is made by the test: `big.fa` (one record of 5,100,000 letters with lower case runs, which the index splits into
chunks, and one of 20,000) and `many.fa` (600 records) from the generator of `makembindex_fixtures.rs`, the databases
`gen/bigdb/big` and `gen/manydb/many` of them, and `gen/mbigdb/big` with the lower case masks of `big.fa`
(`convert2blastmask -masking_algorithm other`, `makeblastdb -mask_id lower`: algorithm id 100). `spill/spill.fa` is
made by the test too: 2,200,000 records of `A`, then 1,900,000, 100,000 and 1,000,000 letters of the generator (row
`old_fa_stride1_rollback`: with `-stride 1` the offsets of the first local id reach into the local id's bits, and the
volume that rolls back cuts its lists where NCBI cuts them). The freeze made the same files with NCBI's programs
(`remaining-gaps-20261112/D/d5/fixture/freeze.py` and `gen_big.py`, probe scripts kept with the session's files, not
in this repository); an index reads only the databases' sequences and masks.

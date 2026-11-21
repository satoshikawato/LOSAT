# makembindex fixtures (line RG item D5a)

`cases.tsv` lists the cases; `ncbi/` holds the official NCBI BLAST+ 2.17.0 `makembindex` output of each: its standard
output, standard error and exit status, and in `<name>.files` the path, size and FNV-1a 64 digest of each file the run
made. `tests/makembindex_fixtures.rs` runs every row in its own copy of `in/` and compares.

`in/`: `s1.fa`, `m.fa`, `empty.fa` and `crlf.fa` written by hand and by the probe of line RG item D5a; `ids.fa` and `mv.fa`
(the first 40 records of the generated `many.fa`) from the probe's generator. The databases were made by NCBI BLAST+
2.17.0: `s1db/s1` of `s1.fa`; `mdb/m` of `m.fa` with the masks of `dustmasker -outfmt maskinfo_asn1_bin` (algorithm id
11); `idsdb/ids` of `ids.fa` with `-parse_seqids`; `mvdb/mv` of `mv.fa` with `-max_file_sz 8KB` (three volumes).
`adb/zz` is a second database of `s1.fa`. `al.nal` is an alias of `mvdb/mv`, `s1db/s1` and `idsdb/ids`; `al2.nal`
of `s1db/s1` and `adb/zz`, whose volumes NCBI indexes in the order of their base names (`s1` before `zz`), not of their
paths.

`gen/` is made by the test: `big.fa` (one record of 5,100,000 letters with lower case runs, which the index splits into
chunks, and one of 20,000) and `many.fa` (600 records) from the generator of `makembindex_fixtures.rs`, the databases
`gen/bigdb/big` and `gen/manydb/many` of them, and `gen/mbigdb/big` with the lower case masks of `big.fa`
(`convert2blastmask -masking_algorithm other`, `makeblastdb -mask_id lower`: algorithm id 100). The freeze made the
same files with NCBI's programs (`remaining-gaps-20261112/D/d5/fixture/freeze.py` and `gen_big.py`, probe scripts kept
with the session's files, not in this repository); an index reads only the databases' sequences and masks.

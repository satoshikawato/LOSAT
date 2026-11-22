# A single-volume `makeblastdb -gi_mask` database whose alias file names its volume (PD-LOSAT-NCBI-DEFECTS exception 22)

Line NP item M4c-1, 2026-11-20. NCBI BLAST+ 2.17.0 (the official Linux x86-64 build, `BLAST_USAGE_REPORT=false`),
`probe.sh` in this directory (inputs: `nucl.fa`, `prot.fa` and their dustmasker/segmasker masks of
`LOSAT/tests/seqdb/make_gi_mask_fixtures.py`). Paths of the run directory are shown as `<W>`.

## What NCBI writes

`makeblastdb -in nucl.fa -dbtype nucl -parse_seqids -title t -out db -mask_data dust_nucl.asnb -gi_mask
-gi_mask_name gm1` writes every file and ends with exit status 2:

```text
BLAST Database error: Illegal configuration: DB alias files are mutually recursive.
```

Its alias file `db.nal` names the database itself by its absolute path:

```text
#
# Alias file created: Nov 20, 2026  7:48 PM
#
TITLE t
DBLIST <W>/v5n/db
MASKLIST gm1 
```

The same for version 4 (`-blastdb_version 4`), for protein (`-dbtype prot`, segmasker masks, version 5 and 4), for
`-out sub/db` and without `-gi_mask_name`. No NCBI program opens the database: `blastdbcmd -db db -info` and
`blastdbcmd -db db -metadata` print the line above, exit status 2.

## With `DBLIST db`

Only that line changed to the database's file name (`DBLIST db`; `DBLIST db` also for `-out sub/db`), NCBI's
programs open each of the six databases, from the database's directory (`-db db`, `-db sub/db`, `-db <W>/v5n/db`):
every run below exits with status 0.

- `blastdbcmd -info`: the database, the algorithm `0 dust window=64; level=20; linker=1` (`0 seg window=12;
  locut=2.2; hicut=2.5` for protein) and the volume `<W>/v5n/db`.
- `blastdbcmd -entry all -outfmt "%o %a %m0"`: the masks of the sequences with a GI; `none` for the sequences
  without one (`AB000005.1`, `seq9`; `AAA00003.1`), whose masks makeblastdb drops. These lines equal `%m11` of the
  database built with the mask column (`-mask_data` without `-gi_mask`) but for those sequences, and equal
  `%m0` of a database of several volumes built with `-gi_mask` (`-max_file_sz 400`, which NCBI writes readable),
  for version 5 and 4, `-out sub/db` and the default name alike; protein: `%m21` of the column likewise.
- `blastdbcmd -entry all -mask_sequence_with 0`: the masked letters in lower case.
- `blastn -query nucl.fa -db db -outfmt 6 -db_soft_mask 0` equals `-db_soft_mask 11` of the column database;
  `-db_hard_mask 0` equals `-db_hard_mask 11` but for the two sequences without a GI (unmasked there: one row of
  400 letters for each, where the column gives two); `blastp -query prot.fa -db db -outfmt 6 -db_soft_mask 0`
  equals `-db_soft_mask 21` of the column database.
- `blastdbcmd -metadata` gives the metadata file that makeblastdb would have written (the same
  `CSeqDB::GetDBMetaData` and JSON writer as makeblastdb's own step, makeblastdb.cpp:1296-1317,
  blastdbcmd.cpp:1206-1221). It lists `db.nal` twice, as the alias file and as a file of the volume, and counts its
  bytes twice:

```text
{
  "version": "1.2",
  "dbname": "db",
  "dbtype": "Nucleotide",
  "db-version": 5,
  "description": "t",
  "number-of-letters": 5020,
  "number-of-sequences": 12,
  "last-updated": "2026-11-20T19:48:00",
  "number-of-volumes": 1,
  "bytes-total": 52376,
  "bytes-to-cache": 1495,
  "files": [
    "db.nal",
    "db.nal",
    "db.ndb",
    "db.nhr",
    "db.nin",
    "db.nnd",
    "db.nni",
    "db.nog",
    "db.nos",
    "db.not",
    "db.nsq",
    "db.ntf",
    "db.nto"
  ]
}
```

From another directory (`-db <W>/v5n/db`) NCBI finds the database but not the GI-mask files, which it looks for
from the current directory and `BLASTDB` (seqdbgimask.cpp, `SeqDB_ResolveDbPath`): `blastdbcmd -info` prints its
head and `BLAST Database error: Could not open gi-mask index files.` (exit status 2), `%m0` gives `Error:
[blastdbcmd] Could not open gi-mask index files.` (exit status 1), `blastn -db_soft_mask 0` exit status 2. The
database of several volumes gives the same three results there. With `BLASTDB=<W>/v5n`, `-db db` reads the masks
from there (exit status 0).

## Rebuilds

`makeblastdb` into the same `-out` over the changed database removes it (`Deleted existing Nucleotide BLAST database
named <W>/re2/db`) together with its alias file: without `-gi_mask` the new build ends with exit status 0 and no
alias file; with `-gi_mask` NCBI writes its self-naming alias file again and ends with exit status 2.

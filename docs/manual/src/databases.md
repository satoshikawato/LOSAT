# BLAST databases

A BLAST database is a set of files that holds sequences in a form the search programs read fast, with their ids, titles, and optionally tax ids and masks. LOSAT uses the database format of NCBI BLAST+. There is no LOSAT-only format: NCBI's programs read the databases that `losat makeblastdb` writes, and LOSAT reads the databases that NCBI's `makeblastdb` writes and the ones NCBI distributes.

You can search without a database. `-query` with `-subject` compares two FASTA files directly. Use a database (`-db`) when you search the same sequences many times, when the sequences are an NCBI database such as swissprot or RefSeq, or when you want ids, taxonomy and masks to come with the sequences.

A search with `-db` gives the output of NCBI BLAST+ 2.17.0 run with `-num_threads 1`, at every `-num_threads` of LOSAT. The same sequences searched with `-subject` and with `-db` differ, as they do in NCBI BLAST+: statistics see the whole database, the `Database:` lines differ, and a database stores a hyphen in a nucleotide sequence as `N` where `-subject` drops it.

## Make a database

```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db
losat blastp -query query.faa -db proteins_db -outfmt 6
```

The database name is the `-out` name without a file extension. See [makeblastdb](programs/makeblastdb.md) for ids, tax ids, masks and volumes.

## Get an NCBI database

NCBI publishes its prebuilt databases at <https://ftp.ncbi.nlm.nih.gov/blast/db/>. Download the archive of a database and unpack it into a directory. This example fetches swissprot:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar -xzf swissprot.tar.gz
```

The taxonomy files that the name fields need (see Taxonomy below) come the same way:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz
```

Large databases come as several archives (`nt.000.tar.gz`, `nt.001.tar.gz`, ...); unpack all of them into the same directory. NCBI's `update_blastdb.pl`, which is part of NCBI BLAST+ and not of LOSAT, downloads and updates them for you.

Then search the database by name. Run the command in the directory that holds the files, or put that directory in `BLASTDB` (see below):

```bash
losat blastp -query query.faa -db swissprot -outfmt 6
```

NCBI's databases such as swissprot, pdbaa, 16S_ribosomal_RNA and the RefSeq nucleotide databases store parsed sequence ids. LOSAT reads them, searches them with every search program and writes their ids as NCBI does.

## Format versions 4 and 5

| | Version 5 | Version 4 |
| --- | --- | --- |
| Written by `losat makeblastdb` | Default | `-blastdb_version 4` |
| Read by LOSAT | Yes | Yes |
| Id and tax id lookup | LMDB files (`.pdb`/`.ndb`, `.ptf`/`.ntf`, `.pot`/`.not`, `.pto`/`.nto`) | ISAM files (`.psi`/`.psd`, `.nsi`/`.nsd`) |
| Current NCBI downloads | Yes | Older databases |

You do not choose a version to search: LOSAT detects it. `blastdbcmd -db DB -info` prints it as `BLASTDB Version`.

## Volumes and alias files

A large database is split into volumes (`name.00`, `name.01`, ...) and an alias file (`.pal` for protein, `.nal` for nucleotide) that lists them. LOSAT searches such a database like any other: give the alias name to `-db`. `losat makeblastdb` writes volumes when a file would pass `-max_file_sz`.

An alias file can also name other databases (`DBLIST`) or select part of a database. LOSAT reads these alias keys as NCBI does: `TITLE`, `DBLIST`, `NSEQ`, `LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`, `FIRST_OID`, `LAST_OID`, `OIDLIST` (an OID list file, such as the `.pxm` and `.nxm` files of `makeblastdb -oid_masks`, applied to each volume below the alias file) and `OID_MASK_TYPE` (the alias files of `blastdb_aliastool -oid_masks 1` write `1`: the definition lines with a predicted RefSeq accession, such as `XP_` or `XM_`, are hidden, and a sequence left with none is searched but not reported, as in NCBI; an RPS profile left with none in a database of several volumes is refused where it would be reported, see [What LOSAT refuses](refusals.md)), `GILIST`, `TILIST` and, over a version 4 database, `SEQIDLIST` (a list file of GIs, TIs or sequence ids, such as `blastdb_aliastool` writes, named relative to the alias file: only the listed sequences are searched, and the definition lines whose ids no list of the alias files names are hidden; a sequence left with none is searched but not reported, as in NCBI. A list in an alias file below one that also names the same volume directly is not read, as in NCBI. NCBI BLAST+ 2.17.0 finds no sequence through a TI list, and LOSAT finds none either). A `TAXIDLIST` over a version 4 database is read nowhere, as in NCBI. `VDBLIST` and `MAX_SEQ_LENGTH` are ignored, as in NCBI. An alias file with `SEQIDLIST` or `TAXIDLIST` over a version 5 database, or `MEMB_BIT`, is refused (see [What LOSAT refuses](refusals.md)). `MASKLIST` names the GI-indexed masks that [`makeblastdb -gi_mask`](programs/makeblastdb.md#gi-indexed-masks) writes, and LOSAT reads them as NCBI does: when `-db` names the database alone and its own alias file has `MASKLIST`, its masks (`blastdbcmd -info`, `-outfmt %m`, `-mask_sequence_with`, and `-db_soft_mask` and `-db_hard_mask` of `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`; `makeblastdb -input_type blastdb` copies them into a mask column) are the GI-indexed masks, with ids 0, 1, ... in the order of `MASKLIST`. The mask files are looked for from the current directory and the `BLASTDB` path, not from the database's directory. In an alias file over that database, or with several names in `-db`, `MASKLIST` is passed over, as NCBI passes it over.

`-db` can name several databases: `-db "swissprot pdbaa"` (the names in double quotes, separated by spaces). They are searched as one database. A name with spaces goes in double quotes inside the value.

## Where LOSAT looks for a database

For a `-db` name without a directory, LOSAT looks, in this order, in

1. the current directory,
2. the directories in the environment variable `BLASTDB` (separated by `:`, or `;` on Windows),
3. the `BLASTDB` entry of `NCBI_CONFIG__BLAST__BLASTDB` or of NCBI's configuration files (`<program>.ini`, then `.ncbirc`, found in the current directory, your home directory and the other places NCBI BLAST+ reads them).

A `-db` name with a directory (`../data/swissprot`, `/data/blast/nt`) is looked up at that path. `losat blastdbcmd -show_blastdb_search_path` prints the directories LOSAT searches for a name.

```ini
# .ncbirc
[BLAST]
BLASTDB=/data/blast
```

The same search path finds the `taxdb` files, the files named by `-seqidlist` and `-taxidlist`, and `taxonomy4blast.sqlite3`. A database that is not found stops the search with NCBI's message, naming the search path.

## Taxonomy

A database can store a tax id for each sequence. `makeblastdb` writes them with `-taxid` or `-taxid_map`; NCBI's databases carry them already.

| You want | LOSAT needs |
| --- | --- |
| The tax id of a hit (`staxid`, `staxids`, outfmt 13 to 16 ids) | The tax ids in the database. |
| Names of a tax id (`ssciname`, `scomname`, `sblastname`, `sskingdom` and their plural forms, `sciname` in outfmt 13 to 16, `blastdbcmd` `%S`, `%L`, `%B`, `%K`) | NCBI's `taxdb.bti` and `taxdb.btd`, in one directory of the search path above (the current directory or `BLASTDB`; not the directory of the `-db` name). |
| The organism report (`-outfmt 18`) | `taxdb.bti` and `taxdb.btd`. Without them the report has its headers only. |
| `-taxids`, `-negative_taxids`, `-taxidlist`, `-negative_taxidlist` | The tax ids in the database. To include the descendants of a tax id, NCBI's `taxonomy4blast.sqlite3`. Without it LOSAT writes NCBI's message (`The -taxids command line option requires additional data files`) and, as NCBI BLAST+ does, does not limit the search by tax id at all: hits of every tax id are reported. |

Without `taxdb`, name fields print `N/A` and LOSAT writes NCBI's warning (`Taxonomy name lookup from taxid requires installation of taxdb database ...`). Download `taxdb.tar.gz` as shown in Get an NCBI database above and unpack it into a directory of `BLASTDB`.

## Masks

A database can hold masks that mark low-complexity or repeat regions, made by `dustmasker`, `segmasker`, `windowmasker` or `convert2blastmask` and added with `makeblastdb -mask_data`. Masks are used in these ways:

- `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` `-db_soft_mask ID` and `-db_hard_mask ID` search with the masks of one algorithm: words inside a soft mask do not start hits, and a hard mask cuts the masked letters out of the search. `ID` is the numeric id or the name shown by `blastdbcmd -info` under "Available filtering algorithms". The alignments in format 0 to 4 of `blastn`, `blastp` and `blastx` show the masked letters of the database sequence in lower case; `tblastn` and `tblastx` show the translated letters in upper case, as NCBI does. `tblastn` and `tblastx` read each mask in every reading frame at the positions divided by 3, and with `-db_hard_mask` they score the masked bases as ambiguous after the first search, as NCBI does (see [What LOSAT refuses](refusals.md)).
- `blastdbcmd -mask_sequence_with ID` writes masked letters in lower case, and `-outfmt "%m<ID>"` prints the mask ranges.

See [dustmasker](programs/dustmasker.md), [segmasker](programs/segmasker.md), [windowmasker](programs/windowmasker.md), [convert2blastmask](programs/convert2blastmask.md) and [makeblastdb](programs/makeblastdb.md).

## Which programs read which databases

| Program | Database it reads |
| --- | --- |
| [`blastn`](programs/blastn.md), [`tblastx`](programs/tblastx.md), [`tblastn`](programs/tblastn.md) | Nucleotide |
| [`blastp`](programs/blastp.md), [`blastx`](programs/blastx.md), [`psiblast`](programs/psiblast.md) | Protein |
| [`rpsblast`](programs/rpsblast.md), [`rpstblastn`](programs/rpstblastn.md) | A profile database, such as NCBI's CDD (files `.rps`, `.loo`, `.aux`, `.freq`; one or several volumes): NCBI's, or one built by [`makeprofiledb`](programs/makeprofiledb.md). |
| [`blastdbcmd`](programs/blastdbcmd.md) | Either; `-dbtype` chooses when a name exists for both. |
| [`dustmasker`](programs/dustmasker.md) (`-infmt blastdb`) | Nucleotide |
| [`windowmasker`](programs/windowmasker.md) (`-infmt blastdb`) | Nucleotide |
| [`segmasker`](programs/segmasker.md) (`-infmt blastdb`) | Protein |
| [`rpsbproc`](programs/rpsbproc.md) | None: it reads the archives that `rpsblast` and `rpstblastn` write, and the CDD data files. |

`-subject` works with `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` and `psiblast`; `rpsblast` and `rpstblastn` take only `-db`, as in NCBI.

A database of the wrong molecule type is reported as not found, with NCBI's message `No alias or index file found for nucleotide database [NAME]` (or `protein`).

## Search options that use the database

- `-dbsize` and `-searchsp` set the database size used for statistics (with `-db` and `-subject`).
- `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist`, `-negative_seqidlist` and `-no_taxid_expansion` select sequences, with `-db`.
- `-gilist`, `-negative_gilist`, `-ipglist` and `-negative_ipglist` need NCBI's GI and IPG databases and are not provided.

## Differences from NCBI BLAST+

LOSAT reads a database as NCBI does, with a few cases that NCBI cannot finish or that depend on memory contents; see [Differences from NCBI BLAST+](differences.md) and the reasons in [the decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md) and [the decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

## See also

- [makeblastdb](programs/makeblastdb.md)
- [blastdbcmd](programs/blastdbcmd.md)
- [Output formats](output-formats.md)
- [What LOSAT refuses](refusals.md)

# makeblastdb

`losat makeblastdb` builds a BLAST database from FASTA files, ASN.1 text or binary `Seq-entry`/`Bioseq-set` files (`-input_type asn1_txt` or `asn1_bin`), or existing BLAST databases (`-input_type blastdb`). Build a database when you search the same sequences many times, when you need NCBI-style sequence ids, taxonomy ids or masks in the search, or when you want to read sequences back with [`blastdbcmd`](blastdbcmd.md).

The files it writes are the files of NCBI's `makeblastdb`. NCBI BLAST+ and LOSAT read each other's databases. The default is format version 5, as in `makeblastdb` 2.17.0. `-blastdb_version 4` writes the older layout.

## Synopsis

```bash
losat makeblastdb -in FASTA -dbtype nucl|prot [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `makeblastdb -in FASTA -dbtype nucl|prot [options]`.

## Input

- `-in` names the FASTA file. `-in -` (the default) reads the standard input.
- Name several files in one `-in` value, separated by spaces: `-in "a.faa b.faa"`. Their sequences go into one database.
- A name ending in `.gz`, `.bz2` or `.zst` is read through gzip, bzip2 or zstd while the database is built. Files made of several gzip members, bzip2 streams or zstd frames (the output of `bgzip`, `pbzip2` and `pzstd`) are read to the end.
- `-dbtype` is required: `nucl` or `prot`.
- `-input_type` is `fasta` (the default), `asn1_txt`, `asn1_bin`, or `blastdb` (see [Copying databases](#copying-databases)).

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-in` | `-` | FASTA or ASN.1 input file or files, or the databases to copy with `-input_type blastdb`. |
| `-dbtype` | required | `nucl` or `prot`. |
| `-out` | the `-in` name | Name of the database to write. |
| `-title` | the `-in` name (with `-input_type blastdb`, the input database's title) | Title stored in the database and shown by the search programs. |
| `-blastdb_version` | `5` | Format version: `5` or `4`. |
| `-parse_seqids` | off | Parse the ids of the definition lines and store them, so that `-entry`, `-seqidlist` and similar lookups find sequences by accession. |
| `-taxid` | none | One tax id for every sequence. |
| `-taxid_map` | none | Text file of lines `SEQUENCE_ID TAXID`. Needs `-parse_seqids`. `-taxid_map -` reads the map from the standard input. |
| `-hash_index` | off | Also write the sequence hash index. |
| `-mask_data` | none | Masks from `dustmasker`, `segmasker`, `windowmasker` or `convert2blastmask` (comma-separated files, ASN.1 text or binary). |
| `-mask_id` | the name in the mask file | Comma-separated names for the masking algorithms. A name cannot be a number. |
| `-mask_desc` | none | Comma-separated free-text descriptions of the algorithms. |
| `-gi_mask` | off | Write the masks as GI-indexed mask files instead of the mask column (see [GI-indexed masks](#gi-indexed-masks)). Needs `-parse_seqids`; not with `-mask_id`. |
| `-gi_mask_name` | the `-mask_data` names | Comma-separated names of the GI-indexed mask files, one for each mask file. Needs `-mask_data` and `-gi_mask`. |
| `-max_file_sz` | `3000000000` | Largest size of one database file, in bytes. A larger database is split into volumes. |
| `-logfile` | standard output and error | Write the log and messages to this file (appended). |
| `-metadata_output_prefix` | none | Path written before each file name in the JSON metadata (`.pjs`, `.njs`), for a database that will be read from another directory. The path separator (`/`, `\` on Windows) is added when the value does not end with one. |
| `-oid_masks` | none | `1` also writes an OID mask file per volume (`.pxm`, `.nxm`) that marks the sequences whose ids are all predicted RefSeq accessions (`XP_`, `XM_`, `XR_`, ...), which NCBI BLAST+ and LOSAT read for an alias file with `OID_MASK_TYPE 1`. |

For every option: `losat makeblastdb -help`.

### Sequence ids with `-parse_seqids`

Without `-parse_seqids` each sequence gets a local id (`BL_ORD_ID:0`, `BL_ORD_ID:1`, ...) and you find it by its position. With `-parse_seqids` LOSAT stores the ids of the definition lines (`sp|P69905.2|...`, `NP_000002.1`, `lcl|name`). Use it whenever you want to fetch sequences by id with `blastdbcmd`, restrict a search with `-seqidlist`, or use `-taxid_map`.

Two sequences with the same id stop the build with `makeblastdb`'s message (`Duplicate seq_ids are found`).

### Taxonomy ids

`-taxid 9606` gives every sequence tax id 9606. `-taxid_map map.txt` gives each sequence its own tax id. The map has one line per sequence, the sequence id and the tax id, separated by white space:

```text
sZ_nz_p 9606
sI_nz_p 562
```

Tax ids are what `-taxids` and `-taxidlist` filters, the tabular `staxid` field and [outfmt 18](../output-formats.md) use. To print scientific names you also need NCBI's `taxdb` files (see [BLAST databases](../databases.md)).

### Masks

To search a database with `-db_soft_mask` or `-db_hard_mask`, or to read masked FASTA with `blastdbcmd -mask_sequence_with`, add masks when you build it:

1. Make a mask file with [`dustmasker`](dustmasker.md) (nucleotides), [`segmasker`](segmasker.md) (proteins), [`windowmasker`](windowmasker.md) (repeats of nucleotides, with a counts file) or [`convert2blastmask`](convert2blastmask.md) (FASTA that another tool has masked). Use `-outfmt maskinfo_asn1_bin` or `maskinfo_asn1_text`, and `-parse_seqids` if the database will have parsed ids.
2. Give it to `makeblastdb` with `-mask_data`. Add `-mask_id` to name the algorithm and `-mask_desc` to describe it.
3. Look up the numeric algorithm id with `blastdbcmd -db DB -info` (the table "Available filtering algorithms").

A mask finds its sequence by the same id as in NCBI BLAST+. Build with `-parse_seqids` when the mask file names accessions. Several masks go in one comma-separated list; give `-mask_id` as many names as there are mask files (a different count stops with `the size of mask_id does not match that of mask_data`).

A mask file that is damaged stops with makeblastdb's own error and exit status. A text in it with a byte outside the visible ASCII range is stored with `#` in the byte's place, and the run reports it, as makeblastdb does.

Locations nested inside each other (`mix` and `equiv`) are read at any depth. NCBI's makeblastdb crashes on nests of about 6,500 to 6,900 levels and more with the usual 8 MiB stack limit; LOSAT writes the database that it writes with no stack limit.

### GI-indexed masks

`-gi_mask` writes the masks of each mask file into six files of their own instead of the mask column: `NAME.gmd` and `NAME.gnd` (the masked ranges), `NAME.gmo` and `NAME.gno` (the GIs of the masked sequences) and `NAME.gmi` and `NAME.gni` (their index). `NAME` is the mask file's entry of `-gi_mask_name`, or the mask file's name as you typed it; it is a path from the current directory, not from the database. The database gets an alias file (`.nal` or `.pal`) whose `MASKLIST` line names these files. Only sequences with a GI get masks: the masks of other sequences are left out without a message, and a mask file that reached no sequence with a GI writes no files and is left out of `MASKLIST`.

The files, messages and exit status are those of NCBI's `makeblastdb` 2.17.0, with one difference:

- A database in one volume opens. NCBI's `makeblastdb` writes an alias file that names the database by its full path, which NCBI's programs take for another alias file: no program can open that database, and the build ends with `BLAST Database error: Illegal configuration: DB alias files are mutually recursive.` and exit status 2. LOSAT's alias file names the volume by its file name (`DBLIST db` for `-out db`), as the alias file of a database in several volumes names its volumes. The build then ends as one of several volumes does, with the metadata file (`.njs`, `.pjs`) and exit status 0, and NCBI's programs read the database and its masks. Every other file is NCBI's.
- An alias file that NCBI's `makeblastdb` left from such a build cannot be read when the old database is removed, by NCBI or by LOSAT. A new build with `-gi_mask` writes its own alias file over it. A build without `-gi_mask` keeps it and ends with the error above, as NCBI's does: delete the `.nal` or `.pal` file first.

NCBI's programs and LOSAT's look for the GI-mask files from the current directory and the `BLASTDB` path, not from the database's directory, so they find them when they run where the build ran (or with `BLASTDB` set to that directory, or with an absolute `-gi_mask_name`). From elsewhere they stop with `Could not open gi-mask index files.`

LOSAT's `blastdbcmd`, `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` read the GI-indexed masks as NCBI's do (`-info`, `-outfmt %m`, `-mask_sequence_with`, `-db_soft_mask`, `-db_hard_mask`), every search program searches such a database, `-input_type blastdb` copies its masks (see [Copying databases](#copying-databases)), and [`makembindex -db_mask`](makembindex.md) leaves them out of an index. A sequence without a GI has no masks of its own: like NCBI's programs, LOSAT's give it the masks of GI 4294967295 when the database has that GI. See [Databases](../databases.md) for when `MASKLIST` applies.

### Hash index

`-hash_index` writes the sequence hash index (`.phi`/`.phd`, `.nhi`/`.nhd`). `blastdbcmd -outfmt %h` prints the hash of each sequence.

### Volumes

`-max_file_sz` is the size, in bytes, at which a database is cut into volumes. A database in several volumes has an alias file (`.pal` or `.nal`) that lists them. The search programs and `blastdbcmd` take the alias name like any database name.

### Copying databases

`-input_type blastdb` writes the sequences of existing databases into a new one. Use it to merge databases, to change the format version (`-blastdb_version`), to cut the volumes at another size, or to add tax ids, masks or a hash index without the original FASTA files.

- `-in` names one database or several, separated by spaces, as `-db` does for the search programs. Alias names work. `-out` is required and must not be the `-in` value.
- Without `-title` the new database takes the title of the input database. With several input databases, give `-title`.
- The sequences keep their definition lines (ids, tax ids, membership and leaf tax ids), their PIGs, the masks of every algorithm and the data of other columns. The ids are always indexed, as with `-parse_seqids`.
- An alias name keeps the sequences its OID ranges and lists keep. Under `OID_MASK_TYPE 1` the definition lines with a predicted RefSeq accession are left out, as in NCBI; a sequence left with no line is refused (NCBI's `makeblastdb` stops with a null-pointer report).
- `-taxid`, `-taxid_map`, `-mask_data`, `-gi_mask`, `-hash_index` and `-oid_masks` work as for FASTA input. `-taxid_map` and `-gi_mask` still need `-parse_seqids` on the command line.
- An input database built with `-gi_mask`, named alone in `-in`, gives its GI-indexed masks as its algorithms: the copy stores them in its mask column, with the algorithm ids that NCBI's `makeblastdb` gives them (11, 12, ... for dust). The mask files are looked for from the current directory and `BLASTDB`, as by every reader: from elsewhere the copy stops with `Could not open gi-mask index files.` With several names in `-in` the GI-indexed masks are passed over and the copy has none, as in NCBI.
- With `-gi_mask` the masks of the input databases go to GI-indexed mask files too, one set for each of their algorithms, named after the algorithm and its id in the input (`dust11` and so on), as NCBI's `makeblastdb` names them. A sequence's masks go under the GIs of all its definition lines.
- `-randomize` (hidden in `makeblastdb`'s `-help`, and in LOSAT's; `-help-full` lists it) copies the sequences in a random order, which depends on the second of the clock. LOSAT reads the clock when the build starts; `makeblastdb` reads it when the copy starts, after reading the input's definition lines. With the same second the order is the one `makeblastdb` 2.17.0 gives on Linux. `makeblastdb` on Windows and macOS shuffles with other C++ libraries, whose order has not been compared. With another `-input_type`, `-randomize` stops with `makeblastdb`'s error.
- With `-title` or several names, an input name that is not a database is reported (`Unable to open input NAME as BLAST db`), and the build stops with the database error. A single name without `-title` that is not a database stops the build earlier, when its title is read, with the database error alone. Both are as in `makeblastdb`.

## Output

`makeblastdb` writes files that start with the `-out` name. The log goes to the standard output (or to `-logfile`).

| Files | Content |
| --- | --- |
| `.pin` `.phr` `.psq` (protein), `.nin` `.nhr` `.nsq` (nucleotide) | Index, headers and sequences. |
| `.pjs` `.njs` | JSON metadata. |
| `.pal` `.nal` | Alias file, for a database in several volumes or built with `-gi_mask`. |
| `.pdb` `.ptf` `.pot` `.pto` (protein), `.ndb` `.ntf` `.not` `.nto` (nucleotide) | Version 5 only: the lookup tables for ids and tax ids. |
| `.psi` `.psd` (protein), `.nsi` `.nsd` (nucleotide) | With `-parse_seqids` in version 4: id indexes. |
| `.pog` (protein), `.nog` (nucleotide) | With `-parse_seqids`, in version 4 and in version 5: id data. |
| `.pos` (protein), `.nos` (nucleotide) | With `-parse_seqids` in version 5: ids of each sequence. |
| `.phi` `.phd` (protein), `.nhi` `.nhd` (nucleotide) | With `-hash_index`. |
| `.naa` `.nab` `.nac` (nucleotide), `.paa` `.pab` `.pac` (protein) | With `-mask_data`: the mask column. |
| `.gmd` `.gnd` `.gmo` `.gno` `.gmi` `.gni` | With `-gi_mask`: the GI-indexed masks, named after `-gi_mask_name` (see [GI-indexed masks](#gi-indexed-masks)). Data files of more than one volume are numbered (`NAME.00.gmd`, ...). |

The sequence and header files have the same bytes as NCBI's. The index and metadata files differ only in the creation date they record. So do the GI-mask index files (`.gmi`, `.gni`) and the alias files. The alias file of a database in one volume built with `-gi_mask` names the volume where NCBI's names the database's full path (see [GI-indexed masks](#gi-indexed-masks)).

## Differences from NCBI BLAST+

- A `-in` file made of several bzip2 streams or zstd frames is read to the end. NCBI's `makeblastdb` stops after the first stream and drops the rest without a message. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- `-gi_mask` on a database in one volume: LOSAT's alias file names the volume by its file name, so the database opens and the build ends normally (exit status 0, with the metadata file). NCBI's names the database by its full path, which no program can open, and ends with exit status 2. See [GI-indexed masks](#gi-indexed-masks) and [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- When the database path (`-out`) contains a space or a double quote, or when an old database of the same name has to be replaced, LOSAT removes and replaces the old files as one database, where `makeblastdb` treats the path as a list of names. Normal names are not affected. See [the decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).
- Database files are compared with `makeblastdb` 2.17.0 on Linux x86-64. On Windows and macOS LOSAT follows NCBI's source where the file layout depends on the platform (the LMDB page size, the size of a `long`); those files have not been compared with NCBI's builds there. The order of `-randomize` is the Linux builds' on every platform: NCBI's Windows and macOS builds shuffle with their own C++ libraries.

The full list is on [Differences from NCBI BLAST+](../differences.md).

## What LOSAT refuses

- `-gi_mask` with two mask files whose GI-mask files are the same files (one `-gi_mask_name` twice, or two names of one path), when both would write one data file at the same time. NCBI's files then depend on how its C++ runtime buffers two streams on one file. LOSAT stops at the sequence where the second would start writing that file (see [What LOSAT refuses](../refusals.md)). Where one of the two writes nothing, or their data files are written one after the other, LOSAT writes NCBI's files.
- `-input_type asn1_txt` and `asn1_bin` untitled nucleotide inputs that trigger `AutoDef` (`AutodefOptions` user object on a non-TLS, non-RefSeq-contig sequence; see [What LOSAT refuses](../refusals.md)). All other `CDeflineGenerator` title paths (`MolInfo`, `Unverified`/`Unreviewed`, `TPA`/`MAG`/`MULTISPECIES`, `BioSource` modifiers, feature-derived protein titles, `NC`/`NM`/`NR`/`WGS`/`TSA`/`TLS`/`PDB`/patent/segmented/map, explicit titles, `ASN1_BlastDefLine`, and `CFastaReader` user objects) work.
- With `-dbtype nucl`, a FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more (such as 4,294,967,293 letters of `A`, `C`, `G` and `T` in a row), where NCBI crashes.
- NCBI settings that change which ids are parsed or which titles are stored (`[BLAST] LONG_SEQID`, a `[BLAST] MAX_SEQID_LENGTH` that `makeblastdb` would use with `-parse_seqids`, `DIAG_*` variables and other `NCBI_CONFIG_*` entries).

For the complete list with the message you see, see [What LOSAT refuses](../refusals.md).

## Examples

Build a protein database and search it:

<!-- example: makeblastdb-1 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db -title "my proteins"
losat blastp -query proteins.faa -db proteins_db -outfmt 6
```

Parse the ids, give each sequence its tax id, and print them back:

<!-- example: makeblastdb-2 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out tax_db -parse_seqids -taxid_map taxids.txt
losat blastdbcmd -db tax_db -entry all -outfmt "%a %l %T"
```

Add DUST masks to a nucleotide database and read them back:

<!-- example: makeblastdb-3 -->
```bash
losat dustmasker -in genome.fna -parse_seqids -outfmt maskinfo_asn1_bin -out genome.dust.asnb
losat makeblastdb -in genome.fna -dbtype nucl -parse_seqids -mask_data genome.dust.asnb -mask_id dust_default -out genome_db
losat blastdbcmd -db genome_db -entry all -outfmt "%a %m100"
```

The algorithm id (100 here) comes from `blastdbcmd -db genome_db -info`.

Write the older format version 4 with a hash index:

<!-- example: makeblastdb-4 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_v4 -blastdb_version 4 -parse_seqids -hash_index
losat blastdbcmd -db proteins_v4 -entry all -outfmt "%a %h"
```

Split a database into volumes of at most 30000 bytes per file (a small size, to show the split), and print the last three sequences:

<!-- example: makeblastdb-5 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out split_db -max_file_sz 30000
losat blastdbcmd -db split_db -entry all -outfmt "%o %l" | tail -n 3
```

## See also

- [BLAST databases](../databases.md)
- [blastdbcmd](blastdbcmd.md): read sequences and information back
- [dustmasker](dustmasker.md), [segmasker](segmasker.md), [windowmasker](windowmasker.md), [convert2blastmask](convert2blastmask.md): make masks
- [What LOSAT refuses](../refusals.md)

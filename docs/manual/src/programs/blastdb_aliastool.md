# blastdb_aliastool

`losat blastdb_aliastool` makes alias files: small text files that name one or more BLAST databases, optionally with a list of the sequences to keep, so that `-db` can search a combination or a subset of databases without copying them. It also converts text lists of GIs and of sequence ids into the binary lists that `-gilist`-style alias keys and `-seqidlist` read fastest, and shows what a binary seqid list holds.

It writes the same files as NCBI's `blastdb_aliastool` 2.17.0, byte for byte apart from the date in them, with the same messages and exit statuses.

## Synopsis

```bash
losat blastdb_aliastool -dblist "DB1 DB2" -dbtype prot -out NAME -title TITLE
losat blastdb_aliastool -db DB -dbtype prot -seqidlist IDS -out NAME
losat blastdb_aliastool -seqid_file_in IDS.txt [-seqid_file_out IDS.bsl] [-seqid_db DB]
losat blastdb_aliastool -gi_file_in GIS.txt [-gi_file_out GIS.bgl]
losat blastdb_aliastool -seqid_file_info IDS.bsl
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastdb_aliastool ...`.

## What it makes

| Command | Writes |
| --- | --- |
| `-dblist "A B"` or `-dblist_file FILE`, with `-dbtype`, `-title`, `-out NAME` | `NAME.pal` (protein) or `NAME.nal` (nucleotide): one database made of `A` and `B`. With `-gilist`, only the GIs of the list. |
| `-db DB -gilist FILE -out NAME` | `NAME.pal`/`.nal` limited to the GIs of `FILE`. A text list is first converted to `NAME.p.gil` (or `.n.gil`). |
| `-db DB -seqidlist FILE -out NAME` | `NAME.pal`/`.nal` limited to the sequence ids of `FILE` (text, or binary from `-seqid_file_in`). |
| `-db DB -taxidlist FILE -out NAME` | `NAME.pal`/`.nal` limited to the taxonomy ids of `FILE`. |
| `-db DB -oid_masks 1 -out NAME` | `NAME.pal`/`.nal` that leaves out the sequences of the database's exclude-model mask (`makeblastdb -oid_masks 1`). Give `-out` in the database's directory. |
| `-num_volumes N -out NAME` | `NAME.pal`/`.nal` over the volumes `NAME.00` to `NAME.N-1`. |
| `-vdblist "SRR1 SRR2"` or `-vdblist_file FILE` | `NAME.nvl`/`.pvl` naming SRA (VDB) runs, or a `VDBLIST` line added to the alias file of `-dblist`. LOSAT writes the file; its searches do not read VDB data. |
| `-gi_file_in FILE` | The binary GI list `-gi_file_out` (default `FILE.bgl`). |
| `-seqid_file_in FILE` | The binary seqid list `-seqid_file_out` (default `FILE.bsl`), titled `-seqid_title` (default: the file's name). With `-seqid_db DB`, only the ids found in `DB` (a version 5 database) are kept, and the list records the database. |
| `-seqid_file_info FILE` | Prints the number of ids, title, date and database of a binary seqid list. |

The alias files hold `NSEQ` and `LENGTH`, the number of sequences and letters the alias selects, as NCBI counts them when the alias file is opened. A list that selects no sequence stops the run with NCBI's error (`No seqs in GILIST were found in BLAST database`), and no alias file is left.

For every option: `losat blastdb_aliastool -help`.

## Differences from NCBI BLAST+

None in the files or the messages, apart from the date written into each file, for the inputs LOSAT runs.

A name or title with bytes that are not UTF-8 (a Latin-1 name, say), in any option or in a line of `-dblist_file` or `-vdblist_file`, is read and written as those bytes, as NCBI's. An alias file of `-seqidlist` or `-taxidlist` over a version 5 database is not yet supported: LOSAT stops with its refusal and leaves no alias file (see [What LOSAT refuses](../refusals.md#blastdb_aliastool)). The same lists over a version 4 database, and every other mode, are supported.

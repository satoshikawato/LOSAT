# blastdbcheck

`losat blastdbcheck` checks BLAST databases for damage: missing or empty files, index entries that point outside the files, sequences without definition lines, alias files whose counts do not add up, and lookups of the sequence ids that fail. It checks one database (`-db`) or every database in a directory (`-dir`, with `-recursive` its subdirectories too), and its exit status is 0 when every test passes and 1 when one fails.

It writes the same report as NCBI's `blastdbcheck` 2.17.0, with the same messages and exit statuses.

## Synopsis

```bash
losat blastdbcheck -db DB [-dbtype prot|nucl] [-full | -stride N | -random N | -ends N] [-verbosity 0-4]
losat blastdbcheck -dir DIR [-recursive] [-verbosity 0-4]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastdbcheck ...`.

## What it tests

| Test | Checks |
| --- | --- |
| Metadata (always) | The database opens, its date and title can be read, the volume's files exist and are not empty, the number and total length of its sequences match its index, and the taxonomy database (`taxdb.bti`, `taxdb.btd`) knows `Homo sapiens` (taxid 9606). With `-legacy`, temporary files of `formatdb` left beside the volume fail the test. |
| Sequences | For each sequence tested: its definition lines and Seq-ids, its letters, its Bioseq, the ids found again through the ISAM indexes (not with `-no_isam`), and with `-must_have_taxids` a taxid. |
| Alias files | Each alias file that `-db` or `-dir` reaches: the files it names exist, and `NSEQ`, `LENGTH`, `FIRST_OID`, `LAST_OID`, `MAXOID` and `OIDLIST` agree with the database. |

Which sequences are tested: every one with `-full`, every `N`-th with `-stride N`, `N` drawn at random with `-random N`, the first and last `N` with `-ends N`; with none of these, 200 drawn at random. The random sample is drawn as NCBI draws it, from a seed taken from the clock, so two runs test different sequences.

`-verbosity` sets how much is written: 0 nothing, 1 only failures, 2 (the default) a summary per volume, 3 each failure's details, 4 the result of each sequence. `-logfile FILE` writes the report to `FILE` instead of the standard output.

For every option: `losat blastdbcheck -help`.

## Differences from NCBI BLAST+

None in the report or the exit status, for the inputs LOSAT runs.

LOSAT refuses `-cdd_delta` (not implemented yet), and the inputs where NCBI crashes, never finishes or reads past the end of a file, or where a result that C leaves undefined decides a test: a cut or empty database file, `-stride 0` or below, a random sample of a database that includes none of its sequences, and an alias file whose counts overflow NCBI's arithmetic where the overflow decides the check (see [What LOSAT refuses](../refusals.md#blastdbcheck)).

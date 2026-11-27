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
| CDD headers (`-cdd_delta`) | For protein DELTA-BLAST profile databases: checks `.aux`, `.loo`, `.rps`, `.freq`, `.wcounts` and `.obsr` headers, magic numbers, profile counts and offset tables. |

Which sequences are tested: every one with `-full`, every `N`-th with `-stride N`, `N` drawn at random with `-random N`, the first and last `N` with `-ends N`; with none of these, 200 drawn at random. The random sample is drawn as NCBI draws it, from a seed taken from the clock, so two runs test different sequences.

`-verbosity` sets how much is written: 0 nothing, 1 only failures, 2 (the default) a summary per volume, 3 each failure's details, 4 the result of each sequence. `-logfile FILE` writes the report to `FILE` instead of the standard output.

For every option: `losat blastdbcheck -help`.

## Differences from NCBI BLAST+

For alias count checks and `-ends` OID bounds that overflow NCBI's 32-bit or 64-bit signed integers, LOSAT computes the values in wider signed integers so overflow does not produce false errors or wrapped OIDs (see [Differences from NCBI BLAST+](../differences.md#blastdbcheck-alias-count-checks-and--ends-near-231)). Every other run writes the same report and exit status as NCBI BLAST+ 2.17.0, including truncated database files whose missing bytes fall within the last 4 KiB memory page.

LOSAT refuses inputs where NCBI crashes, never finishes or reads past its memory mapping: an empty sequence or CDD profile file, a cut file whose offsets pass the last 4 KiB page, a cut string ISAM, a non-index file found by `-dir`, `-stride 0` or below, and a random sample of a database that includes none of its sequences (see [What LOSAT refuses](../refusals.md#blastdbcheck)).

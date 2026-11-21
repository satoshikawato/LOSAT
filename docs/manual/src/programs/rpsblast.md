# rpsblast

`losat rpsblast` searches protein queries against a database of profiles (position-specific scoring matrices). Its best-known use is finding the conserved domains of a protein with NCBI's Conserved Domain Database (CDD). Each hit tells you which profile, and so which domain family, matches which part of your query.

## Synopsis

```bash
losat rpsblast -query QUERY.faa -db PROFILE_DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `rpsblast -query QUERY.faa -db PROFILE_DB [options]`.

`rpsblast` has no `-subject`, as in NCBI BLAST+. The database is always a profile database.

## Input

- **Query.** A protein FASTA file with `-query` (`-query -` reads the standard input). Several queries are searched in one run. `-query_loc START-STOP` searches part of the query. A query may have up to 5,000,000 letters.
- **Database.** `-db` names a profile database made by NCBI's `makeprofiledb`, given without file extension. LOSAT reads the `.rps`, `.loo`, `.aux` and `.freq` files of the database. It does not write profile databases.
- **One volume or several.** A database can be a single volume or many volumes with an alias file (`.pal`) that lists them. Give the alias name, as for NCBI's CDD, which comes in many volumes. The volumes are searched one after the other and the hits are merged for each query.
- **Where to get CDD.** NCBI publishes CDD at ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/. Download the prebuilt RPS-BLAST database there (at the time of writing, the archive `Cdd_LE.tar.gz`), unpack it, and give the path of its `Cdd` database as `-db`. The database has `.freq` files, so the default composition-based statistics work.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-query FILE` | `-` | Protein FASTA query |
| `-db NAME` | | Profile database |
| `-evalue E` | 10 | Report hits with an E-value up to this |
| `-comp_based_stats MODE` | 1 | Composition-based statistics (see below) |
| `-outfmt FORMAT` | 0 | Output format |
| `-out FILE` | standard output | Output file |
| `-max_target_seqs N` | 500 | Number of profiles to keep per query |
| `-max_hsps N` | | Largest number of alignments per profile |
| `-qcov_hsp_perc P` | | Minimum percent of the query covered by an alignment |
| `-seg yes\|no\|"W L H"` | no | Filter low-complexity regions of the query |
| `-lcase_masking` | off | Treat lower-case query letters as masked |
| `-window_size N` | 40 | Two-hit window; 0 selects the one-hit method |
| `-dbsize N`, `-searchsp N` | | Set the effective database length or search space |
| `-num_threads N` | 1 | Number of threads; 0 means one per volume |
| `-mt_mode N` | 0 | Split by volume (0) or by query (1); accepted, the output does not depend on it |

For every option: `losat rpsblast -help`.

### Composition-based statistics

`-comp_based_stats` takes these values:

| Value | Meaning |
| --- | --- |
| `0`, `F`, `f`, `D`, `d` | No composition-based statistics (NCBI's `-help` calls 0 "simplified"; NCBI BLAST+ 2.17.0 and LOSAT run `D` and `d` as 0) |
| `1`, `T`, `t` | Composition-based statistics (the default) |
| `2`, `3` | Run as 1, with a larger preliminary E-value cutoff |

Composition-based statistics need the database's `.freq` files. If the first volume has none, LOSAT prints NCBI's warning and searches without them. If a later volume has none, it stops with NCBI's error.

### Threads

The output is the same at every `-num_threads` value. It is the output of NCBI's `rpsblast -num_threads 1`. NCBI's own multi-thread runs can differ from it.

## Output

Supported `-outfmt` values: 0 to 16, 18 and 20, every format that NCBI's `rpsblast` writes. Formats 6, 7, 10 and 20 take a custom field list. As in NCBI's `rpsblast`, `delim=` is read and checked but the columns keep their tab or comma. Subjects in the report are profile IDs such as `CDD:103`. As in NCBI's `rpsblast`, `positive` and `ppos` count the identities. See [Output formats](../output-formats.md).

Every format is the bytes of NCBI BLAST+ 2.17.0 `rpsblast -num_threads 1`, at any thread count, except that LOSAT names its own version where NCBI names its own (`RPSBLAST (LOSAT 0.2.0)` where NCBI writes `RPSBLAST 2.17.0+`, `LOSAT 0.2.0` in the archive's `ident`) and that a message stored in the archive carries LOSAT's wording where LOSAT words it as its own (the same as `blastp`'s archive). The archive is the input that [rpsbproc](rpsbproc.md) reads: `losat rpsblast -query proteins.faa -db Cdd -outfmt 11 | losat rpsbproc -d cdd_data`.

## Differences from NCBI BLAST+

- Where NCBI's result depends on undefined behavior, LOSAT refuses the input ([the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md)).

## What LOSAT refuses

Not yet supported: `-html`.

Also refused:

- `-sum_stats true`, because NCBI's result then depends on memory contents or crashes.
- `-remote`, `-entrez_query` and `-import_search_strategy`.
- Queries above 5,000,000 letters.
- Profile databases whose files are empty or inconsistent, or whose volumes differ in matrix or gap costs.
- `-subject`, which is not an `rpsblast` option.

NCBI's own errors for a missing `.aux`, `.loo` or `.rps` file come out as in NCBI, with exit status 3. For the full list, see [What LOSAT refuses](../refusals.md).

## Examples

Annotate protein queries against a profile database such as CDD:

<!-- example: rpsblast-1 -->
```bash
losat rpsblast -query proteins.faa -db cdd/Cdd -outfmt 6
```

Without composition-based statistics, with comment lines:

<!-- example: rpsblast-2 -->
```bash
losat rpsblast -query proteins.faa -db cdd/Cdd -comp_based_stats 0 -outfmt 7
```

A database of several volumes, one thread per volume:

<!-- example: rpsblast-3 -->
```bash
losat rpsblast -query proteins.faa -db cdd_volumes/Cdd -num_threads 0 -outfmt 6
```

## See also

- [rpstblastn](rpstblastn.md), for nucleotide queries
- [rpsbproc](rpsbproc.md), to turn archives into domain annotations
- [psiblast](psiblast.md), which makes PSSMs
- [BLAST databases](../databases.md)

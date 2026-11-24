# rpstblastn

`losat rpstblastn` searches nucleotide queries against a database of protein profiles, such as NCBI's Conserved Domain Database (CDD). It translates each query in six frames, so you can find conserved domains in a genome, contig or transcript without first predicting its proteins.

## Synopsis

```bash
losat rpstblastn -query QUERY.fna -db PROFILE_DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `rpstblastn -query QUERY.fna -db PROFILE_DB [options]`.

`rpstblastn` has no `-subject`, as in NCBI BLAST+. The database is always a profile database.

## Input

- **Query.** A nucleotide FASTA file with `-query` (`-query -` reads the standard input). Several queries are searched in one run. `-query_loc START-STOP` searches part of the query. Queries are translated as in [blastx](blastx.md), with the genetic code of `-query_gencode` (default 1).
- **Strand.** `-strand both` (the default), `plus` or `minus` selects which strands of the query are translated.
- **Database.** `-db` names a profile database made by `makeprofiledb` (NCBI's or [LOSAT's](makeprofiledb.md)), given without file extension, in one volume or in several volumes listed by an alias file (`.pal`). It is the kind of database that [rpsblast](rpsblast.md) reads, and CDD is the usual choice. NCBI publishes CDD at ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/. Download the prebuilt RPS-BLAST database there (at the time of writing, the archive `Cdd_LE.tar.gz`), unpack it, and give the path of its `Cdd` database as `-db`.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-query FILE` | `-` | Nucleotide FASTA query |
| `-db NAME` | | Profile database |
| `-strand STRAND` | both | Query strand(s) to translate: `both`, `plus` or `minus` |
| `-query_gencode N` | 1 | Genetic code used to translate the query |
| `-evalue E` | 10 | Report hits with an E-value up to this |
| `-comp_based_stats MODE` | 1 | Composition-based statistics (see below) |
| `-use_sw_tback` | off | Align each profile's matches again with Smith-Waterman (with `-comp_based_stats` 1 to 3) |
| `-outfmt FORMAT` | 0 | Output format |
| `-out FILE` | standard output | Output file |
| `-max_target_seqs N` | 500 | Number of profiles to keep per query |
| `-max_hsps N` | | Largest number of alignments per profile |
| `-qcov_hsp_perc P` | | Minimum percent of the query covered by an alignment |
| `-seg yes\|no\|"W L H"` | no | Filter low-complexity regions of the translated query |
| `-lcase_masking` | off | Treat lower-case query letters as masked |
| `-window_size N` | 40 | Two-hit window; 0 selects the one-hit method |
| `-dbsize N`, `-searchsp N` | | Set the effective database length or search space |
| `-num_threads N` | 1 | Number of threads; 0 means one per volume |
| `-mt_mode N` | 0 | Split by volume (0) or by query (1); the output does not depend on it |

For every option: `losat rpstblastn -help`.

### Composition-based statistics

`-comp_based_stats` takes 0, 1, 2 and 3. The default is 1.

| Value | Meaning |
| --- | --- |
| `0`, `F`, `f`, `D`, `d` | No composition-based statistics (NCBI's `-help` calls 0 "simplified"; NCBI BLAST+ 2.17.0 and LOSAT run `D` and `d` as 0) |
| `1`, `T`, `t` | Composition-based statistics, adjusting the scores of each profile's matches in each frame (the default) |
| `2`, `3` | Run as 1, with a larger preliminary E-value cutoff |

Composition-based statistics need the database's `.freq` files. NCBI's error names the file when it is missing.

### Threads

The output is the same at every `-num_threads` value. It is the output of NCBI's `rpstblastn -num_threads 1`.

## Output

Supported `-outfmt` values: 0 to 16, 18 and 20, every format that NCBI's `rpstblastn` writes. Formats 6, 7, 10 and 20 take a custom field list, including the frame fields such as `qframe`. As in NCBI's `rpstblastn`, `delim=` is read and checked but the columns keep their tab or comma. The reports are those of [blastx](blastx.md), with frames, nucleotide coordinates and NCBI's `RPSTBLASTN` header lines (`RPSTBLASTN (LOSAT 0.2.0)` where NCBI writes `RPSTBLASTN 2.17.0+`). As in NCBI's `rpstblastn`, `positive` and `ppos` count the identities. See [Output formats](../output-formats.md).

Every format is the bytes of NCBI BLAST+ 2.17.0 `rpstblastn -num_threads 1`, at any thread count, except that LOSAT names its own version where NCBI names its own (`RPSTBLASTN (LOSAT 0.2.0)` where NCBI writes `RPSTBLASTN 2.17.0+`, `LOSAT 0.2.0` in the archive's `ident`) and that a message stored in the archive carries LOSAT's wording where LOSAT words it as its own (the same as `blastp`'s archive). The archive is the input that [rpsbproc](rpsbproc.md) reads.

## Differences from NCBI BLAST+

- Where NCBI's result depends on undefined behavior, LOSAT refuses the input ([the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md)).

## What LOSAT refuses

Not yet supported: `-html`; `rpstblastn` in the library API and in LOSAT Web.

Also refused:

- `-query_gencode 32`.
- `-ungapped` without composition-based statistics, where NCBI's `rpstblastn` crashes after the search. With them, you get NCBI's own error.
- `-sum_stats true`, because NCBI's result then depends on memory contents or crashes.
- `-remote`, `-entrez_query` and `-import_search_strategy`.
- Profile databases with the problems that `rpsblast` refuses.
- `-subject`, which is not an `rpstblastn` option (NCBI's argument error, exit status 1).

NCBI's errors for a missing `.aux`, `.loo` or `.rps` file and for an empty `.aux` come as in NCBI, with exit status 3. For the full list, see [What LOSAT refuses](../refusals.md).

## Examples

Find conserved domains in nucleotide sequences:

<!-- example: rpstblastn-1 -->
```bash
losat rpstblastn -query contigs.fna -db cdd/Cdd -outfmt 6
```

Search the forward strand only, and print the frame of each hit:

<!-- example: rpstblastn-2 -->
```bash
losat rpstblastn -query contigs.fna -db cdd/Cdd -strand plus -outfmt "6 qseqid qframe sseqid evalue"
```

## See also

- [rpsblast](rpsblast.md), for protein queries
- [rpsbproc](rpsbproc.md), to turn archives into domain annotations
- [blastx](blastx.md), whose translation and reports `rpstblastn` shares
- [BLAST databases](../databases.md)

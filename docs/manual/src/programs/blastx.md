# blastx

`blastx` translates a nucleotide query in all six reading frames and compares each translation with protein subjects. Use it to find which proteins a transcript, contig or read could encode, and to annotate genes in unannotated DNA.

## Synopsis

```bash
losat blastx -query QUERY.fna -subject PROTEINS.faa [options]
losat blastx -query QUERY.fna -db DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastx -query QUERY.fna -subject PROTEINS.faa [options]`.

## Tasks

Choose a task with `-task`. Options you leave out take the default of the task.

| Task | Use it for | Defaults that differ |
|:---|:---|:---|
| `blastx` (default) | General translated search | Word size 3, threshold 12 |
| `blastx-fast` | A faster search | Word size 5, threshold 20 (uses NCBI's compressed word table) |

## Input

- **Query.** `-query` takes a nucleotide FASTA file with one or more records. The default is the standard input, so `-query -` and no `-query` at all read the standard input. A name ending in `.gz`, `.bz2` or `.zst` is decompressed for you.
- **Genetic code.** `-query_gencode` (default 1) selects the code that translates the query. It takes the NCBI genetic code numbers 1 to 6, 9 to 16, 21 to 31 and 33. `-strand` (`both`, `plus`, `minus`) limits the frames to one strand.
- **Subject.** `-subject` takes a protein FASTA file, in the same forms as `-query`. A subject record can have at most 5,000,000 residues. Protein subjects need no genetic code.
- **Database.** `-db` searches a protein BLAST database written by `makeblastdb` or `losat makeblastdb`, or one of NCBI's downloaded protein databases. Name several databases in one string (`-db "a b"`). `-db` and `-subject` cannot be used together. See [BLAST databases](../databases.md).
- **Database filters.** With `-db` you can restrict the search with `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist` and `-negative_seqidlist` (`-no_taxid_expansion` keeps a taxid from reaching its descendants). `-db_soft_mask` and `-db_hard_mask` use the masks stored in a database that was built with `makeblastdb -mask_data` (see [Masks](../databases.md#masks)).
- **Locations.** `-query_loc START-STOP` searches part of the query (positions in nucleotides), and `-subject_loc START-STOP` part of the subject (with `-subject` only). Positions are 1-based and inclusive.
- **Masking.** SEG is on by default (`12 2.2 2.5`; `-seg no` turns it off). `-lcase_masking` treats lower-case letters of the query as masked.

## Common options

| Option | Default | Meaning |
|:---|:---|:---|
| `-task` | `blastx` | `blastx` or `blastx-fast` |
| `-evalue` | 10 | Report hits with an E-value at or below this value |
| `-outfmt` | 0 | Output format; see [Output](#output) |
| `-out` | standard output | Write the report to a file |
| `-num_threads` | 1 | Number of threads. The output is the same at every thread count |
| `-query_gencode` | 1 | Genetic code for the query |
| `-strand` | `both` | Query strands to translate |
| `-max_target_seqs` | 500 | Keep at most this many subject sequences per query |
| `-max_hsps` | no limit | Keep at most this many alignments per subject |
| `-qcov_hsp_perc` | 0 (off) | Drop alignments that cover less than this percent of the query frame |
| `-matrix` | `BLOSUM62` | Scoring matrix: BLOSUM45, 50, 62, 80, 90 or PAM30, 70, 250 |
| `-gapopen`, `-gapextend` | 11, 1 | Gap costs |
| `-comp_based_stats` | 2 | Composition-based statistics: 0 off, 1, 2, 3 |
| `-use_sw_tback` | off | Align the final results again with Smith-Waterman |
| `-seg` | `12 2.2 2.5` | Low-complexity filter |
| `-sum_stats` | `true` | Use sum statistics for groups of alignments |
| `-max_intron_length` | 0 | Longest intron allowed when linking alignments |
| `-window_size` | 40 | Two-hit window; 0 selects the one-hit method |
| `-culling_limit` | 0 (off) | Delete a hit that at least this many higher-scoring hits envelop |
| `-best_hit_overhang`, `-best_hit_score_edge` | off (0.1 for the one not given) | Keep only the best alignments for each region of the query (the Best Hit algorithm). Values above 0 and below 0.5; not with `-culling_limit` |

For every option: `losat blastx -help`.

## Output

`-outfmt` takes these values:

| Value | Format |
|:---|:---|
| 0 | Pairwise alignments (default) |
| 1 to 4 | Query-anchored alignments (with and without identities, flat) |
| 5 | BLAST XML |
| 6, 7 | Tabular, tabular with comment lines |
| 8, 9 | Seq-align in text and in binary ASN.1 |
| 10 | Comma-separated values |
| 11 | BLAST archive |
| 12 | Seq-align in JSON |
| 13, 14 | BLAST JSON and XML2 as one file per query (need `-out`) |
| 15, 16 | BLAST JSON and XML2 in a single file |
| 18 | Organism report |
| 20 | Comma-separated values with header lines |

Formats 6, 7, 10 and 20 write twelve fields by default: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`. A custom list changes them, for example `-outfmt "6 qseqid sseqid pident qframe evalue"`; `delim=` sets the separator. `qframe` and `frames` give the reading frame. [Output formats](../output-formats.md) describes each format.

## Differences from NCBI BLAST+

- With `-db_hard_mask`, a database sequence whose first mask runs to its end (a masked tail, or a sequence masked entirely) is searched before its mask. NCBI BLAST+ 2.17.0 reads the masked letters there and can crash or report other alignments. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-window_size 0` and a low `-threshold`, NCBI BLAST+ 2.17.0 can read memory past the end of the query batch or of a protein subject and crash. LOSAT reads the letters past the end as the end marker. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-use_sw_tback` in e-value mode (`-sum_stats false`, or a `-max_intron_length` of 1 to 4), NCBI BLAST+ 2.17.0 checks each new Smith-Waterman alignment against a list of kept alignments that it picks by a number: 6 × the query's place in the batch (counted from 0) plus 0 to 5 for frames +1, +2, +3, -1, -2, -3. It keeps one list per query, so a number that is the batch's number of queries or more lies past its lists, and NCBI crashes, fails or goes on with what it reads there. With one query in a batch (for example a search of a single query, or a query that fills a batch of 10,002 letters, BLASTX's batch size, by itself) that is every hit outside frame +1. LOSAT checks such an alignment against an empty list and keeps it, as NCBI does when its lists reach that far (in our tests, NCBI's report of the same queries followed by queries without hits, which keep NCBI inside its lists). See [Differences from NCBI BLAST+](../differences.md#blastx--use_sw_tback-in-e-value-mode).

## What LOSAT refuses

LOSAT stops with an error for these `blastx` settings instead of running something else. The full list is on [What LOSAT refuses](../refusals.md).

- With `-use_sw_tback`, a query batch on which NCBI BLAST+ never finishes or crashes (see [What LOSAT refuses](../refusals.md)).
- `-query_gencode 32`; `blastx` does not accept code 32.
- An `-evalue` that is not a finite number.
- `-matrix IDENTITY` with word sizes 5 to 7 (the default of `-task blastx-fast`), and matrix names that only begin with `IDENTITY`.
- `-threshold` values that overflow NCBI's word-table buffers when the word size is 5 to 7.
- `-remote`, `-entrez_query`, `-gilist`, `-ipglist` and `-import_search_strategy`.
- A protein subject of more than 5,000,000 residues, from `-subject` or from a database.
- `-db_hard_mask` on a database with a sequence whose overlapping or out-of-order masks LOSAT cannot search as NCBI BLAST+ does. One such sequence refuses the whole database. Other databases whose masks overlap, nest or are stored out of order are searched as NCBI BLAST+ searches them, among them the A. thaliana TAIR10 and UniRef50 databases masked by `segmasker`; where NCBI's rows for a sequence depend on the sequences next to it, LOSAT searches it on its own (see [Differences](../differences.md#-db_hard_mask-with-overlapping-masks-that-reach-a-neighbouring-sequence)). See [What LOSAT refuses](../refusals.md#blastx).

## Examples

Translate a transcript file in six frames and search a protein FASTA file.

<!-- example: blastx-1 -->
```bash
losat blastx -query transcripts.fna -subject proteins.faa -query_gencode 1 -outfmt 6
```

Build a protein database, then search it with the `blastx-fast` task. The field list adds the reading frame of each hit.

<!-- example: blastx-2 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db
losat blastx -query transcripts.fna -db proteins_db -task blastx-fast -outfmt "6 qseqid sseqid pident length qframe evalue bitscore"
```

Search only the first 300 nucleotides of each query, on the plus strand.

<!-- example: blastx-3 -->
```bash
losat blastx -query transcripts.fna -subject proteins.faa -query_loc 1-300 -strand plus -outfmt 6
```

## See also

[blastn](blastn.md), [blastp](blastp.md), [tblastn](tblastn.md), [tblastx](tblastx.md), [makeblastdb](makeblastdb.md), [BLAST databases](../databases.md), [Output formats](../output-formats.md), [What LOSAT refuses](../refusals.md)

# blastn

`blastn` compares nucleotide queries with nucleotide subjects. Use it to find where a genome, contig or read set matches another, from near-identical copies (the default `megablast` task) to more distant relatives (`blastn`, `dc-megablast`) and very short queries (`blastn-short`). `rmblastn` scores alignments with the nucleotide matrix alone, and its format 0 report shows raw scores.

## Synopsis

```bash
losat blastn -query QUERY.fa -subject SUBJECT.fa [options]
losat blastn -query QUERY.fa -db DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastn -query QUERY.fa -subject SUBJECT.fa [options]`.

## Tasks

Choose a task with `-task`. Options you leave out take the default of the task.

| Task | Use it for | Defaults that differ |
|:---|:---|:---|
| `megablast` (default) | Highly similar sequences | Word size 28, reward 1, penalty -2, gap costs 0/0 (linear), X-drop for preliminary gapped extension 25 |
| `blastn` | More divergent sequences | Word size 11, reward 2, penalty -3, gap open 5, gap extend 2 |
| `dc-megablast` | Distant, cross-species matches | Discontiguous words (11 of 18 positions, coding template), reward 2, penalty -3, gaps 5/2, two hits within a window of 40 |
| `blastn-short` | Queries of a few dozen bases | Word size 7, reward 1, penalty -3, gaps 5/2, E-value 1000, DUST off |
| `rmblastn` | Scores from the matrix alone, as NCBI's `-task rmblastn` | The defaults of `blastn` with reward 0 and penalty 0 |

With reward 0 and penalty 0 (`-task rmblastn`, or `-reward 0 -penalty 0` with another task), the search scores with NCBI's default nucleotide matrix (reward 1, penalty -3) and reads the X-drops (`-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`) as raw scores. The preliminary search keeps alignments by raw score (at least `-min_raw_gapped_score` when it is given) instead of by E-value; `-evalue` still applies to the final alignments. Format 0 shows raw scores only (` Score = 2722`) and leaves out the list of subjects and the Karlin-Altschul lines; the tabular formats still show bit scores and E-values.

`-template_type` (`coding`, `optimal`, `coding_and_optimal`) and `-template_length` (16, 18 or 21) tune `dc-megablast`.

## Input

- **Query.** `-query` takes a nucleotide FASTA file with one or more records. The default is the standard input, so `-query -` and no `-query` at all read the standard input. A name ending in `.gz`, `.bz2` or `.zst` is decompressed for you.
- **Subject.** `-subject` takes a FASTA file in the same forms as `-query`. LOSAT compares each query with each record of the file, with no database to build first.
- **Database.** `-db` searches a BLAST database written by `makeblastdb` or `losat makeblastdb`, or one of NCBI's downloaded databases. Name several databases in one string (`-db "a b"`). `-db` and `-subject` cannot be used together. See [BLAST databases](../databases.md).
- **Database filters.** With `-db` you can restrict the search with `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist` and `-negative_seqidlist` (`-no_taxid_expansion` keeps a taxid from reaching its descendants). `-db_soft_mask` and `-db_hard_mask` use the masks stored in a database that was built with `makeblastdb -mask_data`.
- **Locations.** `-query_loc START-STOP` searches part of the query, and `-subject_loc START-STOP` part of the subject (with `-subject` only). Positions are 1-based and inclusive.
- **Masking.** DUST masks low-complexity query regions by default (`-dust no` turns it off). `-lcase_masking` treats lower-case letters of the query as masked, and `-soft_masking false` also replaces the masked letters by `N` in the query that is searched. `-filtering_db DB` masks the query regions that match a BLAST database of repeats, such as one made with `makeblastdb` from a repeat library; LOSAT searches it with the fixed settings that NCBI BLAST+ uses for this search.
- **Strands.** `-strand` is `both` (default), `plus` or `minus`.

## Common options

| Option | Default | Meaning |
|:---|:---|:---|
| `-task` | `megablast` | One of the five tasks above |
| `-evalue` | 10 (1000 for `blastn-short`) | Report hits with an E-value at or below this value |
| `-outfmt` | 0 | Output format; see [Output](#output) |
| `-out` | standard output | Write the report to a file |
| `-num_threads` | 1 | Number of threads. The output is the same at every thread count |
| `-max_target_seqs` | 500 | Keep at most this many subject sequences per query |
| `-max_hsps` | no limit | Keep at most this many alignments per subject |
| `-perc_identity` | 0 (off) | Drop alignments below this percent identity |
| `-qcov_hsp_perc` | 0 (off) | Drop alignments that cover less than this percent of the query |
| `-best_hit_overhang`, `-best_hit_score_edge` | off (0.1 for the one not given) | Keep only the best alignments for each region of the query (the Best Hit algorithm). Values above 0 and below 0.5; not with `-culling_limit` |
| `-word_size` | by task | Length of the exact match that starts an alignment |
| `-reward`, `-penalty` | by task | Score of a match and of a mismatch |
| `-gapopen`, `-gapextend` | by task | Gap costs |
| `-dust` | `20 64 1` (`no` for `blastn-short`) | Low-complexity filter: `yes`, `no`, or `LEVEL WINDOW LINKER` |
| `-strand` | `both` | Query strands to search |
| `-dbsize`, `-searchsp` | from the data | Set the effective database size or search space for the E-value |
| `-sorthits`, `-sorthsps` | E-value order | Order of subjects (formats 0 to 4) and of alignments (format 0) |

For every option: `losat blastn -help`.

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
| 17 | SAM (`blastn` only) |
| 18 | Organism report |
| 20 | Comma-separated values with header lines |

Formats 6, 7, 10 and 20 write twelve fields by default: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`. A custom list changes them, for example `-outfmt "6 qseqid sseqid pident evalue qcovs"`; `delim=` sets the separator. `losat blastn -help` shows the field names, and [Output formats](../output-formats.md) describes each format.

## Differences from NCBI BLAST+

- With `-db_hard_mask`, a database sequence whose first mask runs to its end (for example a poly-A tail that is masked) is searched before its mask when the word has 11 or more letters. NCBI BLAST+ 2.17.0 reads the masked letters there and can crash or lose alignments. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- In format 0, a subject title that consists only of punctuation is written; NCBI BLAST+ crashes on it. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

## What LOSAT refuses

LOSAT stops with an error for these `blastn` settings instead of running something else. The full list is on [What LOSAT refuses](../refusals.md).

- `-reward 0` with a penalty other than 0, and `-reward 0 -penalty 0` with the greedy extension of `-task megablast`.
- `-html`, `-remote`, `-entrez_query` and `-gilist`.
- `-use_index true` with a MegaBLAST index (`makembindex`) of only some database volumes, with `-window_size` above 0, or with `-db_hard_mask`; also an index volume whose parts do not fit together, an index not made from the database, a seed that starts at or past the last letter of its sequence, and an old-style index of fewer sequences than the database searched in more than one query batch or query chunk (NCBI reads past the index or the sequence there, or stops with an error). Other searches with an index, and the warning paths (no index found, a first index volume whose header does not read, or a `-word_size` below the index's minimum), run as in NCBI BLAST+.
- `-db_hard_mask` with `-task blastn-short` on a database whose sequences are masked to their end, and on a database whose masks overlap or are stored out of order.
- Extreme `-xdrop_gap`, `-xdrop_gap_final` or `-xdrop_ungap` values for which NCBI BLAST+ overflows its 32-bit arithmetic.

## Examples

Compare two FASTA files with the default `megablast` task and write the tabular report.

<!-- example: blastn-1 -->
```bash
losat blastn -query query.fa -subject subject.fa -outfmt 6
```

Read the query from the standard input and choose the fields of the tabular report.

<!-- example: blastn-2 -->
```bash
cat query.fa | losat blastn -query - -subject subject.fa -task blastn -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovs"
```

Build a database, then search it with commented tabular output.

<!-- example: blastn-3 -->
```bash
losat makeblastdb -in subject.fa -dbtype nucl -out subject_db
losat blastn -query query.fa -db subject_db -task blastn -outfmt 7
```

Search with short sequences, using the `blastn-short` task and a stricter E-value.

<!-- example: blastn-4 -->
```bash
losat blastn -query primers.fa -subject target.fa -task blastn-short -evalue 1e-3 -outfmt 6
```

## See also

[blastp](blastp.md), [blastx](blastx.md), [tblastn](tblastn.md), [tblastx](tblastx.md), [makeblastdb](makeblastdb.md), [BLAST databases](../databases.md), [Output formats](../output-formats.md), [What LOSAT refuses](../refusals.md)

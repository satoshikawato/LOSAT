# tblastx

`tblastx` translates both the query and the subject nucleotide sequences in all six reading frames and compares the translations. Use it to find distant relationships between coding regions in two nucleotide sequences when a DNA comparison finds too little, for example between genomes of related viruses or between a transcript and a genome.

## Synopsis

```bash
losat tblastx -query QUERY.fna -subject SUBJECT.fna [options]
losat tblastx -query QUERY.fna -db DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `tblastx -query QUERY.fna -subject SUBJECT.fna [options]`.

## Input

- **Query.** `-query` takes a nucleotide FASTA file with one or more records. The default is the standard input, so `-query -` and no `-query` at all read the standard input. A name ending in `.gz`, `.bz2` or `.zst` is decompressed for you.
- **Subject.** `-subject` takes a nucleotide FASTA file in the same forms. LOSAT compares each query with each record, with no database to build first.
- **Database.** `-db` searches a nucleotide BLAST database written by `makeblastdb` or `losat makeblastdb`, or one of NCBI's downloaded nucleotide databases. Name several databases in one string (`-db "a b"`). `-db` and `-subject` cannot be used together. See [BLAST databases](../databases.md).
- **Database filters.** With `-db` you can restrict the search with `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist`, `-negative_seqidlist`, `-gilist` and `-negative_gilist` (`-no_taxid_expansion` keeps a taxid from reaching its descendants). `-db_soft_mask` and `-db_hard_mask` use the masks stored in a database that was built with `makeblastdb -mask_data`, or its GI-indexed masks (`makeblastdb -gi_mask`, ids 0, 1, ... or the names of its alias file's `MASKLIST`; see [Masks](../databases.md#masks)).
- **Genetic codes.** `-query_gencode` and `-db_gencode` (both default 1) select the codes that translate the query and the subject. They take the NCBI genetic code numbers 1 to 6, 9 to 16, 21 to 31 and 33. See [Genetic code of the subject](#genetic-code-of-the-subject).
- **Strand.** `-strand` (`both`, `plus`, `minus`) limits the query's frames to one strand. The subject is always searched in all six frames.
- **Locations.** `-query_loc START-STOP` searches part of the query, and `-subject_loc START-STOP` part of the subject (with `-subject` only). Positions are in nucleotides, 1-based and inclusive.
- **Masking.** SEG is on by default (`12 2.2 2.5`; `-seg no` turns it off). `-lcase_masking` treats lower-case letters of the query as masked.

## Genetic code of the subject

What you get depends on how you give the subject:

- **With `-db`, or with the default `-db_gencode 1`:** the same as NCBI BLAST+ 2.17.0. `-query_gencode` always works as in NCBI BLAST+.
- **With `-subject` and a `-db_gencode` other than 1:** LOSAT translates the subject with the code you chose, and uses that translation for the search, the scores, the E-values and the report. NCBI BLAST+ 2.17.0 translates a local `-subject` file with the standard code (code 1) whatever `-db_gencode` says, and only prints the chosen code in the report. Your results therefore differ from NCBI's in that case, and they follow the code you asked for.

This matters when your subject FASTA uses a non-standard code, for example a mitochondrial genome or a Mycoplasma genome, and you set `-db_gencode` for it: the six-frame translation that LOSAT searches is then the one that code gives. If you need NCBI BLAST+'s exact output, leave `-db_gencode` at 1 or search a database with `-db`.

`tblastn` makes the same choice ([tblastn](tblastn.md#genetic-code-of-the-subject)).

## Common options

| Option | Default | Meaning |
|:---|:---|:---|
| `-evalue` | 10 | Report hits with an E-value at or below this value |
| `-outfmt` | 0 | Output format; see [Output](#output) |
| `-out` | standard output | Write the report to a file |
| `-num_threads` | 1 | Number of threads. The output is the same at every thread count |
| `-query_gencode` | 1 | Genetic code for the query |
| `-db_gencode` | 1 | Genetic code for the subject |
| `-strand` | `both` | Query strands to translate |
| `-max_target_seqs` | 500 | Keep at most this many subject sequences per query |
| `-max_hsps` | no limit | Keep at most this many alignments per subject |
| `-qcov_hsp_perc` | 0 (off) | Drop alignments that cover less than this percent of the query frame |
| `-matrix` | `BLOSUM62` | Scoring matrix |
| `-word_size` | 3 | Length of the word that starts an alignment |
| `-threshold` | 13 (or the matrix's suggestion) | Minimum score of a neighbouring word |
| `-window_size` | 40 | Two-hit window; 0 selects the one-hit method |
| `-xdrop_ungap` | 7 | X-drop (in bits) for ungapped extensions |
| `-sum_stats` | `true` | Use sum statistics for groups of alignments |
| `-max_intron_length` | 0 | Only 0 runs: linking with uneven gaps is for blastx and tblastn, and any other value stops with an options error, as in NCBI BLAST+ |
| `-seg` | `12 2.2 2.5` | Low-complexity filter |
| `-culling_limit` | 0 (off) | Delete a hit that at least this many higher-scoring hits envelop |
| `-best_hit_overhang`, `-best_hit_score_edge` | off (0.1 for the one not given) | Keep only the best alignments for each region of the query (the Best Hit algorithm). Values above 0 and below 0.5; not with `-culling_limit` |
| `-subject_besthit` | off | Keep only the best alignment of each subject for each part of the query |
| `-dbsize`, `-searchsp` | from the data | Set the effective database size or search space for the E-value |

For every option: `losat tblastx -help`.

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

Formats 6, 7, 10 and 20 write twelve fields by default: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`. A custom list changes them, for example `-outfmt "6 std qframe sframe"` adds the reading frames to the default fields; `delim=` sets the separator. [Output formats](../output-formats.md) describes each format.

## Differences from NCBI BLAST+

- With `-db_hard_mask`, a reading frame whose first part, from its first letter, ends at a mask with no other part of the frame after it (masks at or near either end of a sequence) is searched over that part alone. NCBI BLAST+ 2.17.0 reads the masked letters of the frame there and can crash or report other alignments, some of them many times. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-db_hard_mask`, a sequence with a mask within its first two bases, or within the one or two bases after its last whole codon, or with an empty mask at its start, is searched: the reading frame that such a mask leaves an empty part at its start is searched whole, and the masked bases stay masked in every reading frame. NCBI BLAST+ 2.17.0 crashes on it in every run. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-subject` and a `-db_gencode` other than 1, LOSAT translates the subject with that code (see [Genetic code of the subject](#genetic-code-of-the-subject)). See [the decision record on genetic codes for local tblastx subjects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md).
- On large genome pairs with the default options, NCBI BLAST+ 2.17.0 turns a cutoff that is out of the range of a double into a word cutoff of 1. In the runs we made (a 1.2 Mb query against a 5.5 Mb genome) it ran out of memory or time and wrote nothing. LOSAT keeps the equal finite cutoff and writes the hits. Where the cutoff stays finite, both give the same bytes. With an extreme `-searchsp` and `-evalue` on small inputs, NCBI BLAST+ reports more rows than LOSAT. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- In format 0, a subject title that consists only of punctuation is written; NCBI BLAST+ crashes on it. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With so many very short queries in one batch that its average translated length is 0, or with a very short query against a very short subject, against one subject sequence, NCBI BLAST+ 2.17.0 crashes; LOSAT writes the linked hits. See [Very many short queries, or a very short query, against one subject](../differences.md#very-many-short-queries-or-a-very-short-query-against-one-subject) and [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

## What LOSAT refuses

LOSAT stops with an error for these `tblastx` settings instead of running something else. The full list is on [What LOSAT refuses](../refusals.md).

- `-db_gencode 32` and `-query_gencode 32`.
- An `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62) when an extension would pass into the next frame of the query.
- A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the ungapped re-evaluation.
- `-remote` and `-entrez_query`.

## Examples

Compare two nucleotide FASTA files in six frames. The field list adds the reading frames of each hit.

<!-- example: tblastx-1 -->
```bash
losat tblastx -query genome_a.fna -subject genome_b.fna -outfmt "6 std qframe sframe"
```

Build a nucleotide database, then search it with genetic code 4 for both sequences. A database search follows NCBI BLAST+ for every code.

<!-- example: tblastx-2 -->
```bash
losat makeblastdb -in genome_b.fna -dbtype nucl -out genome_b_db
losat tblastx -query genome_a.fna -db genome_b_db -query_gencode 4 -db_gencode 4 -outfmt 6
```

Search a FASTA subject with genetic code 4. LOSAT translates the subject with code 4; NCBI BLAST+ would use code 1 here.

<!-- example: tblastx-3 -->
```bash
losat tblastx -query genome_a.fna -subject genome_b.fna -query_gencode 4 -db_gencode 4 -outfmt 6
```

## See also

[blastn](blastn.md), [blastx](blastx.md), [tblastn](tblastn.md), [makeblastdb](makeblastdb.md), [BLAST databases](../databases.md), [Output formats](../output-formats.md), [What LOSAT refuses](../refusals.md)

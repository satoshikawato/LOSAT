# tblastn

`tblastn` compares protein queries with nucleotide subjects translated in all six reading frames. Use it to find where a protein (or a protein profile) is encoded in a genome, a set of contigs or a transcriptome, including DNA that has no gene annotation.

## Synopsis

```bash
losat tblastn -query QUERY.faa -subject GENOME.fna [options]
losat tblastn -query QUERY.faa -db DB [options]
losat tblastn -in_pssm QUERY.pssm -subject GENOME.fna [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `tblastn -query QUERY.faa -subject GENOME.fna [options]`.

## Tasks

Choose a task with `-task`. Options you leave out take the default of the task.

| Task | Use it for | Defaults that differ |
|:---|:---|:---|
| `tblastn` (default) | General translated search | Word size 3, threshold 13 |
| `tblastn-fast` | A faster search | Word size 5, threshold 20 (uses NCBI's compressed word table) |

## Input

- **Query.** `-query` takes a protein FASTA file with one or more records. The default is the standard input, so `-query -` and no `-query` at all read the standard input. A name ending in `.gz`, `.bz2` or `.zst` is decompressed for you. To search with a profile instead, see [Searching with a PSSM](#searching-with-a-pssm).
- **Subject.** `-subject` takes a nucleotide FASTA file in the same forms. LOSAT translates each record in six frames and compares your query with every frame.
- **Database.** `-db` searches a nucleotide BLAST database written by `makeblastdb` or `losat makeblastdb`, or one of NCBI's downloaded nucleotide databases. Name several databases in one string (`-db "a b"`). `-db` and `-subject` cannot be used together. See [BLAST databases](../databases.md).
- **Database filters.** With `-db` you can restrict the search with `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist` and `-negative_seqidlist` (`-no_taxid_expansion` keeps a taxid from reaching its descendants). `-db_soft_mask` and `-db_hard_mask` use the masks stored in a database that was built with `makeblastdb -mask_data` (see [Masks](../databases.md#masks)).
- **Genetic code.** `-db_gencode` (default 1) selects the code that translates the subject. See [Genetic code of the subject](#genetic-code-of-the-subject).
- **Locations.** `-query_loc START-STOP` searches part of the query (in residues), and `-subject_loc START-STOP` part of the subject (in nucleotides, with `-subject` only). Positions are 1-based and inclusive.
- **Masking.** SEG is on by default (`12 2.2 2.5`; `-seg no` turns it off). `-lcase_masking` treats lower-case letters of the query as masked.

## Genetic code of the subject

`-db_gencode` takes the NCBI genetic code numbers 1 to 6, 9 to 16 and 21 to 33. NCBI BLAST+'s own command line does not take 32; `tblastn` does.

What you get depends on how you give the subject:

- **With `-db`, or with the default code 1:** the same as NCBI BLAST+ 2.17.0 (except code 32, which NCBI BLAST+ rejects).
- **With `-subject` and a `-db_gencode` other than 1:** LOSAT translates the subject with the code you chose, and uses that translation for the search, the scores, the E-values and the report. NCBI BLAST+ 2.17.0 translates a local `-subject` file with the standard code (code 1) whatever `-db_gencode` says, and only prints the chosen code in the report. Your results therefore differ from NCBI's in that case, and they follow the code you asked for.

This matters when your subject FASTA uses a non-standard code, for example a mitochondrial genome or a Mycoplasma genome, and you set `-db_gencode` for it: the proteins you search against are then the ones that code gives. If you need NCBI BLAST+'s exact output, leave `-db_gencode` at 1 or search a database with `-db`.

## Searching with a PSSM

`-in_pssm FILE` replaces `-query` and searches with a position-specific scoring matrix, as NCBI's PSI-TBLASTN does. `-in_pssm -` reads the PSSM from the standard input. The report header says `PSITBLASTN`.

- The file is a PSSM in NCBI's ASN.1 text form (`PssmWithParameters`), such as `psiblast -out_pssm` writes. [psiblast](psiblast.md) describes what it accepts. A binary or XML PSSM, or a compressed file, is refused.
- The PSSM's scores replace the matrix, and its Karlin-Altschul values set the cutoffs, E-values and bit scores.
- As in NCBI, `-comp_based_stats` has no effect, alignments are linked only with `-sum_stats true`, and `-task tblastn-fast` changes nothing. The word size is 2 to 4 and the matrix is one of the standard ones.
- The subjects are translated with genetic code 1 whatever `-db_gencode` says. LOSAT still checks the value.
- The report formats are 0 to 4, 6, 7 and 10. Other formats are refused.
- Also refused with `-in_pssm`:
  - `-ungapped`;
  - `-use_sw_tback` on a search that finds an alignment, where NCBI fails, crashes, never finishes or writes rows that change from one run to the next (a search that finds nothing runs);
  - `-db_gencode 32`;
  - `-outfmt 11`;
  - `-mt_mode 1` with more than one thread and `-db`;
  - a PSSM whose query `-seg` masks completely.

## Common options

| Option | Default | Meaning |
|:---|:---|:---|
| `-task` | `tblastn` | `tblastn` or `tblastn-fast` |
| `-evalue` | 10 | Report hits with an E-value at or below this value |
| `-outfmt` | 0 | Output format; see [Output](#output) |
| `-out` | standard output | Write the report to a file |
| `-num_threads` | 1 | Number of threads. The output is the same at every thread count |
| `-db_gencode` | 1 | Genetic code for the subject |
| `-max_target_seqs` | 500 | Keep at most this many subject sequences per query |
| `-max_hsps` | no limit | Keep at most this many alignments per subject |
| `-qcov_hsp_perc` | 0 (off) | Drop alignments that cover less than this percent of the query |
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
| `-subject_besthit` | off | Keep only the best alignment of each subject for each part of the query |

For every option: `losat tblastn -help`.

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

With `-in_pssm`, only 0 to 4, 6, 7 and 10 work.

Formats 6, 7, 10 and 20 write twelve fields by default: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`. A custom list changes them, for example `-outfmt "6 qseqid sseqid pident sframe evalue"`; `delim=` sets the separator. [Output formats](../output-formats.md) describes each format.

## Differences from NCBI BLAST+

- With `-db_hard_mask`, a reading frame whose first part, from its first letter, ends at a mask with no other part of the frame after it (masks at or near either end of a sequence) is searched over that part alone. NCBI BLAST+ 2.17.0 reads the masked letters of the frame there and can crash or report other alignments. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-db_hard_mask`, a sequence with a mask within its first two bases, or within the one or two bases after its last whole codon, or with an empty mask at its start, is searched: the reading frame that such a mask leaves an empty part at its start is searched as outside the mask. NCBI BLAST+ 2.17.0 crashes on it in every run. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-subject` and a `-db_gencode` other than 1, LOSAT translates the subject with that code (see [Genetic code of the subject](#genetic-code-of-the-subject)). See [the decision record on genetic code 32 for local tblastn subjects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TLOSAN-LOCAL-GENCODE-32.md).
- `-db_gencode 32` is accepted (without `-in_pssm`). NCBI BLAST+ 2.17.0 rejects it.
- With `-window_size 0` and a low `-threshold`, NCBI BLAST+ can read memory past its buffers and crash. LOSAT reads the letters past the end as the end marker. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- In format 0, a subject title that consists only of punctuation is written; NCBI BLAST+ crashes on it. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

## What LOSAT refuses

LOSAT stops with an error for these `tblastn` settings instead of running something else. The full list is on [What LOSAT refuses](../refusals.md).

- `-ungapped`.
- With `-use_sw_tback`, the rare query and subject pairs on which NCBI BLAST+ never finishes or crashes (see [What LOSAT refuses](../refusals.md)).
- An `-evalue` of 1.7976931348623157e308 or more.
- `-matrix IDENTITY` with word sizes 5 to 7 (the default of `-task tblastn-fast`).
- `-threshold` values that overflow NCBI's word-table buffers when the word size is 5 to 7.
- `-html`, `-remote`, `-entrez_query`, `-gilist` and `-import_search_strategy`.
- A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback.
- The `-in_pssm` settings listed [above](#searching-with-a-pssm).

## Examples

Search a genome with a protein query.

<!-- example: tblastn-1 -->
```bash
losat tblastn -query proteins.faa -subject genome.fna -outfmt 6
```

Build a nucleotide database, then search it with a protein query and genetic code 4. A database search follows NCBI BLAST+ for every code.

<!-- example: tblastn-2 -->
```bash
losat makeblastdb -in genome.fna -dbtype nucl -out genome_db
losat tblastn -query proteins.faa -db genome_db -db_gencode 4 -outfmt 6
```

Search a FASTA subject with genetic code 4. LOSAT translates the subject with code 4; NCBI BLAST+ would use code 1 here.

<!-- example: tblastn-3 -->
```bash
losat tblastn -query proteins.faa -subject genome.fna -db_gencode 4 -outfmt 6
```

Search with a PSSM instead of a query sequence.

<!-- example: tblastn-4 -->
```bash
losat tblastn -in_pssm query.pssm -subject genome.fna -outfmt 6
```

## See also

[blastp](blastp.md), [blastx](blastx.md), [tblastx](tblastx.md), [psiblast](psiblast.md), [makeblastdb](makeblastdb.md), [BLAST databases](../databases.md), [Output formats](../output-formats.md), [What LOSAT refuses](../refusals.md)

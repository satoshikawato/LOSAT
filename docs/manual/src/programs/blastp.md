# blastp

`blastp` compares protein queries with protein subjects. Use it to find homologs of a protein in a proteome or a protein database, or short peptides in longer proteins.

## Synopsis

```bash
losat blastp -query QUERY.faa -subject SUBJECT.faa [options]
losat blastp -query QUERY.faa -db DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastp -query QUERY.faa -subject SUBJECT.faa [options]`.

## Tasks

Choose a task with `-task`. Options you leave out take the default of the task.

| Task | Use it for | Defaults that differ |
|:---|:---|:---|
| `blastp` (default) | General protein comparison | Word size 3, BLOSUM62, gap open 11, gap extend 1, threshold 11 |
| `blastp-fast` | A faster search | Word size 5, threshold 20 (uses NCBI's compressed word table) |
| `blastp-short` | Very short peptide queries | Word size 2, PAM30, gap open 9, gap extend 1, E-value 20000, threshold 16, window size 15 |

## Input

- **Query.** `-query` takes a protein FASTA file with one or more records. The default is the standard input, so `-query -` and no `-query` at all read the standard input. A name ending in `.gz`, `.bz2` or `.zst` is decompressed for you.
- **Subject.** `-subject` takes a protein FASTA file in the same forms. LOSAT compares each query with each record, with no database to build first.
- **Database.** `-db` searches a protein BLAST database written by `makeblastdb` or `losat makeblastdb`, or one of NCBI's downloaded databases such as `swissprot`. Name several databases in one string (`-db "a b"`). `-db` and `-subject` cannot be used together. See [BLAST databases](../databases.md).
- **Database filters.** With `-db` you can restrict the search with `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist` and `-negative_seqidlist` (`-no_taxid_expansion` keeps a taxid from reaching its descendants). `-db_soft_mask` and `-db_hard_mask` use the masks stored in a database that was built with `makeblastdb -mask_data`, or its GI-indexed masks (`makeblastdb -gi_mask`, ids 0, 1, ... or the names of its alias file's `MASKLIST`; see [Masks](../databases.md#masks)).
- **Locations.** `-query_loc START-STOP` searches part of the query, and `-subject_loc START-STOP` part of the subject (with `-subject` only). Positions are 1-based and inclusive.
- **Masking.** SEG is off by default (`-seg yes` or `-seg "WINDOW LOCUT HICUT"` turns it on). `-lcase_masking` treats lower-case letters of the query as masked.

## Common options

| Option | Default | Meaning |
|:---|:---|:---|
| `-task` | `blastp` | One of the three tasks above |
| `-evalue` | 10 (20000 for `blastp-short`) | Report hits with an E-value at or below this value |
| `-outfmt` | 0 | Output format; see [Output](#output) |
| `-out` | standard output | Write the report to a file |
| `-num_threads` | 1 | Number of threads. The output is the same at every thread count |
| `-max_target_seqs` | 500 | Keep at most this many subject sequences per query |
| `-max_hsps` | no limit | Keep at most this many alignments per subject |
| `-qcov_hsp_perc` | 0 (off) | Drop alignments that cover less than this percent of the query |
| `-matrix` | `BLOSUM62` (`PAM30` for `blastp-short`) | Scoring matrix: BLOSUM45, 50, 62, 80, 90 or PAM30, 70, 250 |
| `-gapopen`, `-gapextend` | 11, 1 (the matrix's best values) | Gap costs; the pairs NCBI allows for the matrix |
| `-word_size` | 3 (2 for `blastp-short`, 5 for `blastp-fast`) | Length of the word that starts an alignment |
| `-threshold` | 11 (by word size and task) | Minimum score of a neighbouring word |
| `-comp_based_stats` | 2 | Composition-based statistics: 0 off, 1, 2, 3. A `u` or `U` after a mode other than 0 (`2u`) also combines each e-value with a P-value of the two sequences' compositions (unified P-values) |
| `-use_sw_tback` | off | Align the final results again with Smith-Waterman |
| `-ungapped` | off | Search without gapped extension: report the ungapped alignments (with `-comp_based_stats 0`) |
| `-seg` | `no` | Low-complexity filter |
| `-window_size` | 40 (15 for `blastp-short`) | Two-hit window; 0 selects the one-hit method |
| `-dbsize`, `-searchsp` | from the data | Set the effective database size or search space for the E-value |
| `-culling_limit` | 0 (off) | Delete a hit that at least this many higher-scoring hits envelop |
| `-best_hit_overhang`, `-best_hit_score_edge` | off (0.1 for the one not given) | Keep only the best alignments for each region of the query (the Best Hit algorithm). Values above 0 and below 0.5; not with `-culling_limit` |
| `-subject_besthit` | off | Keep only the best alignment of each subject for each part of the query |

For every option: `losat blastp -help`.

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

Formats 6, 7, 10 and 20 write twelve fields by default: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`. A custom list changes them, for example `-outfmt "6 qseqid sseqid pident evalue qcovs"`; `delim=` sets the separator. `losat blastp -help` shows the field names, and [Output formats](../output-formats.md) describes each format.

## Differences from NCBI BLAST+

- With `-db_hard_mask`, a database sequence whose first mask runs to its end (a masked tail, or a sequence masked entirely) is searched before its mask. NCBI BLAST+ 2.17.0 reads the masked letters there and can crash or report other alignments. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With `-window_size 0` and a low `-threshold`, NCBI BLAST+ 2.17.0 can read memory past the end of a sequence and crash. LOSAT reads the letters past the end as the end marker, which gives NCBI BLAST+'s own output where its stray bytes are zero. See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).
- With unified P-values (`-comp_based_stats 1u`, `2u` or `3u`), NCBI BLAST+ 2.17.0 adjusts the e-values of a query with the length and the search space of the query that follows it in the same batch of queries (it reads the queries in batches of about 10,000 residues; the environment variable `BATCH_SIZE` sets the size). LOSAT does the same, and writes what NCBI BLAST+ writes with `-num_threads 1` (apart from the version line, which names LOSAT). For the last query of a batch, which includes a search with one query, NCBI BLAST+ reads past the end of its data, and its report can change from run to run. There LOSAT uses the query's own length and search space, which gives the report NCBI BLAST+ writes for that query when a copy of it follows in the same batch and does not change the pieces in which NCBI BLAST+ searches the batch (the environment variable `CHUNK_SIZE`). See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

## What LOSAT refuses

LOSAT stops with an error for these `blastp` settings instead of running something else. The full list is on [What LOSAT refuses](../refusals.md).

- `-use_sw_tback` with `-ungapped` at the first query batch that is searched, and with `-use_sw_tback` the rare query batches on which NCBI BLAST+ never finishes or crashes (see [What LOSAT refuses](../refusals.md)).
- With `-ungapped`: a subject of more than 5,000,000 residues, and a `-matrix` name that begins with `IDENTITY` but is not `IDENTITY`.
- With `-ungapped -window_size 0` and a low `-threshold`: `-outfmt 1` to `4` with a one-hit alignment of length 0 among the alignments shown, where NCBI BLAST+ crashes.
- `-matrix IDENTITY` with word sizes 5 to 7 (the default of `-task blastp-fast`).
- `-threshold` values that overflow NCBI's word-table buffers when the word size is 5 to 7.
- `-remote`, `-entrez_query`, `-gilist`, `-ipglist` and `-import_search_strategy`.
- `-db_hard_mask` on a database with a sequence whose overlapping or out-of-order masks LOSAT cannot search as NCBI BLAST+ does. One such sequence refuses the whole database. Other databases whose masks overlap, nest or are stored out of order are searched as NCBI BLAST+ searches them, among them the A. thaliana TAIR10 and UniRef50 databases masked by `segmasker`; where NCBI's rows for a sequence depend on the sequences next to it, LOSAT searches it on its own (see [Differences](../differences.md#-db_hard_mask-with-overlapping-masks-that-reach-a-neighbouring-sequence)). See [What LOSAT refuses](../refusals.md#blastp).

## Examples

Compare two protein FASTA files and keep hits with an E-value of 1e-5 or lower.

<!-- example: blastp-1 -->
```bash
losat blastp -query query.faa -subject proteins.faa -evalue 1e-5 -outfmt 6
```

Build a protein database, then search it with two threads. The report is the same as with one thread.

<!-- example: blastp-2 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db
losat blastp -query query.faa -db proteins_db -num_threads 2 -outfmt 7
```

Search with a short peptide, using the `blastp-short` task.

<!-- example: blastp-3 -->
```bash
losat blastp -query peptide.faa -subject proteins.faa -task blastp-short -outfmt 6
```

## See also

[blastn](blastn.md), [blastx](blastx.md), [tblastn](tblastn.md), [psiblast](psiblast.md), [makeblastdb](makeblastdb.md), [BLAST databases](../databases.md), [Output formats](../output-formats.md), [What LOSAT refuses](../refusals.md)

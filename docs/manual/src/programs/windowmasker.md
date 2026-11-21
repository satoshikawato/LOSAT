# windowmasker

`losat windowmasker` masks the repeats of nucleotide sequences with the WindowMasker algorithm, from a file of unit counts made for the genome. Use it to see which parts of a genome a nucleotide search would mask, to write a masked FASTA file, or to make masks for [`makeblastdb -mask_data`](makeblastdb.md).

It reads a FASTA file or an existing nucleotide BLAST database, and a counts file in any of the four formats of NCBI's `windowmasker` (`ascii`, `binary`, `oascii`, `obinary`). It writes the same bytes as NCBI's `windowmasker` 1.0.0 in all eight output formats.

LOSAT masks with a counts file; it does not make one yet. Make the counts with NCBI's `windowmasker -mk_counts`, or use the counts files that NCBI publishes for many genomes (the files that `blastn -window_masker_taxid` reads).

## Synopsis

```bash
losat windowmasker -ustat COUNTS -in SEQUENCES.fna [options]
losat windowmasker -ustat COUNTS -infmt blastdb -in DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `windowmasker -ustat COUNTS -in SEQUENCES.fna [options]`.

## Input

- `-ustat` names the counts file. It is required.
- `-in` names the input. `-in -` (the default) reads the standard input.
- `-infmt fasta` (the default) reads FASTA. `-infmt blastdb` reads every sequence of a nucleotide database named by `-in`, in every volume. `-infmt seqids` stops with NCBI's error, as NCBI's does.
- `-parse_seqids` parses the ids of the FASTA definition lines, so the output names sequences by their ids. Use it when the masks are for a database built with `makeblastdb -parse_seqids`.
- `-ids FILE` masks only the sequences that match one of the ids in `FILE` (one per line); `-exclude_ids FILE` masks all others. As in NCBI's, an id matches the first word of a sequence's title, whole or by its `|`-separated parts. The title is NCBI's: HTML entities such as `&amp;` are not decoded and only the dots and spaces at its end are dropped; a third-party id (`tpg`, `tpe`, `tpd`) adds `TPA: ` before it; a database sequence without a title and with a tax id takes the scientific name of its organism from the taxonomy database (`taxdb`), with its first letter in upper case, or `No definition line found` without one. Matching is case-sensitive. The other sequences are left out of the output.
- `-text_match F` (with `-infmt blastdb` only) reads each id of the `-ids` or `-exclude_ids` file as a Seq-id, such as `gb|AB000001.1|` or `AB000001.1`, and matches it against the Seq-ids of each database sequence, as NCBI's does: case does not matter (except in a PDB chain), the versions must agree (`AB000001` does not match `AB000001.1`), and an accession alone does not match an id stored with a name (`P12345` does not match `sp|P12345|ABC_HUMAN`). A database made with `-parse_seqids` stores a FASTA id that is a bare word, such as `rec1`, as `lcl|rec1`; in a database made without it, the second sequence is `gnl|BL_ORD_ID|1`, the third `gnl|BL_ORD_ID|2`, and so on; `gnl|BL_ORD_ID|0` matches no sequence, as in NCBI's. An id that is not a Seq-id, such as a bare word, is written to the standard error as an error and skipped; when no id of the file can be read, every sequence is masked. With `-outfmt seqloc_*` or `maskinfo_*`, a sequence's Seq-id is written as NCBI's writes it once the list has named it: `AB000001.1` in the list is a DDBJ accession, so `gb|AB000001.1|`, and every other GenBank, EMBL or DDBJ accession with the same letters, number of digits and version, is written as a DDBJ id, and a name is written in the list's case.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-ustat` | none | The counts file. |
| `-in` | `-` | Input file or database name. |
| `-infmt` | `fasta` | `fasta`, `blastdb` or `seqids`. |
| `-out` | `-` | Output file. |
| `-outfmt` | `interval` | Output format (below). |
| `-window` | from the counts | Window length. |
| `-t_thres`, `-t_extend`, `-t_high`, `-t_low` | from the counts | Scores that start masking, extend it, and bound the unit scores. |
| `-t_thres_pct`, `-t_extend_pct`, `-t_high_pct`, `-t_low_pct` | none | The same scores as a percentage of the units, from the `##pct:` lines of the counts file. |
| `-set_t_high`, `-set_t_low` | none | Scores used for units above `-t_high` or below `-t_low`. |
| `-dust` | `F` | Also mask low-complexity regions with DUST, as `dustmasker` does. |
| `-dust_level` | `20` | DUST score threshold with `-dust T`. |
| `-parse_seqids` | off | Parse the Seq-ids of FASTA input. |

For every option: `losat windowmasker -help`.

## Output

`-outfmt` is one of:

| Format | What you get |
| --- | --- |
| `interval` | For each sequence, its definition line, then one line `START - END` per masked interval. Positions count from 0 and include both ends. |
| `fasta` | The sequences with masked letters in lower case. |
| `seqloc_asn1_text`, `seqloc_asn1_bin`, `seqloc_xml` | The intervals as Seq-loc objects. A sequence with several Seq-ids is named, here and in `maskinfo_*`, by its gi (not gi 0), else its first accession-type id (GenBank, RefSeq, Swiss-Prot, ...), else its first Seq-id, as NCBI's does. |
| `maskinfo_asn1_text`, `maskinfo_asn1_bin`, `maskinfo_xml` | The intervals as masking information. The `maskinfo_asn1_*` files are the input of `makeblastdb -mask_data`. |

Warnings about the input, such as a counts file without percentages when a `-t_*_pct` option is given, are written to the standard error, as NCBI writes them.

## Differences from NCBI BLAST+

None that change the output of a normal run. When the reader of the output closes a pipe early, LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

## What LOSAT refuses

- `-mk_counts` and `-convert`: LOSAT does not make or convert counts files yet.
- Counts files on which NCBI's result is undefined, such as an `ascii` file with units of more than 16 letters.
- An ids file matched as text (without `-text_match F`) with an empty id (a line that starts with a space or a tab, `>`, `|` or `>|`).

See [What LOSAT refuses](../refusals.md#windowmasker).

## Examples

Print the masked intervals of each sequence:

<!-- example: windowmasker-1 -->
```bash
losat windowmasker -ustat genome.counts -in genome.fna
```

Also mask the low-complexity regions:

<!-- example: windowmasker-2 -->
```bash
losat windowmasker -ustat genome.counts -in genome.fna -dust T
```

Write the masking information as ASN.1 text, as `makeblastdb -mask_data` reads it:

<!-- example: windowmasker-3 -->
```bash
losat windowmasker -ustat genome.counts -in genome.fna -parse_seqids -outfmt maskinfo_asn1_text
```

Mask a database you built earlier:

<!-- example: windowmasker-4 -->
```bash
losat makeblastdb -in genome.fna -dbtype nucl -out genome_db
losat windowmasker -ustat genome.counts -infmt blastdb -in genome_db
```

## See also

- [dustmasker](dustmasker.md)
- [makeblastdb](makeblastdb.md)
- [blastn](blastn.md): `-window_masker_db` masks the query with a counts file during a search

# dustmasker

`losat dustmasker` finds the low-complexity regions of nucleotide sequences with the symmetric DUST algorithm. Use it to see which parts of a genome a nucleotide search would mask, to write a masked FASTA file, or to make masks for [`makeblastdb -mask_data`](makeblastdb.md).

It reads a FASTA file or an existing nucleotide BLAST database. It writes the same bytes as NCBI's `dustmasker` 1.0.0 in all nine output formats.

## Synopsis

```bash
losat dustmasker -in SEQUENCES.fna [options]
losat dustmasker -infmt blastdb -in DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `dustmasker -in SEQUENCES.fna [options]`.

## Input

- `-in` names the input. `-in -` (the default) reads the standard input.
- `-infmt fasta` (the default) reads FASTA. `-infmt blastdb` reads every sequence of a nucleotide database named by `-in`, in every volume.
- `-parse_seqids` parses the ids of the FASTA definition lines, so the output names sequences by their ids (`lcl|name`, accessions). Use it when the masks are for a database built with `makeblastdb -parse_seqids`.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-in` | `-` | Input file or database name. |
| `-infmt` | `fasta` | `fasta` or `blastdb`. |
| `-out` | `-` | Output file. |
| `-outfmt` | `interval` | Output format (below). |
| `-window` | `64` | DUST window length. |
| `-level` | `20` | DUST score threshold for subwindows. |
| `-linker` | `1` | Masked intervals closer than this are merged. |
| `-parse_seqids` | off | Parse the Seq-ids of FASTA input. |
| `-hard_masking` | off | With `-outfmt fasta`, write `N` in place of masked letters instead of lower case. |

For every option: `losat dustmasker -help`.

## Output

`-outfmt` is one of:

| Format | What you get |
| --- | --- |
| `interval` | For each sequence, its definition line, then one line `START - END` per masked interval. Positions count from 0 and include both ends. |
| `acclist` | One line per masked interval: the definition line, the start and the end, separated by tabs. |
| `fasta` | The sequences with masked letters in lower case (`N` with `-hard_masking`). |
| `seqloc_asn1_text`, `seqloc_asn1_bin`, `seqloc_xml` | The intervals as Seq-loc objects. A sequence with several Seq-ids is named, here and in `maskinfo_*`, by its gi (not gi 0), else its first accession-type id (GenBank, RefSeq, Swiss-Prot, ...), else its first Seq-id, as NCBI's does. |
| `maskinfo_asn1_text`, `maskinfo_asn1_bin`, `maskinfo_xml` | The intervals as masking information with the algorithm name and options. The `maskinfo_asn1_*` files are the input of `makeblastdb -mask_data`. |

## Differences from NCBI BLAST+

None that change the output of a normal run. When the reader of the output closes a pipe early, LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

## What LOSAT refuses

`dustmasker` runs every option of `-help`. It stops with an error, where NCBI crashes or reads memory outside its buffers, in these cases:

- A FASTA record of exactly 4,294,967,295 letters, counting the letters of its gap lines. Longer records are read as NCBI reads them.
- A database title in a text encoding that NCBI cannot tell, or a database definition line with bytes that NCBI reads as `#`.
- A database sequence with a tax id whose record is in an empty `taxdb.btd`, or in a damaged `taxdb.btd` that reaches past the end of NCBI's mapping of the file.
- A failed write of a `seqloc_*` or `maskinfo_*` output.

See [What LOSAT refuses](../refusals.md).

## Examples

Print the masked intervals of each sequence:

<!-- example: dustmasker-1 -->
```bash
losat dustmasker -in genome.fna
```

Write the sequences with the low-complexity letters in lower case:

<!-- example: dustmasker-2 -->
```bash
losat dustmasker -in genome.fna -outfmt fasta
```

Use a smaller window and a lower level:

<!-- example: dustmasker-3 -->
```bash
losat dustmasker -in genome.fna -window 32 -level 15
```

Write the masking information as ASN.1 text, as `makeblastdb -mask_data` reads it:

<!-- example: dustmasker-4 -->
```bash
losat dustmasker -in genome.fna -parse_seqids -outfmt maskinfo_asn1_text
```

Mask a database you built earlier:

<!-- example: dustmasker-5 -->
```bash
losat makeblastdb -in genome.fna -dbtype nucl -out genome_db
losat dustmasker -infmt blastdb -in genome_db
```

## See also

- [segmasker](segmasker.md) for proteins
- [convert2blastmask](convert2blastmask.md)
- [makeblastdb](makeblastdb.md)
- [blastn](blastn.md): `-dust` masks the query during a search

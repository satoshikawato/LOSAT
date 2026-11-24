# segmasker

`losat segmasker` finds the low-complexity regions of protein sequences with the SEG algorithm. Use it to see which parts of a protein a search would mask, to write a masked FASTA file, or to make masks for [`makeblastdb -mask_data`](makeblastdb.md).

It reads a FASTA file or an existing protein BLAST database. It writes the same bytes as NCBI's `segmasker` 1.0.0 in all eight output formats.

## Synopsis

```bash
losat segmasker -in SEQUENCES.faa [options]
losat segmasker -infmt blastdb -in DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `segmasker -in SEQUENCES.faa [options]`.

## Input

- `-in` names the input. `-in -` (the default) reads the standard input.
- `-infmt fasta` (the default) reads FASTA. `-infmt blastdb` reads every sequence of a protein database named by `-in`, in every volume.
- `-parse_seqids` parses the ids of the FASTA definition lines, so the output names sequences by their ids. Use it when the masks are for a database built with `makeblastdb -parse_seqids`.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-in` | `-` | Input file or database name. |
| `-infmt` | `fasta` | `fasta` or `blastdb`. |
| `-out` | `-` | Output file. |
| `-outfmt` | `interval` | Output format (below). |
| `-window` | `12` | SEG window length. |
| `-locut` | `2.2` | SEG low cutoff. |
| `-hicut` | `2.5` | SEG high cutoff. |
| `-parse_seqids` | off | Parse the Seq-ids of FASTA input. |

For every option: `losat segmasker -help`.

`segmasker` runs the SEG of the `blastp` and `blastx` query filter, but does not merge overlapping segments, as NCBI's `segmasker`.

## Output

`-outfmt` is one of:

| Format | What you get |
| --- | --- |
| `interval` | For each sequence, its definition line, then one line `START - END` per masked interval. Positions count from 0 and include both ends. |
| `fasta` | The sequences with masked letters in lower case. |
| `seqloc_asn1_text`, `seqloc_asn1_bin`, `seqloc_xml` | The intervals as Seq-loc objects. A sequence with several Seq-ids is named, here and in `maskinfo_*`, by its gi (not gi 0), else its first accession-type id (GenBank, RefSeq, Swiss-Prot, ...), else its first Seq-id, as NCBI's does. |
| `maskinfo_asn1_text`, `maskinfo_asn1_bin`, `maskinfo_xml` | The intervals as masking information with the algorithm name and options. The `maskinfo_asn1_*` files are the input of `makeblastdb -mask_data`. |

## Differences from NCBI BLAST+

None that change the output of a normal run. When the reader of the output closes a pipe early, LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

## What LOSAT refuses

`segmasker` runs every option of `-help`. It stops with an error, where NCBI crashes or reads memory outside its buffers, in these cases:

- Everything that [`dustmasker`](dustmasker.md) refuses: a database title in an encoding NCBI cannot tell, a database sequence with a tax id from an empty or damaged `taxdb.btd`, and a failed write of a `seqloc_*` or `maskinfo_*` output.
- A protein database with residue bytes above 127. `makeblastdb` writes 0 to 27; NCBI's SEG reads larger bytes as negative numbers.
- Consistent with LOSAT's error wording ([Differences from NCBI BLAST+](../differences.md#wording-of-messages)), database errors keep NCBI's message and diagnostic codes without NCBI's C++ source paths and function headers.

See [What LOSAT refuses](../refusals.md).

## Examples

Print the masked intervals of each protein:

<!-- example: segmasker-1 -->
```bash
losat segmasker -in proteins.faa
```

Write the sequences with the low-complexity letters in lower case:

<!-- example: segmasker-2 -->
```bash
losat segmasker -in proteins.faa -outfmt fasta
```

Change the SEG window and cutoffs:

<!-- example: segmasker-3 -->
```bash
losat segmasker -in proteins.faa -window 8 -locut 1.5 -hicut 1.8
```

Write the masking information as ASN.1 text, as `makeblastdb -mask_data` reads it:

<!-- example: segmasker-4 -->
```bash
losat segmasker -in proteins.faa -parse_seqids -outfmt maskinfo_asn1_text
```

Mask a database you built earlier:

<!-- example: segmasker-5 -->
```bash
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db
losat segmasker -infmt blastdb -in proteins_db
```

## See also

- [dustmasker](dustmasker.md) for nucleotides
- [convert2blastmask](convert2blastmask.md)
- [makeblastdb](makeblastdb.md)
- [blastp](blastp.md): `-seg` masks the query during a search

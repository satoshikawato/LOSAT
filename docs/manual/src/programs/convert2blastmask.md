# convert2blastmask

`losat convert2blastmask` turns FASTA that another tool has masked into masking information for [`makeblastdb -mask_data`](makeblastdb.md). Use it when your masking program (RepeatMasker, WindowMasker or your own script) writes lower case letters for the masked regions and you want a BLAST database that remembers those masks.

It writes the same bytes as NCBI's `convert2blastmask` 2.17.0.

## Synopsis

```bash
losat convert2blastmask -in MASKED.fna -masking_algorithm NAME -masking_options TEXT [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `convert2blastmask -in MASKED.fna -masking_algorithm NAME -masking_options TEXT [options]`.

## Input

- `-in` names a FASTA file. `-in -` (the default) reads the standard input.
- The lower case letters of each sequence are the masked intervals.
- Instead of sequence data, a record can hold lines `FROM - TO`, one interval per line. This is the form that `-outfmt interval` writes, so you can edit intervals by hand.
- `-parse_seqids` parses the ids of the definition lines. Use it when the database will be built with `makeblastdb -parse_seqids`, so each mask finds its sequence.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-in` | `-` | Masked FASTA input. |
| `-out` | `-` | Output file. |
| `-outfmt` | `maskinfo_asn1_text` | `maskinfo_asn1_text`, `maskinfo_asn1_bin`, `maskinfo_xml` or `interval`. |
| `-masking_algorithm` | required | Name of the masking program: `dust`, `seg`, `windowmasker`, `repeat`, or `other` for your own. |
| `-masking_options` | required | Free text that describes how the masking was made, for example the command line. |
| `-parse_seqids` | off | Parse the Seq-ids of the FASTA input. |

For every option: `losat convert2blastmask -help`.

## Output

| Format | What you get |
| --- | --- |
| `maskinfo_asn1_text` | Masking information as ASN.1 text. The default. |
| `maskinfo_asn1_bin` | The same as binary ASN.1. |
| `maskinfo_xml` | The same as XML. |
| `interval` | The masked intervals as lines `FROM - TO`. Positions count from 0 and include both ends. |

The three `maskinfo_*` files are the input of `makeblastdb -mask_data`. The algorithm name and the text of `-masking_options` are stored in them. `blastdbcmd -info` shows them for the finished database, as the algorithm name (the `-mask_id` you gave to `makeblastdb`, if any) and its options.

## Differences from NCBI BLAST+

None that change the output of a normal run. A byte outside visible ASCII in `-masking_options` or in an id read with `-parse_seqids` is written as `#`, with NCBI's error message for each byte, as NCBI writes it. When the reader of the output closes a pipe early, LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

## What LOSAT refuses

`convert2blastmask` runs every option of `-help`. It stops with an error in these cases:

- A value that is not UTF-8 of any option but `-out` and `-logfile`, after NCBI's own argument errors and `-dryrun` (NCBI reads the bytes).
- Masking information that LOSAT's writer cannot write as NCBI does.
- A failed write of a `maskinfo_*` output. A failed write of `interval` goes on silently, as in NCBI.

See [What LOSAT refuses](../refusals.md).

## Examples

Write the masked intervals of lower case FASTA:

<!-- example: convert2blastmask-1 -->
```bash
losat convert2blastmask -in masked.fna -masking_algorithm other -masking_options "my masker" -outfmt interval
```

Write the masking information as ASN.1 text:

<!-- example: convert2blastmask-2 -->
```bash
losat convert2blastmask -in masked.fna -masking_algorithm other -masking_options "my masker"
```

Make masks for a database, build the database with them, and read the masks back:

<!-- example: convert2blastmask-3 -->
```bash
losat convert2blastmask -in masked.fna -parse_seqids -masking_algorithm other -masking_options "my masker" -out masked.asnt
losat makeblastdb -in masked.fna -dbtype nucl -parse_seqids -mask_data masked.asnt -mask_id my_masker -out masked_db
losat blastdbcmd -db masked_db -entry all -outfmt "%a %m100"
```

## See also

- [dustmasker](dustmasker.md), [segmasker](segmasker.md)
- [makeblastdb](makeblastdb.md)
- [BLAST databases](../databases.md)

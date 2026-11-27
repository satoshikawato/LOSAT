# windowmasker

`losat windowmasker` masks the repeats of nucleotide sequences with the WindowMasker algorithm, from a file of unit counts made for the genome. Use it to make the counts file of a genome, to see which parts of a genome a nucleotide search would mask, to write a masked FASTA file, or to make masks for [`makeblastdb -mask_data`](makeblastdb.md).

It has NCBI's three modes:

- `-mk_counts` makes a counts file from the sequences of a genome.
- `-convert` writes a counts file in another of the four formats (`ascii`, `binary`, `oascii`, `obinary`).
- `-ustat` masks sequences with a counts file.

It reads FASTA files or nucleotide BLAST databases, and counts files in all four formats. It writes the same bytes as NCBI's `windowmasker` 1.0.0: the same counts files, and the same masks in all eight output formats. So the whole workflow runs without NCBI BLAST+: counts from `-mk_counts` for [`blastn -window_masker_db`](blastn.md), masks for `makeblastdb -mask_data`. You can also use the counts files that NCBI publishes for many genomes (the files that `blastn -window_masker_taxid` reads).

## Synopsis

```bash
losat windowmasker -mk_counts -in GENOME.fna -sformat obinary -out GENOME.counts [options]
losat windowmasker -convert -in COUNTS -sformat ascii -out COUNTS.ascii
losat windowmasker -ustat COUNTS -in SEQUENCES.fna [options]
losat windowmasker -ustat COUNTS -infmt blastdb -in DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `windowmasker -ustat COUNTS -in SEQUENCES.fna [options]`.

## Making a counts file (`-mk_counts`)

- `-in` names the input: FASTA files separated by commas, or with `-infmt blastdb` nucleotide databases. With `-fa_list T`, `-in` is a file that lists the FASTA files, one per line. As in NCBI's, `-in -` is a file named `-`, not the standard input.
- `-unit` sets the length of the units, 1 to 16 letters. Without it, the length follows from `-genome_size`, or from the length of the input, which is then read once more.
- `-sformat` sets the format (`ascii`, the default; `binary`; `oascii` and `obinary`, the optimized formats that `blastn -window_masker_db` reads fastest). `-smem` is the memory in megabytes that the optimized formats may use.
- `-t_low` and `-t_high` bound the counts that the file keeps. `-t_low_pct`, `-t_extend_pct`, `-t_thres_pct` and `-t_high_pct` set the percentages of the four thresholds; as in NCBI's, they apply only when `-out` names a file.
- `-mem` is the memory in megabytes for counting. With less memory the input is read more times; the counts are the same.
- `-ids` and `-exclude_ids` choose the sequences, as in masking. `-meta` adds a note to the file.
- `-checkdup T` first checks the input for duplicated sequences, and writes each possible duplication to the standard error.
- The progress lines (`pass 1`, `pass 2`, ...) go to the standard error, as NCBI writes them.

## Converting a counts file (`-convert`)

`-in` names the counts file, in any of the four formats; `-sformat` and `-smem` set the format written to `-out`. `-meta` replaces the note of the file. The `binary`, `oascii` and `obinary` files that `-convert` writes are the same as those that `-mk_counts` writes from the same counts; an `ascii` file has no statistics lines.

## Masking (`-ustat`)

- `-ustat` names the counts file. It is required.
- `-in` names the input. `-in -` (the default) reads the standard input.
- `-infmt fasta` (the default) reads FASTA. `-infmt blastdb` reads every sequence of a nucleotide database named by `-in`, in every volume. `-infmt seqids` stops with NCBI's error, as NCBI's does.
- `-parse_seqids` parses the ids of the FASTA definition lines, so the output names sequences by their ids. Use it when the masks are for a database built with `makeblastdb -parse_seqids`.
- `-ids FILE` masks only the sequences that match one of the ids in `FILE` (one per line); `-exclude_ids FILE` masks all others. As in NCBI's, an id matches the first word of a sequence's title, whole or by its `|`-separated parts. The title is NCBI's: HTML entities such as `&amp;` are not decoded and only the dots and spaces at its end are dropped; a third-party id (`tpg`, `tpe`, `tpd`) adds `TPA: ` before it; a database sequence without a title and with a tax id takes the scientific name of its organism from the taxonomy database (`taxdb`), with its first letter in upper case, or `No definition line found` without one. Matching is case-sensitive. The other sequences are left out of the output.
- `-text_match F` (with `-infmt blastdb` only) reads each id of the `-ids` or `-exclude_ids` file as a Seq-id, such as `gb|AB000001.1|` or `AB000001.1`, and matches it against the Seq-ids of each database sequence, as NCBI's does: case does not matter (except in a PDB chain), the versions must agree (`AB000001` does not match `AB000001.1`), and an accession alone does not match an id stored with a name (`P12345` does not match `sp|P12345|ABC_HUMAN`). A database made with `-parse_seqids` stores a FASTA id that is a bare word, such as `rec1`, as `lcl|rec1`; in a database made without it, the second sequence is `gnl|BL_ORD_ID|1`, the third `gnl|BL_ORD_ID|2`, and so on; `gnl|BL_ORD_ID|0` matches no sequence, as in NCBI's. An id that is not a Seq-id, such as a bare word, is written to the standard error as an error and skipped; when no id of the file can be read, every sequence is masked. With `-outfmt seqloc_*` or `maskinfo_*`, a sequence's Seq-id is written as NCBI's writes it once the list has named it: `AB000001.1` in the list is a DDBJ accession, so `gb|AB000001.1|`, and every other GenBank, EMBL or DDBJ accession with the same letters, number of digits and version, is written as a DDBJ id, and a name is written in the list's case.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-mk_counts`, `-convert` | off | Make a counts file, or convert one. |
| `-ustat` | none | The counts file to mask with. |
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
| `-sformat` | `ascii` | Format of the counts that `-mk_counts` and `-convert` write. |
| `-unit` | from the genome | Unit length for `-mk_counts`. |
| `-checkdup` | `F` | Check the input of `-mk_counts` for duplicated sequences. |

For every option: `losat windowmasker -help`.

## Output

The counts files are those of NCBI's `windowmasker`, byte for byte, in all four formats. As in NCBI's, `binary` counts written to the standard output start with an extra word, which the readers do not accept: write them to a file.

The masks: `-outfmt` is one of:

| Format | What you get |
| --- | --- |
| `interval` | For each sequence, its definition line, then one line `START - END` per masked interval. Positions count from 0 and include both ends. |
| `fasta` | The sequences with masked letters in lower case. |
| `seqloc_asn1_text`, `seqloc_asn1_bin`, `seqloc_xml` | The intervals as Seq-loc objects. A sequence with several Seq-ids is named, here and in `maskinfo_*`, by its gi (not gi 0), else its first accession-type id (GenBank, RefSeq, Swiss-Prot, ...), else its first Seq-id, as NCBI's does. |
| `maskinfo_asn1_text`, `maskinfo_asn1_bin`, `maskinfo_xml` | The intervals as masking information. The `maskinfo_asn1_*` files are the input of `makeblastdb -mask_data`. |

Warnings about the input, such as a counts file without percentages when a `-t_*_pct` option is given, are written to the standard error, as NCBI writes them.

## Differences from NCBI BLAST+

None that change the output of a normal run. When the reader of the output closes a pipe early, LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

A line of an `-ids` or `-exclude_ids` file (matched as text, without `-text_match F`) that gives an empty id, such as a line that starts with a space or a tab, or a line `>` or `|`, gets NCBI's error `bad id: ...: ignoring` and is skipped; the rest of the output is NCBI's for the file without that line. NCBI writes the same error and then crashes. See [`windowmasker`: an empty id in `-ids` or `-exclude_ids`](../differences.md#windowmasker-an-empty-id-in--ids-or--exclude_ids).

## What LOSAT refuses

- Counts files on which NCBI's result is undefined, such as an `ascii` file with units of more than 16 letters.
- `-mk_counts` of an input without a unit of unambiguous letters (all `N`; or, with `-unit`, no sequence, as when `-ids` keeps none): NCBI writes its percentages as its C runtime prints 0/0.
- `-mk_counts -infmt seqids`: NCBI fetches each sequence from NCBI's servers.
- A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more (such as 4,294,967,293 letters of `A`, `C`, `G` and `T` in a row), where NCBI crashes. Other records of 4,294,967,296 letters or more are read as NCBI reads them.
- `-mk_counts` or `-convert` of units of 1 to 3 letters to `oascii` or `obinary` with an `-smem` that is a multiple of 65536: NCBI's result is undefined.

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

Make the counts file of a genome in the optimized binary format, and mask a search's queries with it:

<!-- example: windowmasker-5 -->
```bash
losat windowmasker -mk_counts -in genome.fna -sformat obinary -out genome.counts
losat blastn -query reads.fna -subject genome.fna -window_masker_db genome.counts
```

Write a counts file as text:

<!-- example: windowmasker-6 -->
```bash
losat windowmasker -convert -in genome.counts -sformat ascii -out genome.counts.txt
```

## See also

- [dustmasker](dustmasker.md)
- [makeblastdb](makeblastdb.md)
- [blastn](blastn.md): `-window_masker_db` masks the query with a counts file during a search

# makembindex

`losat makembindex` builds a MegaBLAST index of a nucleotide BLAST database: a table of where each short word occurs in the database, which `blastn -task megablast -use_index true` reads in place of scanning the database. LOSAT's `blastn` searches with the index and gives NCBI's output, except for an index of only some database volumes, a two-hit search and `-db_hard_mask`, which it refuses (see [What LOSAT refuses](../refusals.md)).

It writes the same files as NCBI's `makembindex` 2.17.0, byte for byte.

## Synopsis

```bash
losat makembindex -input DATABASE [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `makembindex -input DATABASE [options]`.

## Input

- `-input` names a nucleotide BLAST database, as `-db` does for a search. An alias database works too; each of its volumes gets its own index.
- `-db_mask ID` leaves the regions masked by the database's masking algorithm `ID` out of the index. `-show_filters` lists the algorithms the database has, with their ids. Give the number: as in NCBI 2.17.0, a name such as `dust` is not found and the run stops with NCBI's error.
- With `-old_style_index true`, the input can also be FASTA (`-iformat fasta`, the default of that style): `-input FILE`, or the standard input without `-input`. Lower case letters are masked and left out of the index.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-input` | none | The database (or, with `-old_style_index true -iformat fasta`, the FASTA file). |
| `-nmer` | `12` | Length of the indexed words, 8 to 15. |
| `-volsize` | `1536` | Largest size of one index volume, in megabytes. A larger database is split into several volumes. |
| `-db_mask` | none | Id of the database's masking algorithm whose regions are left out. |
| `-show_filters` | off | List the database's masking algorithms and stop. |
| `-legacy` | `true` | The index format that `blastn -use_index` reads. `-legacy false` writes the newer format, with `-stride` and `-ws_hint`. |
| `-stride` | `5` | With `-legacy false`: distance between the database positions that are stored. |
| `-ws_hint` | `28` | With `-legacy false`: the word size you expect to search with. A value below `-nmer` + `-stride` − 1 is raised to it, with a warning. |
| `-idmap` | `false` | Also write `VOLUME.00.idx.map` (one per index volume), the first word of each sequence's title as NCBI makes it, one per line (for a database sequence without a title and with a tax id, the first word of its organism's scientific name from the taxonomy database, or `No`). |
| `-old_style_index` | `false` | One index for the whole input, named by `-output` (deprecated in NCBI). |
| `-output` | none | With `-old_style_index true`: the name of the index files. Ignored otherwise, with a warning. |
| `-iformat` | `blastdb` (`fasta` for the old style) | Input format. |
| `-stat` | none | Write the size of each word's list to this file. |

For every option: `losat makembindex -help`.

## Output

For each volume of the database, `makembindex` writes next to the volume's files:

| File | What it holds |
| --- | --- |
| `VOLUME.00.idx`, `VOLUME.01.idx`, ... | The index, in volumes of at most `-volsize` megabytes. |
| `VOLUME.shd` | The number of sequences and of index volumes. |
| `VOLUME.00.idx.map` | With `-idmap true`: the sequences of that index volume. |

With `-old_style_index true` the files are `OUTPUT.00.idx`, `OUTPUT.01.idx`, ... and there is no `.shd` file.

The standard error shows each file as it is made: `creating NAME...done`. The last line, `creating NAME...removed (empty)`, is the check that no sequence is left; no file is written for it.

The index of the default `-nmer 12` has 64 MB of word table in every volume, whatever the size of the database.

## Differences from NCBI BLAST+

None in the files or the messages of a normal run.

## What LOSAT refuses

`makembindex` runs every option of `-help` (`-xmlhelp` on the command line only). It stops with an error in these cases:

- `-legacy false` with a `-stride` above `-ws_hint` − `-nmer` + 1 (for example `-stride 20` with the default `-ws_hint 28` and `-nmer 12`). NCBI's `makembindex` never finishes such an index; LOSAT stops before writing it.
- `-legacy false` with a `-stride` of 32768 or more, where NCBI's index has no defined result.
- An index volume of more than 4 GiB of packed sequences (a `-volsize` of about 18,000 or more and more than 16 Gbp of sequences in one index volume), where NCBI's index sizes wrap at 32 bits.
- A command-line word that is not UTF-8.
- The database cases that `dustmasker` refuses, such as a definition line with bytes NCBI reads as `#` or a tax id with an empty `taxdb.btd`.
- A database sequence without letters, which only a damaged database has: NCBI's `makembindex` stops with a null-pointer report.

See [What LOSAT refuses](../refusals.md).

## Examples

Index a database with the default settings:

<!-- example: makembindex-1 -->
```bash
losat makeblastdb -in genome.fna -dbtype nucl -out genome
losat makembindex -input genome
```

Leave the low-complexity regions out of the index:

<!-- example: makembindex-2 -->
```bash
losat dustmasker -in genome.fna -outfmt maskinfo_asn1_bin -out genome.dust.asnb
losat makeblastdb -in genome.fna -dbtype nucl -mask_data genome.dust.asnb -out genome
losat makembindex -input genome -show_filters
losat makembindex -input genome -db_mask 11
```

## See also

- [makeblastdb](makeblastdb.md)
- [dustmasker](dustmasker.md)
- [blastn](blastn.md)
- [BLAST databases](../databases.md)

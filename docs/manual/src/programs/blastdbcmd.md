# blastdbcmd

`losat blastdbcmd` reads a BLAST database. Use it to see what a database holds (`-info`), to list the databases in a directory, to write sequences as FASTA or in a format you choose (by id, by position, or all of them), to cut a range or the minus strand, to read taxonomy ids and masks, and to check a database before you search it.

It reads databases made by NCBI's `makeblastdb` or by [`makeblastdb`](makeblastdb.md) of LOSAT, in format versions 4 and 5, including NCBI-distributed databases such as swissprot and pdbaa.

## Synopsis

```bash
losat blastdbcmd -db DB -info
losat blastdbcmd -db DB -entry ID[,ID...] [options]
losat blastdbcmd -db DB -entry all [options]
losat blastdbcmd -db DB -entry_batch FILE [options]
losat blastdbcmd -list DIRECTORY
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blastdbcmd -db DB -info` and so on.

## Input

- `-db` names the database. Without `-db` the name is `nr`, as in NCBI. The name is looked up in the current directory, then in `BLASTDB`, then in `.ncbirc` (see [BLAST databases](../databases.md)). `-db "a b"` (a name with spaces, in double quotes) reads several databases as one.
- `-dbtype` is `nucl`, `prot` or `guess` (the default). Give it when two databases of one name exist.
- `-entry` takes ids separated by commas, or `all`. An id can be an accession, an accession with a version, a name, a local id (`lcl|name`), a GI, or, in a version 4 database, `gnl|BL_ORD_ID|N` (the position `N`). Ids are found only in a database made with `-parse_seqids` (NCBI-distributed databases have them). In a database without parsed ids, use `-entry all`.
- `-entry_batch` reads one entry per line from a file, or from the standard input with `-`. A line can add a range, a strand, or a mask algorithm id after the id, separated by spaces.
- `-range START-STOP` (1-based) and `-strand plus|minus` select a part of one sequence. `-range 10-` runs to the end.
- `-taxids`, `-taxidlist` and `-no_taxid_expansion` select sequences by tax id.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-db` | `nr` | Database name. |
| `-dbtype` | `guess` | `nucl`, `prot` or `guess`. |
| `-info` | off | Print the title, number of sequences, length, date, longest sequence, format version, masks and volumes. |
| `-entry` | none | Sequence ids, comma-separated, or `all`. |
| `-entry_batch` | none | File with one entry per line (id, then optional range, strand, mask id). |
| `-range` | whole sequence | Part of the sequence, `START-STOP`. |
| `-strand` | `plus` | `plus` or `minus` for a nucleotide sequence. |
| `-outfmt` | `%f` | FASTA, or a text with `%` specifiers (below). |
| `-out` | standard output | Output file. |
| `-line_length` | `80` | Letters per line in FASTA. |
| `-target_only` | off | Write only the definition line of the entry you asked for. |
| `-get_dups` | off | Also write sequences with duplicate accessions. |
| `-mask_sequence_with` | none | Write masked letters in lower case, using the numeric id of a masking algorithm. |
| `-tax_info` | off | Print the taxonomy information in the database. |
| `-list` | none | List the databases in a directory. |
| `-show_blastdb_search_path` | off | Print the directories searched for databases. |

For every option: `losat blastdbcmd -help`.

### Output specifiers

`-outfmt` takes a text in which these specifiers are replaced for each sequence (one line per sequence). `%f` (FASTA) is the default; if you give `%f`, other specifiers are ignored.

| Specifier | Meaning |
| --- | --- |
| `%f` | Sequence in FASTA format. |
| `%s` | Sequence letters, without the definition line. |
| `%a` | Accession. |
| `%g` | GI. |
| `%o` | Ordinal id (position in the database, from 0). |
| `%i` | Sequence id. |
| `%t` | Sequence title. |
| `%l` | Sequence length. |
| `%h` | Sequence hash value. |
| `%T`, `%X` | Tax id; tax ids of the leaf nodes. |
| `%S`, `%N` | Scientific name; scientific names of the leaf nodes. |
| `%L`, `%C` | Common name; common names of the leaf nodes. |
| `%B`, `%K` | BLAST name; super kingdom. |
| `%e`, `%P` | Membership integer; PIG. |
| `%m<id>` | Masking data of algorithm `<id>`, as `N-M` ranges separated by `;`, or `none`. |

The names (`%S`, `%L`, `%C`, `%N`, `%B`, `%K`) and `-tax_info` names need NCBI's taxonomy database files; without them LOSAT prints `N/A`, as NCBI does (see [BLAST databases](../databases.md)). `%d` and `%b` (ASN.1 texts) are refused.

`-list DIRECTORY` prints one line per database in the directory. `-list_outfmt` changes the line: `%f` path, `%p` molecule type, `%t` title, `%d` date, `%l` letters, `%n` sequences, `%U` bytes, `%v` format version. `-recursive` goes into subdirectories and `-remove_redundant_dbs` hides databases that an alias file already names.

## Output

`blastdbcmd` writes text: FASTA, the lines of `-outfmt`, the `-info` report, the `-list` lines or the `-tax_info` table. It has no BLAST `-outfmt` numbers. An id that is not found is reported on the standard error, the other entries are still written, and the exit status is 1. With `-entry_batch` the message is `Error: [blastdbcmd] Skipped NOSUCH` (one line for each id not found). With `-entry A,NOSUCH` it is `Error: [blastdbcmd] Entry not found: NOSUCH`. If every id of `-entry` is missing, a second line follows: `Error: [blastdbcmd] Entry or entries not found in BLAST database`.

## Differences from NCBI BLAST+

LOSAT writes six outputs correctly where NCBI's `blastdbcmd` 2.17.0 writes letters or sequences that do not belong to the header or id you asked for. They are all edge cases of masks, ranges and ordinal ids:

- `-long_seqids` with `-range` on the minus strand.
- `gnl|BL_ORD_ID|N` in a version 4 database of several volumes (NCBI writes one sequence per volume).
- `-mask_sequence_with` together with `-strand minus`.
- Two or more `%m` specifiers in one `-outfmt`.
- `-mask_sequence_with` on a protein without a range (NCBI masks only the first third).
- A batch line without masks after a line with masks for the same sequence.

See [the decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

Two more differences are about NCBI failing:

- An accession stored under more than about 500 sequences of a version 5 database stops NCBI with an error; LOSAT looks it up like any other accession.
- `-taxids 1` (the root) makes NCBI search without end; LOSAT expands it and ends.

See [the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md). When the reader of the output closes a pipe early (`blastdbcmd ... | head`), LOSAT stops at once with exit status 0, where NCBI ends with the signal (141). See [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md). The full list is on [Differences from NCBI BLAST+](../differences.md).

## What LOSAT refuses

- `-outfmt` with `%d` or `%b`.
- A `-db` list in which an alias file and another name give the same volume.
- A taxonomy lookup that reads a record of an empty `taxdb.btd`, or a name from a damaged `taxdb.btd` whose record reaches past the end of NCBI's mapping of the file (NCBI crashes or reads memory there). `-info` and formats that look up nothing (`-outfmt "%a %T"`) still work.
- A failed write of `-metadata` output.

All other `blastdbcmd` options work as in NCBI; the full list with messages is on [What LOSAT refuses](../refusals.md).

## Examples

Show what a database holds: the title, the number of sequences, the total length, the date and the format version. The report names the database's path and build date:

<!-- example: blastdbcmd-1 -->
```bash
losat blastdbcmd -db proteins_db -info
```

Print the accession, length and tax id of every sequence, one line each:

<!-- example: blastdbcmd-2 -->
```bash
losat blastdbcmd -db proteins_db -entry all -outfmt "%a %l %T"
```

Write the first 60 letters of the sequence with a given id as FASTA:

<!-- example: blastdbcmd-3 -->
```bash
losat blastdbcmd -db proteins_db -entry XP_000001.1 -range 1-60
```

Write the entries listed in a file, one id per line. If two ids name the same sequence, it is written once per id. An id that the database does not have is reported on the standard error as `Error: [blastdbcmd] Skipped NOSUCH`, and the command exits with status 1, as NCBI's `blastdbcmd` does:

<!-- example: blastdbcmd-4 -->
```bash
losat blastdbcmd -db proteins_db -entry_batch ids.txt
```

Write a masked protein in lower case, and list the mask ranges of every sequence (algorithm id 21 is shown by `-info`). Here the output differs from NCBI BLAST+ 2.17.0 on purpose: without `-range`, NCBI's `blastdbcmd` lowers only the masks in the first third of a protein, while LOSAT lowers all of them (see the list of differences above):

<!-- example: blastdbcmd-5 -->
```bash
losat blastdbcmd -db masked_db -entry P10001.1 -mask_sequence_with 21
losat blastdbcmd -db masked_db -entry all -outfmt "%a %m21"
```

Write a nucleotide range from the minus strand:

<!-- example: blastdbcmd-6 -->
```bash
losat blastdbcmd -db nucl_db -entry NM_100001.1 -strand minus -range 1-50
```

## See also

- [BLAST databases](../databases.md)
- [makeblastdb](makeblastdb.md)
- [What LOSAT refuses](../refusals.md)

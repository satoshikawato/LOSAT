# blast_formatter

`losat blast_formatter` writes the reports of a search again from its BLAST archive, without running the search again. Use it when you saved a search with `-outfmt 11` and want its Seq-aligns later, in ASN.1 or JSON.

It writes the same bytes as NCBI's `blast_formatter` 2.17.0, for archives of every search program with `-subject` or `-db`.

## Synopsis

```bash
losat blast_formatter -archive RESULTS.asn -outfmt 8 [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blast_formatter -archive RESULTS.asn -outfmt 8 [options]`.

## Input

- `-archive` names the BLAST archive: the file that a search wrote with `-outfmt 11`, from LOSAT or from NCBI BLAST+. `-archive -` reads the standard input.
- The archive may be in ASN.1 text (the default of `-outfmt 11`), ASN.1 binary (`ARCHIVE_FORMAT=bin`) or XML (`ARCHIVE_FORMAT=xml`).
- The reports are written for each query batch of the archive, one after another.

## Common options

| Option | Meaning |
| --- | --- |
| `-archive` | The BLAST archive to read. |
| `-outfmt` | `8`, `9` or `12`. |
| `-out` | Output file. |
| `-max_target_seqs`, `-num_descriptions`, `-num_alignments` | How many subject sequences, descriptions and alignments the reports keep, as in NCBI. |
| `-sorthits`, `-sorthsps`, `-line_length`, `-show_gis`, `-html`, `-parse_deflines` | Read as in NCBI. |

`-h`, `-help`, `-help-full`, `-version` (and the full versions), `-xmlhelp` (command line only), `-dryrun`, `-logfile` and `-conffile` work as in NCBI. For every option: `losat blast_formatter -help`.

## Output

| `-outfmt` | What you get |
| --- | --- |
| `8` | The Seq-aligns as ASN.1 text. |
| `9` | The Seq-aligns as binary ASN.1. |
| `12` | The Seq-aligns as JSON. |

## Differences from NCBI BLAST+

None that change the output of a normal run. The options, NCBI's help and usage texts, argument errors, messages and exit statuses are NCBI's, also where the output cannot be written: `-out /dev/full` gives NCBI's messages and exit status, and a standard output on a full device or a closed pipe ends the program as NCBI's does (`terminate called after throwing an instance of 'std::__ios_failure'`, exit status 134), or, after a standard error that could not be written, with exit status 6, as NCBI's. See [Differences from NCBI BLAST+](../differences.md).

## What LOSAT refuses

`blast_formatter` stops with its own error message in these cases:

- `-outfmt` `0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `10`, `11`, `13`, `14`, `15`, `16`, `17`, `18` and `20`.
- A database named in the archive that cannot be opened here as NCBI opens it. NCBI BLAST+ then asks NCBI's servers for its sequences.
- `PRE_FETCH_SEQS_LIMIT=2147483647` with `-outfmt` `8`, `9` or `12`: NCBI BLAST+ then loads every subject sequence ahead, from NCBI's servers where the database lacks one.
- Archives that NCBI's own programs do not write and on which NCBI BLAST+ has no defined result, a stored message of severity Fatal, and archives of the mapper.
- A value in an archive's results that LOSAT does not yet write in the output format, such as a NaN or an infinite real number in `-outfmt` `12`.
- `-rid`: the results of a search that NCBI's servers ran. The search programs refuse `-remote` in the same way.

An empty `-rid`, the archive's other errors and the errors of `-outfmt` come as in NCBI. See [What LOSAT refuses](../refusals.md).

## Example

Write the Seq-aligns of a search again from its archive:

```bash
losat blastn -query query.fna -subject subject.fna -outfmt 11 -out results.asn
losat blast_formatter -archive results.asn -outfmt 8 -out results.seqalign.asn
```

## Library, Python, and LOSAT Web

The library API runs it as `losat::api::run_blast_formatter(argv, input, output, diagnostics)`, like `run_convert2blastmask`. In Python: `pylosat.blast_formatter(*words, input=b"")` returns a `ToolOutput`, like `pylosat.convert2blastmask` (with `input` serving `-archive -`). In the browser, the `losat` npm package runs `blast_formatter` through `engine.runTool(['blast_formatter', ...], { directory })` (supporting both `-subject` archives and local `-db` archives whose database files are provided in `directory`), and LOSAT Web exposes it in the **Format Archive** tab (and via `"Format archive…"` on completed `-outfmt 11` search outputs). See [Using LOSAT as a library](../library.md) and [LOSAT Web](../web.md).

## See also

- [blastn](blastn.md) and the other search programs, which write the archive with `-outfmt 11`
- [Output formats](../output-formats.md)

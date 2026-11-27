# blast_formatter

`losat blast_formatter` writes the reports of a search again from its BLAST archive, without running the search again. Use it when you saved a search with `-outfmt 11` and want its pairwise or query-anchored report, its tabular report or its Seq-aligns later.

It writes the same bytes as NCBI's `blast_formatter` 2.17.0: the Seq-aligns for archives of every search program, the pairwise report, the query-anchored reports (also with `-html`) and the tabular formats for archives of `blastn`, `blastp`, `psiblast`, `deltablast`, PHI-BLAST of a protein query (`psiblast -phi_pattern`) and `rpsblast`, with `-subject` or `-db`.

## Synopsis

```bash
losat blast_formatter -archive RESULTS.asn [-outfmt 0] [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `blast_formatter -archive RESULTS.asn [options]`.

## Input

- `-archive` names the BLAST archive: the file that a search wrote with `-outfmt 11`, from LOSAT or from NCBI BLAST+. `-archive -` reads the standard input.
- The archive may be in ASN.1 text (the default of `-outfmt 11`), ASN.1 binary (`ARCHIVE_FORMAT=bin`) or XML (`ARCHIVE_FORMAT=xml`).
- The reports are written for each query batch of the archive, one after another.

## Common options

| Option | Meaning |
| --- | --- |
| `-archive` | The BLAST archive to read. |
| `-outfmt` | `0` (the default) to `4`, `6`, `7`, `10` or `20` (with any field list and delimiter), `8`, `9` or `12`. |
| `-out` | Output file. |
| `-max_target_seqs`, `-num_descriptions`, `-num_alignments` | How many subject sequences, descriptions and alignments the reports keep, as in NCBI. |
| `-sorthits`, `-sorthsps`, `-line_length`, `-show_gis`, `-html`, `-parse_deflines` | Read as in NCBI. |

`-h`, `-help`, `-help-full`, `-version` (and the full versions), `-xmlhelp` (command line only), `-dryrun`, `-logfile` and `-conffile` work as in NCBI. For every option: `losat blast_formatter -help`.

## Output

| `-outfmt` | What you get |
| --- | --- |
| `0` (the default) | The pairwise report: the description table, each alignment and the statistics that the archive keeps, as NCBI's `blast_formatter` writes them (`-html` gives it in HTML, with NCBI's links). Like NCBI's, it is not always the search program's own `-outfmt 0` report: the archive does not keep everything the search printed, such as PSI-BLAST's sequences found in an earlier round. The subjects of `-db` archives are read from the database the archive names. |
| `1` to `4` | The query-anchored reports (with and without identities, with and without the query anchored), as for `0`. |
| `6`, `7`, `10`, `20` | The tabular report: tab-separated, with comment lines, comma-separated, comma-separated with a header line; any field list (`-outfmt "6 qaccver saccver pident qcovs"`). The subjects of `-db` archives are read from the database the archive names. |
| `8` | The Seq-aligns as ASN.1 text. |
| `9` | The Seq-aligns as binary ASN.1. |
| `12` | The Seq-aligns as JSON. |

## Differences from NCBI BLAST+

None that change the output of a normal run. The options, NCBI's help and usage texts, argument errors, messages and exit statuses are NCBI's, also where the output cannot be written: `-out /dev/full` gives NCBI's messages and exit status, and a standard output on a full device or a closed pipe ends the program as NCBI's does (`terminate called after throwing an instance of 'std::__ios_failure'`, exit status 134). See [Differences from NCBI BLAST+](../differences.md).

## What LOSAT refuses

`blast_formatter` stops with its own error message in these cases:

- `-outfmt` `5`, `11`, `13`, `14`, `15`, `16`, `17` and `18`.
- The tabular formats (`6`, `7`, `10`, `20`) of an archive of a translated search (`blastx`, `tblastn`, `tblastx`, `rpstblastn`) or of queries given as locations, `-subject` sequences with `BL2SEQ_LEGACY` set in `-outfmt 7`, and Seq-aligns whose sequences NCBI BLAST+'s scope may find otherwise; each with the case named.
- A database named in the archive that cannot be opened here as NCBI opens it. NCBI BLAST+ then asks NCBI's servers for its sequences.
- The pairwise and query-anchored reports (`-outfmt` `0` to `4`) of an archive of a translated search (`blastx`, `tblastn`, `tblastx`, `rpstblastn`) or of queries given as locations, of `-subject` sequences with `BL2SEQ_LEGACY` set (NCBI BLAST+ then writes the legacy bl2seq report, which LOSAT's search programs refuse too), of Seq-aligns in a form that NCBI's search programs do not write or whose Seq-ids do not read past those that the description table reads, and in `-outfmt` `1` to `4` of a query Bioseq among its own subjects (NCBI BLAST+ aborts on it); each with the case named.
- `PRE_FETCH_SEQS_LIMIT=2147483647` with `-outfmt` `0` to `4`, `8`, `9` or `12`: NCBI BLAST+ then loads every subject sequence ahead, from NCBI's servers where the database lacks one.
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

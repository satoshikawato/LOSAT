# Output formats

`-outfmt` chooses how a search writes its results. LOSAT reads `-outfmt` as NCBI BLAST+ 2.17.0 does, and writes the same bytes, apart from the places that name the program's version, where LOSAT writes its own name and version: `BLASTN (LOSAT 0.2.0)` where NCBI writes `BLASTN 2.17.0+` (the first line of `0` to `4`, the `# BLASTN` line of `7`, the version of `5` and of `13` to `16`), and `LOSAT 0.2.0` where NCBI writes `2.17.0+` alone (`ident` in `11`, `VN:` in `17`). See [Differences from NCBI BLAST+](differences.md#the-version-line-and-the-error-prefixes-name-losat). The default is `0`.

A value that NCBI rejects gets NCBI's error. A format that NCBI writes but a LOSAT program does not write is refused before the inputs are read, with `output format N is not supported by LOSAT's PROGRAM`. Numbers `19` and `21` belong to NCBI's `igblastn` and `magicblast` and give NCBI's error in every program.

## Which program writes which format

"The five programs" means `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`. Where a program writes a format, it writes it for both `-subject` and `-db` searches.

| `-outfmt` | Format | Written by |
| --- | --- | --- |
| `0` | Pairwise | all eight search programs |
| `1` | Query-anchored, showing identities | all but `psiblast` |
| `2` | Query-anchored, no identities | all but `psiblast` |
| `3` | Flat query-anchored, showing identities | all but `psiblast` |
| `4` | Flat query-anchored, no identities | all but `psiblast` |
| `5` | BLAST XML | all but `psiblast` |
| `6` | Tabular | all eight |
| `7` | Tabular with comment lines | all eight |
| `8` | Seqalign, ASN.1 text | all but `psiblast` |
| `9` | Seqalign, ASN.1 binary | all but `psiblast` |
| `10` | Comma-separated values | all eight |
| `11` | BLAST archive (ASN.1 text, ASN.1 binary or XML) | all but `psiblast` |
| `12` | Seqalign, JSON | all but `psiblast` |
| `13` | Multiple-file BLAST JSON (needs `-out`) | all but `psiblast` |
| `14` | Multiple-file BLAST XML2 (needs `-out`) | all but `psiblast` |
| `15` | Single-file BLAST JSON | all but `psiblast` |
| `16` | Single-file BLAST XML2 | all but `psiblast` |
| `17` | SAM | `blastn` only |
| `18` | Organism report | all but `psiblast` |
| `20` | Comma-separated values with header lines | all eight |

"All eight" is the five programs, `psiblast`, `rpsblast` and `rpstblastn`. Some programs write fewer formats than the table shows:

- `psiblast` writes `0`, `6`, `7`, `10` and `20`; with `-in_pssm` or `-in_msa`, `0`, `6` and `7`.
- `tblastn` with `-in_pssm` writes `0` to `4`, `6`, `7` and `10`.

`rpsbproc` is not a search: it writes its own tab-delimited report (see [rpsbproc](programs/rpsbproc.md)). `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` have their own formats, listed on their pages.

## Pairwise and query-anchored text (0 to 4)

`0` is the classic BLAST report: a header, a table of one-line descriptions, then the alignments. `1` to `4` are NCBI's query-anchored views: the alignments of all subjects are stacked under the query. `1` and `3` print a dot where a subject letter matches the query; `2` and `4` print the letters. `3` and `4` are NCBI's flat variants.

`-num_descriptions`, `-num_alignments`, `-line_length`, `-sorthits` and `-sorthsps` shape these formats. For formats above `4` they are ignored with a warning, as in NCBI.

## Tabular and custom fields (6, 7, 10, 20)

All four write one line per alignment. They differ only in layout:

| Format | Layout |
| --- | --- |
| `6` | Tab-separated, no header. |
| `7` | Tab-separated, with comment lines: the program, the query, the database, the field names, and the number of hits. |
| `10` | Comma-separated, no header. |
| `20` | Comma-separated, with one header line of field names. |

With no field list, the 12 default columns are written:

```text
qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

`std` stands for these 12. To choose columns, put them after the number, separated by spaces:

```bash
losat blastn -query query.fna -subject genome.fna -outfmt "6 qseqid sseqid pident evalue"
```

To change the separator, add `delim=` right after the number (before the field names):

```bash
losat blastn -query query.fna -subject genome.fna -outfmt "10 delim=@ qacc sacc score"
```

As in NCBI BLAST+, `psiblast`, `rpsblast` and `rpstblastn` read and check `delim=` but keep the tab or comma. In every program, a separator that a listed field uses between its entries (`;` with `sallseqid`, `staxids`, `sscinames`, `scomnames`, `sblastnames` or `sskingdoms`, `/` with `frames`, `<>` with `salltitles`) gets NCBI's error.

LOSAT ignores a word that is not a field name; a list with no known field gives the default columns.

### Field names

| Field | Meaning |
| --- | --- |
| `qseqid`, `qgi`, `qacc`, `qaccver`, `qlen` | Query id, GI, accession, accession with version, length. |
| `sseqid`, `sallseqid`, `sgi`, `sallgi`, `sacc`, `saccver`, `sallacc`, `slen` | Subject id, all ids, GI, all GIs, accession, accession with version, all accessions, length. |
| `qstart`, `qend`, `sstart`, `send` | Start and end of the alignment on the query and the subject. |
| `qseq`, `sseq` | Aligned part of the query and the subject. |
| `evalue`, `bitscore`, `score` | E value, bit score, raw score. |
| `length`, `pident`, `nident`, `mismatch`, `positive`, `ppos` | Alignment length, percent identity, identical matches, mismatches, positive-scoring matches, percent positives. |
| `gapopen`, `gaps` | Gap openings, total gaps. |
| `frames`, `qframe`, `sframe` | Query and subject frames (separated by `/`), query frame, subject frame. |
| `btop` | Traceback operations (BTOP). |
| `staxid`, `ssciname`, `scomname`, `sblastname`, `sskingdom` | Subject tax id, scientific name, common name, BLAST name, super kingdom. |
| `staxids`, `sscinames`, `scomnames`, `sblastnames`, `sskingdoms` | The unique values for all tax ids of the subject, separated by `;`. |
| `stitle`, `salltitles` | Subject title, all subject titles (separated by `<>`). |
| `sstrand` | Subject strand. |
| `qcovs`, `qcovhsp`, `qcovus` | Query coverage per subject, per HSP, per unique subject (`qcovus` is for `blastn`). |

That is 50 field names. The GI fields print `0` where a sequence has no GI. The tax id fields come from the database (or print `N/A`); the name fields also need NCBI's `taxdb` files, otherwise they print `N/A` with a warning (see [BLAST databases](databases.md)). All the fields work with `psiblast`, `rpsblast` and `rpstblastn` too.

## XML, JSON and XML2 (5, 12 to 16)

- `5` is the older BLAST XML.
- `15` is BLAST JSON and `16` is BLAST XML2, one document for the whole run.
- `12` is the Seq-align in JSON (the alignments, not the report).
- `13` (JSON) and `14` (XML2) write one file per query and need `-out FILE`. LOSAT names the files as NCBI does: the report of query _n_ goes to `<base>_<n>.json` (13) or `<base>_<n>.xml` (14) next to `FILE`, where `<base>` is `FILE` without a `.json` (13) or `.xml` (14) extension. `FILE` itself is a master file that lists them. As in NCBI, the master file for format 13 names the side files by the `.xml` rule, so `-out o.json` lists `o.json_1.json` while the file written is `o_1.json`. Without `-out`, NCBI's error `Please provide a file name for outfmt 13` is shown.
- With `-db`, formats `13` to `16` write the scientific name of a sequence's tax id (`sciname`) when `taxdb.bti` is found as described in [BLAST databases](databases.md).

## ASN.1 and the archive (8, 9, 11)

- `8` and `9` write the alignments as an ASN.1 Seq-align, in text (`8`) and binary (`9`).
- `11` writes the BLAST archive, which holds the query, the search settings and the results. A message stored in it carries LOSAT's wording where LOSAT words that message as its own (the manual's name, the exception headers), as in the other formats' messages; its `ident` reads `LOSAT 0.2.0` where NCBI writes `2.17.0+`; the rest is NCBI's bytes. `rpsbproc` reads it for the CDD domain report. The environment variable `ARCHIVE_FORMAT` chooses the encoding as in NCBI: unset for ASN.1 text, a value starting with `bin` for ASN.1 binary, `xml` for XML.
- A title, id or database name with a byte that NCBI's stream replaces is written with `#` and NCBI's error message for each byte, as NCBI does.

## SAM (17)

`blastn` writes SAM with `-outfmt 17`. `-outfmt "17 SQ"` adds the target sequence in the SEQ column and `-outfmt "17 SR"` makes the subject the reference; other words are ignored, as in NCBI. With `-ungapped`, each subject has one line that joins all of its alignments with the query, as in NCBI: the CIGAR is the summed length of the alignments as one `M`, the position and clipping span all of them, `AS`, `EV` and `BS` are the best alignment's, and `PI` divides the best alignment's identities by the summed length. The document is NCBI's byte for byte except two things in the `@PG` line: the `CL:` tag holds LOSAT's own command line, and `VN:` reads `LOSAT 0.2.0` where NCBI writes `2.17.0+`. Like NCBI, LOSAT writes the document when the search ends.

## Organism report (18)

`-outfmt 18` is NCBI's "Tax BLAST report". It groups the subjects of each query by the tax id of each definition line, under the names of NCBI's `taxdb.bti` and `taxdb.btd`. A `-subject` record has tax id 0. Without a taxonomy database the report has its headers only. As in NCBI, a subject whose tax id the database does not know is left out of the report but still widens its columns.

## Differences from NCBI BLAST+

- `-outfmt 17`: the `CL:` tag of `@PG` holds LOSAT's command line.
- The version line of formats 0, 5, 7, 11, 13 to 17 names LOSAT, as described above.
- A few reports that NCBI cannot finish (it crashes or never ends) are refused after the part NCBI writes: see [What LOSAT refuses](refusals.md) under Output.

## Examples

<!-- example: output-formats-1 -->
```bash
losat blastn -query query.fna -subject genome.fna -outfmt "6 qseqid sseqid pident evalue"
```

<!-- example: output-formats-2 -->
```bash
losat blastn -query query.fna -subject genome.fna -outfmt 7
```

<!-- example: output-formats-3 -->
```bash
losat blastn -query query.fna -subject genome.fna -outfmt "10 delim=@ qacc sacc score"
```

<!-- example: output-formats-4 -->
```bash
losat blastn -query query.fna -subject genome.fna -outfmt 20
```

## See also

- [BLAST databases](databases.md)
- Each program's page lists the formats it writes.

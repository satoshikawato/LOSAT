# rpsbproc

`losat rpsbproc` turns the output of an RPS-BLAST search into a list of conserved domains. For each query it reports the domain hits, and optionally the functional sites and structural motifs of those domains, using NCBI's Conserved Domain Database (CDD) data files. It writes the same report as NCBI's `rpsbproc` 0.5.1.

## Synopsis

```bash
losat rpsbproc -i ARCHIVE.asn -d DATA_DIR [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `rpsbproc -i ARCHIVE.asn -d DATA_DIR [options]`.

## Input

You need two things.

**1. A BLAST archive.** `rpsbproc` reads the ASN.1 text archive that `rpsblast -outfmt 11` or `rpstblastn -outfmt 11` writes, with one archive per batch of queries. Give it with `-i FILE`, or pipe it to the standard input. It also reads the BLAST XML of `-outfmt 5`, as NCBI's `rpsbproc` does, but prints NCBI's warning that XML support is deprecated.

LOSAT's own `rpsblast` and `rpstblastn` write the archive with NCBI BLAST+ 2.17.0's bytes, so the two can be piped (one archive per batch of queries reaches the standard input as NCBI's do), and NCBI's `rpsblast` and `rpstblastn` work as well. For example:

```bash
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 11 | losat rpsbproc -d cdd_data -t doms -m rep -f
```

**2. The CDD data files.** `rpsbproc` needs six files from the CDD directory of NCBI's FTP site, ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/:

- `cddid.tbl`
- `cdtrack.txt`
- `family_superfamily_links`
- `cddannot.dat`
- `cddannot_generic.dat`
- `bitscore_specific.txt`

Put them in one folder and name it with `-d`. `-d` may also name a folder that has a `data` folder inside. Without `-d`, `rpsbproc` looks in `./data`, in the current folder, and in the folder of the program.

Use the data files of the same CDD release as the database you searched.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-i FILE` | standard input | The archive written by `rpsblast` or `rpstblastn` with `-outfmt 11` (LOSAT's or NCBI's) |
| `-o FILE` | standard output | Output file |
| `-d DIR` | see Input | Folder with the six data files |
| `-t doms\|feats\|both` | both | What to report: domain hits, features (sites and motifs), or both |
| `-m rep\|std\|full` | rep | How many domain hits to report: concise, standard or all |
| `-f` | off | Also list the superfamily of each domain hit |
| `-e E` | 0.01 | E-value shown in the report header; see below |
| `-q` | off | Quiet: no progress on the standard error (the program information is still shown) |

For every option: `losat rpsbproc -help`.

### Which hits you get

- `-t doms` prints the `DOMAINS` section. `-t feats` prints the `SITES` and `MOTIFS` sections. `-t both` prints all three.
- `-m rep` keeps the concise set of domain hits, `-m std` the standard set and `-m full` all hits. Specific and non-specific hits are classified as in NCBI's `rpsbproc`.
- `-f` adds the `SUPERFAMILIES` section, which lists the superfamilies of the domain hits.
- `-e` sets the E-value shown in the report header. It selects no hit of an archive, as in NCBI: the hits in an archive are those that passed the `-evalue` you gave to `rpsblast`. With BLAST XML input, `-e` drops the alignments above it.

## Output

`rpsbproc` writes a tab-delimited report: a header that explains the columns, then one `SESSION` block per search batch with a `QUERY` block for each query. Each query lists its `DOMAINS`, and its `SITES` and `MOTIFS` when you ask for them. The format is NCBI's. It is not one of the `-outfmt` formats; see [Output formats](../output-formats.md) for those.

For nucleotide queries (from `rpstblastn`), each reading frame is classified separately. Its hits are shown with the frame after the query name, such as `Query_1[-2]`, at positions on the plus strand. The `DOMAINS` section is printed even when there are no hits.

Messages and exit codes are NCBI's: 2 for a problem in the data files, 3 for the input file, 4 for the output file and 5 for the input data.

## Differences from NCBI BLAST+

- Where NCBI's `rpsbproc` reads memory it should not, or its result depends on the machine, LOSAT refuses the input or computes the value without the undefined step ([the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md)).
- Hits of one query with equal E-value and bit score keep their input order, where NCBI's order can differ, and so can which of them is `Specific` ([Differences from NCBI BLAST+](../differences.md#rpsbproc-hits-that-tie)).

## What LOSAT refuses

`rpsbproc` refuses inputs on which NCBI's program reads outside its buffers or crashes. Examples are motif checks that depend on bytes past a site, sites outside the query, and archives with parts LOSAT does not read, such as Seq-align segments other than `Dense-seg`, `Std-seg` and `Packed-seg`. It also refuses file names that are not UTF-8 and a `[datafiles]` entry in a registry file. For the full list, see [What LOSAT refuses](../refusals.md).

## Examples

Print the concise domain hits of an archive, with superfamilies:

<!-- example: rpsbproc-1 -->
```bash
losat rpsbproc -i hits.asn -d cdd_data -t doms -m rep -f
```

Read the archive from the standard input and report domains, sites and motifs with the standard set of hits:

<!-- example: rpsbproc-2 -->
```bash
cat hits.asn | losat rpsbproc -d cdd_data -t both -m std
```

## See also

- [rpsblast](rpsblast.md) and [rpstblastn](rpstblastn.md), which search CDD
- [Output formats](../output-formats.md)
- [Coming from NCBI BLAST+](../from-ncbi-blast.md)

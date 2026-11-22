# BLASTX command line

```bash
losat blastx -query query.fna -subject subject.faa -num_threads 1 -outfmt 0
losat blastx -query query.fna -subject subject.faa -outfmt "7 std qframe btop stitle" -out results.txt
```

The FASTA path (`-subject`) writes pairwise format 0,
the query-anchored formats 1 to 4 (with and without identities, and flat), tabular
format 6, commented tabular format 7, CSV formats 10 and 20 (20 with header lines),
BLAST XML format 5, Seq-align formats 8 (text ASN.1), 9 (binary ASN.1) and 12
(JSON), the BLAST archive format 11, BLAST JSON format 15, BLAST XML2 format 16 and
the organism report format 18; `-db` writes the same formats. Format 0 is the
default; format 17 (SAM) is BLASTN's alone.
Formats 15 and 16 take the structure of 5 (one document across the query
batches) and, with `-db`, the taxid and scientific name from `taxdb.bti` when
the taxonomy database is on the search path. Formats 13 and 14 (multiple-file
BLAST JSON and XML2) write the same documents as side files next to `-out FILE`
with a master file at `FILE` (they need `-out`); like 5, 15 and 16 they refuse a
query that is not a regular file, and a run that is refused after its first
batch leaves no side file.
Supported query genetic-code IDs are
`1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33`.
ID 32 is invalid for this BLASTX CLI. Protein subjects have no `db_gencode`.

The search profile supports NCBI's matrices (BLOSUM45, BLOSUM50, BLOSUM62,
BLOSUM80, BLOSUM90, PAM30, PAM70, PAM250) with the gap costs of NCBI's table for
each and their suggested threshold and window size, IDENTITY in an ungapped
search with word sizes 2 to 4, word sizes 2 to 7 (5 to 7 with the compressed
lookup table, as NCBI), composition modes 0 to 3 (`2u` as 2, `5` as 0, as NCBI
reads them), plus/minus/both query strands, SEG and lowercase masks,
soft/hard masking, gapped search or ungapped search with composition mode 0,
sum statistics, intron length, target/HSP limits, culling, and subject-besthit.
Omitted defaults retain their NCBI meanings, including description/alignment
limits. Multiple query and subject records retain input order, with ordinary
NCBI BLASTX batching and split/rejoin boundaries.

Custom formats 6, 7, 10 and 20 take a field list and `delim=` (README, "Output formats"); `std`
expands to the 12 standard columns. Output goes to stdout unless `-out` names a file; `-out -` also
uses stdout.

Each protein subject may contain at most 5,000,000 residues; LOSAT refuses a longer one (NCBI
BLAST+ searches it in overlapping parts). Unsupported options and tasks fail explicitly (README,
"What LOSAT refuses"). `-query -` (or no `-query`) and `-subject -` read standard input as NCBI BLAST+
reads it. Native threads and the serial and threaded command-WASI and Web/reactor entry points are
implemented. BLASTX is part of v0.2.0; the results on the release candidate are in the
[release record](../release/v0.2.0.md). Ordinary native paths such as `/dev/stdin` retain the source
stream behavior, including pipe EOF, warnings, and format-0 header flush before query data.

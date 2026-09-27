# BLASTX native serial CLI on the current branch

```bash
LOSAT blastx -query query.fna -subject subject.faa -num_threads 1 -outfmt 0
LOSAT blastx -query query.fna -subject subject.faa -outfmt "7 std qframe btop stitle" -out results.txt
```

The local FASTA path uses the final Rust runtime results for pairwise format 0,
tabular format 6, and commented tabular format 7. Format 0 is the default.
Supported query genetic-code IDs are
`1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33`.
ID 32 is invalid for this BLASTX CLI. Protein subjects have no `db_gencode`.

The search profile supports BLOSUM62, gap costs 11/1, word sizes 3/5,
composition modes 0/2, plus/minus/both query strands, SEG and lowercase masks,
soft/hard masking, gapped search or ungapped search with composition mode 0,
sum statistics, intron length, target/HSP limits, culling, and subject-besthit.
Omitted defaults retain their NCBI meanings, including description/alignment
limits. Multiple query and subject records retain input order, with ordinary
NCBI BLASTX batching and split/rejoin boundaries.

Custom formats 6/7 accept these 30 fields and `std` in the requested order,
including duplicate fields:

```
qseqid qacc qaccver qlen sseqid sacc saccver slen
qstart qend sstart send qseq sseq evalue bitscore score length
pident nident mismatch positive gapopen gaps ppos
qframe sframe frames btop stitle
```

`std` expands to the standard 12 columns. Unknown fields, an empty entire
format specification, or fields attached to format 0 are refused under the
registered CLI contract. A blank field list following `6 ` or `7 ` uses the
registered standard-field expansion. Output goes to stdout unless `-out` names
a file; `-out -` also uses stdout.

Each protein subject may contain at most 5,000,000 residues. Longer subjects,
database searches, stdin FASTA specified with `-query -` or `-subject -`,
unsupported options/tasks/matrices/gap costs,
and composition modes 1/3 fail explicitly. Native thread counts above one are
unimplemented until Session F. BLASTX WASI/Web, parallel execution, pool-specific
failures, and formal performance measurements remain for F/G. This branch has
not been published as a BLASTX release. Ordinary native paths such as
`/dev/stdin` retain the source stream behavior, including pipe EOF, warnings,
and format-0 header flush before query data.

Session E acceptance requires the public CLI's raw reports, warnings, exits,
file lifecycle, M1–M7 and applicable native M10 comparisons, current B/C/D and
shared-program regressions, quality checks, and independent audit. The evidence
and final acceptance state live in
`docs/evidence/losatx_stage_e_reporting/run-20260926T143040Z-02aa2d/`.

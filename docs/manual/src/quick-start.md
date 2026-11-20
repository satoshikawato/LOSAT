# Quick start

This page takes about five minutes. You run a nucleotide search and a protein search on two FASTA files, build a database, search it, and read the result. If you use NCBI BLAST+, everything here is familiar: LOSAT takes the same options.

You need `losat` on your `PATH` ([Installation](installation.md)) and a few FASTA files. The examples use these names, which stand for your own files:

| File | What it is |
|:---|:---|
| `query.fna`, `subject.fna` | Nucleotide FASTA files |
| `query.faa`, `subject.faa` | Protein FASTA files |

Every program has the form `losat <program> [options]`. `losat <program> -help` lists its options.

## 1. Compare two nucleotide files

`blastn` compares `-query` with `-subject`. There is no database to build first.

<!-- example: quick-start-1 -->
```bash
losat blastn -query query.fna -subject subject.fna -outfmt 6
```

`-outfmt 6` writes one tab-separated line for each alignment. The default task is `megablast`, which finds highly similar sequences quickly. Add `-task blastn` for more distant matches.

## 2. Compare two protein files

`blastp` works the same way. Here `-evalue 1e-5` keeps only the alignments with an E-value of 1e-5 or lower.

<!-- example: quick-start-2 -->
```bash
losat blastp -query query.faa -subject subject.faa -evalue 1e-5 -outfmt 6
```

The other search programs, `blastx`, `tblastn` and `tblastx`, take the same `-query` and `-subject` options; their pages are in the Programs section.

## 3. Build a database and search it

Use a database when you search the same subject many times. `makeblastdb` writes the same files as NCBI's `makeblastdb`.

<!-- example: quick-start-3 -->
```bash
losat makeblastdb -in subject.faa -dbtype prot -out subject_db -logfile makeblastdb.log
losat blastp -query query.faa -db subject_db -evalue 1e-5 -outfmt 6
```

`-dbtype` is `prot` or `nucl`. `-out` names the database; the search then names it with `-db`. Without `-logfile`, `makeblastdb` prints a short report of what it built. The search programs read databases made by LOSAT or by NCBI BLAST+, including NCBI's distributed databases. See [BLAST databases](databases.md).

## 4. Read tabular output

With `-outfmt 6`, each line is one alignment (one HSP) and has 12 columns, in this order:

| # | Column | Meaning |
|--:|:---|:---|
| 1 | `qaccver` | Query sequence id |
| 2 | `saccver` | Subject sequence id |
| 3 | `pident` | Percent identity |
| 4 | `length` | Alignment length |
| 5 | `mismatch` | Number of mismatches |
| 6 | `gapopen` | Number of gap openings |
| 7 | `qstart` | Start of the alignment in the query |
| 8 | `qend` | End of the alignment in the query |
| 9 | `sstart` | Start of the alignment in the subject |
| 10 | `send` | End of the alignment in the subject |
| 11 | `evalue` | Expect value |
| 12 | `bitscore` | Bit score |

For nucleotide searches, `sstart` is larger than `send` when the match is on the minus strand of the subject.

`-outfmt 7` writes the same lines with comment lines that name the program, the query, the database and the fields, which is easier to read by eye.

<!-- example: quick-start-4 -->
```bash
losat blastn -query query.fna -subject subject.fna -outfmt 7
```

To choose your own columns, put the field names in quotes after the number:

<!-- example: quick-start-5 -->
```bash
losat blastn -query query.fna -subject subject.fna -outfmt "6 qseqid sseqid pident length evalue bitscore"
```

The other formats (pairwise text, XML, JSON, SAM and more) are in [Output formats](output-formats.md).

## 5. Read the query from a pipe

Every search program reads its query from the standard input when you give `-query -`. You can also decompress on the way: a file name that ends in `.gz`, `.bz2` or `.zst` is decompressed for you.

<!-- example: quick-start-6 -->
```bash
cat query.fna | losat blastn -query - -subject subject.fna -outfmt 6
```

## 6. Use more threads

`-num_threads` sets the number of threads. The default is 1. The output is the same at every thread count. See [Threads, speed and memory](performance.md).

```bash
losat blastp -query query.faa -db subject_db -num_threads 4 -outfmt 6
```

## Where next

- To call the programs by their NCBI names, see [Coming from NCBI BLAST+](from-ncbi-blast.md).
- For every option of a program, read its page in the Programs section, or run `losat <program> -help`.
- If a command stops with an error, see [Troubleshooting](troubleshooting.md).

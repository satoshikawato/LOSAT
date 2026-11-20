# Compressed `-in` files: what makeblastdb 2.17.0 reads (probe of 2026-11-04)

Line GZ (`losat-baselines/gzinput-20261103/`), NCBI BLAST+ 2.17.0 (the local oracle, `$NCBI_BIN`),
two one-record FASTA files `a.fa` (37 bases) and `b.fa` (32 bases). Only the lines that matter are
quoted; the runs are reproducible with the commands shown.

```text
$ gzip -kc a.fa > a.fa.gz; gzip -kc b.fa > b.fa.gz; cat a.fa.gz b.fa.gz > ab.fa.gz
$ makeblastdb -in ab.fa.gz -dbtype nucl -out db_gz     # exit 0
Adding sequences from FASTA; added 2 sequences in 0.000401974 seconds.
$ blastdbcmd -db db_gz -info | grep sequences
        2 sequences; 69 total bases

$ bzip2 -kc a.fa > a.fa.bz2; bzip2 -kc b.fa > b.fa.bz2; cat a.fa.bz2 b.fa.bz2 > ab.fa.bz2
$ makeblastdb -in ab.fa.bz2 -dbtype nucl -out db_bz2   # exit 0, no warning
Adding sequences from FASTA; added 1 sequences in 0.00503087 seconds.
$ blastdbcmd -db db_bz2 -info | grep sequences
        1 sequences; 37 total bases

$ zstd -q a.fa -o a.fa.zst; zstd -q b.fa -o b.fa.zst; cat a.fa.zst b.fa.zst > ab.fa.zst
$ makeblastdb -in ab.fa.zst -dbtype nucl -out db_zst   # exit 0, no warning
Adding sequences from FASTA; added 1 sequences in 0.000159025 seconds.
$ blastdbcmd -db db_zst -info | grep sequences
        1 sequences; 37 total bases
```

The second bzip2 stream and the second zstd frame are lost silently (`PD-LOSAT-NCBI-DEFECTS`,
makeblastdb exception 7). `losat makeblastdb` reads every stream and frame: 2 sequences for all
three files, the files of `db_parts.*` equal to those of the plain input
(`LOSAT/tests/cli_smoke.rs`, `makeblastdb_reads_compressed_inputs_as_makeblastdb_does`).

Text under a compressed name (`cp a.fa text.fa.bz2`; the same for `.zst` and `.gz`): makeblastdb
warns twice (`CBZip2Decompressor::Process ... 'Magic' leading bytes are missed`,
`CZstdDecompressor::Process ... Unknown frame descriptor`, `CZipDecompressor::Process ... data
error`), says the file does NOT appear to be FASTA (processing anyway), and ends with
`BLAST failed to write output`, exit 6. The NCBI search programs read `a.fa.gz`, `a.fa.bz2` and
`a.fa.zst` as text: `BLAST query error: CFastaReader: Near line 1, there's a line that doesn't look
like plausible data, but it's not marked as defline or comment.`, exit 1.

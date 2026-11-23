# Common workflows

Recipes for everyday tasks: searching a database you build, annotating proteins and contigs, comparing genomes, iterative searches, conserved domains, masks, and running existing tools on LOSAT. Each command uses the options of NCBI BLAST+ 2.17.0; replace the file names with your own. The program pages give every option.

## Search reads against a reference set (`blastn`)

Build a nucleotide database from your reference sequences once, for example a set of 16S rRNA genes, then search your reads or amplicons against it. `blastn` uses its default task, `megablast`, for close matches:

```bash
losat makeblastdb -in 16s_all.fna -dbtype nucl -out 16s_db
losat blastn -query reads.fna -db 16s_db -evalue 1e-5 -max_target_seqs 5 -outfmt 6 -num_threads 16
```

The output is the same at every `-num_threads`, so use as many threads as you have cores.

Compressed reads work directly. A `-query` or `-subject` file whose name ends in `.gz`, `.bz2` or `.zst` is read through its decompressor:

```bash
losat blastn -query reads.fna.gz -db 16s_db -outfmt 6 -num_threads 16
```

LOSAT writes a temporary uncompressed copy of such a file, beside the input or in `TMPDIR` (see [Environment variables](environment.md#temporary-files)). To avoid the copy, decompress into a pipe and read the query from the standard input:

```bash
gzip -dc reads.fna.gz | losat blastn -query - -db 16s_db -outfmt 6 -num_threads 16
```

The WASI command (`losat.wasm`) has no temporary directory: it writes the copy beside the input, or stops with an error if it cannot. Use the pipe there when the input's folder is read-only.

## Annotate proteins against Swiss-Prot (`blastp`, `blastx`)

Download NCBI's `swissprot` database as shown in [BLAST databases](databases.md#get-an-ncbi-database). Then search your proteins against it, with the columns you want in the table:

```bash
losat blastp -query proteins.faa -db swissprot \
  -evalue 1e-5 -max_target_seqs 5 -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"
```

To search nucleotide contigs or transcripts against the same protein database without predicting genes first, use `blastx`. It translates each query in six frames:

```bash
losat blastx -query contigs.fna -db swissprot -evalue 1e-5 -outfmt 6 -num_threads 16
```

## Compare genomes with another genetic code (`tblastx`)

`tblastx` compares two nucleotide sequences through their translations in six frames. For organisms that use another genetic code, give it for both sides, for example code 4 for *Mycoplasma*:

```bash
losat tblastx -query query_genomes.fna -subject subject_genomes.fna \
  -query_gencode 4 -db_gencode 4 -outfmt 6 -num_threads 8
```

With `-subject`, LOSAT translates the subject with the code of `-db_gencode`, as you asked. NCBI BLAST+ translates a `-subject` file with code 1 whatever `-db_gencode` says, so the results of this command differ from NCBI's; see [Differences from NCBI BLAST+](differences.md#a-local--subject-search-of-tblastx-and-tblastn-uses--db_gencode). `tblastn` with `-subject` does the same.

## Iterative homolog search (`psiblast`)

Run up to three rounds of PSI-BLAST against a protein database. Each round searches with a profile (PSSM) built from the hits of the round before; the search stops earlier if no new sequence is found. Keep the tabular hits and the final PSSM for later searches:

```bash
losat psiblast -query seed.faa -db swissprot \
  -num_iterations 3 -evalue 0.001 -outfmt 6 -num_threads 8 \
  -save_pssm_after_last_round -out_pssm seed.asn -out_ascii_pssm seed_pssm.txt
```

Search a nucleotide database with the saved PSSM through `tblastn -in_pssm`, which translates the database sequences in six frames:

```bash
losat makeblastdb -in genome.fna -dbtype nucl -out genome_db
losat tblastn -in_pssm seed.asn -db genome_db -evalue 1e-5 -outfmt 6 -num_threads 8
```

## Conserved domains (`rpsblast`, `rpstblastn`, `rpsbproc`)

Search proteins (`rpsblast`), or nucleotide sequences translated in six frames (`rpstblastn`), against NCBI's Conserved Domain Database (CDD). [rpsblast](programs/rpsblast.md) says where to download the `Cdd` database:

```bash
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 6 -num_threads 8
losat rpstblastn -query contigs.fna -db Cdd -evalue 0.01 -outfmt 6 -num_threads 8
```

Both write the same output at every `-num_threads`: NCBI BLAST+'s output with `-num_threads 1` (apart from the version line, which names LOSAT).

For the domain, site and motif annotation of `losat rpsbproc`, you need a BLAST archive (`-outfmt 11`), which LOSAT's `rpsblast` and `rpstblastn` write as NCBI BLAST+'s do. Annotate it with LOSAT. `cdd_data` is a folder with the six CDD data files listed on the [rpsbproc](programs/rpsbproc.md) page:

```bash
# Step 1: the archive
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 11 -out hits.asn

# Step 2: domains, sites and motifs, with LOSAT
losat rpsbproc -i hits.asn -d cdd_data -t both -m rep -f -o domains.txt
```

## Mask low-complexity regions in a database

Find the low-complexity regions of a genome with `dustmasker`, store them in the database with `makeblastdb -mask_data`, then search with them as soft masks (`-db_soft_mask 11`, the id of DUST masks; `losat blastdbcmd -db genome_db -info` lists the masks a database holds):

```bash
losat dustmasker -in genome.fna -parse_seqids -outfmt maskinfo_asn1_bin -out genome_dust.asnb
losat makeblastdb -in genome.fna -dbtype nucl -mask_data genome_dust.asnb -parse_seqids -out genome_db
losat blastn -query query.fna -db genome_db -db_soft_mask 11 -outfmt 6
```

Words inside a soft mask do not start alignments, but alignments can extend through masked letters. See [BLAST databases](databases.md#masks).

## Run existing tools on LOSAT (`--install-shims`)

Many tools and workflow systems run `blastn`, `blastp`, `makeblastdb` and the other programs by name from your `PATH`. Make the NCBI program names run LOSAT, and put their folder first on your `PATH`:

```bash
losat --install-shims ~/.local/share/losat-shims
export PATH="$HOME/.local/share/losat-shims:$PATH"
blastp -version
```

`blastp -version` now prints `blastp: LOSAT 0.2.0` on its first line and `Package: LOSAT 0.2.0, output compatible with blast 2.17.0` on its second line. The tools need no changes. LOSAT's benchmark ran the BLAST command lines of published tools on their own test data, among them ABRicate, AMRFinderPlus, BUSCO, Kaptive, mlst, OrthoFinder, Prokka, Roary and Trinotate; each output was the same as NCBI BLAST+ 2.17.0's, apart from the line that names the program's version in the formats that have one (for BUSCO, whose NCBI BLAST+ output changes with the thread count, the same as its 1-thread output) ([Threads, speed and memory](performance.md)).

[Coming from NCBI BLAST+](from-ncbi-blast.md) explains the names, how they are made on each platform, and what the command does when a file of the same name is already there.

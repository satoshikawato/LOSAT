# Common workflows

This page collects practical recipes for everyday searches, databases, iterative PSI-BLAST runs, conserved domain annotation, and drop-in pipeline replacement.

## 16S rRNA or amplicon classification (`blastn`)

Build a nucleotide database once, then classify queries using MegaBLAST (default task of `blastn`) across multiple threads:

```bash
losat makeblastdb -in 16s_all.fna -dbtype nucl -out 16s_db
losat blastn -query reads.fna -db 16s_db -evalue 1e-5 -max_target_seqs 5 -outfmt 6 -num_threads 16
```

To stream gzipped reads straight from decompression without creating an extra uncompressed copy in your shell pipeline:

```bash
gzip -dc reads.fna.gz | losat blastn -query - -db 16s_db -outfmt 6 -num_threads 16
```

Native `losat` also reads `.gz`, `.bz2`, and `.zst` inputs directly when given a file path (for instance `-query reads.fna.gz`). On WASI (`losat.wasm`), decompress first or pipe via `-query -`.

## Proteome annotation against Swiss-Prot (`blastp` / `blastx`)

Search a FASTA file of translated CDSs against UniProtKB/Swiss-Prot with custom tabular columns and an organism summary:

```bash
losat blastp -query proteins.faa -db swissprot \
  -evalue 1e-5 -max_target_seqs 5 -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"
```

To search nucleotide contigs directly against a protein database without predicting ORFs first, use `blastx`:

```bash
losat blastx -query contigs.fna -db swissprot -evalue 1e-5 -outfmt 6 -num_threads 16
```

## Pairwise genome comparisons (`tblastx`, `blastn`)

Compare nucleotide assemblies using six-frame translation with non-standard genetic codes (for example Mycoplasma / code 4 for both query and local FASTA subject):

```bash
losat tblastx -query query_genomes.fna -subject subject_genomes.fna \
  -query_gencode 4 -db_gencode 4 -outfmt 6 -num_threads 8
```

Note that LOSAT applies `-db_gencode` to local `-subject` files in `tblastx` and `tblastn`; see [Differences from NCBI BLAST+](differences.md#a-local--subject-search-of-tblastx-and-tblastn-uses--db_gencode).

## Iterative homolog search with `psiblast`

Run three iterations of PSI-BLAST, keeping both the tabular hits and the final PSSM checkpoint for later searches:

```bash
losat psiblast -query seed.faa -db swissprot \
  -num_iterations 3 -evalue 0.001 -outfmt 6 -num_threads 8 \
  -save_pssm_after_last_round -out_pssm seed.asn -out_ascii_pssm seed_pssm.txt
```

Reuse the saved checkpoint against a translated nucleotide genome using `tblastn -in_pssm`:

```bash
losat tblastn -in_pssm seed.asn -db genome_db -evalue 1e-5 -outfmt 6 -num_threads 8
```

## Conserved domain search and annotation (`rpsblast`, `rpstblastn`, `rpsbproc`)

Search protein sequences (`rpsblast`) or six-frame translated DNA (`rpstblastn`) against NCBI's Conserved Domain Database (CDD):

```bash
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 6 -num_threads 8
losat rpstblastn -query contigs.fna -db Cdd -evalue 0.01 -outfmt 6 -num_threads 8
```

`losat rpsblast` and `losat rpstblastn` produce the exact output of NCBI BLAST+ `-num_threads 1` regardless of `-num_threads`.

In v0.2.0, `losat rpsblast` and `losat rpstblastn` write formats `0`, `6`, `7`, and `10` (and `1`–`4` in `rpstblastn`), not `-outfmt 11`. To run domain, site, and motif annotation with `losat rpsbproc`, generate the ASN.1 archive (`-outfmt 11`) with NCBI's `rpsblast` or `rpstblastn`, then run `losat rpsbproc` using the six CDD annotation data files (`cddid.tbl`, `cdtrack.txt`, `family_superfamily_links`, `cddannot.dat`, `cddannot_generic.dat`, `bitscore_specific.txt`):

```bash
# Step 1: produce archive with NCBI rpsblast (-outfmt 11)
ncbi-blast/rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 11 -out hits.asn

# Step 2: annotate domains and features with LOSAT
losat rpsbproc -i hits.asn -d cdd_data -t both -m rep -f -o domains.txt
```

## Masking low-complexity regions before building a database

Mask low-complexity repeats into interval masks, then embed them when creating a v5 database so searches can apply `-db_soft_mask`:

```bash
losat dustmasker -in genome.fna -parse_seqids -outfmt maskinfo_asn1_bin -out genome_dust.asnb
losat makeblastdb -in genome.fna -dbtype nucl -mask_data genome_dust.asnb -parse_seqids -out genome_db
losat blastn -query primers.fna -db genome_db -db_soft_mask 11 -outfmt 6
```

## Drop-in replacement in existing tools (`--install-shims`)

Workflow managers and tools such as Prokka, AMRFinderPlus, ABRicate, Roary, OrthoFinder, Kaptive, MLST, Trinotate, BUSCO, CheckV, and Galaxy call `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `makeblastdb`, and `blastdbcmd` by name on `PATH`. Point them at LOSAT without editing any tool scripts:

```bash
losat --install-shims ~/.local/share/losat-shims
export PATH="$HOME/.local/share/losat-shims:$PATH"
blastp -version   # confirms "Package: LOSAT 0.2.0, output compatible with blast 2.17.0"
```

See [Coming from NCBI BLAST+](from-ncbi-blast.md) for symlink details and safety checks.

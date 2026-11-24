# deltablast

`losat deltablast` finds distant relatives of a protein. It first looks for conserved domains in your query by searching a domain database, builds a position-specific scoring matrix (PSSM) from the domains it finds, and searches a protein database with that matrix. Further rounds work as in [psiblast](psiblast.md).

## Synopsis

```bash
losat deltablast -query QUERY.faa -db DB [options]
losat deltablast -query QUERY.faa -subject SUBJECT.faa [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `deltablast -query QUERY.faa -db DB [options]`.

## Input

- **Query.** A protein FASTA file with `-query` (`-query -` reads the standard input). It can hold several queries.
- **Database or subject.** `-db` names a protein BLAST database (see [BLAST databases](../databases.md)). `-subject` searches a protein FASTA file instead, and `-subject_loc START-STOP` searches part of the subject.
- **Domain database.** `-rpsdb NAME` names the database of conserved domains (default `cdd_delta`). LOSAT finds it as it finds `-db`: in the current folder, then in the folders of `BLASTDB`. It must be a DELTA-BLAST domain database: NCBI's `cdd_delta`, which you can download from NCBI's FTP site (`blast/db/cdd_delta.tar.gz`), or one made by NCBI's `makeprofiledb -dbtype delta`. Besides the files of an RPS database it needs the `.freq`, `.wcounts` and `.obsr` files; where one is missing you get NCBI's error. See [BLAST databases](../databases.md).
- **Database filters.** `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-no_taxid_expansion`, `-seqidlist` and `-negative_seqidlist` keep every round to the sequences they select.

## Rounds and PSSMs

Round 1 searches the PSSM of the conserved domains that the domain database finds in your query. A domain counts when its E-value is below `-domain_inclusion_ethresh` (default 0.05). A query without domain hits starts from the scoring matrix's scores.

`-num_iterations N` sets the largest number of rounds (default 1). `0` runs until the search converges. After round 1, the rounds are those of `psiblast`: LOSAT builds a PSSM from the alignments with an E-value below `-inclusion_ethresh` (default 0.002) and uses it in the next round.

The PSSM files are written as in `psiblast`; the PSSM of round 1 is the PSSM of the domains:

- `-out_pssm FILE` and `-out_ascii_pssm FILE` get the PSSM at the start of each later round.
- `-save_pssm_after_last_round` also saves the PSSM built from the last round's hits.
- `-save_each_pssm` saves the PSSM of every round as `FILE.1`, `FILE.2` and so on.

`-comp_based_stats` takes 0 to 3 (default 1; `psiblast`'s default is 2) and NCBI's letters (`D`, `d`, `T` and `t` mean 1, `F` and `f` mean 0). Later rounds run mode 2 and mode 3 as mode 1, and LOSAT prints NCBI's warning, as `psiblast` does.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-query FILE` | `-` | Protein FASTA query |
| `-db NAME` | | Protein database |
| `-subject FILE` | | Protein FASTA file to search instead of a database |
| `-rpsdb NAME` | `cdd_delta` | Domain database |
| `-num_iterations N` | 1 | Number of rounds; 0 means until convergence |
| `-domain_inclusion_ethresh E` | 0.05 | E-value a domain must beat to enter the PSSM of round 1 |
| `-inclusion_ethresh E` | 0.002 | E-value a hit must beat to enter the PSSM of a later round |
| `-pseudocount N` | 0 | Pseudocount used when building the PSSM |
| `-evalue E` | 10 | Report hits with an E-value up to this |
| `-out_pssm FILE` | | Write the PSSM (ASN.1 text) |
| `-out_ascii_pssm FILE` | | Write the PSSM as a readable matrix |
| `-save_pssm_after_last_round` | off | Also save the PSSM built after the last round |
| `-save_each_pssm` | off | Save the PSSM of every round |
| `-comp_based_stats N` | 1 | Composition-based statistics, 0 to 3 |
| `-matrix NAME` | BLOSUM62 | Scoring matrix |
| `-gapopen N`, `-gapextend N` | 11, 1 | Gap costs (the matrix's own values for other matrices) |
| `-max_target_seqs N` | 500 | Number of subject sequences to keep |
| `-seg yes\|no\|"W L H"` | no | Filter low-complexity regions of the query |
| `-lcase_masking`, `-soft_masking` | off | Use lower-case letters of the query and subjects as masks; apply filters as soft masks |
| `-num_threads N` | 1 | Number of threads; the output is the same at every value |
| `-outfmt FORMAT` | 0 | Output format |
| `-out FILE` | standard output | Output file |

The other options of NCBI's `deltablast` that `psiblast` also has work as they do in LOSAT's `psiblast` (`-word_size`, `-threshold`, `-gap_trigger`, the X-drop values, `-qcov_hsp_perc`, `-max_hsps` and so on). For every option: `losat deltablast -help`.

## Output

`losat deltablast` writes `-outfmt 0`, `6` and `7`. Formats 6 and 7 take a custom field list, for example `-outfmt "6 qseqid sseqid evalue"`. With `-outfmt 0`, the report shows each round. With `-outfmt 7`, the comment lines count rounds. See [Output formats](../output-formats.md).

The output is the same as NCBI BLAST+ 2.17.0's `deltablast -num_threads 1` (standard output, error messages, exit status and PSSM files), at every thread count. The program line says `DELTABLAST (LOSAT 0.2.0)`, and error messages start with `LOSAT ...`, as in the other programs.

## Differences from NCBI BLAST+

- `deltablast` is not in LOSAT Web (the browser refuses it, as it refuses `psiblast`), and not in the library API or `pylosat`.

## What LOSAT refuses

Not yet supported:

- `-show_domain_hits` and `-query_loc`.
- `-outfmt` other than 0, 6 and 7. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors.
- `-html` and `-export_search_strategy`.

Also refused:

- `-remote`, `-entrez_query`, `-import_search_strategy`, `-gilist` and `-negative_gilist`.
- `-use_sw_tback`.
- A domain database of several volumes.
- A `-subject` file with more records than `-max_target_seqs`.
- Some runs of several queries, where NCBI crashes or runs only the first of them through the later rounds.

See [What LOSAT refuses](../refusals.md#deltablast) for the full list.

## Examples

Search a database for three rounds and print a table:

```bash
losat deltablast -query query.faa -db proteins -num_iterations 3 -outfmt 6
```

Use a domain database in another folder and keep the PSSM of the domains:

```bash
losat deltablast -query query.faa -db proteins -rpsdb /data/cdd_delta -out_ascii_pssm query.pssm.txt
```

## See also

- [psiblast](psiblast.md), whose later rounds `deltablast` runs
- [rpsblast](rpsblast.md), the search of the domain database
- [makeblastdb](makeblastdb.md), to build a protein database
- [Output formats](../output-formats.md)

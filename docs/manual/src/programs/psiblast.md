# psiblast

`losat psiblast` finds distant relatives of a protein. It searches a protein database with your query, builds a position-specific scoring matrix (PSSM) from the hits it trusts, and searches again with that matrix. You can repeat this for several rounds, or until no new sequences appear.

## Synopsis

```bash
losat psiblast -query QUERY.faa -db DB [options]
losat psiblast -query QUERY.faa -subject SUBJECT.faa [options]
losat psiblast -in_pssm QUERY.pssm -db DB [options]
losat psiblast -in_msa ALIGNMENT.aln -db DB [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `psiblast -query QUERY.faa -db DB [options]`.

## Input

- **Query.** A protein FASTA file with `-query` (`-query -` reads the standard input), a saved PSSM with `-in_pssm`, or an alignment with `-in_msa`, in place of `-query`. `-in_pssm -` reads the PSSM from the standard input, and `-in_msa -` the alignment. `-query_loc START-STOP` searches part of the query.
- **Several queries.** Each query runs its own rounds, one query after another.
- **Database or subject.** `-db` names a protein BLAST database, with or without parsed sequence IDs (see [BLAST databases](../databases.md)). `-db "a b"` searches several databases at once. `-subject` searches a protein FASTA file instead. `-subject_loc START-STOP` searches part of the subject.
- **Database filters.** `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-no_taxid_expansion`, `-seqidlist` and `-negative_seqidlist` keep every round to the sequences they select.
- **Alignment.** `-in_msa FILE` starts from a multiple sequence alignment in Clustal, FASTA (with gaps), PHYLIP, NEXUS, Sequin or MultAlin format; LOSAT finds the format as NCBI BLAST+ does, from the first lines of the file. LOSAT builds the PSSM from the alignment as NCBI BLAST+ does, and searches it as it searches an `-in_pssm` PSSM, with the same rounds and PSSM files. The first sequence is the query. `-msa_master_idx N` takes the N-th sequence instead (1 is the first), and `-ignore_msa_master` builds the PSSM without counting the query sequence. A lower-case letter of the query sequence is left out of the alignment. The PSSM from the alignment is built with a pseudocount of 0; `-pseudocount` applies to the PSSMs of later rounds. `-matrix` and the gap costs name the scoring of the PSSM. The text after the sequence name in a FASTA, PHYLIP or NEXUS file is the sequence's title: it appears in the `# Query:` line, in the PSSM file and in the messages, without its `[name=value]` modifiers. The modifiers are read as NCBI reads them (a warning for one it does not know or that is discouraged, an error for a value it rejects), and what they set (organism, molecule type, comments, publications, protein names, strand and topology) is written with the query in the PSSM file and the search strategy as NCBI writes it. With a `molecule` or `mol-type` modifier that sets the query's molecule to `other`, the search is as for a protein, `-outfmt 0` writes the whole report with a line in place of each alignment (`Sequence with id <id> no longer exists in database...alignment skipped`, one for each HSP, no heading and no score lines, the id being the subject's: its own id, `BL_ORD_ID:<oid>` for a database sequence without ids, `Subject_N` for the N-th `-subject` record), as NCBI does, and a tabular report (`-outfmt 6` or `7`) whose fields read the query's letters (the default fields do) ends at its first row with NCBI's `CSeqVectorException`, exit status 255, after the rows and comment lines before it, as NCBI does. A master sequence that is the start of a longer sequence above it in the file is refused: put the master above the longer one (see [Refusals](../refusals.md)). The library and pylosat do not take `-in_msa`: make the PSSM with `-out_pssm` and read it there.
- **PSSM file.** `-in_pssm` reads a PSSM in NCBI's ASN.1 text form, the form `-out_pssm` writes, or in its ASN.1 binary form. It is searched with the BLOSUM62 matrix and gap costs 11/1.
  - The reports show the PSSM's query as the file holds it: lower-case letters stay lower case and match no subject letter, and a gap residue (`-`) is shown as `-`. NCBI BLAST+ shows the same.
  - A query residue that its coding does not have (such as `*` in `iupacaa`) stops the run with `LOSAT engine error: Warning: Invalid residues found at positions ...` and exit status 3, as in NCBI BLAST+ (which starts the message with `BLAST engine error:`). [tblastn](tblastn.md) `-in_pssm` does the same.

## Rounds and PSSMs

Round 1 is a protein search of your query, as in `blastp`. After each round, LOSAT builds a PSSM from the alignments with an E-value below `-inclusion_ethresh` and uses it in the next round.

`-num_iterations N` sets the largest number of rounds (default 1). `0` runs until the search converges. A run ends when a round finds no hits, when no sequence has newly entered the model (the report says `Search has CONVERGED!`), or when N rounds are done.

The PSSM files are written like this:

- `-out_pssm FILE` and `-out_ascii_pssm FILE` get the PSSM at the start of each later round. With `-in_pssm` or `-in_msa`, they get the input PSSM first.
- `-save_pssm_after_last_round` also saves the PSSM built from the last round's hits.
- `-save_each_pssm` saves the PSSM of every round as `FILE.1`, `FILE.2` and so on, where `FILE` is the name given to `-out_pssm` or `-out_ascii_pssm`.

`-comp_based_stats` takes 0 to 3 (default 2). Round 1 runs in the mode you give. Later rounds, and any `-in_pssm` or `-in_msa` search, run mode 2 and mode 3 as mode 1, and LOSAT prints NCBI's warning once per query.

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-query FILE` | `-` | Protein FASTA query |
| `-in_pssm FILE` | | Start from a saved PSSM instead of `-query` |
| `-in_msa FILE` | | Start from a multiple sequence alignment (Clustal, FASTA, PHYLIP, NEXUS, Sequin or MultAlin) instead of `-query` |
| `-msa_master_idx N` | 1 | With `-in_msa`: the sequence that is the query (1 is the first) |
| `-ignore_msa_master` | off | With `-in_msa`: build the PSSM without counting the query sequence |
| `-db NAME` | | Protein database |
| `-subject FILE` | | Protein FASTA file to search instead of a database |
| `-num_iterations N` | 1 | Number of rounds; 0 means until convergence |
| `-inclusion_ethresh E` | 0.002 | E-value a hit must beat to enter the PSSM |
| `-pseudocount N` | 0 | Pseudocount used when building the PSSM |
| `-evalue E` | 10 | Report hits with an E-value up to this |
| `-out_pssm FILE` | | Write the PSSM (ASN.1 text) |
| `-out_ascii_pssm FILE` | | Write the PSSM as a readable matrix |
| `-save_pssm_after_last_round` | off | Also save the PSSM built after the last round |
| `-save_each_pssm` | off | Save the PSSM of every round |
| `-comp_based_stats N` | 2 | Composition-based statistics, 0 to 3 |
| `-matrix NAME` | BLOSUM62 | Scoring matrix |
| `-gapopen N`, `-gapextend N` | 11, 1 | Gap costs (the matrix's own values for other matrices) |
| `-max_target_seqs N` | 500 | Number of subject sequences to keep |
| `-seg yes\|no\|"W L H"` | no | Filter low-complexity regions of the query |
| `-num_threads N` | 1 | Number of threads; the output is the same at every value |
| `-outfmt FORMAT` | 0 | Output format |
| `-out FILE` | standard output | Output file |

For every option: `losat psiblast -help`.

## Output

| Search | Supported `-outfmt` |
| --- | --- |
| FASTA query (`-query`) | 0, 6, 7, 10, 20 |
| PSSM or alignment query (`-in_pssm`, `-in_msa`) | 0, 6, 7 |

Formats 6, 7, 10 and 20 take a custom field list, for example `-outfmt "6 qseqid sseqid evalue"`. As in NCBI's `psiblast`, `delim=` is read and checked but the columns keep their tab or comma. With `-outfmt 0`, the report shows each round, and the description table of a later round has NCBI's two parts: sequences used in the model and found again, and sequences not found before or not below the threshold before. With `-outfmt 7`, the comment lines count rounds. See [Output formats](../output-formats.md).

## Differences from NCBI BLAST+

- Where NCBI's result depends on undefined behavior, LOSAT computes the value without it or refuses the input ([the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md)).

## What LOSAT refuses

Not yet supported:

- `-phi_pattern`.
- `-outfmt` other than 0, 6, 7, 10 and 20 (with `-in_pssm` or `-in_msa`, other than 0, 6 and 7). `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors.

Also refused:

- `-remote`, `-entrez_query`, `-html`, `-import_search_strategy`, `-gilist` and `-ipglist`.
- `-use_sw_tback`.
- With `-in_pssm` or `-in_msa`: any matrix or gap cost other than BLOSUM62 with 11/1. With `-in_pssm`: a PSSM in XML form.
- With `-in_msa`, where NCBI BLAST+ crashes: a query sequence of X or O alone, and `-matrix IDENTITY`.
- With `-in_msa`, where NCBI BLAST+ crashes or reads memory that the file does not give: a PHYLIP header with 0 sequences, an empty command (a lone `;`) in a NEXUS file, a FASTA alignment whose rows have different lengths where the shorter row is read past its end, a `secondary-accession` range that runs backwards or has a number of more than 18 digits. See [Refusals](../refusals.md).
- A `-subject` file whose sequences all have no letters.

`psiblast` has no `-task` and no `-ungapped`. For the full list, see [What LOSAT refuses](../refusals.md).

## Examples

Search a database for three rounds and print a table:

<!-- example: psiblast-1 -->
```bash
losat psiblast -query query.faa -db proteins -num_iterations 3 -outfmt 6
```

Iterate until convergence, keep the final PSSM, then search again from that PSSM:

<!-- example: psiblast-2 -->
```bash
losat psiblast -query query.faa -db proteins -num_iterations 0 -out_pssm query.pssm -out_ascii_pssm query.pssm.txt -save_pssm_after_last_round -outfmt 7
losat psiblast -in_pssm query.pssm -db proteins -comp_based_stats 1 -outfmt 6
```

Search a FASTA file instead of a database:

<!-- example: psiblast-3 -->
```bash
losat psiblast -query query.faa -subject subjects.faa -outfmt 6
```

## See also

- [blastp](blastp.md), the search that round 1 runs
- [makeblastdb](makeblastdb.md), to build a protein database
- [rpsblast](rpsblast.md), to search a database of profiles with a protein query
- [Output formats](../output-formats.md)

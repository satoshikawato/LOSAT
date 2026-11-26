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
- **Database filters.** `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-no_taxid_expansion`, `-seqidlist`, `-negative_seqidlist`, `-gilist`, `-negative_gilist`, `-ipglist` and `-negative_ipglist` keep every round to the sequences they select.
- **Alignment.** `-in_msa FILE` starts from a multiple sequence alignment in Clustal, FASTA (with gaps), PHYLIP, NEXUS, Sequin or MultAlin format; LOSAT finds the format as NCBI BLAST+ does, from the first lines of the file. LOSAT builds the PSSM from the alignment as NCBI BLAST+ does, and searches it as it searches an `-in_pssm` PSSM, with the same rounds and PSSM files. The first sequence is the query. `-msa_master_idx N` takes the N-th sequence instead (1 is the first), and `-ignore_msa_master` builds the PSSM without counting the query sequence. A lower-case letter of the query sequence is left out of the alignment. The PSSM from the alignment is built with a pseudocount of 0; `-pseudocount` applies to the PSSMs of later rounds. `-matrix` and the gap costs name the scoring of the PSSM. The text after the sequence name in a FASTA, PHYLIP or NEXUS file is the sequence's title: it appears in the `# Query:` line, in the PSSM file and in the messages, without its `[name=value]` modifiers. The modifiers are read as NCBI reads them (a warning for one it does not know or that is discouraged, an error for a value it rejects), and what they set (organism, molecule type, comments, publications, protein names, strand and topology) is written with the query in the PSSM file and the search strategy as NCBI writes it. With a `molecule` or `mol-type` modifier that sets the query's molecule to `other`, the search is as for a protein, `-outfmt 0` writes the whole report with a line in place of each alignment (`Sequence with id <id> no longer exists in database...alignment skipped`, one for each HSP, no heading and no score lines, the id being the subject's: its own id, `BL_ORD_ID:<oid>` for a database sequence without ids, `Subject_N` for the N-th `-subject` record), as NCBI does, and a tabular report (`-outfmt 6` or `7`) whose fields read the query's letters (the default fields do) ends at its first row with NCBI's `CSeqVectorException`, exit status 255, after the rows and comment lines before it, as NCBI does. A master sequence that is the start of a longer sequence above it in the file is refused: put the master above the longer one (see [Refusals](../refusals.md)). In the Rust library (`Inputs::msa` / `Inputs::msa_database`) and `pylosat` (`Alignment.from_file`), pass the alignment as the query; see [Using LOSAT as a library](../library.md).
- **PSSM file.** `-in_pssm` reads a PSSM in NCBI's ASN.1 text form, the form `-out_pssm` writes, in its ASN.1 binary form, or in NCBI's XML form (as `datatool -px` writes it), with NCBI's messages and errors for each form. It is searched with the BLOSUM62 matrix and gap costs 11/1.
  - The file's Karlin-Altschul values are used as they are, also when they are negative, infinite (`PLUS-INFINITY`) or not a number (`NOT-A-NUMBER`): E-values then print `0.0` and bit scores `inf` or `-1.0`, as in NCBI BLAST+. A query whose Seq-inst omits `length` stops with NCBI BLAST+'s query set-up error (`Warning: Cannot get sequence length `, exit status 3), as in NCBI BLAST+. A query whose Seq-inst `length` is one less than its letters is searched as its first `length` letters, as in NCBI BLAST+. A PSSM of frequency ratios only whose `length` is one less than its letters or more stops with NCBI BLAST+'s PSSM engine error (exit status 3), as in NCBI BLAST+. Other lengths are refused. With a gapped Lambda so small that NCBI BLAST+'s search for the hit-saving cutoff passes 32-bit integers (such as 1e-9), or a gapped Lambda below 0, no score reaches the cutoff and the report has no hits; NCBI BLAST+ never ends there, or writes the same report (see [Differences from NCBI BLAST+](../differences.md#-in_pssm-with-a-very-small-or-negative-gapped-lambda)).
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
| FASTA query (`-query`) | 0 to 12, 18, 20 |
| PSSM or alignment query (`-in_pssm`, `-in_msa`) | 0 to 12, 18, 20 |

Formats 6, 7, 10 and 20 take a custom field list, for example `-outfmt "6 qseqid sseqid evalue"`. As in NCBI's `psiblast`, `delim=` is read and checked but the columns keep their tab or comma. With `-outfmt 0`, the report shows each round, and the description table of a later round has NCBI's two parts: sequences used in the model and found again, and sequences not found before or not below the threshold before. With `-outfmt 7`, the comment lines count rounds. With `-outfmt 5`, one XML document holds every round, one `Iteration` per round. With `-outfmt 1` to `4`, each round's report has its query-anchored section in place of the alignments of `-outfmt 0`, and the report ends with the footer of `-outfmt 0`. With `-outfmt 8`, `9` and `12`, each round writes its own Seq-annot or Seq-align set. With `-outfmt 18`, each round writes its organism report, from every hit of the round: as in NCBI, that report is not cut to `-max_target_seqs`, so a round can show more subjects than the other formats. A PSSM whose query has the molecule `other` is searched as a protein, as in NCBI: the alignments of `-outfmt 0` are NCBI's line `Sequence with id ... no longer exists in database...alignment skipped`, and the query-anchored section of `-outfmt 1` to `4` is an empty line. NCBI reads such a query as a nucleotide sequence when it merges a round's alignments for that section; at a round whose merge throws NCBI's `CAlnException`, LOSAT stops with NCBI's exception (exit status 255) after the report before the section, and at a round where NCBI reads past the end of a row's starts, crashes, reads memory it has freed or computes a length below zero, LOSAT refuses that round, after the report before it. Formats 13 to 16 are not yet supported. With `-outfmt 11`, each round writes its own archive, and the archive of a later round holds the PSSM that the round searched. With `-subject`, `-outfmt 11` stops at the first round that searches a PSSM (round 2 for a FASTA query, round 1 for `-in_pssm` or `-in_msa`), after the archives of the rounds before it, with NCBI's error `Empty reference for CSearchDatabase.` and exit status 3, as NCBI BLAST+ does. See [Output formats](../output-formats.md).

## Differences from NCBI BLAST+

- Where NCBI's result depends on undefined behavior, LOSAT computes the value without it or refuses the input ([the decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md)).

## What LOSAT refuses

Not yet supported:

- `-phi_pattern`.
- `-outfmt` 13 to 16. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors.

Also refused:

- `-remote`, `-entrez_query`, `-html` and `-import_search_strategy`.
- `-use_sw_tback`.
- With `-in_pssm` or `-in_msa`: any matrix or gap cost other than BLOSUM62 with 11/1.
- With `-in_msa`, where NCBI BLAST+ crashes: a query sequence of X or O alone, and `-matrix IDENTITY`.
- With `-in_msa`, where NCBI BLAST+ crashes or reads memory that the file does not give: a PHYLIP header with 0 sequences, an empty command (a lone `;`) in a NEXUS file, a FASTA alignment whose rows have different lengths where the shorter row is read past its end, a `secondary-accession` range that runs backwards or has a number of more than 18 digits. See [Refusals](../refusals.md).
- A `-subject` file whose sequences all have no letters.
- `-outfmt` 1 to 4 with a PSSM whose query has the molecule `other`, at a round whose merge reads past the end of a row's starts, crashes, reads memory it has freed or computes a length below zero in NCBI BLAST+. The rounds before it are written.

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

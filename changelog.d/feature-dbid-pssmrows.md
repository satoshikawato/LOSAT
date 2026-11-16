- psiblast `-db` with parsed Seq-ids: the PSSM (`-out_pssm`, `-out_ascii_pssm`, and each round of
  `-num_iterations` other than 1) now builds its multiple alignment as NCBI BLAST+ 2.17.0 does, with
  one row for each Seq-align id and following hits whose ids NCBI takes for one sequence in one row.
  Copies of one record that follow each other (the same sequence in two databases of `-db "a b"`)
  gave LOSAT another PSSM and other later rounds than NCBI's; hits such as `gb|X|` and `gb|X.1|`
  that follow each other are no longer refused. Refused instead, where psiblast compares them: ids
  that NCBI can take for one sequence though they differ (in letter case, or one accession under
  `gb`, `emb` and `dbj`), and two copies of one record with another hit between them (NCBI crashes
  there).
- `-db` with parsed Seq-ids (all six programs): a sequence whose Seq-align id is a `gim|` id is now
  refused where a report reads the sequence, as for the other ids NCBI BLAST+ 2.17.0 does not find
  again: NCBI's id handles never take a `gim|` id for one with itself, so it writes `Unknown` (and in
  BLASTN, on the minus strand, `Error pre-fetching sequence data`), where LOSAT wrote the id.
- psiblast `-db` with parsed Seq-ids, in rounds (`-num_iterations` other than 1) or with a PSSM file:
  a reported sequence that no PSSM reads (its alignments above `-inclusion_ethresh`) and that the
  formats do not read (`6 qseqid evalue`) is now written as NCBI BLAST+ 2.17.0 writes it, also where
  its Seq-align id is one that a report reading the sequence refuses (a `gim|` id, an id the lookup
  finds as another sequence). LOSAT refused every such run.

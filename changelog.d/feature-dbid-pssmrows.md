- psiblast `-db` with parsed Seq-ids: the PSSM (`-out_pssm`, `-out_ascii_pssm`, and each round of
  `-num_iterations` other than 1) now builds its multiple alignment as NCBI BLAST+ 2.17.0 does, with
  one row for each Seq-align id and following hits whose ids NCBI takes for one sequence in one row.
  Copies of one record that follow each other (the same sequence in two databases of `-db "a b"`)
  gave LOSAT another PSSM and other later rounds than NCBI's; hits such as `gb|X|` and `gb|X.1|`
  that follow each other are no longer refused. Refused instead: ids that NCBI can take for one
  sequence though they differ in letter case or are one accession under `gb`, `emb` and `dbj`, and
  two copies of one record with another hit between them (NCBI crashes there).

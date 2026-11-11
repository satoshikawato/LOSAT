- psiblast `-db` with parsed Seq-ids: the PSSM (`-out_pssm`, `-out_ascii_pssm`, and each round of
  `-num_iterations` other than 1, also with `-in_pssm`) now reads the letters of each included hit
  by its Seq-align id as NCBI BLAST+ 2.17.0 does. Where the database's lookup of that id gives
  another sequence first (`gb|X|` of one sequence, `gb|X.1|` of another), the PSSM takes that
  sequence's letters; where that sequence is too short for the hit, or none is found, NCBI's warning
  `Failed to retrieve sequence` is written for the hit and its row is left out of the PSSM. LOSAT
  refused every such run.

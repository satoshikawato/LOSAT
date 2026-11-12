- psiblast `-db` with parsed Seq-ids: the PSSM (`-out_pssm`, `-out_ascii_pssm`, and each round of
  `-num_iterations` other than 1, also with `-in_pssm`) now reads the letters of each included hit
  by its Seq-align id as NCBI BLAST+ 2.17.0 does. Where the database's lookup of that id gives
  another sequence first (`gb|X|` of one sequence, `gb|X.1|` of another), the PSSM takes that
  sequence's letters; where that sequence is too short for the hit, or none is found, NCBI's warning
  `Failed to retrieve sequence` is written for the hit and its row is left out of the PSSM. LOSAT
  refused every such run.
- BLASTN, TBLASTX and psiblast `-db` with parsed Seq-ids: the prefetch reads by each Seq-align id
  the sequence that NCBI's scope gets, also for a hit whose definition line the reports do not find:
  BLASTN's and TBLASTX's `Error pre-fetching sequence data` for a minus-strand alignment, and the
  taxonomy database's texts of the sequence read, are now NCBI's.
- psiblast `-db` where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off and the taxonomy
  database is damaged: the taxonomy database's texts of the lookups that the PSSM's reads make
  (also for a database without parsed Seq-ids, read by `gnl|BL_ORD_ID|n`, and where a PSSM error
  stops the run) are now written where NCBI BLAST+ 2.17.0 writes them, after the query's report and
  among the PSSM's warnings. The offset error of the last entry comes as often as NCBI's 256-slot
  taxonomy cache of a volume asks for it: twice for one sequence whose definition lines hold its
  taxid before and after another taxid of the same slot, again for a later sequence after such a
  taxid. LOSAT wrote nothing there, or refused the run. Still refused: a later query or round that
  reads such a sequence again where NCBI's count depends on whether it loads it again.

- `-db` with databases built with `-parse_seqids`: copies of one record (the same definition lines but for their
  tax ids, in two volumes or in two databases of an alias or of `-db "a b"`) that follow each other in a query's
  hits are no longer refused in outfmt 0, 5 and 13 to 16 (all six programs). They are shown in one hit, as NCBI
  BLAST+ 2.17.0 shows them: one description line and one heading, the names and length of the first copy, its
  HSPs and then those of the next (`Hsp_num` going on in outfmt 5 and 13 to 16). `-num_descriptions` and
  `-num_alignments` count such a hit as one subject, and so does `psiblast`'s cut of its longer list to
  `-max_target_seqs` in every format (it keeps the copy that follows the last subject kept; outfmt 6 writes a row
  for each copy). Copies with other tax ids are still refused in outfmt 13 to 16. Other sequences whose ids NCBI
  takes for one (`pdb|1ABC|` and `pdb|1ABC|AA`, `gb|X|` and `gb|X.1|`) are still refused in those formats.

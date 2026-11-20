- `losat windowmasker -text_match F` with `-infmt blastdb`: the `-ids` and `-exclude_ids` lists are
  matched as Seq-ids against each database sequence's Seq-ids, as NCBI BLAST+ 2.17.0's
  `windowmasker` does (it was refused); an id that is not a Seq-id is written as NCBI's error and
  skipped.

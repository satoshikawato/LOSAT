- `losat windowmasker -text_match F` with `-infmt blastdb`: the `-ids` and `-exclude_ids` lists are
  matched as Seq-ids against each database sequence's Seq-ids, as NCBI BLAST+ 2.17.0's
  `windowmasker` does (it was refused); an id that is not a Seq-id is written as an error
  (`can't understand id`) and skipped. With `-outfmt seqloc_*` or `maskinfo_*`, the Seq-ids of the
  sequences are written as NCBI writes them once a list has named them (a DDBJ accession of the
  list makes GenBank ids of the same prefix come out as DDBJ ids).
- `losat dustmasker`, `losat segmasker` and `losat windowmasker` with `-outfmt seqloc_*` or
  `maskinfo_*`: a sequence with several Seq-ids is named by its gi (not gi 0), else its first accession-type
  id (GenBank, RefSeq, ...), else its first Seq-id, as NCBI BLAST+ 2.17.0 does (LOSAT wrote the
  first Seq-id).

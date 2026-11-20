- `-db` with a database built with `-parse_seqids`: a sequence whose only Seq-id is a stored `gnl|BL_ORD_ID|n`
  with n other than its number in its volume (makeblastdb keeps such a FASTA id as written) is now read as a
  sequence with parsed Seq-ids. NCBI BLAST+'s formatters look it up by n plus the first ordinal id of its
  volume and find no sequence or another one, so LOSAT refuses it in the formats that read the subject's
  sequence (outfmt 0, 5, 6, 7, 10, 13 to 17, 20) and writes n plus the first ordinal id of the volume in the
  Seq-annots and the archive (8, 9, 11, 12), as NCBI BLAST+ 2.17.0 does. Before, LOSAT wrote the sequence's
  own number and title in every format, where NCBI BLAST+ writes `Unknown` or another sequence's names, or
  stops.
  A sequence whose `gnl|BL_ORD_ID|n` is its own number is read as before.

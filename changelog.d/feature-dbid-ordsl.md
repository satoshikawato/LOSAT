- `-seqidlist` and `-negative_seqidlist` naming `gnl|BL_ORD_ID|n`, where a sequence with such an id lies in a volume whose
  first OID is not 0 (a later database of a `-db` list, or a later volume of a database), most often a database made with
  `makeblastdb -input_type blastdb` from one without `-parse_seqids`: NCBI BLAST+ finds the sequence by the number stored in
  the accession table, then tests its definition lines with the volume's first OID added. A positive list can thus select a
  sequence whose every line it then removes, and NCBI leaves its hits out of the report without a word; LOSAT wrote them.
  LOSAT now refuses such a search, naming the sequence, as for other seqid lists that keep a sequence without a line: a
  positive list before the search, whether or not that sequence would be hit, a negative list when the sequence is to be
  reported. Where such a sequence keeps another line, LOSAT now writes the definition lines NCBI writes (it tested the stored
  number and kept or dropped the other line). The sequences of the first volume are searched as before.

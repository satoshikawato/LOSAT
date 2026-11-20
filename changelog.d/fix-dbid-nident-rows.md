- BLASTP, PSIBLAST, TBLASTN, BLASTX and TBLASTX `-db` with `PRE_FETCH_SEQS_LIMIT` turning NCBI's
  prefetch off and a tabular field list whose only field that reads the subject is `nident` (such
  as `6 qseqid nident` or, in the translated programs, `6 qseqid nident sframe`) now look a subject
  up only where one of its shown HSPs has identities, as NCBI BLAST+ 2.17.0's formatter reads the
  subject for `nident` only on such a row: a subject whose shown HSPs have no identity (possible
  with a low `-threshold`) no longer gets a taxonomy database message that NCBI does not write.
  LOSAT wrote it there with exit status 0.

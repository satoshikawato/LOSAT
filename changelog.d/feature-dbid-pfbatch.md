- `-db` with parsed Seq-ids (BLASTN, TBLASTN, BLASTX, TBLASTX): with a positive
  `PRE_FETCH_SEQS_LIMIT`, a tabular field list that reads nothing of a reported sequence (for
  example `-outfmt "6 qseqid evalue"`) is no longer refused for what NCBI BLAST+'s prefetch cannot
  read (a Seq-align id whose lookup fails; a minus-strand alignment of a sequence it gets no
  Bioseq of) where every batch that reports the sequence has more Seq-aligns than the limit:
  NCBI's prefetch reads nothing in such a batch, and LOSAT writes NCBI BLAST+ 2.17.0's rows. Each
  batch counts its own Seq-aligns, as BLASTP and psiblast already did.

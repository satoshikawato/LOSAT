- `-db` with parsed Seq-ids (all six programs): where NCBI BLAST+'s prefetch of a query batch
  fails and no report reads the sequence (a tabular field list such as `-outfmt "6 qseqid evalue"`,
  outfmt 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647`, psiblast's sequences beyond
  `-max_target_seqs`), LOSAT writes NCBI's `Warning: [<program>] Error pre-fetching sequence data`
  for that batch and the rows, as NCBI BLAST+ 2.17.0 does, instead of refusing the run. The
  prefetch fails at a Seq-align id whose lookup fails (`gnl|BL_ORD_ID|0` that makeblastdb stores
  with a text tag) and, in BLASTN and TBLASTX, at a minus-strand alignment of a sequence it gets
  no Bioseq of; BLASTN `-ungapped` now takes that strand from the sequence's first HSP, as NCBI's
  one Seq-align of an ungapped search does. A failed lookup in a batch whose other sequences make
  the taxonomy database write a message is still refused: NCBI loads them in an order that LOSAT
  does not reproduce.

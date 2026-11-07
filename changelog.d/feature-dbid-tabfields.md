- `-db` with parsed Seq-ids (all programs): a tabular field list that reads nothing of a reported
  sequence (for example `-outfmt "6 qseqid evalue"`, the coordinates, `length`, `gaps`, `qcovhsp`,
  `qcovs`; in BLASTN also `pident`, `nident` and `mismatch`) is no longer refused for a sequence
  that NCBI BLAST+'s formatters find as another sequence or not at all, nor for a stored title
  whose id it does not print: it writes NCBI BLAST+ 2.17.0's rows. LOSAT refused such runs, where
  NCBI exits 0 with the same rows for any lookup. What NCBI's prefetch of the batch reads for such
  a sequence is refused as for outfmt 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647`; psiblast
  with a PSSM file, whose PSSM NCBI makes from the letters it reads by the Seq-align id, is refused
  as before.

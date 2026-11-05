- TBLASTN and TBLASTX `-db` (TBLASTX where it writes every report after its batch loop: one query
  batch, a custom tabular field list, outfmt 10 and 20, the XML and Seq-annot formats) no longer refuse
  a run for the taxonomy database whose damaged last entry shares NCBI's taxonomy cache slot with
  another taxid of a volume when NCBI BLAST+ 2.17.0 fetches no subject at all: fields that read none
  (such as `6 qseqid evalue`) with `PRE_FETCH_SEQS_LIMIT` turning the prefetch off (a positive limit
  below a batch's alignments), a query batch without hits, and a query without hits in a batch whose
  prefetch is off. LOSAT refused these runs, which NCBI completes without writing any taxonomy
  message, with exit status 1. TBLASTX's per-batch report stream (outfmt 0, 6 and 7 with two or more
  batches) still refuses such a run at its first batch.

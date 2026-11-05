- TBLASTN and TBLASTX `-db` no longer refuse a run for a taxonomy database whose damaged last entry
  shares NCBI's taxonomy cache slot with another taxid of a database volume when NCBI BLAST+ 2.17.0
  fetches no subject: fields that read none (such as `6 qseqid evalue`) with `PRE_FETCH_SEQS_LIMIT`
  turning the prefetch off (also a positive limit below a batch's alignments), and query batches or
  queries without hits. LOSAT refused these runs, which NCBI completes without any taxonomy message,
  with exit status 1. TBLASTX writes the reports of a run with such a taxonomy database after its
  batch loop (not batch by batch), so a later batch's refusal still comes before any report.

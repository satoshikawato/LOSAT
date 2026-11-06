- BLASTN `-db` no longer refuses a run for a taxonomy database whose damaged last entry shares NCBI's
  taxonomy cache slot with another taxid of a database volume when NCBI BLAST+ 2.17.0 fetches no
  subject: query batches without hits, queries without hits whose reports would fetch with
  `PRE_FETCH_SEQS_LIMIT` turning the prefetch off, and fields that read none (such as
  `6 qseqid evalue`) with a positive `PRE_FETCH_SEQS_LIMIT` below a batch's alignments. LOSAT refused
  these runs, which NCBI completes without any taxonomy message, with exit status 1. A run that
  fetches a subject is still refused, at the first query batch that fetches one; where BLASTN writes
  each batch's reports as it takes the batch, that refusal comes after the reports of the batches
  before.

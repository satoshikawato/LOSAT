- PSIBLAST `-db` with `taxdb.bti` on the search path now writes NCBI BLAST+ 2.17.0's messages about a
  damaged taxonomy database (a wrong magic number, the count of entries, an offset error) in outfmt
  0, 6, 7, 10 and 20, where NCBI's prefetch of each query's subjects writes them: also for subjects
  of its traceback hit list beyond `-max_target_seqs`, which no report shows. LOSAT wrote none, with
  exit status 0; an empty `taxdb.btd`, where NCBI psiblast dumps core, is refused as BLASTP refuses
  it. With `PRE_FETCH_SEQS_LIMIT` turning NCBI's prefetch off for a query batch (0, or below the
  batch's number of alignments), BLASTP and PSIBLAST write these messages after each query's own
  warnings, for the subjects its report shows, as NCBI does; BLASTP wrote them before.

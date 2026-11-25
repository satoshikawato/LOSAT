- When a `blastp` search of several query batches (`BATCH_SIZE`, or queries whose letters fill a batch) with the
  compressed lookup table (`-task blastp-fast`, or a word size of 5 to 7) meets a batch whose neighboring words fill
  NCBI BLAST+'s 1024 overflow banks (a very low `-threshold`), LOSAT now writes the results of the query batches before
  that batch and then stops with its error, as NCBI BLAST+ 2.17.0 prints each batch's results before it builds the next
  batch's table (and then corrupts its heap); before, it stopped before any output. A split query's piece whose table
  fills them stops its batch in the same place. Searches that do not fill the banks are unchanged.

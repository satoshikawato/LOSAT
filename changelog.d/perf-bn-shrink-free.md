- BLASTN frees each query batch's lookup-table CSR and BLASTNA query buffer without raising glibc's mmap threshold
  (`LOSAT_X_BNSHRINKFREE`), so the heap no longer keeps freed space resident between batches (to be measured:
  dc-megablast of 31,168 cDNAs at one thread).

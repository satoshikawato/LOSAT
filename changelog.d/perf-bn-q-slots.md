- BLASTN `-subject` searches whose query batches are searched ahead in slots (`LOSAT_X_QBATCHPAR`, 2+ threads with a
  query file NCBI splits) no longer hold the whole query file: the query store of `LOSAT_X_BNQDROP` / `LOSAT_X_BNQPACK`
  now serves the slots too (`LOSAT_X_BNQSLOTS`; to be measured: dc-megablast of 31,168 cDNAs at 32 threads).

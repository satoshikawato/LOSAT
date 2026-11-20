- `blastn` with a query batch split into chunks (`CHUNK_SIZE`) over a database whose average
  sequence length is too short now stops with NCBI BLAST+ 2.17.0's error (`BLAST engine error: The
  average subject length is too short`, exit status 3) also with megablast gap costs above 32767
  or reward 0 and penalty 0; LOSAT refused those options for the whole batch first. As in a batch
  that is not split, the subjects are checked before the greedy extension's set-up.

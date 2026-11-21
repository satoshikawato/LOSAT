- `blastn` with a query batch split into chunks (`CHUNK_SIZE`) over a database whose average
  sequence length is too short now stops with NCBI BLAST+ 2.17.0's error (`BLAST engine error: The
  average subject length is too short`, exit status 3) also with megablast gap costs above 32767
  or reward 0 and penalty 0; LOSAT refused those options for the whole batch first. As in a batch
  that is not split, the subjects are checked before the greedy extension's set-up.
- `blastn` refuses matrix-only scoring with the greedy extension in a later query batch after
  writing the reports of the batches before, the output NCBI BLAST+ writes before it crashes
  there. The refusal of megablast gap costs above 32767 is now listed in the manual.

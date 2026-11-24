- blastdbcmd: a sequence whose definition lines include one without a Seq-id (a hand-edited database;
  makeblastdb writes none) is written as NCBI BLAST+ 2.17.0 writes it instead of being refused: the FASTA
  header is left out when the first line has no Seq-id, and `-range`, `-mask_sequence_with`, `%a` and `%i`
  of such a sequence end the run with NCBI's null pointer message (exit 1) after the sequences before it.

- BLASTX's composition-based statistics solve their Newton problems four at a time across the frames of a match
  (`LOSAT_X_BXNEWTONLANES`, on by default; `=0` solves each one where it arises): blastx of AP027131 against the
  AP027133 proteins takes 0.933 of the cycles at one thread and 0.944 at eight threads on the workstation, and 0.977
  of the CPU in the Wasm command. A thread keeps the work of about 32 problems at a time (16 threads: 73.0 -> 76.6 MiB).
  Outputs do not change.

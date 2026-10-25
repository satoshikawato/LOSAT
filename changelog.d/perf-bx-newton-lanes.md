- BLASTX's composition-based statistics solve their Newton problems four at a time across the frames of a match
  (`LOSAT_X_BXNEWTONLANES`, on by default; `=0` solves each one where it arises): blastx of AP027131 against the
  AP027133 proteins takes 0.946 of the cycles at one thread on the workstation (3.66 -> 3.46 s) and 0.97 in the Wasm
  command. At eight threads 0.946 (0.740 -> 0.708 s). Outputs do not change.

- LOSAT Web gives one file where two options of a psiblast run name the same file, with the command line's bytes:
  `-logfile L -out_pssm L` gives `L` as the PSSM file (39,579 bytes in the test run, as NCBI BLAST+ 2.17.0 and
  LOSAT's command line make it), labeled as the `-out_pssm` file; `-logfile E.1 -out_pssm E -save_each_pssm` gives
  `E.1` as the first round's PSSM; `-logfile L -out_ascii_pssm L` gives the ASCII PSSM. Before, the Web GUI showed
  the file as the `-logfile` file, empty or holding the log's messages in place of the PSSM.

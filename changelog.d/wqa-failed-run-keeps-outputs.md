- LOSAT Web: a run that fails keeps the output that the command line writes before its error, byte for byte, and the
  results page shows it with the error, with Export, as it shows the output of a completed run. Before, a failed run kept
  no output, where the command line leaves the batches written before the error in its output file (the outfmt 0
  prolog, the hits of the batches before the failing one). The other outputs of that run (outfmt 0, 6 and 7 besides the
  command line's own) are not kept. Engine ABI version 7.

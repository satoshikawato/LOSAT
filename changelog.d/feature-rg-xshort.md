- `tblastn` with an alignment on a subject of fewer than 12 bases that reaches the traceback (in practice with
  `-window_size 0`) now stops with an error that names the cause instead of crashing. NCBI BLAST+ 2.17.0 crashes
  there too: it reads a translation of the subject's frame that it never makes for such a short subject. Searches
  where no alignment on such a subject reaches the traceback give NCBI BLAST+'s output as before.

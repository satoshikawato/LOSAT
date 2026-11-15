- Engine: a BLASTP call's BLAST archive (outfmt 11), its tabular reports (outfmt 6, 7, 10
  and 20) and its organism report (outfmt 18) are written by their own functions, moved
  out of `run_call` into their own module. No output changes.

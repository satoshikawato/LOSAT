- Engine: a BLASTP call's Seq-annots (outfmt 8, 9 and 12), its BLAST XML, XML2 and JSON
  reports (outfmt 5 and 13 to 16) and its pairwise and query-anchored reports (outfmt 0 to
  4) are written by their own functions, moved out of `run_call` into the module of the
  other reports. No output changes.

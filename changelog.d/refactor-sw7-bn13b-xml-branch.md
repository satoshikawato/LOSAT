- Engine: a BLASTN call's BLAST XML, XML2 and JSON reports (outfmt 5 and 13 to 16) are written
  by their own function, moved out of `post_process_hits_and_write` into the module of the
  call's other report writers. No output changes.

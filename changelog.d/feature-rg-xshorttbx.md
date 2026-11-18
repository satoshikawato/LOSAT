- A `tblastx` search that finds an alignment on a subject (or a `-subject_loc` range) shorter than 12 bases, in
  practice only with `-window_size 0`, now stops with an error that names the cause instead of writing rows. NCBI
  BLAST+ 2.17.0 crashes on the same search: it reads a translation of the subject that it never makes for so short a
  subject. The reports of the query batches before are written in full first (NCBI BLAST+ writes them before it
  crashes, except what its output buffer still holds: with `-outfmt 0`, the end of the last one); searches where
  no such alignment is found give NCBI BLAST+'s output as before.
- A `tblastx` search of a database under `-negative_seqidlist` that stops at a sequence header it cannot read now
  writes the reports of the query batches before it first, as NCBI BLAST+ does before its engine error; before, it
  wrote none of them.

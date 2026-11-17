- A `tblastx` search that finds an alignment on a subject (or a `-subject_loc` range) shorter than 12 bases, in
  practice only with `-window_size 0`, now stops with an error that names the cause instead of writing rows. NCBI
  BLAST+ 2.17.0 crashes on the same search: it reads a translation of the subject that it never makes for so short a
  subject. The reports of the query batches before are written, as NCBI BLAST+ writes them before it crashes;
  searches where no such alignment is found give NCBI BLAST+'s output as before.

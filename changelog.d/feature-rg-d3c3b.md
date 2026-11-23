- psiblast and tblastn `-in_pssm` read NCBI's pointer forms in a PSSM, in ASN.1 text and binary,
  as NCBI BLAST+ 2.17.0 does: `NULL` for an optional member (`params`, `identifier`,
  `intermediateData`, `finalData`), `:Type value` of the declared type (`params :PssmParameters`, a
  query Seq-id `:Seq-id`), `@n` with NCBI's `Unrecognized format for PSSM`. A PSSM without a
  query, or with `query NULL`, ends with NCBI's exception of its empty Seq-entry (exit status 255),
  where LOSAT refused it. A value stated as `:Type` of another type than its place declares ends
  with NCBI's `Unrecognized format for PSSM` (exit status 1), as in NCBI BLAST+ 2.17.0. A list
  element stated as `NULL` is read as NCBI reads it: in the query's descriptors it ends the run
  with NCBI's `Attempt to access NULL pointer.` engine error (exit status 3) on every output
  format; where NCBI's program does not use the list the search runs, and `-out_pssm` and
  `-export_search_strategy` end at the element with NCBI's write error (exit status 255).
  Still refused: a `:Type` value whose reading prints a message, meets another such value or
  stops with an error that NCBI does not report as an unrecognized format, or that sits in a
  skipped part, and a `NULL` element of a list that NCBI's object manager dereferences (a Seq-id
  of the query, an annotation, a feature, a location part and the like: NCBI's report of memory
  addresses or a crash).
- makeblastdb `-mask_data` reads its mask files with the table-driven reader of `-in_pssm` (one
  ASN.1 reader for text and binary); outputs unchanged.

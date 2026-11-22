- psiblast and tblastn `-in_pssm` read NCBI's pointer forms in a PSSM, in ASN.1 text and binary,
  as NCBI BLAST+ 2.17.0 does: `NULL` for an optional member (`params`, `identifier`,
  `intermediateData`, `finalData`), `:Type value` of the declared type (`params :PssmParameters`, a
  query Seq-id `:Seq-id`), `@n` with NCBI's `Unrecognized format for PSSM`. A PSSM without a
  query, or with `query NULL`, ends with NCBI's exception of its empty Seq-entry (exit status 255),
  where LOSAT refused it. Still refused: a value stated as `:Type` of another type than its place
  declares, and a list element stated as `NULL` that NCBI's program dereferences later.
- makeblastdb `-mask_data` reads its mask files with the reader of `-in_pssm` and rpsbproc (one
  ASN.1 reader for text and binary); outputs unchanged.

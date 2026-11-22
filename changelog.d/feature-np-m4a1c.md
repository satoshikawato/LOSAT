- psiblast and tblastn `-in_pssm`: a PSSM query residue that its coding does not have (`*` or `-` in iupacaa, a
  lower-case ncbieaa letter, an ncbistdaa code of 28 or more, a byte NCBI's reader repairs as `#`) now stops the run
  as NCBI BLAST+ 2.17.0 does: `LOSAT engine error: Warning: Invalid residues found at positions <list> ` (the first
  20), exit status 3, after the outfmt 0 prolog, the search strategy and psiblast's round-1 `-out_pssm`. LOSAT refused
  such a PSSM, and searched a `*` of iupacaa as the stop residue. `-out_ascii_pssm` of such a query is refused (NCBI
  prints the residue's letter from memory it never set).
- psiblast and tblastn `-in_pssm`: a PSSM query in lower-case iupacaa letters or with the gap residue (`-` of
  ncbieaa, 0 of ncbistdaa) is now searched and reported as NCBI BLAST+ 2.17.0 does: the search reads lower-case
  letters as their capitals, and the reports show them in lower case, identical or positive to no subject letter
  (pident, nident, positive, mismatch, btop, the middle line, Identities and Positives, tblastn's outfmt 1 to 4); the
  gap residue shows as `-`, identical to a gap of the subject and no end gap. `-out_pssm` keeps the lower case and
  `-out_ascii_pssm` prints the capitals; psiblast's round 2 stops on the gap residue's column with NCBI's error.
  LOSAT refused both forms.
- makeblastdb `-mask_data`: a mask file with a string byte outside the visible ASCII range is read as NCBI BLAST+
  2.17.0 reads it (`#` for the byte, with NCBI's error message naming where it is), and a `UTF8String` where a
  `VisibleString` belongs with NCBI's warning, once per run; LOSAT refused both. A damaged mask file ends the run with
  makeblastdb's own report and exit status (the exception report, "BLAST ran out of memory" for a string no host can
  hold), where LOSAT wrote its own text.
- makeblastdb `-mask_data`: NCBI's pointer syntax in a mask file is read as NCBI BLAST+ 2.17.0 reads it: `NULL` for
  an optional member such as `fuzz-from`, `@n` with makeblastdb's error report, `:Type value` of the declared type.
  A `:Type` of another type, and a list element stated as `NULL` where makeblastdb dereferences it (its report
  prints the run's memory addresses), are refused. A binary length of 8 bytes whose end wraps past 2^64 now reads
  as makeblastdb reads it.

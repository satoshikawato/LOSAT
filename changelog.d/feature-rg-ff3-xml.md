- psiblast and tblastn `-in_pssm` read a PSSM in NCBI's XML form (as NCBI's `datatool -px` writes
  it), as NCBI BLAST+ 2.17.0 does: the same values, messages and errors as its ASN.1 text and binary
  twins. That is `Unrecognized format for PSSM` for any error of the reader (also for a PSSM in the
  standard XML form without scope prefixes, which NCBI does not read either), NCBI's exception with
  exit status 255 for an INTEGER that is not a number or does not fit and for a `nil` attribute
  whose value is not a boolean, namespace prefixes and declarations read as NCBI reads them, NCBI's
  message (803.7) with the line for a control byte in a string, and an OPTIONAL or DEFAULT
  ENUMERATED without its `value` left unset. LOSAT refused every XML PSSM.
- psiblast `-in_pssm` with `-out_pssm` writes the author names of a PSSM in ASN.1 binary or XML
  (`Auth-list` `names`, which NCBI reads again only when it writes them in ASN.1 text) whose read
  replaces a byte outside the visible ASCII range: NCBI's message (803.7) of that read, then `#` in
  the file, as NCBI does. LOSAT refused them. Still refused: author names whose read fails (a number
  too large for its type), where NCBI writes part of the file and stops with an exception.
- psiblast `-outfmt 11` writes the XML and ASN.1 binary archives (`ARCHIVE_FORMAT=xml`, `bin`) of the
  PSSM files that LOSAT refused there: unnamed SEQUENCE and CHOICE values (a User-object descriptor,
  annotations), NULL, UTF8String (in XML with NCBI's message for a control byte), NaN and infinite
  real numbers in ASN.1 binary, and author names, copied as they were read when the file has the
  archive's format and read again otherwise, with NCBI's message (803.7) for a byte outside the
  visible ASCII range, as NCBI BLAST+ 2.17.0 does. After `-out_pssm` or `-export_search_strategy`,
  which write the PSSM in ASN.1 text first, the archive and a second `-out_pssm` file write the names
  as that first write read them, and the message comes once, as in NCBI. Still refused: author names
  in another format than the archive whose read fails (a number too large for its type).

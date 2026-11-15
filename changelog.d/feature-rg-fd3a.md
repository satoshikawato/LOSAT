- psiblast and tblastn `-in_pssm` read the PSSM's ASN.1 text as NCBI BLAST+ 2.17.0 reads it, every member of every
  type it reaches, with a type table generated from NCBI's ASN.1 modules (`scripts/asn_spec_gen.py`, checked against
  NCBI's generated code by `scripts/asn_spec_check.py`): the PSSMs NCBI refuses (an unknown member, a value of the
  wrong type or without a name, a missing member, an upper-case member name) get NCBI's "Unrecognized format for PSSM"
  where LOSAT searched them; a number outside its member's range and a member tag NCBI cannot convert get NCBI's C++
  exception report (exit status 255); a byte outside the visible ASCII range is read as `#` with NCBI's message (the
  search goes on, as NCBI's), and so are the messages of a `Variation-ref`'s deprecated members. The members that
  LOSAT skipped are kept (line RG item D3-1).
  A query's `annot` (feature tables, alignments, graphs, locations, Seq-tables) is indexed as NCBI's object manager
  indexes it: the PSSM is searched, with NCBI's messages and errors where NCBI writes them; the few cases whose
  NCBI output names a memory address or depends on the machine are refused, each with its reason. The Seq-ids of
  the annotations are one sequence where NCBI's Seq-id index takes them for one (accessions, names, local and
  general ids and patent numbers in any letter case; versions and releases exactly; `giim` ids as NCBI's index of
  them takes them), and their messages come in NCBI's order (by gi, accession and general id number). A Seq-table
  bit-bvector column whose size alone is read (an int-delta over it, an extra fuzz column) is no longer refused. A PSSM file that
  starts with a UTF-8 byte order mark gets NCBI's "Unrecognized format for PSSM".

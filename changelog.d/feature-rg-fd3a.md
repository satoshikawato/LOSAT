- psiblast and tblastn `-in_pssm` read the PSSM's ASN.1 text as NCBI BLAST+ 2.17.0 reads it, every member of every
  type it reaches, with a type table generated from NCBI's ASN.1 modules (`scripts/asn_spec_gen.py`, checked against
  NCBI's generated code by `scripts/asn_spec_check.py`): the PSSMs NCBI refuses (an unknown member, a value of the
  wrong type or without a name, a missing member, an upper-case member name) get NCBI's "Unrecognized format for PSSM"
  where LOSAT searched them; a number outside its member's range and a member tag NCBI cannot convert get NCBI's C++
  exception report (exit status 255); a byte outside the visible ASCII range is read as `#` with NCBI's message (the
  search goes on, as NCBI's), and so are the messages of a `Variation-ref`'s deprecated members. The members that
  LOSAT skipped are kept (line RG item D3-1).
  A query's `annot` (feature tables, alignments, graphs, locations, Seq-tables) is indexed as NCBI's object manager
  indexes it: the PSSM is searched, with NCBI's messages and errors where NCBI writes them. Refused, each with its
  reason: the cases whose NCBI output names a memory address, depends on the machine or reads outside NCBI's data,
  well-formed BitMagic data of a Seq-table (until the Owner decides on its decode), and the annotation's text that
  NCBI prints after its `Failed to parse location of` message. The Seq-ids of the annotations are one sequence
  where NCBI's Seq-id index takes them for one (accessions, names, local and general ids and patent numbers in any
  letter case; versions and releases exactly; `giim` ids by the `giim` handles alive in the run, across the
  annotations), and their messages come in NCBI's order (by gi, accession and general id number). A Seq-table
  bit-bvector column whose size alone is read (an int-delta over it, an extra fuzz column) is no longer refused. A PSSM file that
  starts with a UTF-8 byte order mark gets NCBI's "Unrecognized format for PSSM". With tblastn `-outfmt 11`, a PSSM
  that NCBI cannot read (`Error: vector::reserve` too) or whose annotations stop NCBI's object manager gets NCBI's
  error archive, each of the reader's and the object manager's messages saved whole with its severity; such a stop
  is no longer refused in tblastn's formats, in psiblast's `-outfmt 10` and `20`, for the query's titles, with
  `-out_pssm` (no file written, an existing one left as it was) and `-export_search_strategy` (an empty file), and
  for the PSSM's own values that NCBI reads after the stop (its sizes, scores, statistics, letters and ids). The
  library API's PSSM query (`Inputs::pssm`) ends as the command line's run: a PSSM file that NCBI ends without a word
  gives a call without an error, and an error after the object manager's `Failed to parse location of` message has
  the command line's text (its line feed).
  Two library calls in one process (one after the other, at once on several threads, or two `run_batches`
  iterators advanced in turn) each keep their own diagnostic stream: a call whose stream the object manager fails
  no longer mutes another call's warnings and errors.

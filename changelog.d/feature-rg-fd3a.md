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
  reason: the cases whose NCBI output names a memory address, depends on the machine or reads outside NCBI's data
  (for BitMagic data of a Seq-table, see line RG item Q-BM), and the annotation's text that
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
  for everything that NCBI reads after the stop: the PSSM's own values (its sizes, scores, statistics, letters, their
  alphabet and ids), the Seq-ids of `-subject` records under `-parse_deflines`, the PSSM files, the search strategy and
  LOSAT's limits of the search (`-evalue`, `-use_sw_tback`, gap costs, tblastn's `-xdrop_gap`, `-dbsize` with a list
  that keeps no sequence, report sums). A query Seq-entry set is no longer refused: NCBI's object manager takes its
  annotations (the members' first as it makes them, the set's own first as it indexes them), and NCBI then stops at
  the set (`CSeq_entry::GetSeq()`). A Seq-id twice in the query's entry gets NCBI's `duplicate Bioseq id` exception,
  before the index of the annotations (it was the index's message or LOSAT's refusal), for every two Seq-ids that
  NCBI's Seq-id handles take as one (equal ones, and ones that differ only in what NCBI's index does not look at: the
  case of letters, GenBank, EMBL or DDBJ for an accession alone, PDB's `chain` and `chain-id`, a `giim` id's other
  database) and with NCBI's spelling of each (the first Seq-id of a handle, with the case of the letters of each);
  with `-subject` it is refused only where a subject record's spelling may show in the message. Values nested in
  themselves (Bioseq-sets, Seq-loc mixes, User-objects and the like) are read to any depth on every target, as LOSAT
  searched such PSSMs before (now with NCBI's errors inside them; nested Bioseq-sets are a query Seq-entry set, above):
  NCBI BLAST+ may crash on a nest of thousands of levels with its default 8 MiB stack (no output), where LOSAT writes
  what NCBI writes with a larger stack. The
  library API's PSSM query (`Inputs::pssm`) ends as the command line's run: a PSSM file that NCBI ends without a word
  gives a call without an error, and an error after the object manager's `Failed to parse location of` message has
  the command line's text (its line feed).
  Two library calls in one process (one after the other, at once on several threads, or two `run_batches`
  iterators advanced in turn) each keep their own diagnostic stream: a call whose stream the object manager fails
  no longer mutes another call's warnings and errors.

- Library internals (line LIB, P-V5b): the crate's private items document without rustdoc warnings (49 -> 0); the
  tabular id and title columns have one form, the one a search writes (the test-only readers of `seqid` are gone and
  the ledger test `tabular_id_and_title_fields_are_ncbis` runs `ParsedSubject` and the tabular writer); BLASTX's
  search no longer takes a `diagnostics` flag that only a test set; blastp's query-block search keeps declining
  several outputs it would open itself, now with a test of why (two outputs to one file, a file that cannot be
  created). The command line's output and help texts are unchanged.

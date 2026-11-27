- psiblast `-db` with parsed Seq-ids: where the PSSM (`-out_pssm`, `-out_ascii_pssm`, the archive
  of `-outfmt 11`, and each round of `-num_iterations` other than 1) reads an included hit by its
  Seq-align id after the round's report and NCBI BLAST+ 2.17.0 gives another database sequence for
  that id (a sequence loaded before whose definition line holds it, by the data loader's registry or
  the scope, or that the scope matches weakly: `ref|NP_100001|` for `ref|NP_100001.1|`), the PSSM now
  takes that sequence's letters, or, where they end before the hit's range does, writes NCBI's
  warning `Failed to retrieve sequence` and the query's letters for the hit, as NCBI BLAST+ does.
  LOSAT refused every such run. Still refused: such a read where LOSAT does not know which sequence
  NCBI BLAST+ gives (an order of its loads that their Seq-ids do not decide, the 1,000-entry registry
  at its limit), and one whose load of the other sequence can write texts of a damaged taxonomy
  database or the warning of a UTF8String header.

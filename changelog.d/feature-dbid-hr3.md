- psiblast `-db` with parsed Seq-ids: where the PSSM (`-out_pssm`, `-out_ascii_pssm`, the archive
  of `-outfmt 11`, and each round of `-num_iterations` other than 1) reads an included hit by its
  Seq-align id after the round's report and NCBI BLAST+ 2.17.0 gives another database sequence for
  that id (a sequence loaded before whose definition line holds it, by the data loader's registry or
  the scope, also by an earlier query or round of the run, or that the scope matches weakly:
  `ref|NP_100001.1|` for the Seq-align id `ref|NP_100001|`), the PSSM now takes that sequence's
  letters, or, where they end before the hit's range does, writes NCBI's warning `Failed to retrieve
  sequence` and the query's letters for the hit, as NCBI BLAST+ does. LOSAT refused every such run.
  Such a read loads the other sequence for the rest of the run, as NCBI BLAST+'s does. Still
  refused: such a read where LOSAT does not know which sequence NCBI BLAST+ gives (an order of its
  loads that their Seq-ids do not decide, the 1,000-entry registry at its limit), and one whose
  load of the other sequence can write texts of a damaged taxonomy database or the warning of a
  UTF8String header.
- blastn and tblastx `-db` with parsed Seq-ids: a minus-strand alignment's prefetch where the
  1,000-entry registry may or may not still hold another sequence for the subject's Seq-align id,
  and the lookup gives the subject, is no longer refused (either gives a sequence, as in NCBI
  BLAST+).
- psiblast `-db` with parsed Seq-ids: a PSSM read of a later query whose Seq-align id an earlier
  query's prefetch registered for two sequences, in an order that LOSAT does not follow, is now
  refused; LOSAT took the letters of the sequence that the lookup gives, where NCBI BLAST+ may give
  the other's (a later round of the earlier query finds the ids in its scope and leaves the
  registry as it was).

- `-db` searches where NCBI BLAST+ gives another database sequence (the holder) for a reported subject's Seq-id (its
  data loader or its scope holds a sequence loaded before under that id, or matches it weakly, or its prefetch resolves
  it to that sequence first) now write NCBI BLAST+'s record of the holder instead of LOSAT's refusal in the tabular
  formats (`-outfmt` 6, 7, 10 and 20), BLAST XML (5), BLAST XML2 and JSON (13 to 16) and the organism report (18) of
  BLASTP, psiblast, BLASTN (every task), BLASTX, TBLASTN and TBLASTX: the ids and title of the holder's definition line
  that NCBI BLAST+ loaded it by, its other definition lines, tax ids and length, and its letters, with the subject's
  alignment coordinates and scores. For example, after a first query's report that shows `gi|5001|ref|NP_100001.1|`, a
  second query's hit on another database's `ref|NP_100001.1|` is written with the first sequence's names, as NCBI BLAST+
  writes it. LOSAT still refuses, naming the holder, where it does not write that record as NCBI BLAST+ does: the
  pairwise reports (outfmt 0 to 4) and outfmt 17, a holder that the call does not report, a holder shorter than the
  subject's aligned range in BLASTX (or, in TBLASTN and TBLASTX, of another length), a holder of other letters in a
  TBLASTN table that prints identities without the aligned letters, a table whose only field of the subject is `slen`
  where a query batch's prefetch is off, a
  subject whose loads give it in one place and itself or a failing lookup in another, a positive id list, a header with a
  UTF8String, taxonomy database texts, a holder whose lines do not hold the subject's id or hold it only as an id
  of another type of the same accession, and an id that differs from the holder's in letter case only.

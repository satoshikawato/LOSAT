- `blastx -ungapped -window_size 0` with a low `-threshold` no longer crashes (a panic) where NCBI's
  one-hit word finder finds ungapped alignments of negative length, and gives NCBI BLAST+ 2.17.0's
  output there, with `-subject` and `-db`: such alignments take part in the linking of the others and
  are then deleted, as in NCBI BLAST+.
- The alignments of length 0 that such searches report are written as NCBI BLAST+ writes them, in
  every format that `blastx` supports: tabular rows of empty ranges with `-db`, NCBI's exception
  (exit status 255) after the rows before it where NCBI's formatters stop (outfmt 0, `qcovs`, the
  tabular rows with `-subject`), the cut document and the error message of BLAST XML, XML2 and JSON,
  and the Seq-aligns of outfmt 8, 9, 11 and 12. LOSAT no longer crashes when such an alignment
  touches a lower-case or SEG mask. With `-subject`, an alignment of length 0 at the subject's first
  letter is refused (an error, exit status 1): NCBI BLAST+'s Seq-loc mapper cuts it into two
  alignments without scores, which NCBI BLAST+ writes as rows of the query's Seq-align id, stops
  on with an exception or runs out of memory on, depending on the frame, the fields and its memory
  limit.
- A subject whose every reported alignment has length 0: with `-db`, NCBI BLAST+'s warning `Error
  pre-fetching sequence data` on standard error, as NCBI writes it (every format but 8, 9, 11 and
  12; 8, 9 and 12 too with `PRE_FETCH_SEQS_LIMIT=2147483647`; the archive never). With `-subject`,
  NCBI BLAST+'s Seq-loc mapper finds no sequence type in such a subject and makes alignments
  without scores of query nucleotides against the rest of the subject's letters: the tabular
  reports (outfmt 6, 7, 10 and 20) of a subject's only such alignment now show NCBI BLAST+'s row
  of it (its coordinates, letters and counts, an E-value of 0.0 and a bit score of -1.0; where
  its query row would start before the query's first nucleotide, the row of the query's Seq-align
  id or the exception that NCBI BLAST+ gives), where LOSAT wrote the alignment of length 0. LOSAT refuses (exit status 1) the other formats of such
  a subject, two or more such alignments (NCBI BLAST+ writes more rows, stops with an exception
  or runs out of memory), and a table that reads the letters of a query nucleotide outside its
  `-query_loc` interval.

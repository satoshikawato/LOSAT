- `-db` with databases built with `-parse_seqids`: the same record in two volumes or two databases with other tax
  ids, whose title NCBI can make from the tax ids (an empty title with `taxdb.bti` on the search path, a protein
  title ending in `]`), is no longer refused in outfmt 5, 6, 7, 10, 17 and 20 and the tabular field lists without
  tax id fields (all six programs): NCBI BLAST+ 2.17.0 writes those bytes whichever copy it shows. It is still
  refused in outfmt 0, which shows the generated title, in the tabular formats where the first word of that title
  replaces an `lcl|Subject_N` id (for a nucleotide, or a protein title ending in `]`), and in outfmt 13 to 16.
- `-db` with databases built with `-parse_seqids`: a nucleotide sequence with an `lcl|Subject_N` id, an empty title
  and a tax id is now refused in the tabular formats when `taxdb.bti` is on the search path (BLASTN, TBLASTN,
  TBLASTX). NCBI BLAST+ 2.17.0 replaces the id by the first word of the title it makes from the tax ids (`Homo`);
  LOSAT wrote `Subject_N`, with exit status 0.

- `-db` with databases built with `-parse_seqids`: the same record in two volumes or in two databases of an
  alias (or of `-db "a b"`) that differs only in its tax ids is no longer refused where no report writes tax
  ids or their names (outfmt 6, 7, 10 and 20 without `staxid`, `staxids` or a name field, and SAM; all six
  programs; outfmt 0 and 5 still refuse the two copies as following hits that NCBI shows as one, as for two
  copies with the same tax ids). NCBI BLAST+ 2.17.0 shows the copy that its lookup gives first for both, the same bytes; LOSAT
  refused the sequence. Outfmt 13 to 16 and the tabular tax id and name fields still refuse it (NCBI writes
  that copy's tax ids for both), as does every format where a title can depend on the tax ids (an empty
  title with `taxdb.bti` on the search path, a protein title ending in `]` that a report rewrites) or where
  the taxonomy database's messages depend on which copy NCBI fetches (the messages of opening it where only
  one copy has a tax id, the offset error of a tax id that a copy holds).
- `-db` with databases built with `-parse_seqids` and a taxonomy database whose lookup of a tax id writes the
  offset error: the same record with that tax id in two volumes is now refused. NCBI BLAST+ 2.17.0 writes the
  error once, for the volume of the copy that its lookup gives first; LOSAT wrote it once for each volume,
  with exit status 0. The messages of opening the taxonomy database, which come once per run, are unchanged.

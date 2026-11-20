- `-db` with databases built without `-parse_seqids`: a sequence with an empty title and a tax id
  (makeblastdb `-taxid`) that the taxonomy database names is now refused in outfmt 0 where the report
  shows it (all six programs). NCBI BLAST+ 2.17.0 makes its title from the tax id
  (`unnamed protein product [Homo sapiens]`, for a nucleotide `Homo sapiens`); LOSAT wrote the title
  without the organism, with exit status 0. The other formats and a tax id that the taxonomy database
  does not name are written as before.

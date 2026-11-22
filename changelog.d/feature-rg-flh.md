- `psiblast -in_pssm` and `tblastn -in_pssm`: when a PSSM query repeats a Bioseq id, the error message writes the
  repeated accession under the database name of the first accession of its prefix, digit count and version
  (for example `emb|AAA00001.1|` after `emb|AAA00002.1|`), as NCBI BLAST+ does, instead of under its own.

- `makeblastdb -oid_masks 1` writes the OID mask file of each volume (`.pxm`, `.nxm`) that marks the
  sequences whose ids are all predicted RefSeq accessions, as NCBI BLAST+ 2.17.0 does; LOSAT refused
  the option.

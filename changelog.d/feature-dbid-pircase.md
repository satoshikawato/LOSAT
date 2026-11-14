- `-db` with parsed Seq-ids in version 5 databases, in every search program: a `pir` or `prf` id, or a text id
  with a name and no accession, that another sequence holds in another letter case in the same alternative and
  release (`pir||PIRNAME` and `pir||pirname`, `gb||locusy` and `gb||LOCUSY` in two databases of a `-db` list or an
  alias) is now refused where a report reads the sequence, as in version 4. NCBI BLAST+ 2.17.0 takes the two ids
  for one and writes the names, and in outfmt 0 the letters, of the sequence it reads first for both; LOSAT wrote
  each sequence's own names (version 5 keys these ids in their own letter case, so LOSAT's lookup did not find
  the other). Ids of other alternatives or releases (`gb||locusy` and `emb||LOCUSY`, `sp||NAMEY` and `tr||namey`)
  stay apart, as NCBI's.
- Version 5 lookups of a Seq-align id open the `acc2oid` table once per search instead of once per reported
  sequence.

- LOSAT Web: an index scan keeps only the residue counts that occur in each record, not a table of 256 counts
  (2 KiB per record). `register` of 28,000 16S sequences takes 57 MiB less of the Data worker's linear memory, which does not
  shrink.

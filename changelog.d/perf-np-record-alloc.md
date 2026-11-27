- makeblastdb: builds databases of many records faster: about 1.2 times as fast on a proteome
  (1.1 times with `-mask_data`), 1.05 to 1.07 times on cDNA, and 1.3 to 1.4 times on files of
  50,000 primers or peptides (1.1 times with `-parse_seqids`); genomes and other inputs of few
  long records take the same time. The databases and messages do not change.

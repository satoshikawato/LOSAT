- makeblastdb: builds databases of many records faster, with the same memory: about 1.2 times as
  fast on a proteome (1.1 to 1.2 times with `-parse_seqids` or `-mask_data`), 1.04 to 1.09 times
  on cDNA, and 1.3 to 1.5 times on files of 50,000 primers or peptides (1.2 times with
  `-parse_seqids`); genomes and other inputs of few long records take the same time. The
  databases and messages do not change.

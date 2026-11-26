- makeblastdb: builds databases from FASTA about 1.3 to 1.6 times as fast on a genome, cDNA and a
  proteome (1.1 to 1.5 times with `-parse_seqids` or `-mask_data`), and stays faster on short-line,
  short-record, spaced and ambiguous inputs, with the same memory; segmasker and dustmasker are
  2 to 6 % faster on the same inputs. The databases and masks do not change.

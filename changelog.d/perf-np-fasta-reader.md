- makeblastdb: FASTA input is read about 1.4 to 1.7 times as fast with the same memory (the 34 MB A. thaliana
  proteome: 0.16 s to 0.11 s; 10 MB of cDNA: 33 ms to 20 ms; the 146 Mb D. melanogaster genome: 0.65 s to 0.47 s);
  segmasker and dustmasker, which read FASTA the same way, are 2.5-6.5 % faster. The reader looks for line ends,
  residues and the end of a title 32 or 64 bytes at a time instead of one byte at a time. The databases and masks do not
  change.

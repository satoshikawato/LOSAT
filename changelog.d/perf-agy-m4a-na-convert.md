- makeblastdb: a nucleotide database with ambiguous bases (a genome FASTA with runs of N, for example) builds about
  1.3 times as fast with the same memory (a 146 Mb genome: 0.65 s to 0.49 s): the bases of such a sequence are converted
  64 at a time instead of a table lookup and a branch per base. The databases do not change.

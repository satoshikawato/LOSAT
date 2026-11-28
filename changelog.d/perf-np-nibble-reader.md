- makeblastdb: building a nucleotide database from FASTA takes less memory and time on genomes and other long
  sequences, because the reader keeps each sequence's bases at half a byte each instead of a whole byte. The databases do
  not change, and protein input, dustmasker, windowmasker, segmasker, makembindex and convert2blastmask are unchanged.

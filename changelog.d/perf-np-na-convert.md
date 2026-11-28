- makeblastdb: a nucleotide database builds faster when its sequences hold bases other than A, C, G and T (a genome with
  runs of N, for example) or gaps; the databases do not change. The Drosophila melanogaster genome (BDGP6, 144 Mbp) takes
  0.70 of the time at one thread (421-430 ms to 294-298 ms), with the same memory.

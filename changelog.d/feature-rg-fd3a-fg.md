- psiblast and tblastn `-in_pssm` guess the file's format as NCBI BLAST+ 2.17.0 does, with every test NCBI runs before
  it decides on binary ASN.1: a file that NCBI takes for another format (a VCF or GenBank line, a NEXUS tree, a
  compressed or SRA file) gets "Unsupported format for PSSM", a file with a byte that is not text (a byte order mark,
  UTF-16, a NUL) gets NCBI's "Unrecognized format for PSSM" or its silent end, a PSSM after more than 8 KB of
  comments or spaces is searched, and a Glimmer3 line whose frame passes the range of `int` gets NCBI's C++ exception
  report. Only a binary ASN.1 PSSM, an XML PSSM and the few files whose format NCBI decides by reading memory past
  them are refused. makeblastdb's input and mask-file checks use the same port.

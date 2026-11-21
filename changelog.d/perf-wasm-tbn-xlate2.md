- Web (WebAssembly): TBLASTN against a nucleotide database takes less time in Node, with the same
  reports. At one thread, one human protein against the human genome (GRCh38) runs in about 19 %
  less time (19.4 s to 15.7 s) and against 31,168 fly cDNAs in about 10 % less; searches of many
  proteins, whose time goes to the scan and the extensions, change by about 1 %. x86-64 CPUs
  without AVX2 run the same faster translation, with 1 to 3 % fewer instructions per search;
  native builds with AVX2 run as before.

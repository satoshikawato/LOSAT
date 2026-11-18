- Web (WebAssembly): TBLASTN against a nucleotide database takes less time in Node, with the same
  reports. At one thread, one human protein against three bacterial genomes runs in about 19 % less
  time and against 31,168 fly cDNAs in about 9 % less; searches of many proteins, whose time goes to
  the scan and the extensions, change by about 1 %. The six frames of the database are read with
  about a third of the instructions (the 128-bit loop starts each phase on a whole byte and makes
  its constants once). x86-64 CPUs without AVX2 run the same loop, with 1 to 3 % fewer instructions
  per search; native builds with AVX2 run as before.

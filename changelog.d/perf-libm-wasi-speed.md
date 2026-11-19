- WebAssembly: the port of glibc's `log` computes the same bits with less work (the split of each
  table value is precomputed, and one check of the last rounding replaces two), about 40 % less time
  per call. blastp, tblastn and blastx with composition-based statistics (the default), which call
  it some 12 to 20 million times in a small search, took 6 to 9 % less time in the serial WASI
  command; the output does not change.

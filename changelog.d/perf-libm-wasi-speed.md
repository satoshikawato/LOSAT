- WebAssembly: LOSAT's `log`, the code of Arm Optimized Routines (MIT) built to give glibc 2.35's
  bits, computes the same bits with less work (the split of each table value is precomputed, and
  one check of the last rounding replaces two), about 40 % less time per call. blastp, tblastn and
  blastx with composition-based statistics (the default), which call it some 12 to 20 million times
  in a small search, took about 5 % less time in the serial WASI command (1 to 9 % across cases and
  runs); the output does not change.

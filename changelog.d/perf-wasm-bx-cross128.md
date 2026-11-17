- Web (WebAssembly): BLASTX runs faster on queries with many lookup hits per cell, with the same
  reports (a 50-kb piece of the fly genome against the fly proteome at one thread: NUMBERS_TBD of the
  time before). Its two-hit seeds are now tested in batches with simd128, as the native build does
  with AVX2 (`LOSAT_X_BXCROSS128`, on by default; also used by x86-64 CPUs without AVX2, not timed
  there). Native builds with AVX2 give the same reports and run the same instructions, one function
  laid out otherwise (cycles NUMBERS_TBD).

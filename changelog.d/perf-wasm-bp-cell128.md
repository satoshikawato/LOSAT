- Web (WebAssembly): BLASTP runs faster on query sets with many lookup hits per cell, with the same
  reports (1,000 proteins against a 6-species proteome at one thread: about 18 % fewer instructions;
  100 proteins 16 %, 10 proteins 11 %). Its two-hit seeds are now tested in batches of eight with
  simd128, as the native build does with AVX2 (`LOSAT_X_BPCELL128`, on by default; also used by
  x86-64 CPUs without AVX2, not timed there). Short single queries run as before. Native builds with
  AVX2 give the same reports with the same machine code.

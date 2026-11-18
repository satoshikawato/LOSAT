- Web (WebAssembly): BLASTP takes less time in Node, with the same reports. At one thread against a
  6-species proteome, 1,000 query proteins run in about 12 % less time, 10 proteins about 7 %, and one
  protein of 880 to 1,800 residues about 3 to 5 %. Its two-hit seeds are now tested in batches of eight
  with simd128, as the native build does with AVX2 (`LOSAT_X_BPCELL128`, on by default; x86-64 CPUs
  without AVX2 use it too, not timed there). Native builds with AVX2 run as before.

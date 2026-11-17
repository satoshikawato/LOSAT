- Web (WebAssembly): BLASTX runs about a quarter faster on queries with many lookup hits per cell (a
  50-kb piece of the fly genome against the fly proteome at one thread: x0.75 of the time before), as
  it now tests its two-hit seeds in batches across the lookup cells with simd128, as the native build
  does with AVX2 (`LOSAT_X_BXCROSS128`, on by default; also used by x86-64 CPUs without AVX2). The
  reports do not change, and neither do the native AVX2 builds.

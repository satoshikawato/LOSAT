- Speed: `LOSAT_WASI_THREADS_DEBUG` is read once per process instead of at every stage report. BLASTN
  reported a stage once per subject chunk, so each database sequence paid a read of the environment
  (about 340 instructions); a megablast search of 20,000 short subjects runs about 3 % fewer
  instructions. Output is unchanged.

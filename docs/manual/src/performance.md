# Threads, speed and memory

## Threads

Every search program takes `-num_threads N` (default 1). LOSAT uses the threads for both kinds of
search:

- **`-db`**: the database is searched in parallel, as in NCBI BLAST+.
- **`-subject`**: LOSAT also searches FASTA subjects with `N` threads. NCBI BLAST+ runs a
  `-subject` search on one thread and prints
  `'num_threads' is currently ignored when 'subject' is specified.`, so with LOSAT a pairwise
  comparison of two genomes uses every thread you give it.

The output does not depend on the thread count: the same command writes the same bytes with any
`-num_threads`, and on the release test set those are the bytes of NCBI BLAST+ 2.17.0 run with one
thread (the release tests compare 1, 2, 4, 8 and 16 threads). You can raise `-num_threads` without
checking the results again.

`-mt_mode` is accepted for compatibility, and the output does not depend on it either.
`blastn`, `blastp`, `blastx` and `tblastn` take 0 (auto), 1 (split by queries) or 2 (split by
database). `rpsblast` and `rpstblastn` take 0 (split by database volumes) or 1 (split by queries).
`tblastx` and `psiblast` have no `-mt_mode`, as in NCBI BLAST+ 2.17.0.

LOSAT does not lower `-num_threads` to the number of CPUs, and it prints no warning about it.
NCBI BLAST+ does both. Asking for more threads than the machine has cores rarely helps; a value
up to the number of physical cores is a good start.

## Memory

Memory use grows with the size of the query batch, the database volume being searched and the
number of threads. If a large search runs out of memory:

- lower `-num_threads`: memory use rises with the thread count;
- split a very large query file into several runs.

When an allocation fails, LOSAT stops at once. NCBI BLAST+ prints `BLAST ran out of memory`
and exits with code 4 (see [Differences from NCBI BLAST+](differences.md)).

## WebAssembly

The WASI builds (`losat.wasm`) give the same output as the native program but run slower, because
WebAssembly's vector instructions are narrower than those of current CPUs. The threaded WASI build
(`wasm32-wasip1-threads`) takes `-num_threads` like the native program. See [Platforms and WebAssembly](platforms.md).

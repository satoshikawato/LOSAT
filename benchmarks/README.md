# LOSAT benchmarks

Each directory holds one set of measurements of LOSAT against NCBI BLAST+ 2.17.0: the inputs, the
measured values and the scripts that make its tables and figures.

| Set | Measured | Machine |
| --- | --- | --- |
| [v0.2.0](v0.2.0/README.md) | Twelve search modes of `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast` and `rpstblastn` at 1 and 16 threads; `makeblastdb` and `rpsbproc`; scaling with threads (1 to 16), query count and database size (up to NCBI's `core_nt`, 319.4 GB); the BLAST steps of published tools; the WebAssembly build. LOSAT `93ca2f652` (v0.2.0 release candidate) and development builds, named per figure. | Two Intel Xeon Gold 6338 CPUs (64 cores), Linux |
| [db6-small](db6-small/README.md) | `-db` searches (megablast, blastn, blastp, blastx) on small public databases at 1, 8 and 16 threads, and database construction, during v0.2.0 development | Two Intel Xeon Gold 6338 CPUs |
| [v0.1.0](v0.1.0/README.md) | LOSAT v0.1.0: output comparison and execution time of `blastn`, `blastp` and `tblastx` with `-subject`, native at 1 and 8 threads and WebAssembly | Intel Core i9-14900HX, Ubuntu 24.04 under WSL2 |

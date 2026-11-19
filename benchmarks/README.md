# LOSAT Benchmarks

Each release benchmark suite is stored in a versioned directory with frozen input metadata, raw measurements, deterministic snapshot tables, and publication figure rendering scripts.

## Suites

| Version | Status | Scope | Notes |
| :--- | :--- | :--- | :--- |
| **[v0.2.0](v0.2.0/README.md)** | **Current release** | All 12 search modes (`blastn`/`blastp`/`blastx`/`tblastn`/`tblastx`/`psiblast`/`rpsblast`/`rpstblastn`), utilities (`makeblastdb`, `rpsbproc`), thread scaling (`1`–`16`), large DB scaling (`340` genes to `core_nt` 298 GiB), downstream bioinformatics tools, and WebAssembly | Measured on dual Intel Xeon Gold 6338 (`ws`) against NCBI BLAST+ 2.17.0 & `rpsbproc` 0.5.1 |
| [v0.1.0](v0.1.0/README.md) | Historical archive | Initial pairwise `-subject` benchmarks (`blastn`, `blastp`, `tblastx`) | Intel Core i9-14900HX (WSL2) |

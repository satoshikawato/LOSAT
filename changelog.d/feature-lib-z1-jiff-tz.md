- makeblastdb: on Linux, macOS and WASI the local time of the creation date and the log is read
  with the crate jiff (MIT) from the zone that glibc finds for NCBI's Linux build (`TZ`, `TZDIR`,
  `/etc/localtime`), in place of LOSAT's own reader (line LIB; Owner 2026-11-20 「TZ=A」).
  Zone files with leap seconds now give glibc's time; a `TZ` that glibc reads outside POSIX's
  forms gives UTC (docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md). Windows is unchanged.

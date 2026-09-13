# Invalid isolated build

The initial isolated source copy omitted unchanged, HEAD-identical `LOSAT/build.rs`. Consequently `losat_wasi_threads` was not set, and n8 exited explicitly as unsupported. This is an experiment setup failure, not a candidate trap or speed regression. All rows in this run are excluded from adoption, including the n1 rows that completed.

The original threaded artifact is preserved by SHA-256 under `../invalid-variant-artifacts/`. The complete source is `../baseline-source-v2.tar.gz` with `../snapshot-v2.json`; subsequent isolated builds include `build.rs`. The original baseline binaries were built in the repository, include the build script, and are unaffected.

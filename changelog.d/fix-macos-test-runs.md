- Test-only: the two macOS `cargo test` jobs of the platform checks no longer fail for reasons of the
  test environment. The blastp query-block test that links a side file to `/dev/full` runs on Linux only
  (macOS has no `/dev/full`); the blastn test that a stalled stage lets the helpers leave tries up to 40
  stages or 20 s instead of five, and skips a pool larger than the CPU count; the outfmt 11 archive
  checks name the number of differing bytes, the first offset and the bytes around it when a case
  differs from the frozen NCBI bytes. No program output changes.

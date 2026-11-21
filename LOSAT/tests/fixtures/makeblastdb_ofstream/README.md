# Observation tables of the C++ output file stream (O1a, line LIB-O)

What they are: the behaviour of a `std::ofstream` (as makeblastdb 2.17.0's C++ runtime has it) under scripted
operations, with every system call made on the file and the stream state after each operation. They are the
specification that `LOSAT/src/writedb/ofstream.rs` is tested against
(`src/writedb/ofstream_fixture_tests.rs` replays every row).

## Origin of the rows (honest provenance)
- Program: `harness/ofprobe.cpp` (our own) built with g++ 13.3.0 (Ubuntu 24.04), dynamic libstdc++.so.6 from the same
  g++ 13 series; run on the local development host (WSL2), 2026-11-20 (rows_boundaries, rows_failures, rows_states) and
  2026-11-21 (rows_random, rows_close). Command: `LD_PRELOAD=shim.so SHIM_PATH=<out> SHIM_LOG=<log> ./ofprobe <script> <out>`,
  scripts made by `harness/gen2.py` (boundaries), `gen3.py` (failures), `genstates.py` (states), `genrand.py 11 800`
  (rows_random), `seq_close.txt` = the 7 hand-written sequences then `genrand.py 21 800` (rows_close, one run), then `parse.py`.
- Shim (`harness/shim.c`): logs `write`, `writev`, `lseek`, `close`/`fclose` of the output file; fails `write` past a byte
  limit (mode 0: the call reaching the limit is short, later calls fail with ENOSPC; mode 1: a call that would pass the limit fails
  and writes nothing); makes `close(2)` fail with a chosen errno (`K errno` rows; the descriptor is released either way).
- NOT regenerated on ws: ws's libstdc++ is the GCC 11.4 one and behaves differently (it writes a 4095-byte block at once), so the
  ws tables are not these tables. makeblastdb 2.17.0 (official Linux x64 release) carries GCC 13.2.0 (`strings`), the 13 series.
- Link to the NCBI binary: end-to-end runs of the official makeblastdb 2.17.0 and `losat makeblastdb` under the same shim, each
  output suffix failing at byte limits around the file sizes and the 8191-byte boundaries (and, for close, at `close(2)`):
  1,592 write-failure cases, 1,334 byte-identical (exit status, messages, calls, files); the other 258 are not stream behaviour
  and are the same before this rewrite (90 `.pjs`/`.njs` message wordings, 16 `.pal` message prefixes, 152 where NCBI re-opens
  a truncated `.pin`/`.nin` and exits -11 or 2 while LOSAT goes on; followed up separately). 198 cases with `close(2)` failing
  (EIO, ENOSPC, EDQUOT) on every output suffix of protein, nucleotide, masked and large databases: exit status, messages and
  files identical to the run without the failure, so a failed close is not visible outside the stream. Run on ws
  (`kawato-ProLiant-DL360-Gen10-Plus`, 2026-11-20/21).

## Files
`rows_boundaries.txt` (held bytes b x write size n, 1274 sequences), `rows_failures.txt` (108), `rows_states.txt` (every state bit
set before every operation, 416), `rows_random.txt` (800, with close failures), `rows_close.txt` (807: close failing with EIO,
ENOSPC, EDQUOT, with and without held bytes, failed writes, repeated close, seekp, setstate). Format: header of
`ofstream_fixture_tests.rs`.

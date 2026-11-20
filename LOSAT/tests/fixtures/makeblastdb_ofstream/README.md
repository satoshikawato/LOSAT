# Observation tables of the C++ output file stream (O1a, line LIB-O)

What they are: the behaviour of a `std::ofstream` (as makeblastdb 2.17.0's C++ runtime has it) under scripted
operations, with every system call made on the file and the stream state after each operation. They are the
specification that `LOSAT/src/writedb/ofstream.rs` is tested against
(`src/writedb/ofstream_fixture_tests.rs` replays every row).

How they were made (own programs, in `harness/`): `ofprobe.cpp` runs the operations and prints `rdstate()`;
`shim.c`, preloaded, logs `write`/`writev`/`lseek`/`close`/`fopen` of the output file and fails `write` past a byte
limit (mode 0: the call reaching the limit is short, later calls fail with ENOSPC; mode 1: a call that would pass
the limit fails and writes nothing); `parse.py` prints one line per operation.
Binary: `ofprobe` built with g++ 13.3.0 (Ubuntu 24.04), dynamic libstdc++.so.6, on the local development host,
2026-11-20. makeblastdb 2.17.0 (official Linux x64) is statically linked with GCC 13.2.0's libstdc++ (`strings`);
under the same shim it writes every database file as plain `write`/`close` calls of such a stream
(e.g. `.pjs`: `w 481`, `w 1`; `.pos`: `w 32`, `w 28`, `l 8`, `w 24`). The tables are not yet regenerated on ws.

Files: `rows_boundaries.txt` (held bytes b x write size n, 1274 sequences: `gen2.py`), `rows_failures.txt`
(108: `gen3.py`), `rows_states.txt` (every state bit set before every operation: `genstates.py`),
`rows_random.txt` (800 sequences: `genrand.py 11 800`). Format: see the header of the test file.

# rpsbproc growth observations (line LIB-O, item O1c, 2026-11-20)

What the string and the vector of NCBI's rpsbproc 0.5.1 motif check do as they grow, observed from real
objects. `src/rpsbproc/prosite_growth.rs` is written from these tables (not from the C++ library's source);
`prosite_growth_tests.rs` and `sites_tests.rs` check every row.

## Binary and compiler

- NCBI rpsbproc 0.5.1 Linux x64 (`~/losat-bench/tools/env_rpsbproc-0.5.1/bin/rpsbproc` on ws): libstdc++ linked
  statically; its `.comment` reads `GCC: (conda-forge gcc 13.4.0-16) 13.4.0`. The existing frozen cases
  (`../rpsbproc/oracle/rd13..rd18`) show a string `max_size()` of 2^63 - 1: the build is C++20.
- The observation programs are built on ws with conda-forge `gxx_linux-64=13.4.0` (gxx_impl 13.4.0 build 20,
  libstdcxx-devel 13.4.0 build 120), `-O2 -std=c++20 -static-libstdc++` (`max_size()`: 2^62 - 1 with `-std=c++17`,
  2^63 - 1 with `-std=c++20`; only C++20 matches the binary).

## Files

- `string_ops.tsv`, `string_ops_large.tsv`, `string_ops_host.tsv`: `programs/strprobe.cpp` (list of operations
  `R n` reserve, `A n` append(n, c), `P` push_back, `C` clear on one `std::string`; per operation `len/capacity/peak`,
  or `!what|len/capacity/peak` for an exception; `peak` = the most heap bytes live at once during the operation,
  counted by an allocator on a `basic_string` that runs in step with a plain `std::string`, whose size and capacity
  must equal). Operation lists: `programs/gen_strings.py` (seed 20261120). `string_ops_host.tsv` holds the rows
  that hold 2^30 bytes or more at once (run on ws, 500 GB of memory): the model refuses them as host dependent.
- `vector_ops.tsv`: `programs/vecprobe.cpp` (a `std::vector<size_t>` and a stack of copies of it: `G n` clear and
  push n, `K` clear, `Z` pop_back, `S` push a copy, `O` take the top back (`v = top; pop`), `N n` replace the vector by
  one of n elements; per operation `capacity/size/stack entries/stack capacity sum/peak records`), lists from
  `programs/gen_vectors.py`.
- `binary_allocs.tsv`: the allocations of NCBI's rpsbproc binary itself for edited motifs (`programs/runmotif.sh`,
  the `LD_PRELOAD` shim `programs/mallog.c`, `programs/runall_motifs.sh`, `programs/make_binary_allocs.py`, motifs in
  `programs/motifs.tsv`), ties the tables to the binary.

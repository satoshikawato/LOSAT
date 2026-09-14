# Independent C-NA64 build and Wasm identity audit

The read-only ncbi_parity_auditor supports the build-validation and exact-Wasm
reuse boundary. All four planned stages exit0: fmt,621 tests passed /3 ignored /
0 failed, Clippy with -D warnings, and threaded command/reactor build.

The two new Wasm binaries and four JS files are byte-identical to the old bundle,
including the historical integrated-threaded-default metadata bindings:

- command: `5b4f42831efa61ca73ead78f7675e53e5384ce713b9fa5ed5d481aeb3c052295`
- reactor: `d58a40248c1a3de5d33b5ed1e30cc6ea9ee5926cd37a9c385fcc16e4dc39762f`

Actual Wasm fingerprints contain only `["-C","target-feature=+simd128"]`. Bundle
metadata differs only in argv/cwd, not target, imports, exports, shared memory,
features or Rust/Node versions. Native x86_64 alignment flags do not enter Wasm.
The qualified native release binary remains `a1ab2845...`.

The auditor verified all589 source files at each of three locations,87 fixed
launch files and14 driver hashes. Monotonic endpoints show79.528 seconds between
the last native N/P measurement and build start, and99.677 seconds between build
completion and the first remaining-group process. Heavy work did not overlap.

The launch binds nine groups, including the new protected LC738874/LC738875
threshold check at E-values10/100/10000 with the new native n1 and threaded n4
artifacts: three fresh oracles and six comparisons are planned, not yet claimed
as results by this audit.

Only prior Wasm observations with identical input, option and runner conditions
may be reused. This does not reuse old native qualification, qualify unmeasured
X allocator live bytes, or establish frozen release certification. Original NP
reuse and fixed16-call FAIL records and hashes remain unchanged. Remaining gates
and capacity must finish before final adoption.

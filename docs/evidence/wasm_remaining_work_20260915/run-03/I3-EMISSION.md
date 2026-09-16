# I3 serial emitted-code check

Actual serial command artifact: `3da152983cb98317d176cfa2759c74bd7a6b496dac0498758a1a09f0f6450094`.
The name section and full disassembly contain no scan helper function. The
actual `link_hsp_group_ncbi` function 1054 contains the indexed scan and no
call to `scan_large_gap_predecessors`. Its saved body is 16,961 WAT lines,
133 locals. Size itself is not a performance result.

See `work/I3/serial-emission/function-indices.json`, `function-1054.wat`,
and `emitted.wat`. The 802-case serial visit/state/trace differential is
byte-exact. Fresh-oracle raw checks and the new fixed performance set are
separate requirements; inlining alone does not prove non-regression.

## Threaded code

The actual I3 threaded artifact is
`95e58d009c060c196065bede76bbdaed998e1bd491f6fc4b615535c1bfb2ea22`.
Its full group function 1340 and scan helper 1368 WAT are byte-identical to
I2b. Their hashes are in `I3-threaded-code-equivalence.json`. The helper call
remains in the original group location. This is code equivalence for those
two functions, not complete-artifact identity or a new performance result.

# I2b: corrected target selector for I2

The I2 scan body and prefix rationale in I2-POLICY.md are unchanged. An
independent audit found that Rust 1.92 does not expose the `atomics` cfg for
this target, so I2 did not enable the intended cursor. I2b uses the existing
`all(losat_wasi_threads, feature="wasm-threads")` boundary from build.rs and
utils/threading.rs. I2 sources and artifacts remain preserved and unmeasured.

I2b has a separate source manifest and build directory. The corrected
standalone differential requires the frozen source hash, explicitly mirrors
build.rs's cfg, and retains the real threaded ABI. It passed all 800 generated
cases plus two empty-array boundaries on Native, serial and threaded Wasm,
trace off/on, with byte-exact complete visit/state and trace transcripts.

After emitted-code verification and the six-case threaded raw/thread gate,
the first performance condition is Mje n8 compiled-module reuse: two fixed
AB/BA sessions, one warmup plus three samples per version, ordinary Node,
exclusive execution, full lifecycle, unchanged 5%/50 ms and memory guards.
Stop on a failed session without extending or replacing samples. Only a
passing candidate proceeds to the remaining acceptance work.

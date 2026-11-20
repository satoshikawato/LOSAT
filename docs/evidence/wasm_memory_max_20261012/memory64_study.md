# Memory64 (64-bit WebAssembly) feasibility for LOSAT Web — source-backed study

Date: 2026-10-12. Scope: read-only research, no builds, no installs. All web sources were fetched in this session (2026-10-12) unless a row says otherwise.
Marks: **V** = read by me in a primary source (or run locally) in this session; **S** = secondary source; **I** = my inference, not verified.
Fetched copies of the sources: `$TASK_DIR/tmp/m64src/`; probe scripts: `.../tmp/m64tools/probe1.js`, `probe2.js`.

## (a) Verdict

1. **Not for now** as a LOSAT build: the threaded build rests on `wasm32-wasip1-threads`, and there is no wasm64 WASI target in rustc or wasi-libc, with upstream stating "no plans"; the only Rust wasm64 target is tier 3, nightly + `build-std`, without thread spawn.
2. Platform side is two-thirds ready: Chrome 133+, Firefox 134+ and Node 24+ ship Memory64 (16 GiB cap, shared memories work), but Safari 27.0 (stable, 2026-09-14) does not; it is only in Safari Technology Preview 251. A wasm64 build would be a third build next to the two wasm32 ones, not a replacement.
3. Cost besides the port: explicit bounds checks (Firefox/Safari always; V8 a cheaper constant compare), reported as 10% to over 100% slower. Going past 4 GiB is cheaper on wasm32 by keeping the database outside linear memory (route C). Revisit when Safari stable ships Memory64 and a Rust wasm64 route with threads exists upstream.

The earlier answer from memory ("browsers ship Memory64; the Rust threaded WASI target for wasm64 is missing") is half right: Chrome and Firefox ship it, Safari stable does not; the Rust statement is right and stronger than stated (no wasm64 WASI target of any kind, and none planned).

## (b) Findings per question

### Q1. Standard status

| # | Claim | Source | Date / version | Quote | Mark |
|---|---|---|---|---|---|
| 1.1 | Memory64 is a finished proposal, part of Wasm 3.0 | https://github.com/WebAssembly/proposals/blob/main/finished-proposals.md | main, WG vote 2025-07-23 | "\| [Memory64][memory64] \| Sam Clegg \| [WG 2025-07-23][WG-2025-07-23] \| core \| 3.0 \|" | V |
| 1.2 | Wasm 3.0 is the live standard and has 64-bit address space | https://webassembly.org/news/2025-09-17-wasm-3.0/ | 2025-09-17 | "Memories and tables can now be declared to use i64 as their address type instead of just i32." | V |
| 1.3 | The web limits a 64-bit memory to 16 GiB | same | 2025-09-17 | "on the web, a 64-bit memory is limited to 16 gigabytes" | V |
| 1.4 | JS API: `address: "i64"` in the memory descriptor | https://webassembly.github.io/spec/js-api/ | Editor's Draft 9 October 2026 | "enum AddressType { "i32", "i64", };" and "dictionary MemoryDescriptor { required AddressValue initial; AddressValue maximum; AddressType address; };" | V |
| 1.5 | Addresses and page counts of an i64 memory are BigInt in JS | same | same | "If addrtype is i64, Let n be ? ToBigInt(v). If n < 0 or n > 2^64 − 1, throw a TypeError." | V |
| 1.6 | Declared limit (validation) | same, §8 | same | "The maximum min or max field of a 64-bit memory is 2^37-1 pages (2^53 - 2^16 bytes)." | V |
| 1.7 | Runtime limit: 262,144 pages = 16 GiB | same, §8 | same | "The maximum size of a 64-bit memory is 262,144 pages (16 GiB)." | V |
| 1.8 | Threads (shared memory) is still phase 4, not in 3.0; shared memory64 is a combination of 3.0 and a phase-4 proposal | https://github.com/WebAssembly/proposals/blob/main/README.md | main | "### Phase 4 - Standardize the Feature (WG) ... \| [Threads][threads] \| Conrad Watt \|" | V |
| 1.9 | V8 14.6 is stricter than 1.6: a declared maximum above 262,144 pages is rejected (constructor and module validation) | local probe (below) | Node 26.8.2 / V8 14.6.202.34 | "Property 'maximum': value 262145 is above the upper bound 262144" | V |

### Q2. Engine support and limits

| # | Claim | Source | Date / version | Quote | Mark |
|---|---|---|---|---|---|
| 2.1 | Chrome: enabled by default in 133, desktop, Android and WebView | https://chromestatus.com/feature/5070065734516736 (read through `/api/v0/features?q=Memory64`) | entry updated 2025-08-22 | `"status": {"text": "Enabled by default", ... "milestone_str": "133"}, ... "desktop": 133, "android": 133, "webview": 133` | V |
| 2.2 | V8 maximum: 16 GiB on 64-bit, 2 GiB − 64 KiB on 32-bit builds | https://github.com/v8/v8/blob/main/src/wasm/wasm-limits.h | main, file last changed 2026-05-12 | "constexpr size_t kV8MaxWasmMemory64Pages = kSystemPointerSize == 4 ? 32'767 // = 2 GiB - 64Kib : 262'144; // = 16 GiB" | V |
| 2.3 | V8: a shared memory64 works, grows past 4 GiB up to 16 GiB | local probe (below) | Node 26.8.2 / V8 14.6 | see probe output | V |
| 2.4 | Firefox: shipped in 134 | https://spidermonkey.dev/blog/2025/01/15/is-memory64-actually-worth-using.html | 2025-01-15 | "the Memory64 proposal for WebAssembly has finally been released in both Firefox 134 and Chrome 133" | V |
| 2.5 | Firefox maximum: 16 GB (ArrayBuffer byte-length limit) | https://github.com/mozilla-firefox/firefox/blob/main/js/src/vm/ArrayBufferObject.h and `js/src/wasm/WasmMemory.cpp` (`MaxMemoryPages`) | main (WasmMemory.cpp last changed 2026-08-11) | "static constexpr size_t ByteLengthLimit = size_t(16) * 1024 * 1024 * 1024;  // 16 GB." | V |
| 2.6 | Firefox: shared i64 memories and atomics on them are in the engine test suite (source read, not run) | https://github.com/mozilla-firefox/firefox/blob/main/js/src/jit-test/tests/wasm/memory64/basic.js | main | "(memory i64 ${initial} ${max ?? ''} ${shared ? `shared` : ''}))" | V |
| 2.7 | Safari: added in Technology Preview 251 only | https://webkit.org/blog/18194/release-notes-for-safari-technology-preview-251/ | 2026-08-26 | "Added support for WebAssembly Memory64. (318257@main)" | V |
| 2.8 | Safari 27.0 stable does not list Memory64 (its only new WebAssembly feature is JSPI); 27.2 beta notes do not mention it either | https://developer.apple.com/documentation/safari-release-notes/safari-27-release-notes | "Released September 14, 2026 — 27.0 (20625.1.29)" | "Added support for WebAssembly JavaScript Promise Integration (JSPI)." (no occurrence of "Memory64" in the page data) | V (absence) |
| 2.9 | WebKit trunk has it on by default | https://github.com/WebKit/WebKit/blob/main/Source/WTF/Scripts/Preferences/UnifiedWebPreferences.yaml | main, 2026-10-10 | "WasmMemory64Enabled: type: bool status: stable ... jscOptionName: useWasmMemory64 ... defaultValue: true" | V |
| 2.10 | WebKit tracking bug still open | https://bugs.webkit.org/show_bug.cgi?id=300538 | reported 2025-10-10, modified 2026-06-05 | "NEW 300538 Add support for Memory64 in WebAssembly" | V |
| 2.11 | JSC: shared memory64 is tested; the ceiling is between 4 and 16 GiB and depends on the platform's address-space budget (relevant to iOS) | https://github.com/WebKit/WebKit/blob/main/JSTests/wasm/stress/memory64-maximum-limits.js ; `Source/JavaScriptCore/runtime/BufferMemoryHandle.h` | main | "What its buffer reports is bounded twice over: by the largest representable ArrayBuffer byte length, and by the address space one reservation may claim on this platform"; "maxGrowableBufferReservationBytes = Gigacage::primitiveAddressSpaceBudget / 4" | V |
| 2.12 | Compat tables agree: Chrome 133, Firefox 134, Safari preview only, mobile mirrors desktop | https://github.com/mdn/browser-compat-data/blob/main/webassembly/memory64.json ; https://github.com/WebAssembly/website/blob/main/features.json | BCD 2026-09-01; features.json 2026-09-11 | BCD: `"safari": {"version_added": "preview"}`, `"safari_ios": "mirror"`; features.json Safari: `["flag", "Supported in Safari Technical Preview 251"]` | S |
| 2.13 | Node: no flag from 24.0.0 (V8 13.6). CI Node 24.21.0 is therefore enough | https://github.com/nodejs/node/blob/main/doc/changelogs/CHANGELOG_V24.md | 2025-05-06, 24.0.0 | "The V8 engine is updated to version 13.6, which includes several new JavaScript features: ... [WebAssembly Memory64]" | V |
| 2.14 | iOS Safari: no primary statement found; BCD only mirrors desktop | — | — | — | I |

### Q3. Performance

| # | Claim | Source | Date / version | Quote | Mark |
|---|---|---|---|---|---|
| 3.1 | Slowdown range reported by SpiderMonkey | https://spidermonkey.dev/blog/2025/01/15/is-memory64-actually-worth-using.html | 2025-01-15 (21 months old) | "This performance penalty depends on the workload, but it can range from just 10% to over 100%—a 2x slowdown just from changing your pointer size." | V |
| 3.2 | Cause: no 4 GiB guard-region trick | same | same | "The size of the WebAssembly address space is the same as the size of the host address space. Therefore, we must pay the cost of bounds checks on every access, and as a result, Memory64 is slower." | V |
| 3.3 | Advice | same | same | "The only reason to use Memory64 is if you actually need more than 4GB of memory." | V |
| 3.4 | Firefox still always uses explicit checks | https://github.com/mozilla-firefox/firefox/blob/main/js/src/wasm/WasmMemory.cpp | main, 2026-08-11 | "On all 64-bit platforms, we may use explicit bounds checking or huge memories for memory32, but must always use explicit bounds checking for memory64." | V |
| 3.5 | V8 reduces the cost: it reserves 16 GiB per memory64 and uses the trap handler | https://github.com/v8/v8/commit/cbe471028e (and 8b4bcf30bd, 2023-11-21) | 2025-01-17 | "when trap handling is enabled we always reserve kMaxMemory64Size (=16GB) for all 64-bit memories and when such an allocation fails, we fail without retrying with smaller allocations." | V |
| 3.6 | V8 still emits one compare per access, against a constant, without loading the memory size | https://github.com/v8/v8/blob/main/src/wasm/turboshaft-graph-interface.cc | main | "Bounds check `index` against `kMaxMemory64Size - end_offset`, such that at runtime `index + end_offset` will be within `kMaxMemory64Size`, where the trap handler can handle out-of-bound accesses." | V |
| 3.7 | The flag is on by default in the local Node | `node --v8-options` | Node 26.8.2 | "--wasm-memory64-trap-handling (Use trap handling for Wasm memory64 bounds checks) type: bool default: --wasm-memory64-trap-handling" | V |
| 3.8 | JSC: fast (signal-handler) memory is only for 32-bit addresses, so memory64 is bounds-checked | https://github.com/WebKit/WebKit/blob/main/Source/JavaScriptCore/wasm/WasmMemory.cpp | main, 2026-09-30 | "if (Options::useWasmFastMemory() && desiredMemoryMode.value_or(MemoryMode::Signaling) == MemoryMode::Signaling && !addressType.is64Bit()) {" | V |
| 3.9 | Paper: software checks cost over 100%, virtual-memory tricks bring it to 12.7% | https://2024.splashcon.org/details/vmil-2024-papers/5/Performant-Bounds-Checking-for-64-Bit-WebAssembly (DOI 10.1145/3689490.3690400) | VMIL, 2024-10-20 (abstract read, not the PDF) | "we can reduce the bounds checking overhead from more than 100% when using software bounds checks to only 12.7% using a combination of shadow memory and guard pages" | V |
| 3.10 | Counter-example: 64-bit integer workloads can gain (multiprecision arithmetic; not LOSAT's profile) | https://pari.math.u-bordeaux.fr/archives/pari-users-2601/msg00039.html | 2026-01-23, Chromium 143 / Firefox 147 | "Only for Chromium the performance doubled (down to 108ms from 219ms)." | S |
| 3.11 | No slowdown figures published by V8 or in the Emscripten docs were found; no figure exists for LOSAT kernels | — | — | — | I |

### Q4. Rust toolchain and tooling

| # | Claim | Source | Date / version | Quote | Mark |
|---|---|---|---|---|---|
| 4.1 | `wasm64-unknown-unknown` is tier 3 with std "?"; identical in the 1.92.0, stable and nightly docs | https://doc.rust-lang.org/nightly/rustc/platform-support.html (and `/1.92.0/`, `/stable/`) | nightly 2026-10-10 | "Tier 3 targets are those which the Rust codebase has support for, but which the Rust project does not build or test automatically, so they may or may not work. Official builds are not available." | V |
| 4.2 | No prebuilt std; needs a self-built std | https://github.com/rust-lang/rust/blob/master/src/doc/rustc/src/platform-support/wasm64-unknown-unknown.md | master (page text from 2025-04-20, partly outdated: it still calls memory64 "in-progress") | "Rust does not yet ship pre-compiled artifacts for this target. To compile for this target, you will either need to build Rust with the target enabled ..., or build your own copy of `std` by using `build-std` or similar." | V |
| 4.3 | rustup ships no `rust-std` for any wasm64 target | https://static.rust-lang.org/dist/channel-rust-nightly.toml ; `channel-rust-1.92.0.toml` | manifests 2026-10-11 and 2025-12-11 | `[pkg.rust-std.target.wasm32-*]` entries only (7 on nightly, 6 on 1.92.0); none for wasm64 | V |
| 4.4 | `build-std` is nightly-only and unstable | https://doc.rust-lang.org/nightly/cargo/reference/unstable.html#build-std | nightly | "This feature is still in very early stages of development" | V |
| 4.5 | std on this target has no I/O, like `wasm32-unknown-unknown` | 4.2 page | same | "many I/O functions such as `std::fs` and `std::net` will simply return error" | V |
| 4.6 | The target page warns about compiler bugs and missing tests | 4.2 page | same | "you may encounter LLVM bugs"; "the Rust project doesn't run any tests for this target" | V |
| 4.7 | rustc has exactly one wasm64 target; no wasm64 WASI or Emscripten target | https://github.com/rust-lang/rust/blob/master/compiler/rustc_target/src/spec/mod.rs ; local `rustc +1.92.0 --print target-list` | master 2026-10-07; 1.92.0 | `("wasm64-unknown-unknown", wasm64_unknown_unknown),` is the only wasm64 entry | V |
| 4.8 | No wasm64 WASI targets are planned; WASIp1 has no 64-bit ABI | https://github.com/WebAssembly/wasi-sdk/issues/212 (alexcrichton) | 2026-02-14 | "There are no plans to release `wasm64-wasi{,-threads,p1,p1-threads}` targets at this time. Development of WASIp1 has ceased and adding 64-bit support would require redefining the existing ABI, so it's not "just change this thing there"." | V |
| 4.9 | Same from the Rust side | https://github.com/rust-lang/rust/issues/146944 (alexcrichton) | 2025-11-01 | "For example WASI has no definition on 64-bit linear memories at this time." | V |
| 4.10 | wasi-libc wasm64 support is an unmerged PR since 2023; the component-model ABI issue that blocks `wasm64-wasip2` is open | https://github.com/WebAssembly/wasi-libc/pull/444 ; https://github.com/WebAssembly/component-model/issues/22 | PR opened 2023-10-26, last activity 2025-10-11; issue last activity 2025-12-08 | PR title "Support wasm64" (state OPEN); issue "Interaction with the memory64 proposal" (OPEN) | V |
| 4.11 | wasi-threads itself is a legacy proposal | https://github.com/WebAssembly/wasi-threads (README) | 2024-07-17 | "this proposal is considered a legacy proposal, retained for engines that can only support WASI v0.1 (`preview1`)" | V |
| 4.12 | On non-WASI wasm with atomics, std locks work on wasm64 but `std::thread::spawn` is unsupported | https://github.com/rust-lang/rust/blob/master/library/std/src/sys/thread/mod.rs , `sys/thread/unsupported.rs`, `sys/sync/futex/wasm.rs` | master | "all(target_family = "wasm", target_feature = "atomics") => { mod wasm; pub use wasm::sleep; ... pub use unsupported::{ DEFAULT_MIN_STACK_SIZE, Thread, ..."; "Err(io::Error::UNSUPPORTED_PLATFORM)"; "#[cfg(target_arch = "wasm64")] use core::arch::wasm64 as wasm;" | V |
| 4.13 | SIMD intrinsics for wasm64 are nightly-only | https://doc.rust-lang.org/nightly/core/arch/wasm64/index.html | core 1.101.0-nightly (2026-10-10) | "This is a nightly-only experimental API. (simd_wasm64 #90599)" | V |
| 4.14 | wasm-bindgen supports wasm64 since 0.2.119 | https://github.com/wasm-bindgen/wasm-bindgen/pull/5004 | merged 2026-04-27 | "Add full `wasm64-unknown-unknown` (memory64) support to wasm-bindgen" | V |
| 4.15 | but its threads transform rejects memory64, so `wasm-bindgen-rayon` cannot be used on wasm64 | https://github.com/wasm-bindgen/wasm-bindgen/blob/main/crates/cli-support/src/transforms/threads/mod.rs ; issue https://github.com/wasm-bindgen/wasm-bindgen/issues/5330 (open, 2026-09-10, no replies) | main (file last changed 2026-07-08; latest release 0.2.129, 2026-09-25) | "// TODO: support memory64. This transform currently emits i32 addresses for thread counters, locks, stacks, and allocator calls. if mem.memory64 { bail!("threading does not currently support memory64"); }" | V |
| 4.16 | Rayon on the web needs an adapter; the adapter needs nightly + build-std even on wasm32 | https://github.com/rayon-rs/rayon/blob/main/README.md ; https://github.com/RReverser/wasm-bindgen-rayon (README, repo last pushed 2025-11-21) | main | "By default, when building to WebAssembly, Rayon will treat it as any other platform without multithreading support and will fall back to sequential iteration."; "you'll need to use the nightly compiler toolchain and pass some flags to rebuild the standard library"; "This crate was tested with `nightly-2025-11-15`." | V |
| 4.17 | Emscripten wasm64 is mature on the C side | https://github.com/emscripten-core/emscripten/blob/main/ChangeLog.md | 3.1.72 (2024-11-19); 6.0.1 (2026-06-22); current 6.0.12 (2026-10-08) | "The `MEMORY64` setting is no longer experimental."; "The `-sUSE_PTHREADS` and `-sMEMORY64` flags have been deprecated in favor of the more standard `-pthread` and `-m64` (or `--target=wasm64`) flags." | V |
| 4.18 | Rust has no `wasm64-unknown-emscripten` target yet: the PR adding it as tier 3 was opened a week ago and is waiting on review | https://github.com/rust-lang/rust/pull/163808 ; https://github.com/rust-lang/libc/pull/5611 | opened 2026-10-05, both OPEN | "Add a new Emscripten target that uses 64-bit WebAssembly memories."; "There's currently 113 test failures for the existing 32-bit Emscripten target" | V |
| 4.19 | Emscripten ABI is not stable across versions; the Rust docs recommend build-std (nightly) | https://github.com/rust-lang/rust/blob/master/src/doc/rustc/src/platform-support/wasm32-unknown-emscripten.md | master | "The Emscripten compiler toolchain does not follow a semantic versioning scheme that clearly indicates when breaking changes to the ABI can be made." | V |
| 4.20 | Multi-memory is standard (Chrome 120, Firefox 125 per features.json), but Rust cannot place data in a second memory; rustc only exposes the `multimemory` feature switch | https://webassembly.org/news/2025-09-17-wasm-3.0/ ; local `rustc --print target-features` | 2025-09-17; 1.92.0 and 1.97.0 | "a single module can now declare (define or import) multiple memories and directly access them"; "multimemory - Enable multiple memories." | V (standard), I (no Rust-level use) |
| 4.21 | Data can stay outside linear memory: OPFS sync handles read into a (shared) buffer inside a worker | https://fs.spec.whatwg.org/ | Living Standard, 8 October 2026 | "[Exposed=DedicatedWorker, SecureContext] interface FileSystemSyncAccessHandle { unsigned long long read(AllowSharedBufferSource buffer, optional FileSystemReadWriteOptions options = {});" | V |
| 4.22 | wasm-ld supports wasm64 together with a shared, imported memory; its wasm64 maximum is 16 GiB | https://github.com/llvm/llvm-project/blob/main/lld/wasm/Writer.cpp ; `lld/wasm/SyntheticSections.cpp` | main, 2026-10-08 | "if (ctx.arg.is64.value_or(false)) { // TODO: Update once we decide on a reasonable limit here ... maxMemorySetting = 1ULL << 34; }"; "if (ctx.arg.sharedMemory) import.Memory.Flags \|= WASM_LIMITS_FLAG_IS_SHARED; if (is64) import.Memory.Flags \|= WASM_LIMITS_FLAG_IS_64;" | V |
| 4.23 | Binaryen (`wasm-opt`) and wasm-tools handle memory64 | https://github.com/WebAssembly/binaryen/blob/main/src/tools/tool-options.h ; https://github.com/bytecodealliance/wasm-tools/blob/main/crates/wasmparser/src/features.rs | main | ".addFeature(FeatureSet::Memory64, "memory64")"; "pub memory64: MEMORY64(1 << 14) = true;" | V |

### Q5. What a port would force to be redone in LOSAT (quick look at `/home/kawato/losat-work`, commit b28e42e1)

| # | Item | Evidence | Mark |
|---|---|---|---|
| 5.1 | Target-gated code: 98 `target_arch = "wasm32"` sites in 17 files (`LOSAT/src`, `web/adapter/src`) and 6 uses of `core::arch::wasm32` in 4 files must also accept wasm64 (the SIMD ones need the nightly `simd_wasm64` feature, 4.13) | local `grep` | V |
| 5.2 | Pointer width: `usize` and pointers become 8 bytes (`rustc +1.92.0 --print cfg --target wasm64-unknown-unknown` gives `target_pointer_width="64"`); every table of `usize`/pointers doubles, which raises memory use and cache pressure before any data above 4 GiB is loaded | local command | V, effect I |
| 5.3 | ABI v2 is 32-bit: `docs/web/abi_v2.md` has "`emit(stream: u32, ptr: u32, len: u32)`" and "All functions return `i32`". A wasm64 ABI passes pointers and lengths as i64, which JS sees as BigInt (1.5) | repo doc; spec | V |
| 5.4 | The hand-written JS WASI shim cannot be widened: preview1 has no 64-bit ABI (4.8), so either a private 64-bit variant of each import is defined, or WASI is dropped and all host calls become `losat_host` imports | 4.8, 4.5 | I |
| 5.5 | Threads: `wasi_thread_start` / `thread-spawn` from wasi-libc are not available; per-worker stack and TLS set-up and a rayon `spawn_handler` must be written (what wasm-bindgen's transform does, and refuses for memory64, 4.15) | 4.12, 4.15 | I |
| 5.6 | Build identity: the pinned stable 1.92.0 with rustup's prebuilt std cannot produce wasm64; a pinned nightly, `rust-src` and a self-built std become part of the identity (`web/adapter/tools/check_build_identity.py`, `.cargo/config.toml` tables are keyed by wasm32 targets). `getrandom 0.2` with the `js` feature would not apply as is | repo files; 4.2–4.4 | V (facts), I (consequence) |
| 5.7 | All byte-parity checks against native and all performance baselines must be rerun for the new artifact: other target, other std build, other allocator (dlmalloc on `*-unknown-unknown`, 4.2), other integer widths | — | I |

### Local probes (Node 26.8.2, V8 14.6.202.34-node.28, x64, 31 GiB RAM; each ran in well under a second; pages were not touched)

`probe1.js` (constructor and validation limits), code in essence:
```js
new WebAssembly.Memory({address:'i64', initial:1n})                                  // A
new WebAssembly.Memory({address:'i64', initial:1})                                   // B
new WebAssembly.Memory({address:'i64', initial:1n, maximum:262144n, shared:true})    // C
new WebAssembly.Memory({address:'i64', initial:1n, maximum:262145n, shared:true})    // D
new WebAssembly.Memory({address:'i64', initial:1n, maximum:2n**37n-1n})              // E
WebAssembly.validate(module with memory limits flags 0x07 (max|shared|i64), min 1, max 262144 / 262145)  // I, J
```
```
A => 65536 ArrayBuffer
B => THROW TypeError: Cannot convert 1 to a BigInt
C => 65536 SharedArrayBuffer
D => THROW RangeError: WebAssembly.Memory(): Property 'maximum': value 262145 is above the upper bound 262144
E => THROW RangeError: WebAssembly.Memory(): Property 'maximum': value 137438953471 is above the upper bound 262144
I validate (memory i64 1 262144 shared) => true      J validate (memory i64 1 262145 shared) => false
```
`probe2.js` (growth, for `shared` false and true; `m = new WebAssembly.Memory({address:'i64', initial:1n, maximum:262144n, shared})`):
```
m.grow(65551n)  => prev=1 (bigint) byteLength=4296015872 = 4.0010 GiB (0.1 ms)
m.grow(65520n)  => prev=65552  byteLength=8 GiB      m.grow(131072n) => prev=131072 byteLength=16 GiB
m.grow(1n)      => THROW RangeError: WebAssembly.Memory.grow(): Maximum memory size exceeded
new Uint8Array(m.buffer).length => 17179869184        (same results for shared=true, SharedArrayBuffer)   process rss: 44 MiB
```
Toolchain: `rustup toolchain list` = stable, 1.92.0, 1.97.0 (no nightly, no `rust-src`). `rustc +1.92.0 --print target-list | grep wasm` = wasm32-unknown-emscripten, wasm32-unknown-unknown, wasm32-wali-linux-musl, wasm32-wasip1, wasm32-wasip1-threads, wasm32-wasip2, wasm32-wasip3, wasm32v1-none, wasm64-unknown-unknown.

## (c) Routes and rough effort (all estimates are I; they assume the current agent-driven workflow and have wide error bars)

**Route A — nightly `wasm64-unknown-unknown` + `-Zbuild-std` + own thread glue. About 4–7 weeks, high risk. Chrome, Firefox, Node only.**
- Pin a nightly and `rust-src`; build std with `+atomics,+bulk-memory,+simd128`; link with `--shared-memory --import-memory --max-memory` of at most 16 GiB (4.22, 1.9).
- Replace WASI: std has no files, stdio, clock, arguments or randomness on this target (4.5); every import the reactors use becomes a `losat_host` import (a new ABI version).
- Write the thread runtime: worker instantiation on the shared memory, stack and TLS per thread, rayon `spawn_handler`, shutdown (5.5). No upstream code to lean on (4.15).
- JS: `address:'i64'` memory, BigInt pointers and lengths on every call, views above 4 GiB, feature detection with fallback to the wasm32 builds for Safari.
- Rust source: the cfg sites of 5.1, a `usize`/`u32` audit, allocator behaviour above 4 GiB.
- A third artifact in build identity and CI; full parity rerun; A/B timing against wasm32 in V8 and SpiderMonkey before committing to it.

**Route B — Emscripten (`-m64 -pthread`) through Rust. Blocked today; about 5–8 weeks after rust-lang/rust#163808 lands, medium to high risk.**
- No in-tree Rust target yet (4.18); until then a custom target JSON on nightly would be needed. Also nightly + build-std (4.19).
- `std::thread` works on pthreads, so rayon needs no glue, and Emscripten already handles BigInt pointers and wasm64 with threads (4.17).
- The hand-written WASI shim, worker protocol and Node runner are replaced by the Emscripten runtime (pthread worker pool, its file system layer); emsdk becomes a build dependency whose ABI is not versioned (4.19).
- Same parity and performance reruns as route A; same Safari gap.

**Route C — stay on wasm32, lift the cap to 4 GiB, keep the database outside linear memory. About 2–4 weeks, low risk, works in all three browsers and Node on the pinned stable toolchain.**
- Lift the 1 GiB cap of the threaded build up to the wasm32 ceiling of 4 GiB (the cap is Owner decision 1 of this task).
- Store the database in volumes (OPFS or Blob), read one volume at a time into linear memory with `FileSystemSyncAccessHandle.read` (4.21), search volume by volume with whole-database statistics for E-values, stream results out through `emit`.
- No bounds-check penalty, no new target; parity work is limited to the chunked path. It matches the DB-search direction recorded on 2026-10-12 (large-DB search, E-values by whole-DB size).

Recommendation (I): take route C now. Re-open Memory64 when both hold: Safari stable ships it, and either `wasm64-unknown-emscripten` is in rustc or wasm-bindgen's threads transform accepts memory64 (issue #5330). A serial "hello" on `wasm64-unknown-unknown` would be a cheap first probe (under a day) but needs a nightly toolchain installed, which this study was not allowed to do.

## (d) Open points (not verified)

- iOS Safari: no primary statement. The JSC ceiling depends on the platform's address-space budget (2.11); the real limit on an iPhone or iPad is unknown.
- Which Safari release will carry Memory64. WebKit trunk has it on by default (2.9) but Safari 27.0 does not list it (2.8); the 27.1 notes page was not found (HTTP 404).
- Chrome on Android: chromestatus says 133 (2.1). Whether V8's 16 GiB reservation per memory64 (3.5) succeeds on phones, and the usable size there, were not checked.
- Firefox: shared memory64 growth past 4 GiB was read in test sources (2.6), not run; no Firefox binary was used. Firefox for Android limits not checked.
- Performance: 3.1 is from January 2025 and is SpiderMonkey's own statement; no V8 or JSC figures were found, and nothing was measured on LOSAT kernels. The V8 mitigation (3.5–3.7) suggests a smaller penalty in Chrome and Node than in Firefox, which is an inference.
- Whether LOSAT builds for `wasm64-unknown-unknown` with `+atomics,+simd128` at all was not tried (no builds allowed). Whether `RUSTC_BOOTSTRAP=1` on the pinned 1.92.0 could replace a nightly for build-std is unverified.
- Behaviour of `TextDecoder`, `postMessage` and typed-array views on a SharedArrayBuffer above 4 GiB in each browser was not checked.
- An unofficial wasm64 WASI sysroot exists (mentioned in wasi-sdk issue #212, 2025-12-17: "Just using your sysroot seems to be enough"); it is C-only, unsupported, and was not evaluated for Rust.
- Safari multi-memory: features.json lists none, a WebKit bug comment (bug 300538, 2026-02-28) says it landed; not verified.

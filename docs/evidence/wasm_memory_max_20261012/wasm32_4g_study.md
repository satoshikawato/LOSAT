# Study: raising the wasm32 shared-memory maximum of LOSAT from 1 GiB toward 4 GiB

Written 2026-10-12. Read-only research; no repository was edited and nothing was built.
Confidence marks: **V** = read in a primary source (or observed in a local probe) in this session,
**S** = secondary source, **I** = inference. Ages are relative to 2026-10-12.

URL prefixes used below:
RUST = `https://github.com/rust-lang/rust/blob/1.92.0` (tag commit 2025-12-08, 10 months) ·
LLVM = `https://github.com/llvm/llvm-project/blob/llvmorg-21.1.3` (LLD 21.1.3, the version inside rust-lld 1.92.0, about 12 months) ·
V8 = `https://github.com/v8/v8/blob/14.6.202.34` (commit 2026-04-23, 6 months; the V8 inside Node 26.8.2) ·
FX = `https://github.com/mozilla-firefox/firefox/blob/main` (160.0a1, fetched today) ·
WK = `https://github.com/WebKit/WebKit/blob/main` (fetched today) · WKB = `https://bugs.webkit.org/show_bug.cgi?id=` ·
EM = `https://github.com/emscripten-core/emscripten/blob/main` (6.0.12, 2026-10-08, 4 days) ·
WL = `https://github.com/WebAssembly/wasi-libc/blob/3f7eb4c7d6ede4dde3c4bffa6ed14e8d656fe93f` (2025-07-26, 14 months; the commit in wasi-sdk-27, which Rust 1.92.0 CI downloads: RUST`/src/ci/docker/host-x86_64/dist-various-2/Dockerfile` L88) ·
NODE = `https://github.com/nodejs/node/blob/v26.8.2`.

## (a) Summary

1. The 1 GiB is a default link argument of the rustc target (`--max-memory=1073741824`), the same in 1.92.0 and master. A later `-C link-arg=--max-memory=N` overrides it, because wasm-ld uses the last occurrence (checked with the bundled rust-lld).
2. Declare **65,535 pages** (`--max-memory=4294901760`), not 65,536. At exactly 4 GiB, wasi-libc `sbrk(0)` and the dlmalloc segment end wrap to 0, and Rust forbids an allocation that reaches the last byte of the address space. Emscripten stops one page short for the same reason. The cost is 64 KiB.
3. Safe to assume: on 64-bit desktop engines the address-space reservation does not depend on the declared maximum (V8 8 GiB, Firefox 4 GiB + 32 MiB, JavaScriptCore fast memory 4 GiB + 8 MiB). Going from 1 GiB to 4 GiB costs nothing there.
4. Not safe to assume: that a successful `new WebAssembly.Memory` can grow to its maximum. V8 without the trap handler (every Android build, 32-bit builds, Node under an address-space limit) and SpiderMonkey without huge memory reserve the declared maximum and silently fall back to a smaller reservation; `grow` then fails early.
5. Not safe to assume: iOS/iPadOS. JavaScriptCore reserves the maximum up front when it cannot get a fast memory, and throws `RangeError: Out of memory`. WebKit sources from before the iOS 18 line give small-address-space builds a 2 GiB cage and no fast memory, which matches the reported failure of a 2 GiB maximum on iOS 16.4 (the mapping to devices is inferred). Leaks across reloads were reported up to iOS 18.
6. 32-bit V8 and SpiderMonkey cap a wasm32 memory at 32,767 pages and reserve at most 1 GiB up front.
7. The bounds-check strategy does not depend on the maximum in V8 or JavaScriptCore. SpiderMonkey without huge memory keeps a 32-bit limit only when the maximum is below 65,536 pages (one more reason for 65,535).
8. JS side: pointers at 2 GiB or above arrive as negative numbers (`>>> 0` needed), `subarray` with a negative start silently reads from the end, and the Chromium and Firefox IDL for `TextDecoder.decode` does not accept shared views.
9. Node 26.8.2: a 4 GiB shared memory works, but `node:wasi` returns `EINVAL` for pointers at 2 GiB or above in WASI calls that take an i64 parameter (`clock_time_get`, `fd_seek`, ...). Open bug nodejs/node#62671; the fix PR #66636 (2 days old) is not merged.
10. Not verified on real hardware: iOS/iPadOS, Chrome on Android (64-bit and 32-bit), Firefox on Android, Windows commit behaviour. See (d).

## (b) Evidence by question

### Q1. Where the 1 GiB comes from and how to override it

| # | Claim | Source | Quote | Conf. |
|---|---|---|---|---|
| 1.1 | The target spec sets the maximum | RUST`/compiler/rustc_target/src/spec/targets/wasm32_wasip1_threads.rs` L20-23; local `rustc +1.92.0 --print target-spec-json` prints the same list | `&["--import-memory", "--export-memory", "--shared-memory", "--max-memory=1073741824"],` | V |
| 1.2 | master has the same link args | same file on `master` (last change 2026-06-09, 4 months); `diff` against 1.92.0 only touches os/env/arch enums, `main_needs_argc_argv`, `llvm_target` | `"-Wl,--max-memory=1073741824",` | V |
| 1.3 | User link args come after target args | RUST`/compiler/rustc_codegen_ssa/src/back/link.rs` L2245, L2452, L1905 | `add_pre_link_args(cmd, sess, flavor);` ... `add_user_defined_link_args(cmd, sess);` → `cmd.verbatim_args(&sess.opts.cg.link_args);` | V |
| 1.4 | wasm-ld takes the last `--max-memory` | LLVM`/lld/Common/Args.cpp` L28; LLVM`/lld/wasm/Driver.cpp` L639. Local run of the bundled `rust-lld -flavor wasm` (LLD 21.1.3) on `crt1-reactor.o`: `--max-memory=1073741824 --max-memory=4294901760` gives `max=65535 pages`; reversed order gives 16384 | `auto *a = args.getLastArg(key);` | V |
| 1.5 | Page multiple, 4 GiB ceiling | LLVM`/lld/wasm/Writer.cpp` L441, L482-487. Local: `4294967297` → "maximum memory must be aligned to the page size (65536 bytes)"; `4295032832` → "maximum memory too large, cannot be greater than 4294967296"; `4294967296` → 65536 pages | `uint64_t maxMemorySetting = 1ULL << 32;` | V |
| 1.6 | Shared memory always gets a maximum | LLVM`/lld/wasm/Writer.cpp` L496-500; `SyntheticSections.cpp` L367; threads proposal `https://github.com/WebAssembly/threads/blob/main/proposals/threads/Overview.md` L159 | `if (ctx.arg.sharedMemory && maxMemory == 0) {` ... `maxMemory = memoryPtr;` / "If linear memory is marked as shared, the maximum memory size must be specified." | V |
| 1.7 | The host may pass a smaller maximum, never a larger one | `https://github.com/WebAssembly/spec/blob/main/specification/wasm-3.0/2.2-validation.subtyping.spectec` L215-218 (repo main, JS API file last changed 2026-10-02); local Node probe `tmp/probe_link.js` | `-- if n_1 >= n_2` / `-- (if m_1 <= m_2)?`; Node: "memory import has a larger maximum size 65536 than the module's declared maximum 65535" | V |
| 1.8 | Spec ceiling | `https://github.com/WebAssembly/spec/blob/main/document/js-api/index.bs` L2238 | "The maximum `min` or `max` field of a 32-bit memory is 65,536 pages (4 GiB)." | V |

How to override: `RUSTFLAGS="-C link-arg=--max-memory=4294901760"` (or `rustflags` under `[target.wasm32-wasip1-threads]`). The target's linker flavor is `wasm-lld` (local target-spec JSON), so the flag goes to the linker as written.

### Q2. Reservation of a shared memory per engine

| # | Engine and configuration | 1 GiB maximum | 4 GiB maximum | Source and quote | Conf. |
|---|---|---|---|---|---|
| 2.1 | V8, 64-bit, trap handler on (desktop Linux, Windows, macOS; x64 and arm64) | 8 GiB | 8 GiB | V8`/src/objects/backing-store.cc` L35, L66: `constexpr size_t kFullGuardSize32 = uint64_t{8} * GB;`. Up to V8 12.0 it was `kFullGuardSize = uint64_t{10} * GB` (branch `12.0-lkgr`). Local Node 26.8.2: VmSize +8,388,608 kB for both maximums; VmRSS unchanged after `grow` to the maximum | V |
| 2.2 | V8 on Android (any CPU) | no trap handler | no trap handler | V8`/src/trap-handler/trap-handler.h` L19-26, L62-68: `((V8_OS_LINUX && !V8_OS_ANDROID) \|\| V8_OS_WIN \|\| V8_OS_DARWIN)` and `#error "The V8 trap handler should not be enabled on Android"` | V |
| 2.3 | V8 without guard regions (Android 64-bit; Node with `--disable-wasm-trap-handler` or a small `RLIMIT_AS`) | 1 GiB | 4 GiB | same file L73 `return byte_capacity;`, L454-466: "Retry with smaller maximum pages at each retry." `size_t sizes[] = {maximum_pages - delta, maximum_pages - 2 * delta, maximum_pages - 3 * delta, initial_pages};`. Local: VmSize +1.00 GiB and +4.00 GiB; see appendix for the silent fallback | V |
| 2.4 | V8, 32-bit (32-bit Chrome on Android, 32-bit desktop) | min(max, 1 GiB) | 1 GiB; engine limit 32,767 pages | V8`/src/wasm/wasm-limits.h` L42-44: `? 32'767   // = 2 GiB - 64Kib`; V8`/src/wasm/wasm-objects.cc` L946-963: "We try to reserve the maximum, but at most the allocation_maximum to avoid OOMs." | V |
| 2.5 | V8 bounds-check strategy | by trap handler only | same | V8`/src/wasm/wasm-module.cc` L47-51: `} else if (trap_handler::IsTrapHandlerEnabled()) {` ... `memory->bounds_checks = kTrapHandler;` ... "If the trap handler is not enabled, fall back to explicit bounds checks." Flags: `wasm_max_mem_pages` (default `kMaxUInt32`), `wasm_enforce_bounds_checks` | V |
| 2.6 | Chrome 4 GiB support | | since M83 | `https://v8.dev/blog/4gb-wasm-memory` (2020-05-14, 6.4 years): "you can now use up to 4GB of memory in WebAssembly applications" ... "Make sure to test on Chrome M83" | V |
| 2.7 | Chrome on Android in practice | | reported to work on a 5 GB phone | WKB`269777` (2024-02-20, by a Unity engineer): "Xiaomi Redmi Note 10 5G ... has 5GB of RAM, and reliably acquires 4GB in Chrome 121". 64-bit Chrome was first limited to Android 10+ with 8 GB RAM (2021 press reports); today's split is unknown | S |
| 2.8 | SpiderMonkey, 64-bit with huge memory | 4 GiB + 32 MiB + 64 KiB | same | FX`/js/src/wasm/WasmMemory.cpp` comment: "we reserve 4GiB + R bytes of memory regardless of the memory's byteLength"; FX`/js/src/wasm/WasmMemory.h` L233-244 `HugeMappedSize = HugeIndexRange + HugeOffsetGuardLimit + HugeUnalignedGuardPage;` | V |
| 2.9 | When huge memory is on | address width of 38 bits or more | | FX`/js/src/wasm/WasmProcess.cpp` L256-266: "conservatively chosen to continue using huge memory on our smallest address space system, Android on ARM64 (39 bits)"; `MinAddressBitsForHugeMemory = 38`. Platforms: "We have only tested huge memory on x64, arm64 and riscv64" | V |
| 2.10 | SpiderMonkey without huge memory | max + 64 KiB | max + 64 KiB, halving on failure | FX`/js/src/vm/ArrayBufferObject.cpp` L1738-1757: "try to reserve the biggest chunk in the range [initialPages, clampedMaxPages) using log backoff" | V |
| 2.11 | SpiderMonkey, 32-bit | 1 GiB | 1 GiB; limit 32,767 pages | FX`/js/src/wasm/WasmMemory.cpp` L337-343, L363-374: "the maximum heap size as observed by wasm code is one wasm page less than 2GB"; "prevent applications specifying a large max ... from unintentially OOMing the browser" | V |
| 2.12 | Firefox above 2 GiB | | since Firefox 89 (2021-06-01) | MDN Firefox 89 notes (`https://developer.mozilla.org/en-US/Firefox/Releases/89`): "ArrayBuffers can now be created with a length greater than 2GB-1 (up to 8GB) on 64-bit systems (Firefox bug 1703505)". That wasm memories above 2 GiB came with it is from a search summary; Bugzilla was not opened. Today: `ByteLengthLimit = size_t(16) * 1024 * 1024 * 1024;` | V / S |
| 2.13 | SpiderMonkey 32-bit bounds-check limit | yes | only below 65,536 pages | FX`/js/src/wasm/WasmModuleTypes.h` L925-928: `limits.maximum.value() < (0x100000000 / PageSizeInBytes(pageSize()))` | V |
| 2.14 | JavaScriptCore fast memory (signal based) | 4 GiB + 128 pages | same | WK`/Source/JavaScriptCore/runtime/BufferMemoryHandle.cpp` L73 `return PageCount::maxMemory32Bytes + fastMappedRedzoneBytes();`; WK`/Source/JavaScriptCore/runtime/OptionsList.h` L600 `maxNumWasmFastMemories, hasCapacityToUseLargeGigacage() ? 8 : 3`; WK`/Source/bmalloc/bmalloc/Gigacage.h`: `#if BPLATFORM(IOS_FAMILY) \|\| BCPU(ADDRESS32)` → `hasCapacityToUseLargeGigacage = false`, budget `(hasCapacityToUseLargeGigacage ? 64 : 16)` GiB | V |
| 2.15 | JavaScriptCore fallback for a shared memory (no fast memory left) | 1 GiB | min(max, budget / 4) = 4 GiB on iOS | WK`/Source/JavaScriptCore/wasm/WasmMemory.cpp` L169, L251: `tryAllocateGrowableBoundsCheckingMemory(reservedMaximumBytes)`; on failure `return nullptr;`. No step-down. Two tries with `vm.heap.collectSync(CollectionScope::Full)` between them | V |
| 2.16 | JavaScriptCore history | | | WebKit snapshot 2023-03-01 (commit 3c0664a5): `primitiveGigacageSize = (hasCapacityToUseLargeGigacage ? 32 : 2) * bmalloc::Sizes::GB;`, fast memory only if `BOS_EFFECTIVE_ADDRESS_WIDTH > 36`. Commit 416a9ccb (2024-03-15): "Remove JSValue Gigacage and donate region to Primitive Gigacage on iOS". Snapshot 2025-06-01 (commit 0a0952a2): `(hasCapacityToUseLargeGigacage ? 64 : 16)`. WKB`222097` c19 (WebKit developer, 2024-08-06): "believed to have been resolved by https://commits.webkit.org/277115@main, which only went into iOS 18, not 17.6" | V |
| 2.17 | iOS failures reported | | | WKB`255103` (2023-04-06, still NEW; iPhone 14 Pro Max, iOS 16.4): `"maximum": 2147483648 / 65536, "shared": true` → "You will get "Out Of memory"". WKB`260017` (2023-08, WONTFIX; growable SharedArrayBuffer with a 4 GiB `maxByteLength` fails on iPhone SE 2020): "designed to allocate virtual memory region with maxByteLength ... So large size can fail, and expected." WKB`250569` c6 (2023-03): `new WebAssembly.Memory({initial:10, maximum:65535, shared:true})` throws after leaks; "if I set the `maximum` to something much lower like `10`, it works" | V (text), S (device facts) |
| 2.18 | iOS sizes reported to work, and leaks | about 1.9 GiB unshared on iPhone Xs | | WKB`222097` c5 (WebKit developer, 2023-01): "On iOS 16.4, iPhone Xs, I similarly always get 1.9062gb if we ask for it from the beginning." WKB`281657` (2024-10-17, fixed 285841@main on 2024-10-29): "Safari will retain the memory between page reloads". WKB`256023` c5 (iPadOS 18.0, 2024-09): error "on approximately the 17th reload". emscripten-core/emscripten#19374 (2023-05): "The number of reloads before getting the error seems to depend on the `MAXIMUM_MEMORY` value" | S |
| 2.19 | Emscripten guidance for iOS with pthreads | none found | | EM`/src/settings.js` and EM`/site/source/docs/porting/pthreads.rst` contain no iOS note on `MAXIMUM_MEMORY` (grep for iOS/Safari) | V (absence in these two files) |

Which shipped iOS versions take the fast-memory path was not determined: the cage size and `hasCapacityToUseLargeGigacage` depended on the SDK constant `MACH_VM_MAX_ADDRESS` until commit 70592813 (2025-09-09, "iOS build should not use a large Gigacage capacity regardless of MACH_VM_MAX_ADDRESS"). **I**.

### Q3. The 4 GiB edge in the toolchain

| # | Claim | Source | Quote | Conf. |
|---|---|---|---|---|
| 3.1 | Emscripten stops one page short | EM`/src/lib/libcore.js` L185-189 | "Stay one Wasm page short of 4GB: while e.g. Chrome is able to allocate full 4GB Wasm memories, the size will wrap back to 0 bytes in Wasm side for any code that deals with heap sizes" / `Math.min(MAXIMUM_MEMORY, FOUR_GB - WASM_PAGE_SIZE)` | V |
| 3.2 | wasi-libc `sbrk(0)` wraps at 65,536 pages | WL`/libc-bottom-half/sources/sbrk.c` (main differs only by formatting, 2025-12-19). Local hand-assembled wasm in Node: `memory.size * 65536` as i32 is `0` at 65,536 pages and `-65536` (0xFFFF0000) at 65,535 | `return (void *)(__builtin_wasm_memory_size(0) * PAGESIZE);` | V (source and arithmetic); I (libc itself not executed) |
| 3.3 | `memory.grow` result is unambiguous | same file; local: `grow(1)` at 65,535 returns 65535, at 65,536 returns -1 | `if (old == SIZE_MAX) {` ... `return (void *)(old * PAGESIZE);` | V |
| 3.4 | dlmalloc segment end wraps | WL`/dlmalloc/src/malloc.c` L2702-2703, L4129-4130; file identical on main (last change 2024-05-17, so no upstream fix) | `((char*)(A) >= S->base && (char*)(A) < S->base + S->size)` | V (code); I (effect at 65,536 pages not executed) |
| 3.5 | The same bug class is confirmed in Rust's dlmalloc port | `https://github.com/rust-lang/rust/issues/144199` (2025-07-20, closed 2025-07-26; wasm32-unknown-unknown) | "`Segment::top` can overflow back to 0, so the infinite loop in `segment_holding` might be avoided by" ... | V |
| 3.6 | Rust: no allocation may reach the last byte | same issue, comment 2025-07-22; RUST`/library/core/src/ptr/mod.rs` L128-129 | "at least the very last byte of the address space must be unallocatable." / "`base + size <= usize::MAX`", "`size <= isize::MAX`" | V |
| 3.7 | Workaround named in that issue | same issue | "limit the memory size with the LLD flag `--max-memory` to exclude the last page of address space" | V |
| 3.8 | Single allocation limit | RUST`/library/core/src/alloc/layout.rs` L28-30; WL malloc.c L4129 `ssize < HALF_MAX_SIZE_T` | "The size, when rounded up to the nearest multiple of `align`, does not overflow `isize`" | V |
| 3.9 | Pointers above 2 GiB are valid Rust addresses | follows from 3.6 (addresses are `usize`; only the size is bounded by `isize::MAX`) | | I |
| 3.10 | LOSAT's target uses wasi-libc malloc, so the dlmalloc-rs fix does not apply | RUST`/library/std/src/sys/alloc/mod.rs` L72-78 (`target_os = "wasi"` → `mod unix;` → `libc::malloc`); `https://github.com/alexcrichton/dlmalloc-rs/pull/54` (2025-07-22) | "Don't allocate the final addresses on wasm" | V / I |
| 3.11 | wasi-libc tracker | GitHub issue search (4GB, 2GB, sbrk, malloc): only #220 "sbrk only works for limited inputs" (2020) and #377 (2023, `__heap_end`); nothing on the 4 GiB edge | | V (search; absence is not proof) |

### Q4. JS-side hazards above 2 GiB

| # | Claim | Source | Quote or observation | Conf. |
|---|---|---|---|---|
| 4.1 | i32 values reach JS signed | v8.dev article (2020); local Node: an export returning 0x90000000 gives `-1879048192`, `>>> 0` gives `2415919104` | "doing something like manually comparing a signed pointer to an unsigned one will (on address 2GB and above) return false" | V |
| 4.2 | Negative pointers fail silently in some APIs | local Node (`tmp/probe_js.out`) | `subarray(-4,-2)` returns bytes counted from the end, no error; `new DataView(buf,-4,2)` throws RangeError; `u8[-n]` is `undefined` | V |
| 4.3 | Length limits, V8 | branch `8.3-lkgr` `include/v8.h`: `: 0xFFFFFFFF;` (one byte short of a full 4 GiB view); `9.0`-`11.0`: `uint64_t{1} << 32`; V8`/include/v8-array-buffer.h` L401-410: `kMaxSafeBufferSizeForSandbox` (`32ULL * GB - 1`) or 2^53-1; 32-bit `std::numeric_limits<int>::max()`. Local Node: `new Uint8Array(sab).length === 4294967296` works | | V |
| 4.4 | Length limits, others | FX`/js/src/vm/ArrayBufferObject.h` L231-236 (16 GB; `INT32_MAX` on 32-bit). WebKit `runtime/PageCount.h`: `MAX_ARRAY_BUFFER_SIZE (1ull << 32)` in the 2023-03 and 2025-06 snapshots, `(1ull << 34)` on main (since 2026-08-07), `0x7fffffffu` without large typed arrays | a 65,536-page buffer sits exactly at the limit of shipped Safari | V |
| 4.5 | 4 GiB does not fit in 32 bits | v8.dev; local: `4294967296\|0` and `4294967296>>>0` are both 0 | "4GB itself is a 33-bit number" | V |
| 4.6 | `TextDecoder` on shared views | spec `https://github.com/whatwg/encoding/blob/main/encoding.bs` L1462: `decode(optional AllowSharedBufferSource input, ...)`. Chromium `third_party/blink/renderer/modules/encoding/text_decoder.idl` L41: `decode(optional [PassAsSpan] BufferSource input, ...)`. Firefox `dom/webidl/TextDecoder.webidl` L25: `decode(optional BufferSource input, ...)`. WebKit `Source/WebCore/dom/TextDecoder.idl` L45: `AllowSharedBufferSource?`. Node 26.8.2 decodes a shared view (local) | copy with `.slice()` before decoding | V (IDL); I (TypeError in Chrome and Firefox not run) |
| 4.7 | Views after growth | local: after `grow(1)`, `memory.buffer` is a new object; the old view keeps length 1,507,328. EM pthreads.rst L158 | "Pthreads + memory growth (`ALLOW_MEMORY_GROWTH`) is especially tricky" | V |
| 4.8 | Emscripten history | EM`/ChangeLog.md` entry 1.39.14 (2020-05-01); EM`/src/settings.js` L203-205 | "Optionally support 2GB+ heap sizes. To do this we make the JS code have unsigned pointers" / "If the maximum heap size is over 2GB, then pointers must be unsigned in JavaScript" | V |

### Q5. Node.js

| # | Claim | Source | Quote or observation | Conf. |
|---|---|---|---|---|
| 5.1 | No Node-specific size cap | local `node --v8-options`: `--wasm-max-mem-pages ... default: 4294967295`; `maximum:65536` and `65535` both create; `65537` → "value 65537 is above the upper bound 65536" | | V |
| 5.2 | Reservation and commit | local `tmp/probe_mem.js`: VmSize 1,416,240 → 9,805,104 kB at create for 16384, 65535 and 65536; `grow` to the maximum returns 23; VmRSS stays at 46-50 MB; last byte readable and writable | | V |
| 5.3 | Node turns the trap handler off under a small `RLIMIT_AS` | NODE`/src/node.cc` L1064-1082; flag `--disable-wasm-trap-handler` | "if it's not enough to allocate at least one cage for WASM, then the trap-handler-based bound checks cannot be used." | V |
| 5.4 | `node:wasi` rejects high pointers on the slow path | NODE`/src/node_wasi.cc` L300-301; `https://github.com/nodejs/node/issues/62671` (2026-04-10, open); `https://github.com/nodejs/node/pull/66636` (2026-10-10, open) | `return value->IsUint32();` / "all WASI functions taking `uint32_t` pointers rejected them with `EINVAL` ... The fast API path already received raw `uint32_t` and was unaffected." | V |
| 5.5 | Which calls are affected on 26.8.2 | local `tmp/probe_wasi64.js`, `tmp/probe_wasi.js` | pointer 0x80000000: `clock_time_get` errno 28, `fd_seek` errno 28, `args_sizes_get` 0, `fd_write` 0 (also on the last page, 200,000 calls). With `--no-turbo-fast-api-calls`, `fd_write` also returns 28 | V |
| 5.6 | Reach in LOSAT | WL`/libc-bottom-half/cloudlibc/src/libc/time/clock_gettime.c` L14 `__wasi_clock_time_get(clock_id->id, 1, &ts)` (a stack address); WL`/libc-top-half/musl/src/thread/pthread_create.c` L330 "As we allocate stack with malloc() instead of mmap/mprotect"; LOSAT uses `node:wasi` in `web/tools/wasi-test-runner.mjs` and `web/adapter/tests/v_abi.js` | a worker thread whose stack was allocated above 2 GiB would fail `Instant::now()` under `node:wasi`; the main thread stack is low (`--stack-first`) | I (each step V, chain not executed) |
| 5.7 | uvwasi bounds checks use `size_t` | NODE`/deps/uvwasi/src/wasi_serdes.c` L254-256 | `return end > offset && size <= (end - offset);` | V |

## (c) Recommendation

**Declared maximum: 65,535 pages** (`-C link-arg=--max-memory=4294901760`). Support: 3.1 (Emscripten), 3.2 and 3.4 (wrap to 0 in the bundled wasi-libc), 3.5 to 3.7 (confirmed failure of the same dlmalloc logic at the top of the address space, Rust's rule, and the `--max-memory` workaround named there), 2.13 (SpiderMonkey 32-bit limit), 4.4 and 4.5 (buffer length at an engine limit and not representable in 32 bits). Nothing found argues for 65,536.

**Host step-down ladder (pages):** 65535 → 49152 → 32768 → 24576 → 16384 → 12288 → 8192 → 4096, then fail with a clear message.

1. The first rung equals the module's declared maximum. A larger host maximum is a `LinkError` (1.7); a rung below the module's initial 23 pages is also a `LinkError`.
2. Step down only when the constructor throws `RangeError`. This is the JavaScriptCore fallback path (2.15, 2.17) and hard address-space limits. JavaScriptCore already runs a synchronous full GC before it gives up (2.15), so no delay between rungs is needed.
3. Starting at the top is harmless elsewhere: desktop V8, 64-bit Firefox and JavaScriptCore fast memory reserve a fixed size (2.1, 2.8, 2.14), and V8 without guard regions and SpiderMonkey without huge memory fall back by themselves (2.3, 2.10).
4. Keep 16384: it is the value shipped so far, and a 2 GiB cage (older WebKit, 2.16) cannot hold a 2 GiB reservation; a 2 GiB shared maximum is reported to fail on iOS 16.4 (2.17).
5. Constructor success is not capacity. Treat a failed `memory.grow` (allocation failure in wasm) as a normal "out of memory at N MiB" result and report `memory.buffer.byteLength`. Do not pre-grow the real memory to test it: wasi-libc starts its heap at `__heap_end` and extends it with `memory.grow`, so pre-grown pages would not be used (WebAssembly/wasi-libc#377, 2023-01). **I** for the consequence.
6. Optional, from the local probe only (**I**): V8's internal fallback is coarse. Under a 1.8 GB address-space limit, a declared 65536 gave 23 growable pages while a declared 16384 gave 4114. If a run fails at a size below the next rung, releasing the memory and retrying one rung lower can yield more real capacity. Release depends on GC.
7. Log the rung used and the size at failure.

Host JS must use `>>> 0` on every pointer and length from wasm, must not decode shared views directly in Chromium or Firefox, and must re-read `memory.buffer` after growth (Q4). For Node test runners, either keep test heaps below 2 GiB or avoid `node:wasi` for the above-2-GiB cases until nodejs/node#66636 ships (5.4 to 5.6).

## (d) Not verified; needs a real device or a test that was out of scope

- **iPhone and iPad** (Safari, and Chrome on iOS, which uses WebKit): which rung succeeds, whether the fast-memory path is taken, and behaviour after reloads. Classes: 4 GB RAM iPhone (SE 2020/2022, older base models), 8 GB or more iPhone Pro, iPad; iOS 17.x if still supported by LOSAT, iOS 18.x, and the current release.
- **Chrome on Android**: 64-bit build on a mid-range phone (4 to 6 GB RAM), 32-bit Chrome or WebView on low-RAM devices (expected cap 32,767 pages and a 1 GiB reservation, 2.4), and whether a 4 GiB reservation succeeds in practice. Chromium bug tracker pages were not opened (JS-rendered).
- **Firefox on Android** (arm64): whether huge memory is active (2.9).
- **Safari on macOS** (Apple silicon and Intel): limit of 8 fast memories and the fallback after that.
- **Windows** (Chrome, Edge, Firefox): commit charge when a shared memory grows; no source was read for this.
- **wasi-libc at exactly 65,536 pages**: the wrap was read in source but not executed, because builds were out of scope. It does not matter if 65,535 is adopted.
- **`Instant::now()` on a thread with a stack above 2 GiB under `node:wasi`** (5.6): chain of verified steps, not executed.
- The file that defines SpiderMonkey's `WASM_SUPPORTS_HUGE_MEMORY` was not found; Bugzilla 1703505 was not opened.
- Unity, Godot, Squoosh and SQLite-wasm trackers were not searched. From ffmpeg.wasm only #299 (2022, iPhone 12, iOS 15.2, "RangeError: Out of Memory" at load) was read.

## Appendix: local probes (Node v26.8.2, V8 14.6.202.34-node.28, x64 Linux under WSL2)

Scripts and raw outputs are in `$TASK_DIR/tmp/`: `probe_mem.js` (.out), `probe_cap.js` (.out), `probe_link.js` (.out), `probe_wasi.js` (`probe_wasi_65536.out`, `probe_wasi_slowpath.out`, `probe_wasi_tiers.out`), `probe_wasi64.js` (.out), `probe_js.out`, `memsec.py`. The wasm modules are assembled by hand inside the scripts; no toolchain was used. The linker check ran the bundled `rust-lld -flavor wasm` on the sysroot's `crt1-reactor.o` and wrote only to `tmp/` (output deleted).

`probe_cap.js` (complete):

```js
const fs = require('fs');
const vm = () => { const s = fs.readFileSync('/proc/self/status','utf8'); return (s.match(/VmSize:\s+(\d+) kB/)||[])[1]; };
const max = Number(process.argv[2]);
const before = vm();
let m; try { m = new WebAssembly.Memory({ initial: 23, maximum: max, shared: true }); } catch (e) { console.log(`max ${max}: CREATE FAIL ${e.constructor.name}: ${e.message}`); process.exit(0); }
const after = vm();
let pages = 23, step = max - 23, firstErr = '';
while (step > 0 && pages < max) { const d = Math.min(step, max - pages); try { m.grow(d); pages += d; } catch (e) { firstErr = firstErr || `${e.constructor.name}: ${e.message}`; step = Math.floor(step / 2); } }
console.log(`node ${process.version} execArgv ${JSON.stringify(process.execArgv)} declared max ${max}: VmSize delta at create = ${after - before} kB (${((after - before) / 1048576).toFixed(2)} GiB); growable to ${pages} pages (${(pages / 16384).toFixed(3)} GiB); buffer.byteLength ${m.buffer.byteLength}${firstErr ? '; first grow error: ' + firstErr : ''}`);
```

| Run | Declared max | VmSize added at create | Growable to | Note |
|---|---|---|---|---|
| default | 16384 | 8,388,608 kB | 16384 pages | guard regions |
| default | 65536 | 8,388,608 kB | 65536 pages | byteLength 4294967296 |
| `--disable-wasm-trap-handler` | 16384 | 1,048,576 kB | 16384 pages | |
| `--disable-wasm-trap-handler` | 65535 | 4,194,240 kB | 65535 pages | |
| `ulimit -v 3000000` | 32768 | 1,572,736 kB | 24582 pages | then "RangeError: WebAssembly.Memory.grow(): Unable to grow instance memory" |
| `ulimit -v 3000000` | 65535 | 1,115,200 kB | 16401 pages | constructor did not throw |
| `ulimit -v 3000000` | 65536 | 1,049,728 kB | 16402 pages | constructor did not throw |
| `ulimit -v 1800000` | 16384 | 263,296 kB | 4114 pages | |
| `ulimit -v 1800000` | 65536 | 1,472 kB | 23 pages | no growth at all |

`probe_wasi64.js` output (errno 28 = EINVAL; 76 is the rights error for seeking stdout, which shows the call reached uvwasi):

```
ptr 0x00010000: clock_time_get(i32,i64,i32) errno=0  fd_seek(i32,i64,i32,i32) errno=76  args_sizes_get(i32,i32) errno=0  | after 100k calls: clock_time_get failures=0
ptr 0x7ffffff0: clock_time_get(i32,i64,i32) errno=0  fd_seek(i32,i64,i32,i32) errno=76  args_sizes_get(i32,i32) errno=0  | after 100k calls: clock_time_get failures=0
ptr 0x80000000: clock_time_get(i32,i64,i32) errno=28  fd_seek(i32,i64,i32,i32) errno=28  args_sizes_get(i32,i32) errno=0  | after 100k calls: clock_time_get failures=100000
ptr 0x9c3f0000: clock_time_get(i32,i64,i32) errno=28  fd_seek(i32,i64,i32,i32) errno=28  args_sizes_get(i32,i32) errno=0  | after 100k calls: clock_time_get failures=100000
```

`probe_wasi.js` (declared max 65536, default flags): `grow(1) at 65535 -> 65535  size() = 65536  bytes() = memory.size*65536 in i32 = 0`; `grow(1) at 65536 -> -1`; `fd_write` from 0x80001100, 0xfffe0100 and 0xffff0100 all `errno 0`.

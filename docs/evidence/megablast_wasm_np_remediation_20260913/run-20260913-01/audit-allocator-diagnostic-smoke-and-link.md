# Independent audit: allocator diagnostic ABI and link boundary

The required read-only `ncbi_parity_auditor` supports both r2 artifacts for
limited public malloc-family accounting in the fixed untimed search diagnostic.
Actual source, artifact, link map, Node, libc and LLD hashes were independently
checked. Baseline has158 original/159 instrumented files; candidate159/160.

| Version | Artifact SHA-256 | Link map SHA-256 |
|---|---|---|
| baseline | `6709385d1efce204678df453e5f8f70e7ab6fff34b1075b67c34e48803aff57b` | `80eb1967ba8bb77c70dc7c753aae30b2570e99117b9956c3f11c795890f3a252` |
| candidate | `08faf9c4a7e67969389fe0397d8380c26fe8b8f61936091118b48ddde58048b2` | `9c7106d1f0417d48572b4d5494a1a967db93030f56a900b67b6bcbabc6edd707` |

Both versions pass two disposable selftests. Each has call deltas malloc9,
__libc_malloc1, free12, __libc_free4, calloc2, __libc_calloc1, realloc5,
posix_memalign2 and aligned_alloc1; expected failed allocations4. Live usable
bytes, live blocks and accounting errors return to their initial values.
The Rust Vec is checked while alive and immediately after drop.

Each version completes26 std-thread workers/104 lifecycle events across
1,8,8,8,1 cycles, with exact IDs, per-ID ordering and exit0. Realloc counts
increase64 times per worker. Post-wait live bytes are
baseline[1055036,1055036,5616,5616,5616] and
candidate[1055036,5616,5616,5616,5616]; live blocks decrease5→4.
A pre-wait/post-wait88-byte difference in the first baseline cycle confirms why
sampling after actual worker exits matters. No zero-live assumption is made.

The installed archive also has reallocarray/__reallocarray, but those objects
are not linked in these artifacts; their implementation references realloc.
The claim is limited to the nine linked public malloc-family entry/alias names,
not that every installed libc public allocator name is among those nine.
There is no observed double counting. Counts represent live usable bytes and
blocks, not requested payload, chunk headers, free chunks, RSS or linear memory.
They do not establish equivalence of instrumented scheduling to uninstrumented
artifacts or unlimited memory stability.

The first baseline smoke failed at code6 before realloc/alignment/Vec/thread
checks. Missing malloc(MAX) and free(NULL) calls were consistent with compiler
builtin optimization; the exact eliminating LLVM pass was not identified.
A direct exported-wrapper probe returned NULL and counted one expected failure.
The failed source/artifact/smoke are preserved. r2 makes original public/alias
callee pointers opaque, preserving actual --wrap function-table resolution.
The separate pre-compilation source-binding guard failure is also retained.

The search driver and host were independently reviewed after correcting the
nonexistent B1-artifacts/check_wasi_reactor.js path before launch. They use the
original frozen root helper hash. Added host operations only read counters
before/after, read retained result length, enforce untimed jobs and bind imports.
The fixed4 cases×6 calls×ABBA=96 scope, per-case instances, interleaving,
runPair, all five input-buffer frees, worker exit wait and retained LAST_RESULT
are unchanged. No selftest or extra clear runs in a search instance.

This audit supports diagnostic execution only. Actual search outputs and live
series require separate review. Original reuse and fixed16-call memory failures
remain FAILED; these smokes do not waive them.

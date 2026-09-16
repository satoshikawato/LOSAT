# I2: retain I1's real entries, remove repeated helper addressing

I1 is not adopted: its first fixed Mje compiled-module n8 session failed the
time guard (9.681366 → 10.311214 s, +6.5058%). Retain that failure; do not extend
its samples or subtract host startup to accept it. The normal-runtime cold n8
improvements remain real measurements of I1, not proof of acceptable reuse.

Independent code inspection found a concrete difference: the I1 TurboFan scan
recomputes the 28-byte helper address on b0=false; the optimized baseline keeps
that address. Both also have a per-visit Rust helper bounds check. This does not
numerically attribute the entire regression: one candidate sample also had slow
worker readiness, and another had an adjustable-clock disagreement. The chosen
metric remains monotonic full-lifecycle time.

I2 keeps I1's per-HSP Wasm function entry and changes the threaded-Wasm cursor to
an immutable slice prefix. A single macro owns the original selection/trace
statements; Native and plain-WASI use the indexed driver. The prefix excludes
sentinels 0/1. `split_last` chooses the old predecessor, assigning the tail is
the old decrement, and a b0 jump retains `next_larger.saturating_sub(1)` elements.
NCBI constructs strictly backward links at link_hsps.c:675-684,876-885. No-visit
ranges return before borrowing a slice; no allocation or search-state cache is
introduced. The source remains isolated from the production tree.

First check actual emitted code and exact differential visits/selection/trace,
including real jumps to 0/1 and a nonadjacent interior target. Run the focused
threaded raw/thread gate. Then screen the failed Mje n8 compiled-module condition
first: the same fixed AB/BA sessions, warmup plus three samples, full lifecycle,
normal Node, same 5%/50 ms and RSS guards, no concurrent jobs and no extension of
failed samples. Only a passing candidate proceeds to the remaining original
normal-runtime, Native/serial, group-state, reuse, oracle and browser checks.
Historical unavailable C8d timings supply no acceptance evidence.

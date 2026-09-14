# Independent read-only source finding: delayed epoch retirement

The actual dependencies are rayon-core1.13.0, crossbeam-deque0.8.6,
crossbeam-epoch0.9.18 and crossbeam-utils0.8.21. The required parity auditor
compared13 inspected source files against cached crate archive bytes and
Cargo.lock checksums. The allocator-reference snapshot retains the five epoch
files needed for the paths below. This dependency source explains ownership,
not NCBI alignment behavior.

- `src/default.rs:17–35`: a static default Collector and per-thread LocalHandle.
- `src/internal.rs:320–330`: registering a Local allocates it.
- `src/collector.rs:97–102` → `internal.rs:503–532`: handle drop finalizes the
  worker, pushes its bag to the global queue even when empty, then marks Local
  deleted without immediately freeing it.
- `src/sync/queue.rs:98–103`: enqueue allocates a node.
- `src/internal.rs:148–151,199–215`: retirement requires two epochs; one collect
  processes at most eight bags. Collection occurs on the first and each128th
  pin of a Local (`internal.rs:428–435`).
- `src/sync/list.rs:238–260` and `src/sync/queue.rs:155–167`: unlinking Local and
  retiring an old queue head themselves defer destruction, creating additional
  retirement stages. These functions exist in the actual linked artifact.

A worker's Local and queued bag node provide a source-grounded hypothesis for
small post-exit retained allocations. The observed +17,600 bytes/+16 blocks at
n8 is not yet assigned to those types by type-specific counters, so that numeric
attribution is not claimed. Joining/dropping the Rayon pool does not destroy the
static default collector or require its retired records to become immediately
zero. A standalone std-thread smoke is therefore not a full search-memory oracle.

The first completed baseline P132 diagnostic series has constant live bytes and
blocks after call2, while linear memory grows again at call5. This is a finite,
instrumented observation. It does not rewrite the original uninstrumented
six-call or fixed16-call plateau failures, and it does not prove an unlimited
memory bound. No production cleanup, extra epoch flush, longer blind repeat or
new memory exception has been added.

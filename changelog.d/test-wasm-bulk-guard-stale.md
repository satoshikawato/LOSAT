- Test tooling: the stale-bound test of the shared-memory guard's `memory.fill` and `memory.copy` runs in a child
  process that stays in V8's baseline compiler with stack checks off, and checks the bytes written. The three tests it
  replaces passed a guard without its bounds compare. No product change.

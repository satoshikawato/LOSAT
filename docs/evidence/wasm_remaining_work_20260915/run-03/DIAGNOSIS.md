# I1: actual hot-loop entries in normal V8

I1 failed the fixed compiled-module reuse guard and is not adopted. Its normal
Node cold n8 measurements improved substantially, but do not override that
failure. I2b retains the real scan entries and changes the threaded-WASI cursor
to avoid repeated helper addressing. Production Rust has not been changed.

## What the new evidence establishes

The emitted threaded Wasm retains the group function at index 1340 and adds the
actual predecessor-scan function at index 1368. The group calls that function
inside its original scan-eligible HSP branch. The existing host guard appends
two memory helpers; it does not renumber these existing function indices.

| Emitted function | Baseline WAT lines / locals | I1 WAT lines / locals |
|---|---:|---:|
| Group function | 16,796 / 132 | 15,883 / 125 |
| Extracted scan | inline in group | 1,264 / 20 |

These are emitted representation sizes, not speed measurements. Complete WAT,
selected function bodies and artifact-bound indices are in `work/*/`.

In a separate Mje n8 diagnostic, top instruction PCs were matched to printed
code ranges and required a preceding matching V8 code-registration record in
each isolate. `tier-pc-evidence.json` retains the per-isolate counts/ranges:

| Target sampled | Liftoff PCs | TurboFan PCs |
|---|---:|---:|
| Baseline group function | 23,058 | 0 |
| I1 extracted scan | 4,519 | 7,225 |

The I1 TurboFan samples occur in child isolates (7,224 in two children and one
in another). The main-script versus worker-entry compilation records identify
their roles. This demonstrates actual optimized scan execution in the diagnostic,
not merely allocation of TurboFan code. It does not count all calls, attribute
all elapsed time, or prove identical tier scheduling in unprofiled execution.
The profiler/code-print runs were substantially slower than normal execution;
their elapsed times are excluded from every speed comparison.

[V8's compilation documentation](https://v8.dev/docs/wasm-compilation-pipeline)
explains that new calls can use completed optimized code while an already-active
Wasm call continues with its original code. The observed PCs support this
mechanism as a material contributor. The normal-runtime paired measurements in
`I1-primary-normal/` independently test whether the candidate is useful.

## Correctness before the screen

- Fresh NCBI raw/thread gates: baseline 24/24 and I1 24/24, covering both primary
  inputs, thresholds 10/100/10000, and valid-query no-hit; Native/threaded n1/n8.
- Two meaningful Native unit tests passed: actual 0/1 jumps, unchanged initial
  state, nonidentity mapping, strict score/trim ties and copied f64 bits.
- `scan-differential-v2`: 800 cases per target/trace mode, Native and serial Wasm,
  trace off/on. Full visit/state stdout and diagnostic stderr match exactly.
  The original `scan-differential` attempt also passed but represented the four
  display-only hit coordinates as i32. Version 2 corrects them to production
  usize; the earlier attempt is preserved and is not the accepted type evidence.
- The independent actual-source audit found the extracted scan tokens unchanged;
  reinserting them recreates every token of the original group function.

## Later checks and the next candidate

- Normal Node search-body medians: Mje n8 22.630 → 12.360 s; AP n8
  40.763 → 21.125 s. I1 n1 was 16.723 and 30.838 s respectively, so n8
  was faster than n1 on both inputs. These are candidate measurements only.
- Native Mje n1 changed by +0.42%, within the fixed guard.
- Mje n8 compiled-module reuse failed the first AB session: 9.681366 →
  10.311214 s (+6.5058%). The planned stop rule ended the screen; BA was
  not run. Samples were retained without extension or startup subtraction.
- `I1-group-buffered`: all 352 groups, 1,927 rounds and 11,494 visits match
  exactly on Native, serial and threaded Wasm, with trace off/on. The first
  Native attempt timed out while writing every line to the mounted filesystem.
  The replacement uses buffered pipe capture with the same 120 s timeout;
  these diagnostic runs are excluded from performance evidence.
- Independent machine-code inspection found an extra helper-address
  recomputation in I1's optimized b0=false path. I2 uses a borrowed prefix
  cursor to address this, preserving the statement body and backward jumps.
- The initial I2 cfg incorrectly depended on `target_feature="atomics"`,
  which this Rust target does not expose. I2 is superseded without timing.
  I2b uses the existing `all(losat_wasi_threads, feature="wasm-threads")`
  selector. The independent audit confirmed empty scans, 0/1 and interior
  backward jumps, macro hygiene and unchanged selection/trace statements.

Diagnostic failures are retained: the first scan generator depended on an old
comment; the next stopped on the missing atomics cfg. An I2b diagnostic was
started before source preparation finished and had a dead-stripped thread ABI;
it is invalid evidence. The replacement requires a matching frozen source
manifest and retains the real thread ABI with one diagnostic spawn/join.

I2b must pass the failed reuse condition first, then the remaining target,
state, reuse, oracle and browser checks. No result is imported from lost
round-02 files.

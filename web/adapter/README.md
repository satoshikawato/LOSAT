# web/adapter — LOSAT Web engine adapter

A Rust crate that exports LOSAT Web ABI v2 ([`docs/web/abi_v2.md`](../../docs/web/abi_v2.md))
over the engine crate `LOSAT/`, for the two WASI reactors that the application runs:
`losat-web-serial.wasm` and `losat-web-threads.wasm`. ABI v1 (`losat_web_*`, used by
gbdraw) is linked in unchanged. [`web/AGENTS.md`](../AGENTS.md) governs this crate; engine
changes follow the root [`AGENTS.md`](../../AGENTS.md).

| Path | Content |
|---|---|
| `src/lib.rs` | the exports |
| `src/run.rs` | argv parsing (`validate`) and `run`: one `run_local` call writes every supported format, the HSP records and the diagnostics |
| `src/store.rs` | registered inputs (`register`, `release`) |
| `src/scan.rs` | the index scan (`scan_*`), reproducing `bio::io::fasta` |
| `src/describe.rs` | `describe`, generated from the engine's clap definitions |
| `src/emit.rs`, `src/json.rs` | output streams to the host, JSON |
| `tests/scan_properties.rs` | property tests of the scan against the parser |
| `tests/v_abi.js` | V-ABI: both reactors under Node, compared with the native CLI |
| `tools/check_build_identity.py` | build identity with the certified engine builds (plan TD-6) |
| `tools/build_reactors.py` | builds both reactors and records their identity |
| `tools/v_abi_cases.py` | the V-ABI searches (quick for CI, full for gate records) |

The release profile and the Wasm rustflags are copies of LOSAT's; the identity check
compares them, and every dependency version shared by the two `Cargo.lock` files. Build
outputs go outside the worktree (`--target-dir`).

```bash
python3 web/adapter/tools/check_build_identity.py
(cd web/adapter && cargo +1.92.0 test --locked)
RUSTUP_TOOLCHAIN=1.92.0 python3 web/adapter/tools/build_reactors.py --target-dir <dir> --output-dir <reactors>
(cd LOSAT && cargo +1.92.0 build --release --locked)
python3 web/adapter/tools/v_abi_cases.py --suite quick --out <cases.json>
node web/adapter/tests/v_abi.js --native LOSAT/target/release/LOSAT \
  --serial <reactors>/losat-web-serial.wasm --threads <reactors>/losat-web-threads.wasm \
  --cases <cases.json> [--out <dir>]
```

The full V-ABI suite (`--suite full`) uses the recorded input spellings of the
regression cases and needs the Gate A lexical root described in
[`docs/evidence/losat_web_e1a/README.md`](../../docs/evidence/losat_web_e1a/README.md).

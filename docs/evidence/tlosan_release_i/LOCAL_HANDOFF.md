# TLOSAN v0.2.0 Session I local archive handoff

Source/package candidate: `005e3d4b6cba6b5808334088fe9595c89efe01f8` on `feature/tlosan-tblastn-v0.2.0`. Session H's final evidence commit is `2976bd5f427cc5448315b6605ab85ee0787545a5`. This session changes release tooling and evidence outside the Cargo package; the assembler checks that `LOSAT/` and `README.md` still match the exact source candidate. The `.crate` contains `.cargo_vcs_info.json` with the candidate SHA. No output-affecting source or package input changed.

The independent Stage G certification covers only its declared local TBLASTN fixtures, not arbitrary inputs or options. Stage G's 213-entry and Session H's 30-entry SHA-256 manifests passed in this session. Session H's native 396/396 genetic-code matrix, native real/no-hit 72/72, command-WASI real/no-hit 90/90, code-32 9/9, negative 14/14, quality gates, package inspection and clean offline install are reused because their source, package inputs, candidate binaries, and acceptance conditions are unchanged. The separate Stage G focused TBLASTX 12-case gate is not the 20-case supported-profile audit.

## Exact artifacts and extracted execution

The v0.2.0-specific [contract](../../release/v0.2.0_rc_contract.json) pins the exact candidate, three binary hashes, three archive hashes, source crate hash, code-32 smoke fixture inputs and expected outfmt-6 output. It also declares the 27-code local TBLASTN outfmt 0/6/7 scope, both narrow subject-code exceptions, host and unsupported boundaries, and the ban on cross-contract speedup claims. The [assembler](../../release/assemble_v020.py) checks every input hash, creates reproducible `tar.gz` archives with fixed member metadata, extracts each archive, runs `losat 0.2.0`, then runs local code-32 TBLASTN from the extracted Native, serial command-WASI and threaded command-WASI files. It writes [handoff.json](handoff.json) and [artifact checksums](artifacts.sha256). The script SHA-256 is `d90ad2e236f69aafeefef1265e87e6675d9731c7af47769c26d7ea38bf7e8219`; contract SHA-256 is `2bf98d5fa461a6138f5fa6681ca9a50cba407d2c5f5fa35c8d963976e79d9bbf`.

| Artifact | Bytes | SHA-256 | Extracted check |
| --- | ---: | --- | --- |
| `LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz` | 1,513,369 | `5feb69c35e7a5aafea26bebdecd20f7b0fa2dfc3e362355e139958f0787e85b1` | PASS |
| `LOSAT-0.2.0-wasm32-wasip1-command.tar.gz` | 909,369 | `29c1a78afdadab1f32d1436abb09d0dff6ad0aa571cafa373561c052649ea36c` | PASS |
| `LOSAT-0.2.0-wasm32-wasip1-threads-command.tar.gz` | 1,011,689 | `c826cb99ac7abd7f918c3be5fbb975463dbc0215ae09ce86579c3c41df49f6ef` | PASS |
| `LOSAT-0.2.0.crate` | 2,813,168 | `50009387f331797ad380c1d54df0b10ded3d9c93b7055b083010406c2e6cddb8` | Session H clean install PASS; unchanged bytes |

All three extracted executable checks exited zero, printed `losat 0.2.0`, emitted no stderr, and produced the registered selected-code NCBI API outfmt-6 SHA-256 `0835fffadec3753e875f250c31ca075083f67a9ad91373d71bd5c197913010cf`. The threaded WASI check requested four threads. The archive members are only the program, required WASI host files where applicable, and MIT LICENSE. The final local files reside in `/tmp/tlosan-v020-session-i-final-scope-a-20260926/`; its `SHA256SUMS` verified all five listed outputs. A second assembly under `/tmp/tlosan-v020-session-i-final-scope-b-20260926/` produced an identical checksum file, including `handoff.json`.

Native execution was checked on Ubuntu 24.04 x86-64, glibc 2.39; the binary references symbols up to `GLIBC_2.34` and uses `libgcc_s.so.1`, `libm.so.6`, and `libc.so.6`. Serial command-WASI requires a WASI preview1 command host and one thread. Threaded command-WASI requires a compatible preview1 host with shared memory and `wasi.thread-spawn`. Both were executed with the included host scripts under Node 26.8.2 and `--experimental-wasi-unstable-preview1`. These host checks do not certify additional Native OS or architectures, browser/reactor modules, or arbitrary runtimes.

## Reproduce

From this branch, with the retained Session H binaries and exact source crate still at their recorded paths:

```bash
python3 docs/release/assemble_v020.py \
  --crate /tmp/tlosan-v020-cargo-package-target-d/package/LOSAT-0.2.0.crate \
  --output-dir /tmp/tlosan-v020-replay-new
(cd /tmp/tlosan-v020-replay-new && sha256sum --check SHA256SUMS)
```

For a fresh source rebuild, check out commit `005e3d4b6cba6b5808334088fe9595c89efe01f8` in an isolated worktree; use Rust/Cargo 1.92.0, the two installed command-WASI targets, and the three Session H locked release build commands in [Session H](../tlosan_release_h/LOCAL_ARTIFACTS.md). Build the `.crate` there with `cargo package --locked --offline` using a temporary target directory, then run this branch's assembler with `--repo-root <candidate-worktree> --crate <candidate-crate> --output-dir <fresh-/tmp-dir>`. Every binary, crate and archive hash must match the contract; a mismatch is a failed exact-candidate replay and requires a newly reviewed candidate. The assembler never calls NCBI executables. NCBI binaries and the code-32 API oracle remain comparison-only test tools.

The candidate serial command-WASI module also replayed the v0.1.0 TBLASTX supported-profile probes `p03`, `p12`, and `d06` from the committed manifest. All three exited zero with empty stderr and matched their frozen native LOSAT output hashes byte-for-byte; `d06` retains only the approved local-subject code-4 classification. Commands, fixture hashes, Node version, and output hashes are in [tblastx_serial_wasi.json](tblastx_serial_wasi.json), generated by [run_tblastx_serial_wasi.py](run_tblastx_serial_wasi.py). This is separate from the 20-case Native audit and Stage G focused 12-case regression.

The separate Native TBLASTX 20-case supported-profile gate completed on the exact candidate: 14 NCBI raw-byte matches and six manifest-designated selected-subject-code classifications, all frozen LOSAT hashes, required three-run repeatability, empty stderr, IDs and implicit/explicit code-1 control passed. [Audit method and replay](TBLASTX_AUDIT.md), [strict aggregate](tblastx_20_summary.json), and [independent audit](INDEPENDENT_AUDIT.md) provide the distinct Session I result. Session H's 13/20 partial run is superseded only for this gate; Stage G's focused 12 cases remain a separate regression record.

Final decision: **GO for the bounded exact-SHA local handoff** described in [v0.2.0 readiness](../../release/v0.2.0.md). No unresolved unexpected NCBI difference remains in its declared gates. Additional targets/profiles and arbitrary inputs/options are unproven. The final [SHA-256 evidence manifest](evidence.sha256) includes retained Session I records, scripts, release contract and decision files; verify it from this directory.

No tag, publication, registry upload, distribution, signing, or deployment occurred. The Stage G absolute times measure its former fixed binaries; LOSAT `-subject` and NCBI `-db` timings do not support a speedup ratio.

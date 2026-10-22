# docs/evidence/losat_web_e1d (archived)

LOSAT Web E1d (session S05) gate record: adapter and Wasm ABI v2, the full V-ABI run (`run-20260928T221453Z/v-abi-full/`). Verdict in the archived README: completion criteria met; independent audit supported. Session: [session_s05_e1d_adapter_abi_v2.md](../../losat_web_gui_sessions/session_s05_e1d_adapter_abi_v2.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `losat_web_e1d.tar.gz` (111,991 bytes) |
| Asset sha256 | `79807ba2b8df26b073b6184c24828ab88d92f4bc1cecb9c9533aac5b2d84bed3` |
| Files | 175 (1,426,258 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/losat_web_e1d/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p losat_web_e1d.tar.gz
echo "79807ba2b8df26b073b6184c24828ab88d92f4bc1cecb9c9533aac5b2d84bed3  losat_web_e1d.tar.gz" | sha256sum -c
tar -xzf losat_web_e1d.tar.gz
```

One file at a time: `git show 249c7454:docs/evidence/losat_web_e1d/<path>`.

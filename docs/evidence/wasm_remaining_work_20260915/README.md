# docs/evidence/wasm_remaining_work_20260915 (archived)

Wasm remaining work, rounds 01-03 (2026-09-15): H1a serial WASI module reuse, I4 TBLASTX predecessor scanning, adoption, recovery, diagnosis and policy records. Plan: [wasm_remaining_work_implementation_plan_20260915.md](../../archive/wasm_remaining_work_implementation_plan_20260915.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `wasm_remaining_work_20260915.tar.gz` (962,667 bytes) |
| Asset sha256 | `21ef7e7c1e71ac7dd33898100bbeb2dd82539579070ddab7abcfe9beb5746090` |
| Files | 1,244 (11,211,994 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/wasm_remaining_work_20260915/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p wasm_remaining_work_20260915.tar.gz
echo "21ef7e7c1e71ac7dd33898100bbeb2dd82539579070ddab7abcfe9beb5746090  wasm_remaining_work_20260915.tar.gz" | sha256sum -c
tar -xzf wasm_remaining_work_20260915.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

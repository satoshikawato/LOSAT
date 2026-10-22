# docs/evidence/wasm_threading_20260913 (archived)

Wasm threading remediation evidence (2026-09-13): cold-process and reuse performance tables and the packed evidence tarball. Report: [wasm_threading_remediation_results_20260913.md](../../wasm_threading_remediation_results_20260913.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `wasm_threading_20260913.tar.gz` (49,849,484 bytes) |
| Asset sha256 | `60bf37de19c650528c85c58722608e898631425b852ee498d4d91ce7d9802a47` |
| Files | 7 (49,992,081 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/wasm_threading_20260913/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p wasm_threading_20260913.tar.gz
echo "60bf37de19c650528c85c58722608e898631425b852ee498d4d91ce7d9802a47  wasm_threading_20260913.tar.gz" | sha256sum -c
tar -xzf wasm_threading_20260913.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

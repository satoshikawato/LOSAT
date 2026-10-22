# docs/evidence/wasm_four_program_20260915 (archived)

Four-program Wasm optimization, run-01 (2026-09-15): adoption decision (ADOPTED 2026-09-15), report, three-pair results, and the packed execution evidence. Plan: [wasm_four_program_optimization_plan_20260915.md](../../wasm_four_program_optimization_plan_20260915.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `wasm_four_program_20260915.tar.gz` (122,254,143 bytes) |
| Asset sha256 | `78c16ff02f2b826131ee405808588c364ee8a0464c65e17920d4d8e7878d19c5` |
| Files | 93 (136,310,079 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/wasm_four_program_20260915/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p wasm_four_program_20260915.tar.gz
echo "78c16ff02f2b826131ee405808588c364ee8a0464c65e17920d4d8e7878d19c5  wasm_four_program_20260915.tar.gz" | sha256sum -c
tar -xzf wasm_four_program_20260915.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

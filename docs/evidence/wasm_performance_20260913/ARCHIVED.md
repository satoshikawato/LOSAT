# docs/evidence/wasm_performance_20260913: archived files

2,942 of the 2,964 files in this directory (17,054,501 of 17,161,786 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c): no code, test, CI workflow, web build, skill, gate script,
or session prompt reads them (a script of the same run that stays here may need a
restore before it runs). Files that code reads stay here, and so do the summaries (`*.md`),
scripts, manifests, and small top-level records (22 files, 107,285 bytes). Links in the
summaries that point at archived files resolve again after a restore. The runs' own `SHA256SUMS` and
`run-20260913-02/checksum-summary.json` went with the files they list; `ARCHIVED_MANIFEST.tsv`
carries the sha256 of every archived file.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `wasm_performance_20260913.tar.gz` (2,730,052 bytes) |
| Asset sha256 | `9c270cdf63c1a2c42d20cbd653ee9d220a1c8a6721ac9d47ee010a78b465d51c` |
| Files | 2,942 (17,054,501 bytes) |
| Taken from | dev `9d80c770e32c6d045caee7979f7b21e307a16452` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/wasm_performance_20260913/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p wasm_performance_20260913.tar.gz
echo "9c270cdf63c1a2c42d20cbd653ee9d220a1c8a6721ac9d47ee010a78b465d51c  wasm_performance_20260913.tar.gz" | sha256sum -c
tar -xzf wasm_performance_20260913.tar.gz
```

One file at a time: `git show 9d80c770:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `commit-scope.json` | 1 | 1,615,460 | none |
| `run-20260913-01/` | 2,344 | 13,730,128 | 15 stay |
| `run-20260913-02/` | 597 | 1,708,913 | 3 stay |

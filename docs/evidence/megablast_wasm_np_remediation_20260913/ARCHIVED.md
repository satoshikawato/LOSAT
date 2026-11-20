# docs/evidence/megablast_wasm_np_remediation_20260913: archived files

129 of the 196 files in this directory (3,983,447 of 4,185,959 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (67 files, 202,512 bytes). Links in the summaries that point at archived files
resolve again after a restore. The run's own `run-20260913-01/SHA256SUMS` went with the files it lists;
`ARCHIVED_MANIFEST.tsv` carries the sha256 of every archived file.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `megablast_wasm_np_remediation_20260913.tar.gz` (771,823 bytes) |
| Asset sha256 | `9125215439e55a0311cedb331255706ddcce823ec8454ef4dbafa533ed036bd9` |
| Files | 129 (3,983,447 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/megablast_wasm_np_remediation_20260913/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p megablast_wasm_np_remediation_20260913.tar.gz
echo "9125215439e55a0311cedb331255706ddcce823ec8454ef4dbafa533ed036bd9  megablast_wasm_np_remediation_20260913.tar.gz" | sha256sum -c
tar -xzf megablast_wasm_np_remediation_20260913.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20260913-01/` | 129 | 3,983,447 | 67 stay |

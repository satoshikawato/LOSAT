# docs/evidence/tlosan_stage_a: archived files

139 of the 157 files in this directory (195,610 of 241,407 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (18 files, 45,797 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `tlosan_stage_a.tar.gz` (25,124 bytes) |
| Asset sha256 | `f69d059fe92be3adbcadb41052af956082cbb7b7d15c28b3434451db9c173c3b` |
| Files | 139 (195,610 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/tlosan_stage_a/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p tlosan_stage_a.tar.gz
echo "f69d059fe92be3adbcadb41052af956082cbb7b7d15c28b3434451db9c173c3b  tlosan_stage_a.tar.gz" | sha256sum -c
tar -xzf tlosan_stage_a.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `api_20260923_verified/` | 29 | 140,776 | 7 stay |
| `cli_20260923_verified/` | 110 | 54,834 | 2 stay |

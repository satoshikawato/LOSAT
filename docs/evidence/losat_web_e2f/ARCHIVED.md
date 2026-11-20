# docs/evidence/losat_web_e2f: archived files

222 of the 232 files in this directory (2,281,491 of 2,414,236 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (10 files, 132,745 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2f.tar.gz` (179,163 bytes) |
| Asset sha256 | `73324cea54feba84fe422bd1e320ccd405f4ac8bd23462017b9477f707307876` |
| Files | 222 (2,281,491 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2f/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2f.tar.gz
echo "73324cea54feba84fe422bd1e320ccd405f4ac8bd23462017b9477f707307876  losat_web_e2f.tar.gz" | sha256sum -c
tar -xzf losat_web_e2f.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `audit_round2/` | 2 | 91,933 | 3 stay |
| `run-20260930T163727Z/` | 217 | 2,187,725 | 1 stay |
| `run-20261001T130537Z/` | 3 | 1,833 | none |

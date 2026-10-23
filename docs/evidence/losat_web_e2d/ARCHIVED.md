# docs/evidence/losat_web_e2d: archived files

151 of the 186 files in this directory (3,501,917 of 3,781,787 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (35 files, 279,870 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2d.tar.gz` (360,302 bytes) |
| Asset sha256 | `4de9758a0a4fca4accbd29d5940113f01337206094541262e5f3959ae78444ed` |
| Files | 151 (3,501,917 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2d/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2d.tar.gz
echo "4de9758a0a4fca4accbd29d5940113f01337206094541262e5f3959ae78444ed  losat_web_e2d.tar.gz" | sha256sum -c
tar -xzf losat_web_e2d.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `INVENTORY.tsv` | 1 | 211,762 | none |
| `inventory/` | 5 | 212,210 | 11 stay |
| `run-20261005T133711Z/` | 145 | 3,077,945 | 3 stay |

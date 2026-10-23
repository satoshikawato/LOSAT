# docs/evidence/losat_web_w4: archived files

88 of the 94 files in this directory (1,477,210 of 1,577,611 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (6 files, 100,401 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_w4.tar.gz` (138,228 bytes) |
| Asset sha256 | `42f5d529b91e291a36f643fa6b3a98129903d9e6e8e8027be82a1972f1522c00` |
| Files | 88 (1,477,210 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_w4/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_w4.tar.gz
echo "42f5d529b91e291a36f643fa6b3a98129903d9e6e8e8027be82a1972f1522c00  losat_web_w4.tar.gz" | sha256sum -c
tar -xzf losat_web_w4.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20261009T122850Z/` | 30 | 537,284 | 1 stay |
| `run-20261009T135520Z-after-review/` | 20 | 324,634 | 1 stay |
| `run-20261009T140235Z-measure/` | 9 | 139,658 | none |
| `run-20261009T141300Z-measure/` | 9 | 151,002 | none |
| `run-20261009T142446Z-after-review/` | 20 | 324,632 | 1 stay |

# docs/evidence/losat_web_e2i: archived files

432 of the 464 files in this directory (5,160,368 of 5,416,489 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (32 files, 256,121 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2i.tar.gz` (498,613 bytes) |
| Asset sha256 | `0855bc2e7b14ee8df9adcfbc313f2034ea52aaf82bac7a513f355fd5460f7d24` |
| Files | 432 (5,160,368 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2i/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2i.tar.gz
echo "0855bc2e7b14ee8df9adcfbc313f2034ea52aaf82bac7a513f355fd5460f7d24  losat_web_e2i.tar.gz" | sha256sum -c
tar -xzf losat_web_e2i.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `INVENTORY.tsv` | 1 | 174,668 | none |
| `fixtures-20261003T103816Z/` | 5 | 45,013 | none |
| `fixtures-20261003T164528Z/` | 5 | 60,932 | none |
| `inventory/` | 12 | 336,584 | 8 stay |
| `mutants/` | 3 | 9,401 | 3 stay |
| `run-20261003T141051Z/` | 332 | 3,077,607 | 1 stay |
| `run-20261003T165813Z/` | 60 | 955,475 | 1 stay |
| `sweeps/` | 14 | 500,688 | 1 stay |

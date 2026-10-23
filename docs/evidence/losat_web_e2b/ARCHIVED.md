# docs/evidence/losat_web_e2b: archived files

340 of the 386 files in this directory (4,354,901 of 4,848,717 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (46 files, 493,816 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2b.tar.gz` (654,564 bytes) |
| Asset sha256 | `f5317d6619f1787683598393a1219ac0e0f37d5566fc6ae725c6be7e74423842` |
| Files | 340 (4,354,901 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2b/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2b.tar.gz
echo "f5317d6619f1787683598393a1219ac0e0f37d5566fc6ae725c6be7e74423842  losat_web_e2b.tar.gz" | sha256sum -c
tar -xzf losat_web_e2b.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `INVENTORY.tsv` | 1 | 235,279 | none |
| `audit/` | 2 | 244,767 | 7 stay |
| `closed_stdout/` | 1 | 2,766 | 1 stay |
| `inventory/` | 14 | 555,797 | 16 stay |
| `investigations/` | 4 | 10,639 | 4 stay |
| `run-20261002T084101Z/` | 2 | 20,440 | none |
| `run-20261003T021918Z/` | 316 | 3,285,213 | 1 stay |

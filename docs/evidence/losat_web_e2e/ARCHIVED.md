# docs/evidence/losat_web_e2e: archived files

377 of the 465 files in this directory (6,935,705 of 8,296,142 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c): no code, test, CI workflow, web build, skill, gate script,
or session prompt reads them (a script of the same run that stays here may need a
restore before it runs). Files that code reads stay here, and so do the summaries (`*.md`),
scripts, manifests, and small top-level records (88 files, 1,360,437 bytes). Links in the
summaries that point at archived files resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2e.tar.gz` (1,171,043 bytes) |
| Asset sha256 | `452be2d7e26ba628de6af6cf73bf19b252a3d3e4bae70a3dcb3712ec30aa495a` |
| Files | 377 (6,935,705 bytes) |
| Taken from | dev `9d80c770e32c6d045caee7979f7b21e307a16452` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2e/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2e.tar.gz
echo "452be2d7e26ba628de6af6cf73bf19b252a3d3e4bae70a3dcb3712ec30aa495a  losat_web_e2e.tar.gz" | sha256sum -c
tar -xzf losat_web_e2e.tar.gz
```

One file at a time: `git show 9d80c770:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `INVENTORY.tsv` | 1 | 785,936 | none |
| `audit/` | 6 | 219,541 | 21 stay |
| `gencode_api/` | 2 | 1,280 | 2 stay |
| `inventory/` | 13 | 787,280 | 27 stay |
| `run-20261004T163746Z/` | 314 | 3,535,331 | 5 stay |
| `s08pa/` | 23 | 708,654 | 11 stay |
| `s08pb/` | 15 | 40,546 | 6 stay |
| `sweeps/` | 3 | 857,137 | none |

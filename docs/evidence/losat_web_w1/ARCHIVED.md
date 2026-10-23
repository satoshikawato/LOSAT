# docs/evidence/losat_web_w1: archived files

76 of the 82 files in this directory (1,694,699 of 1,758,792 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (6 files, 64,093 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_w1.tar.gz` (265,897 bytes) |
| Asset sha256 | `4bbe16d597444749f083a8c03bd14c176fabdbe53fb77edbcc3a3e3a4d2a49de` |
| Files | 76 (1,694,699 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_w1/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_w1.tar.gz
echo "4bbe16d597444749f083a8c03bd14c176fabdbe53fb77edbcc3a3e3a4d2a49de  losat_web_w1.tar.gz" | sha256sum -c
tar -xzf losat_web_w1.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20261003T043424Z/` | 28 | 403,511 | none |
| `run-20261003T053627Z/` | 34 | 435,340 | none |
| `tblastn_profile/` | 14 | 855,848 | 2 stay |

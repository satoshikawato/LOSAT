# docs/evidence/losat_web_e2a: archived files

194 of the 203 files in this directory (1,828,776 of 1,987,909 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (9 files, 159,133 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2a.tar.gz` (142,461 bytes) |
| Asset sha256 | `7e94ef27a3554b81d0264ab322227fc8c5b6ee8c81896833f7e00a1b1132f8a6` |
| Files | 194 (1,828,776 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2a/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2a.tar.gz
echo "7e94ef27a3554b81d0264ab322227fc8c5b6ee8c81896833f7e00a1b1132f8a6  losat_web_e2a.tar.gz" | sha256sum -c
tar -xzf losat_web_e2a.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20260929T000046Z/` | 6 | 26,990 | none |
| `run-20260929T012339Z/` | 188 | 1,801,786 | 1 stay |

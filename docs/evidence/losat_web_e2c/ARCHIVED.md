# docs/evidence/losat_web_e2c: archived files

225 of the 244 files in this directory (2,393,426 of 2,590,490 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (19 files, 197,064 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2c.tar.gz` (175,173 bytes) |
| Asset sha256 | `dd4fcfaaccfca128f269ad3ab11b5aa38112b8ee1815083b841e67dc927ce314` |
| Files | 225 (2,393,426 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2c/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2c.tar.gz
echo "dd4fcfaaccfca128f269ad3ab11b5aa38112b8ee1815083b841e67dc927ce314  losat_web_e2c.tar.gz" | sha256sum -c
tar -xzf losat_web_e2c.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20260930T022908Z/` | 225 | 2,393,426 | 2 stay |

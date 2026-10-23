# docs/evidence/tlosan_stage_d: archived files

109 of the 931 files in this directory (119,220 of 2,011,690 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (822 files, 1,892,470 bytes). Links in the summaries that point at archived files
resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `tlosan_stage_d.tar.gz` (10,832 bytes) |
| Asset sha256 | `b3ce83ec8d3e11e28bc78517cf3df0e7531622e1a77a1fefe03683bf932be654` |
| Files | 109 (119,220 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/tlosan_stage_d/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p tlosan_stage_d.tar.gz
echo "b3ce83ec8d3e11e28bc78517cf3df0e7531622e1a77a1fefe03683bf932be654  tlosan_stage_d.tar.gz" | sha256sum -c
tar -xzf tlosan_stage_d.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `all_codes_20260925/` | 81 | 91,267 | 410 stay |
| `code32_local_api_20260925_full/` | 28 | 27,953 | 4 stay |

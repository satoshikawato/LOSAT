# docs/evidence/losat_web_e2g: archived files

250 of the 290 files in this directory (3,833,693 of 4,083,236 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c, batch 2): no code, test, CI workflow, web build, skill, gate
script, or session prompt reads them (a script of the same run that stays here may need a restore before it
runs). Files that code reads stay here, and so do the summaries (`*.md`), scripts, manifests, and small
top-level records (40 files, 249,543 bytes). Links in the summaries that point at archived files
resolve again after a restore. `stage2/PROMPT_COMMON.md`, the brief of the stage 2 agents,
names `stage1_functions.tsv` and `stage1_refs.tsv`; restore them before reusing it.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losat_web_e2g.tar.gz` (498,828 bytes) |
| Asset sha256 | `e9adb49f6d1c93ef0352c33336e3e98747a8f4444abb37d36054998f33b456ba` |
| Files | 250 (3,833,693 bytes) |
| Taken from | dev `cd8e0f89e7abaf96e87a55d3cf618dccc4642006` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losat_web_e2g/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losat_web_e2g.tar.gz
echo "e9adb49f6d1c93ef0352c33336e3e98747a8f4444abb37d36054998f33b456ba  losat_web_e2g.tar.gz" | sha256sum -c
tar -xzf losat_web_e2g.tar.gz
```

One file at a time: `git show cd8e0f89:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `INVENTORY.tsv` | 1 | 452,670 | none |
| `ci_wasi_timeout/` | 7 | 32,621 | 4 stay |
| `order/` | 3 | 179,863 | 3 stay |
| `run-20261002T042123Z/` | 230 | 2,236,746 | 1 stay |
| `stage1_functions.tsv` | 1 | 149,433 | none |
| `stage1_refs.tsv` | 1 | 350,251 | none |
| `stage2/` | 7 | 432,109 | 2 stay |

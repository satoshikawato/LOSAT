# docs/evidence/losatx_stage_g_authority_v3: archived files

15 of the 27 files in this directory (18,734,498 of 18,780,982 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c): no code, test, CI workflow, web build, skill, gate script,
or session prompt reads them (a script of the same run that stays here may need a
restore before it runs). Files that code reads stay here, and so do the summaries (`*.md`),
scripts, manifests, and small top-level records (12 files, 46,484 bytes). Links in the
summaries that point at archived files resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `losatx_stage_g_authority_v3.tar.gz` (8,594,001 bytes) |
| Asset sha256 | `df97f4722d65f6682c9b729c6cbc2e39d061105cb266f81834c664aa3ad2e15c` |
| Files | 15 (18,734,498 bytes) |
| Taken from | dev `9d80c770e32c6d045caee7979f7b21e307a16452` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/losatx_stage_g_authority_v3/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p losatx_stage_g_authority_v3.tar.gz
echo "df97f4722d65f6682c9b729c6cbc2e39d061105cb266f81834c664aa3ad2e15c  losatx_stage_g_authority_v3.tar.gz" | sha256sum -c
tar -xzf losatx_stage_g_authority_v3.tar.gz
```

One file at a time: `git show 9d80c770:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `run-20260928T121516Z/` | 15 | 18,734,498 | 12 stay |

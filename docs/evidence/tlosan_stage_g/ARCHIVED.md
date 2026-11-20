# docs/evidence/tlosan_stage_g: archived files

142 of the 214 files in this directory (10,127,645 of 19,041,108 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c): no code, test, CI workflow, web build, skill, gate script,
or session prompt reads them (a script of the same run that stays here may need a
restore before it runs). Files that code reads stay here, and so do the summaries (`*.md`),
scripts, manifests, and small top-level records (72 files, 8,913,463 bytes). Links in the
summaries that point at archived files resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `tlosan_stage_g.tar.gz` (2,614,997 bytes) |
| Asset sha256 | `da6a87c86cb473816c948bbb07a5c48d5502f53b2f3431faf90d145d04fe0cf3` |
| Files | 142 (10,127,645 bytes) |
| Taken from | dev `9d80c770e32c6d045caee7979f7b21e307a16452` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/tlosan_stage_g/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p tlosan_stage_g.tar.gz
echo "da6a87c86cb473816c948bbb07a5c48d5502f53b2f3431faf90d145d04fe0cf3  tlosan_stage_g.tar.gz" | sha256sum -c
tar -xzf tlosan_stage_g.tar.gz
```

One file at a time: `git show 9d80c770:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `benchmark/` | 59 | 142,737 | 4 stay |
| `cargo_test.log` | 1 | 69,637 | none |
| `code32_samples/` | 13 | 15,995 | none |
| `full_real/` | 13 | 9,777,679 | 2 stay |
| `full_real_code1_cross/` | 1 | 13,617 | none |
| `native_codes/` | 1 | 248 | 2 stay |
| `native_codes_multi/` | 28 | 81,630 | 2 stay |
| `native_physical/` | 1 | 255 | 2 stay |
| `native_real/` | 22 | 24,481 | 6 stay |
| `wasm/` | 2 | 1,327 | 2 stay |
| `wasm_real/` | 1 | 39 | 2 stay |

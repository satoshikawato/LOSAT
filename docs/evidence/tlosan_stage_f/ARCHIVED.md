# docs/evidence/tlosan_stage_f: archived files

60 of the 78 files in this directory (8,129,972 of 8,195,088 bytes) moved to a GitHub
release asset on 2026-10-23 (slim-20261022 P2c): no code, test, CI workflow, web build, skill, gate script,
or session prompt reads them (a script of the same run that stays here may need a
restore before it runs). Files that code reads stay here, and so do the summaries (`*.md`),
scripts, manifests, and small top-level records (18 files, 65,116 bytes). Links in the
summaries that point at archived files resolve again after a restore.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261023](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261023) |
| Asset | `tlosan_stage_f.tar.gz` (155,539 bytes) |
| Asset sha256 | `cf33595b25562c7b74ebd8f3cb29feae10576488f6b937c5e4536f6be3486d48` |
| Files | 60 (8,129,972 bytes) |
| Taken from | dev `9d80c770e32c6d045caee7979f7b21e307a16452` |

[`ARCHIVED_MANIFEST.tsv`](ARCHIVED_MANIFEST.tsv) lists the path, bytes, and sha256 of every archived file.
To restore them at the repository root (the tar writes `docs/evidence/tlosan_stage_f/...` next to the files that
stayed):

```bash
gh release download evidence-archive-20261023 -R satoshikawato/LOSAT -p tlosan_stage_f.tar.gz
echo "cf33595b25562c7b74ebd8f3cb29feae10576488f6b937c5e4536f6be3486d48  tlosan_stage_f.tar.gz" | sha256sum -c
tar -xzf tlosan_stage_f.tar.gz
```

One file at a time: `git show 9d80c770:<path>` (`<path>` as in `ARCHIVED_MANIFEST.tsv`).

## Where the archived files were

| Path | Archived files | Archived bytes | Files that stay |
| --- | ---: | ---: | --- |
| `native_codes/` | 3 | 988,316 | none |
| `native_codes_multi/` | 30 | 2,161,079 | none |
| `native_physical/` | 3 | 3,897,969 | none |
| `negative/` | 1 | 8,937 | none |
| `output_file/` | 17 | 33,798 | none |
| `stage_e_generated/` | 2 | 44,874 | none |
| `wasm/` | 4 | 994,999 | none |

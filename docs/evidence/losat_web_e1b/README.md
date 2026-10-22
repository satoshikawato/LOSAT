# docs/evidence/losat_web_e1b (archived)

LOSAT Web E1b (session S03) gate record: core entry for TBLASTN, output hashes after the change (`after/hashes.tsv`). Verdict in the archived README: completion criteria met; V-PERF over +5 % on some runs, left to the maintainer. Session: [session_s03_e1b_core_entry_tblastn.md](../../losat_web_gui_sessions/session_s03_e1b_core_entry_tblastn.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `losat_web_e1b.tar.gz` (52,250 bytes) |
| Asset sha256 | `cd4af1f218cb365113e68d122387643c3eb24f532c2c1569dd6a4bcc2132edba` |
| Files | 27 (634,463 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/losat_web_e1b/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p losat_web_e1b.tar.gz
echo "cd4af1f218cb365113e68d122387643c3eb24f532c2c1569dd6a4bcc2132edba  losat_web_e1b.tar.gz" | sha256sum -c
tar -xzf losat_web_e1b.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

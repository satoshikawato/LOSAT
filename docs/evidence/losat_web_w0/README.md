# docs/evidence/losat_web_w0 (archived)

LOSAT Web W0 (session S01) gate record: contract and skeleton. Verdict in the archived README: completion criteria met; `PD-LOSAT-WEB-APP-BOUNDARY` approved by the maintainer on 2026-09-29. Session: [session_s01_w0_contract_skeleton.md](../../losat_web_gui_sessions/session_s01_w0_contract_skeleton.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `losat_web_w0.tar.gz` (6,754 bytes) |
| Asset sha256 | `0ccc39ba4baf4b852bf13b3943151710a45ff3b4382f6c904909497b8e0ccfa8` |
| Files | 9 (14,059 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/losat_web_w0/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p losat_web_w0.tar.gz
echo "0ccc39ba4baf4b852bf13b3943151710a45ff3b4382f6c904909497b8e0ccfa8  losat_web_w0.tar.gz" | sha256sum -c
tar -xzf losat_web_w0.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

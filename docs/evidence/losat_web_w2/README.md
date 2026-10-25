# docs/evidence/losat_web_w2 (archived)

LOSAT Web W2 (session S10) gate record: data layer. Verdict in the archived README: completion criteria met. Session: [session_s10_w2_data_layer.md](../../losat_web_gui_sessions/session_s10_w2_data_layer.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `losat_web_w2.tar.gz` (14,887 bytes) |
| Asset sha256 | `a409d6da4ebbc9d5af0a55196ff5399beffd609690b39a35ed73f4884ebd0323` |
| Files | 9 (58,085 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/losat_web_w2/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p losat_web_w2.tar.gz
echo "a409d6da4ebbc9d5af0a55196ff5399beffd609690b39a35ed73f4884ebd0323  losat_web_w2.tar.gz" | sha256sum -c
tar -xzf losat_web_w2.tar.gz
```

One file at a time: `git show 249c7454:<path>` (`<path>` as in `MANIFEST.tsv`).

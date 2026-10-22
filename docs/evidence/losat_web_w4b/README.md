# docs/evidence/losat_web_w4b (archived)

LOSAT Web W4b (session S13b) gate run: NCBI-style search and results screens, E2E and screen logs, `run_gate.sh`. Session: [session_s13b_w4b_ncbi_style_ui.md](../../losat_web_gui_sessions/session_s13b_w4b_ncbi_style_ui.md); mapping: [ncbi_ui_mapping.md](../../web/ncbi_ui_mapping.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `losat_web_w4b.tar.gz` (61,187 bytes) |
| Asset sha256 | `af647fc31382b659c98c757f960edcf2d427d70c2b68fc26729584d153465e3c` |
| Files | 32 (587,027 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/losat_web_w4b/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p losat_web_w4b.tar.gz
echo "af647fc31382b659c98c757f960edcf2d427d70c2b68fc26729584d153465e3c  losat_web_w4b.tar.gz" | sha256sum -c
tar -xzf losat_web_w4b.tar.gz
```

One file at a time: `git show 249c7454:docs/evidence/losat_web_w4b/<path>`.

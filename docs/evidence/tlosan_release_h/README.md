# docs/evidence/tlosan_release_h (archived)

TLOSAN Session H (v0.2.0 release readiness, decision NO-GO): local artifact and release-readiness evidence (`LOCAL_ARTIFACTS.md`) and the independent read-only audit (`INDEPENDENT_AUDIT.md`). Cited by [v0.2.0.md](../../release/v0.2.0.md).

The files moved to a GitHub release asset on 2026-10-22 (slim-20261022 P2a): no code, test,
script, CI workflow, or skill read them.

| Item | Value |
| --- | --- |
| Release | [evidence-archive-20261022](https://github.com/satoshikawato/LOSAT/releases/tag/evidence-archive-20261022) |
| Asset | `tlosan_release_h.tar.gz` (55,924 bytes) |
| Asset sha256 | `f14d4e9c7bc06c8623d725f1a3ba44eebd72e9b2bf8fd4bd224af0acf4b16e43` |
| Files | 31 (1,347,344 bytes) |
| Taken from | dev `249c7454a82ac9e44ccfc35883d590c1ee395b26` |

[`MANIFEST.tsv`](MANIFEST.tsv) lists the path, bytes, and sha256 of every file. To restore the
files at the repository root (the tar writes `docs/evidence/tlosan_release_h/...`):

```bash
gh release download evidence-archive-20261022 -R satoshikawato/LOSAT -p tlosan_release_h.tar.gz
echo "f14d4e9c7bc06c8623d725f1a3ba44eebd72e9b2bf8fd4bd224af0acf4b16e43  tlosan_release_h.tar.gz" | sha256sum -c
tar -xzf tlosan_release_h.tar.gz
```

One file at a time: `git show 249c7454:docs/evidence/tlosan_release_h/<path>`.

# Security Policy

## Supported Versions

LOSAT v0.1.0 was published on 2026-09-17. Security fixes go to the newest minor
release line, as a new patch release, and to the `dev` branch.

| Version | Supported |
| --- | --- |
| Latest patch release of the newest minor line (`v0.2.x` once v0.2.0 is published; `v0.1.0` until then) | Yes |
| `dev` branch (unreleased) | Best effort |
| Older releases and snapshots | No |

## Reporting A Vulnerability

If the repository has GitHub private vulnerability reporting enabled, use that
channel first.

If private reporting is not available, open a GitHub issue with a high-level
summary and avoid posting exploit details, sensitive inputs, or private data.
Ask the maintainers for a private disclosure channel in that issue.

Please include:

- Affected LOSAT version or commit.
- Platform and target, such as native CLI, `wasm32-wasip1`, or threaded Wasm.
- Reproduction steps using minimal public input when possible.
- Whether the issue affects runtime behavior, generated output, crashes, or
  artifact integrity.

## Security Scope

Security fixes must preserve the repository's BLAST parity requirements. NCBI
BLAST C/C++ source remains the behavioral authority for search semantics, and
NCBI BLAST+ must not be added as a runtime dependency, build dependency, feature
fallback, or implementation component while addressing a vulnerability.

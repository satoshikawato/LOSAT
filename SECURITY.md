# Security Policy

## Supported Versions

LOSAT v0.1.0 is published. Security fixes go to the latest patch release of the
latest minor version (`v0.2.x` once v0.2.0 is released); until then, `v0.1.0` is
supported on a best-effort basis. Older releases are not supported.

| Version | Supported |
| --- | --- |
| Latest patch of the latest minor release (`v0.2.x` once released) | Yes |
| `v0.1.0` until v0.2.0 is released | Best effort |
| `dev` branch snapshots | Best effort |
| Older releases and snapshots | No |

## Reporting A Vulnerability

If the repository has GitHub private vulnerability reporting enabled (Security
tab, "Report a vulnerability"), use that channel first.

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

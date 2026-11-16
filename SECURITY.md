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

Report a vulnerability privately through GitHub's private vulnerability
reporting: the repository's Security tab, "Report a vulnerability". Please do
not open a public issue, pull request, or discussion for a vulnerability.

Only if private reporting is not available to you, open a GitHub issue with a
high-level summary and no exploit details, sensitive inputs, or private data,
and ask the maintainer for a private disclosure channel in that issue.

Please include:

- Affected LOSAT version or commit.
- Platform and target, such as native CLI, `wasm32-wasip1`, or threaded Wasm.
- Reproduction steps using minimal public input when possible.
- Whether the issue affects runtime behavior, generated output, crashes, or
  artifact integrity.

## Disclosure

The report stays private while the fix is prepared. After the fix is released,
the maintainer publishes a GitHub security advisory (GHSA) for the
vulnerability, naming the affected versions and the version that fixes it.

## Security Scope

Security fixes must preserve the repository's BLAST parity requirements. NCBI
BLAST C/C++ source remains the behavioral authority for search semantics, and
NCBI BLAST+ must not be added as a runtime dependency, build dependency, feature
fallback, or implementation component while addressing a vulnerability.

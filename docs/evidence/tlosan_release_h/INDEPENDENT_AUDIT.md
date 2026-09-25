# Session H independent read-only audit

Candidate: `005e3d4b6cba6b5808334088fe9595c89efe01f8`.

Verdict: **SUPPORTED for the documented NO-GO release-readiness decision and bounded local artifact/parity evidence.** This is not a release GO, a whole-input/whole-option parity certification, or an LOSAT-versus-NCBI speedup finding.

The independent `ncbi_parity_auditor` checked the 29-file evidence manifest that existed before this audit note, the actual size and SHA-256 of all four candidate artifacts, the embedded candidate commit and 178 members of the `.crate`, raw comparison rows (396/396 native code matrix, 72/72 native real/no-hit, 90/90 command-WASI real/no-hit, code 32 at 9/9, and unsupported paths at 14/14), and the all-features test log (no failures). The auditor also reviewed the qualified README/Cargo/CLI wording and the release decision's fixture scope, genetic-code exceptions, unfinished TBLASTX 13/20 audit, mixed `-db`/`-subject` timing boundary, and absent v0.2.0 archive/target contract. No further correction was requested.

Stage G's independent verdict is separate and remains limited to its declared TBLASTN fixtures and fixed-build absolute performance. The full LOSAT v0.2.0 release blockers in [the decision](../../release/v0.2.0.md) remain open. This note was added after the read-only audit; the final Session H evidence manifest includes it and is verified separately.

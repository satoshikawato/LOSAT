# Session I independent release-facing audit

Date: 2026-09-26. Reviewer: `ncbi_parity_auditor`, read-only agent `/root/release_audit`. No implementation or evidence edits were delegated to the reviewer.

## Verdict

**SUPPORTED for the bounded exact-SHA local v0.2.0 handoff** of source/package candidate `005e3d4b6cba6b5808334088fe9595c89efe01f8` on Linux x64 Native, serial command-WASI and threaded command-WASI, plus the unchanged source crate. No proven defect or remaining parity blocker exists within the documented scope. This verdict does not authorize tagging, publication, registry upload, distribution or deployment.

## Independent checks

- Reclassified and rehashed all 20 unique TBLASTX raw pairs for the committed manifest independently from the aggregate: 14 exact NCBI matches and only the six designated local-subject genetic-code classifications. Every LOSAT raw hash equals its frozen canonical output. Manifest commands, input hashes, executable identities, output IDs, empty stderr, all 15 required three-run checks, p12 implicit/explicit code-1 control and p14 query-code-4/default-subject control passed. No sorting or tolerance accepted a parity difference.
- Verified the exact candidate serial command-WASI p03/p12/d06 outputs, module hash, commands and empty stderr against the frozen Native hashes: 3/3.
- Checked the narrow subject-code source route against pinned NCBI C/C++ commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. Referenced worktree differences from Git blobs were only CRLF. The TBLASTX exception is distinct from the TBLASTN-only product decision; neither permits another behavioral difference.
- Reviewed the strict aggregate predicates and final decision/contract wording. Stage G fixture certification, Session H's historical 13/20 partial run, the focused TBLASTX 12-case gate and Session I's completed 20-case audit are separate records.
- Verified final contract SHA-256 `2bf98d5fa461a6138f5fa6681ca9a50cba407d2c5f5fa35c8d963976e79d9bbf` and final handoff SHA-256 `c9b52062f3f6c025e78583f932df38ae8b04542427bb53fdb66a6f5ba3eb943d` against two final assemblies. Both five-entry checksum sets passed and were identical; all three extracted version/code-32 execution checks passed. Archive, binary and source-crate hashes remain pinned.
- Rechecked Stage G's 213/213 and final Session H's 30/30 SHA-256 manifest entries. Unchanged candidate source/package/runtime evidence remains applicable; Session I does not change `LOSAT/` or `README.md`.
- Confirmed host requirements, unsupported paths, additional-platform limits and the absence of a LOSAT `-subject` versus NCBI `-db` speedup claim. Stage G absolute timings are not attributed to the metadata-updated candidate.

## Completion boundary

The final Session I SHA-256 manifest is generated after this audit note and covers retained evidence, scripts, contract, readiness decision and release entry point. Local manifest/link/syntax checks complete the handoff record. Additional targets/profiles, arbitrary inputs/options and every publication action remain outside this verdict. A changed source/package/archive input requires updated candidate provenance, hashes, affected gates and independent review.

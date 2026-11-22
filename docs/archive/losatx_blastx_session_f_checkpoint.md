# LOSATX BLASTX — accepted Session F source checkpoint

This commit records the Rust BLASTX implementation completed through Sessions
A–F: query preparation, preliminary search, linking and composition-adjusted
traceback, formats 0/6/7, deterministic native subject parallelism, command WASI,
and both reactor ABI entrances. It includes the fixed oracle fixtures required
by the source unit tests and the shared caller corrections.

[The checkpoint record](losatx_blastx_session_f_checkpoint.json) maps all 485
accepted source and fixture files to SHA-256, the frozen candidate, NCBI
authority, and twelve accepted binary fingerprints. Git commit identity is
additional provenance; it does not replace those source and binary identities.

Before staging, all 485 working files matched the accepted closure, and both
the full F implementation validator and the artifact integrity validator
were rerun successfully (actual exit status 0). Session F
has an actual independent `ncbi_parity_auditor` PASS with no unresolved findings,
a full implementation acceptance receipt, and a separate artifact integrity
receipt. The frozen records bind candidate SHA-256
`2757f2fb21dc7902d28e3ee9d1415a8d313b223ed1bb15fa23d4cb3e191a9a3f`.

The full raw evidence, generated outputs, binaries, and source archive remain
local at the evidence directory named in the checkpoint record. They are not
included wholesale in this implementation commit. This document is not a
standalone release certification package. Small canonical test fixtures are
tracked at their existing source-referenced paths.

The frozen F implementation validator also checks pre-commit HEAD
`eb7cb676f4b401169928a7e85dde42d9739738f9` and the pre-commit index. After this
commit its HEAD/index checks must refuse the new Git metadata. Do not modify
frozen evidence or reset Git state to make them pass. Session G must verify the
frozen receipts and raw evidence, match both the new Git blobs and working files
to the accepted closure, and bind a new precondition to the current commit.

Plain `wasm32-wasip1` remains serial. Real threading requires
`wasm32-wasip1-threads` and `wasm-threads`. Worker overlap evidence establishes
concurrency and lifecycle, not throughput. BLASTX accepts the pinned CLI's 26
query genetic codes, including 33 and excluding 32; no TBLASTN or TBLASTX
genetic-code exception extends to BLASTX.

Session G, formal benchmarks, distribution-platform certification, and release
publication remain unrun. The existing native serial CLI document records the
Session E scope; Session G will align README and current capability documents
with its final verified candidate.

[Next-session prompt](losatx_blastx_v0.2.0_sessions/session_g_next_prompt.md).

# Evidence storage

The run directory keeps concise reports, top-level JSON/TSV indices, patches and
checksums eligible for Git. Its local `.gitignore` excludes execution subtrees,
large raw/profile data, measured binaries, the archive and its detailed file map.
The user authorized a source/evidence commit and push including all tested-build
prerequisites; `commit-build-input-binding.json` records that dependency check.

The verified local handoff archive preserves actual source snapshots, measured
LOSAT artifacts, runtime files, input FASTA files, scripts, logs and raw evidence.
`handoff-archive.json` records its path, SHA-256 and byte size; the local
`handoff-file-map.json` records complete member hashes. The verification record
checks every member. NCBI executables, Rust/Node installations and Cargo dependency
caches remain external, versioned prerequisites.

The archive captures the pre-commit handoff state. Later commit bookkeeping in
STATUS, REPORT, COMMIT_HANDOFF and this file, plus the source-to-index binding
and independent commit review, is recorded in Git and in the final SHA256SUMS. Those documentation updates do
not change the measured source or raw outputs. Exact source snapshots retain
their original line endings; Git source text follows the existing LF convention.

A Git clone contains the implementation and build inputs, but does not include
the local execution subtrees, measured binaries or this archive. Keep the verified
archive with any future transfer of the complete evidence. The archive has not
been uploaded by this handoff.

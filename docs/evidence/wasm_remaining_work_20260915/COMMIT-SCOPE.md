# Committed implementation and evidence

The implementation includes H1a serial WASI module reuse from round 01 and I4
TBLASTX predecessor scanning from round 03. See
[the final results](run-03/I4-RESULTS.md) and
[the implementation plan](../../wasm_remaining_work_implementation_plan_20260915.md).

This commit retains the decision records, measurement metadata, source/artifact
hashes, environment samples, comparison results, focused-check logs and replay
scripts. Rejected candidates and failed attempts remain identified in those
records. The long-code4 oracle command is retained in
[its manifest](run-03/recovered-run01/long-code4/manifest.json): NCBI `-db` with
`-query_gencode 4 -db_gencode 4`, using AP027133 to build the database.

The approximately 5.4 GB local evidence tree also contains frozen source copies,
compiled artifacts, raw result bytes, complete state transcripts, generated
machine code and round-01 archives. These generated files remain at their
recorded local paths and are not included in this commit. Metadata references to
those files require the local evidence tree; a Git checkout alone cannot rerun
the archived-byte integrity checks. No separate evidence bundle was published.

The implementation files are committed with the existing Git LF line endings.
Pre-existing CRLF working-tree differences and unrelated changes are preserved
locally. The recorded test/source hashes describe the exact tested working-tree
bytes; `commit-source-binding.json` also records their LF-normalized committed
equivalents. No algorithm or host behavior changed during commit preparation.

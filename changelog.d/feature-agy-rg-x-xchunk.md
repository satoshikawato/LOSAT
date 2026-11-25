- `blastn`: matched NCBI BLAST+ 2.17.0's failing step order and output prefixes for unmakable
  query chunks (`CHUNK_SIZE` / `OVERLAP_CHUNK_SIZE`) via `ErrorAfterBatchReports`, kept 64-bit
  wrapping `size_t` chunk-range arithmetic in `query_split.rs` without debug-build underflow,
  handled traceback setup when no chunk of a split batch has a valid context, and checked split
  chunks during duplicate-id batch setup (`PD-LOSAT-NCBI-DEFECTS` 1.107).

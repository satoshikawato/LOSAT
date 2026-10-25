- TBLASTX's application and report code moved from `algorithm/tblastx/blast_engine/run_impl.rs` to
  `algorithm/tblastx/tblastx_app.rs` (NCBI `app/blast/tblastx_app.cpp`) and
  `algorithm/tblastx/blast_format.rs` (NCBI `format/blast_format.cpp`); `run_impl.rs` keeps the search;
  outputs do not change.

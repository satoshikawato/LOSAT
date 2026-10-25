- BLASTN's application and report code moved from `algorithm/blastn/blast_engine/run.rs` to
  `algorithm/blastn/blastn_app.rs` (NCBI `app/blast/blastn_app.cpp`) and
  `algorithm/blastn/blast_engine/blast_format.rs` (NCBI `format/blast_format.cpp`); `run.rs` keeps the
  search; outputs do not change.

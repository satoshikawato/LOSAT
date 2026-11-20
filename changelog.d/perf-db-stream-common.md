- The core of TBLASTN's `-db` stream (the check, the lengths, the reads, the block of each thread) moved to
  `algorithm/common/x_db_stream.rs`, with the program that its errors name as an argument, so that TBLASTX can use the same
  stream; TBLASTN's stages over it stay in `tblastn/x_db_stream.rs`. Outputs and errors do not change.

- Three tests that failed only on macOS in the five-platform run no longer fail there for reasons of
  the test: the side-file write-error test of the query blocks needs `/dev/full` and now runs on
  Linux only (macOS has none: the open fails instead of the write); the test that a stalled BLASTN
  seed stage lets the helpers leave lengthens the owner's sleep with each of its five attempts (60 to
  960 ms), so a helper job that starts late on a loaded machine is still covered; and the outfmt 11
  archive tests (BLASTN, BLASTP, BLASTX, TBLASTN) name the first differing byte, the number of
  differing bytes and the bytes around it when stdout or stderr differs from NCBI's. No program
  output changes.

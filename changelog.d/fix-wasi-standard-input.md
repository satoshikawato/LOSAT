- WASI command builds (`wasm32-wasip1`, `wasm32-wasip1-threads`): every search program (blastn,
  blastp, blastx, tblastn, tblastx, psiblast, rpsblast, rpstblastn) now reads standard input as
  natively: `-query -`, a run without `-query` (its default is `-`), and `-subject -` where the
  program reads it. Before, each such run stopped with "reading the query from standard input ('-')
  on this platform is not supported by LOSAT's <PROGRAM>" (exit status 1). LOSAT now reads the
  process's standard input in place, as NCBI reads `cin`, instead of a copy of its descriptor, which
  WASI cannot make; native builds give the same bytes as before. The WASI command checks
  (`LOSAT/tests/check_wasm_threading.py`) compare standard-input runs of every search program with
  the native command.

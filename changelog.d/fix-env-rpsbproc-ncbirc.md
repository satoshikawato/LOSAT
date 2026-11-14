- `losat rpsbproc` with a `.ncbirc` that has a syntax error now reports it as NCBI's rpsbproc 0.5.1 does
  (`Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: ` and the exception, naming
  that rpsbproc's own build) and goes on with the run, instead of LOSAT's refusal (exit status 1).

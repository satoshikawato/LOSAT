- CI: the fast and nightly output regressions no longer link the checkout to
  `/mnt/c/Users/genom/GitHub/LOSAT` with sudo. The TBLASTN Stage G cases, whose frozen outfmt 0/7
  bytes name that path, run from the repository root with a relative run root of the same length
  (`LOSAT/tests/frozen_paths.py`; outfmt 0 wraps the path at 68 columns), and the run root in their
  output and stderr is replaced by the frozen prefix before hashing, in `capture_outputs.py` and in
  the full V-ABI suite (Owner 2026-11-05: replace the prefix when comparing). It is the one
  replacement of a checkout-path prefix; this change adds no other rule, and the existing ones (the
  SAM @PG CL tag, the database lines of the approved TBLASTX rows, the drop-in corpus's run
  directory and makeblastdb report lines) stay. The frozen files are unchanged, and the
  comparisons pass from any checkout location.

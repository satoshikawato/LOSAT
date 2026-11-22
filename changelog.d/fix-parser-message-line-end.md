- LOSAT's own argument refusals and syntax errors (`-html`, `-remote`, an unknown option, a missing
  value) end with a line end, as clap's own errors do, so the shell prompt no longer lands on the
  message's line. LOSAT Web returns the same text from its ABI v2 `validate`; ABI v1 is unchanged.

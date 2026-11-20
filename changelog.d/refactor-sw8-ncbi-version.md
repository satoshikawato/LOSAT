- Engine: the NCBI BLAST+ version that LOSAT writes into its reports, archives, SAM `@PG` line and
  `-version` text ("2.17.0+") is defined once, at the crate root, instead of in 22 places (four
  constants, fifteen inline literals and three "0.1.0" fallbacks). No output changes.

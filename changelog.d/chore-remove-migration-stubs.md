- The library no longer has 24 modules that re-exported other modules under NCBI file names (two
  also held a few unused constants and presence-vector helpers) and that nothing used: `LOSAT::format`,
  `api::{blast_nucl_options, blast_options_handle, blast_results, tblastx_options}`,
  `blastinput::{blast_args, blastn_args, blastp_args, tblastx_args}` and `core::{aa_ungapped,
  blast_aalookup, blast_aascan, blast_diagnostics, blast_engine, blast_extend, blast_gapalign,
  blast_lookup, blast_nalookup, blast_options, blast_parameters, gencode_singleton, greedy_align,
  link_hsps, na_ungapped}`; `core::blast_hits` drops its glob re-export of TBLASTX's chaining module.
  Several pointed at another program's code (`core::blast_gapalign` re-exported TBLASTX's). The
  `losat` binary is byte-identical.

- `losat makeprofiledb` now supports `-dbtype delta` (row PS14) on both ASN.1 text and `-binary`
  scoremat inputs, writing `.wcounts` (fixed-point weighted residue frequencies) and `.obsr`
  (run-length compressed effective observations) alongside `.rps`, `.loo`, `.aux`, `.freq`, and
  version 5 / version 4 protein database volumes (including multi-volume `.pal` alias files,
  `-obsr_threshold`, `-exclude_invalid`, `-index`, `-taxid`, `-taxid_map`, and `-logfile`
  buffering), matching NCBI BLAST+ 2.17.0 byte for byte across 32 frozen oracle cases (`D01`-`D32`)
  and `rpsblast` checks (`D01`, `D04`). `-dbtype cobalt` remains refused, and `-dbtype delta`
  inputs where NCBI has undefined behaviour (`byRow TRUE`, `numRows < 28` with weighted
  frequencies, empty `numIndeptObsr`, partial trailing frequency columns, out-of-range or `NaN`
  `Uint4` conversions, and 32-bit offset overflows; `PD-LOSAT-NCBI-DEFECTS` 1.109) are refused
  explicitly.

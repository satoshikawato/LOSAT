- LOSAT Web: `-parse_deflines` in "Additional options" runs as on the command line in BLASTN,
  BLASTP, BLASTX, TBLASTN and TBLASTX, where it was refused ("-parse_deflines is not supported in
  the browser yet"). The outputs and warnings are the CLI's for the same command line, an input
  that the CLI refuses under the option fails the run with the CLI's error, and the results name
  the records by the IDs that the outputs show (`AB000001.1` for `>gb|AB000001.1| ...`). A record
  that is both a query and a subject and has masks gets the other outputs, and outfmt 0 is refused
  as the CLI refuses it. In the library, BLASTX's `run_local` no longer refuses `-parse_deflines`,
  and two empty inputs end with NCBI's `Empty CBlastQueryVector` as on the command line. ABI v2 is
  version 4 (`losat_web2_register_inputs`).

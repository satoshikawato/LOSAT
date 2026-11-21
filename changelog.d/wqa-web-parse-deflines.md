- LOSAT Web: `-parse_deflines` in "Additional options" runs as on the command line in BLASTN,
  BLASTP, BLASTX, TBLASTN and TBLASTX, where it was refused ("-parse_deflines is not supported in
  the browser yet"). The outputs and warnings are the CLI's for the same command line, an input
  that the CLI refuses under the option fails the run with the CLI's error, and the results name
  the records by the IDs that the outputs show (`AB000001.1` for `>gb|AB000001.1| ...`). A record
  that is both a query and a subject and has masks gets the other outputs, and outfmt 0 is refused
  as the CLI refuses it. In the library, BLASTX's `run_local` no longer refuses `-parse_deflines`,
  and two empty inputs end with NCBI's `Empty CBlastQueryVector` as on the command line. ABI v2 is
  version 4 (`losat_web2_register_inputs`).
- LOSAT Web: a run that wrote its other outputs but could not write outfmt 0, or outfmt 6 and 7,
  now says so in its results and in the Outputs view, and says where the error is: under Warnings
  in Run details, or that it was not kept (a BLAST+ error in TBLASTX with `-logfile /dev/null`).
  LOSAT's own refusal of an output is under Warnings whatever `-logfile` is, as the command line
  writes it to its standard error. Before, TBLASTX with `-logfile` wrote these errors to that file,
  where the command line's run of the same options writes none of them, and the results said that
  outfmt 0 showed the alignments of "the first 0 subjects", or that the query had no hits.
- LOSAT Web: Run details says where a run's warnings are when TBLASTX's `-logfile` takes them: in
  that file under Outputs, or not kept (`-logfile /dev/null`). Before, it said "The engine wrote no
  warnings." for such a run, whose warnings were in the file; it still says so when the file is
  empty and the other outputs wrote no warnings.

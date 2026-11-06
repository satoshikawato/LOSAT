- LOSAT Web: "Additional options" takes `-outfmt 13` and `14` (multiple-file BLAST JSON and XML2)
  and BLASTX's `-outfmt 11` (BLAST archive). For 13 and 14 the run's output is one ZIP archive of
  the files that the CLI writes with `-out results`: a side file per query (`results_1.json`, ...)
  and the master file that lists them (`results.json` or `results.xml`), byte for byte the CLI's.
  The Outputs offer it for Export (`…outfmt13.zip`), "Download All" names it as BLAST+'s
  `-outfmt` help does ("Multiple-file BLAST JSON (outfmt 13)"), and its command adds
  `-out results`. The engine's in-memory search
  (`run_local`) hands the side files of 13 and 14 to a receiver of the caller
  (`OutputSink::Files`) in all five programs; the command line's `-out` is unchanged.

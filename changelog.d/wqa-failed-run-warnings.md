- Web: the warnings of a run are the command line's for the run's options, then what the command
  line writes for the other formats that the run writes (outfmt 0, 6 and 7). A run that fails shows
  only the warnings that the command line writes before its error: no formatting warning of a format
  that its options do not ask for, and no warning or error of the run's other search. A run of two
  searches (`-num_descriptions`, `-num_alignments`) completes or fails as the search of its own
  `-outfmt`, so `blastp -num_alignments 0 -outfmt 6` is refused with "No hits are being saved", as on
  the command line. TBLASTX's file of `-logfile` holds what the command line's file holds; the
  warnings and BLAST+ errors that the run's other outputs add are under Warnings (not kept with
  `-logfile /dev/null`, which keeps none on the command line).

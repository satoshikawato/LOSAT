- Web search form: the BLASTP selects offer only what LOSAT's BLASTP reads (BLOSUM62, gap costs
  11/1, word sizes 3 and 5, compositional adjustments 0 and 2), TBLASTX's word size select offers 3,
  and TBLASTN's matrix select offers BLOSUM62 alone at word sizes 5 to 7; "Other" still takes any
  value, which the engine refuses with the CLI's message. dc-megablast's Template length and
  Template type show 18 and coding, and changing one writes both, as BLAST+ requires them together.
  Additional options refuses `-h` / `-help`, `-db` and the database filters, and `-mt_mode` in a
  short sentence instead of the CLI's help or its message about `-subject` / `-num_threads`. A
  negative value such as `-penalty -3` is no longer shown as an option name in the run's options and
  commands (WA-FM-02, -03, -04, -05, -07).

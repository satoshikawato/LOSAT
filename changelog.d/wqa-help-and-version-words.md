- LOSAT Web: `-h`, `-help`, `-help-full`, `--help`, `-xmlhelp`, `-version`, `-version-full`, `-version-full-xml`
  and `-version-full-json` in "Additional options" print what the command line prints for them in BLASTN, BLASTP,
  BLASTX, TBLASTN and TBLASTX, where they were refused ("-help is not available in the browser", "the option
  -version runs on the command line of LOSAT's BLASTN only: ..."). Run LOSAT needs no query or
  subject for them; the run searches nothing, its results say so, and the Outputs show the text under "Help" or
  "Version" and save it. ABI v2 is version 10 (stream 11). In the library API and pylosat these words stay refused;
  the refusal now names LOSAT Web among the places that run them.
- LOSAT Web's ABI v2: an argv without `-subject` (for example `blastp -query -h`, whose `-h` is the query's name)
  fails with the command line's error ("Either a BLAST database or subject sequence(s) must be specified") where it
  searched the registered subject. With `-logfile`, a help or version word makes the file of `-logfile` and then
  prints, as the command line does.

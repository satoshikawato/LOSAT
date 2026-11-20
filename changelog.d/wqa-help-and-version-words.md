- LOSAT Web: `-h`, `-help`, `-help-full`, `--help`, `-version`, `-version-full`, `-version-full-xml` and
  `-version-full-json` in "Additional options" print what the command line prints for them in BLASTN, BLASTP,
  BLASTX, TBLASTN and TBLASTX, where they were refused ("-help is not available in the browser", "the option
  -version runs on the command line of LOSAT's BLASTN only: ..."). Run LOSAT needs no query or
  subject for them; the run searches nothing, its results say so, and the Outputs show the text under "Help" or
  "Version" and save it. ABI v2 is version 5 (stream 11). `-xmlhelp` stays refused in LOSAT Web ("-xmlhelp is not
  available in the browser"). In the library API and pylosat these words stay refused; the refusal of the version
  words now names LOSAT Web among the places that run them.

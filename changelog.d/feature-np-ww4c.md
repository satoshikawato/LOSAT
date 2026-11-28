- LOSAT Web: a **Masking** tab runs `dustmasker`, `segmasker` and `convert2blastmask` in the browser
  with every option of their command lines (the form is read from each tool's `-xmlhelp`), on FASTA
  that is pasted or opened (`.gz`, `.bz2` and `.zst` are inflated), or for `dustmasker` and
  `segmasker` on a BLAST database (`-infmt blastdb`). The outputs are previewed and downloaded, and an
  output in `maskinfo_asn1_bin` or `maskinfo_asn1_text` has "Use as -mask_data in makeblastdb", which
  adds it to the makeblastdb tab's `-mask_data`.
- Web engine (npm package `losat`, ABI version 17): `runTool`, `validateTool` and `describeTool` take
  `dustmasker`, `segmasker` and `convert2blastmask` (`WasmEngine.runTool(['dustmasker', ...],
  { directory })`), with the command line's exit status, standard output, standard error and files
  (`-out`, `-logfile`) byte for byte.
- LOSAT Web's `windowmasker`: the check of a `-logfile` or `-out` that names a file of an
  `-infmt blastdb` input finds the database where the run finds it, through the `BLASTDB` of a
  `windowmasker.ini` in the run's directory too (the check read a `WINDOWMASKER.ini`).
- Web engine: `validateTool` refuses an `-out` that `runTool` refuses (a file of the run's
  directory, or the standard error, a device or a file of procfs such as `/dev/stderr`) where the
  run would refuse it, after the same steps: a run that first fails on its input (a missing `-in`
  or `-archive`) gives the command line's error, which `validateTool` leaves to the run instead of
  refusing the `-out`; `blast_formatter` and windowmasker's `-convert`, which open `-out` only
  after reading their input, refuse such an `-out` in the run.
- LOSAT Web's Supported mark of a tool run grades a value typed into a text box by the code path it
  selects, one rule for every tool: convert2blastmask's `-masking_options` by whether its text holds
  a byte outside visible ASCII (which NCBI writes as `#` with an error message) and
  `-masking_algorithm` by the program it names (`other` for any other name), windowmasker's `-meta`
  as a note written as given, where they were graded by their literal text, so most typed values
  read "Not yet supported". New NCBI BLAST+ 2.17.0 fixtures (91 runs) cover every combination of the
  Masking tab's choices that NCBI completes.
- LOSAT Web's tool tabs write a text that the tool writes or stores as given as it is typed, spaces
  at its ends included, where they were cut: makeblastdb's `-title`, `-mask_id` and `-mask_desc`,
  convert2blastmask's `-masking_algorithm` and `-masking_options`, and windowmasker's `-meta`. An
  empty value is given in Additional options (`-title ''`), as on the command line; an empty box
  writes no option. The other text boxes (names, numbers, keywords) and the search forms are as before.
- LOSAT Web: `-help-full` in Additional options prints the tool's help in every tool tab
  (windowmasker, Format Archive, makeblastdb and the Masking tools), as `-help` does and as the
  command line does, where only makeblastdb's tab read it as a help run.

- LOSAT Web: a **Masking** tab runs `dustmasker`, `segmasker` and `convert2blastmask` in the browser
  with every option of their command lines (the form is read from each tool's `-xmlhelp`), on FASTA
  that is pasted or opened (`.gz`, `.bz2` and `.zst` are inflated), or for `dustmasker` and
  `segmasker` on a BLAST database (`-infmt blastdb`). The outputs are previewed and downloaded, and an
  output in a `maskinfo_*` format is marked as the masking data that `makeblastdb -mask_data` reads.
- Web engine (npm package `losat`, ABI version 17): `runTool`, `validateTool` and `describeTool` take
  `dustmasker`, `segmasker` and `convert2blastmask` (`WasmEngine.runTool(['dustmasker', ...],
  { directory })`), with the command line's exit status, standard output, standard error and files
  (`-out`, `-logfile`) byte for byte.
- LOSAT Web's `windowmasker`: the check of a `-logfile` or `-out` that names a file of an
  `-infmt blastdb` input finds the database where the run finds it, through the `BLASTDB` of a
  `windowmasker.ini` in the run's directory too (the check read a `WINDOWMASKER.ini`).
- Web engine: `validateTool` refuses an `-out` that is the standard error, a device or a file of
  procfs (`/dev/stderr`, `/proc/self/comm`) for every tool, as `runTool` refuses it, where it
  accepted the words and the run then refused them.

- `losat makeblastdb -gi_mask` and `-gi_mask_name` write the masks of `-mask_data` as NCBI's GI-indexed mask files
  (`.gmd`, `.gnd`, `.gmo`, `.gno`, `.gmi`, `.gni`) and the alias file whose `MASKLIST` names them, with the files, log,
  messages and exit status of makeblastdb 2.17.0. One difference: for a database of one volume NCBI's alias file names
  the database itself, so that no program can open it (the build ends with `DB alias files are mutually recursive`,
  exit status 2); LOSAT's names the volume by its file name, so the database opens, the build ends normally, and
  NCBI's programs read it and its masks. LOSAT's own programs do not read GI masks yet. Two mask files whose GI-mask
  files are the same files are refused only where both would write one data file at the same time.

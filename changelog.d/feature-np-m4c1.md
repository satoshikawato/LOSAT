- `losat makeblastdb -gi_mask` and `-gi_mask_name` write the masks of `-mask_data` as NCBI's GI-indexed mask files
  (`.gmd`, `.gnd`, `.gmo`, `.gno`, `.gmi`, `.gni`) and the alias file whose `MASKLIST` names them, with the files, log,
  messages and exit status of makeblastdb 2.17.0. As in NCBI, a database of one volume built this way cannot be opened
  (the build ends with `DB alias files are mutually recursive`, exit status 2); one of several volumes can, by NCBI's
  programs. LOSAT's own programs do not read GI masks yet. Two mask files whose GI-mask files would be the same files
  are refused.

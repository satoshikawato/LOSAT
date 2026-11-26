- All programs (command line and library API): `-logfile` pointing to an input file read by the
  run (such as `-query`, `-subject`, `.ncbirc` / `ncbi.ini`, `<program>.ini`, `-conffile`, `-in`,
  `-ustat`, `-in_pssm`, `-in_msa`, `-gilist`, `-ipglist`, `-seqidlist`, `-taxidlist`,
  `-window_masker_db`, `-window_masker_taxid` counts, `-import_search_strategy`, `-mask_data`,
  `-taxid_map`, `-entry_batch`, `-ipg_batch`, `-ids`, `-exclude_ids`, `-i`, `rpsbproc` data files,
  or a BLAST database alias, alias list, GI-mask, MegaBLAST index, or component file of `-db` /
  `-filtering_db`, including relative paths, symlinks, hard links, and missing inputs that
  `-logfile` would create) is now refused with exit status 1 before opening or creating the file,
  leaving existing inputs unchanged, matching LOSAT Web's rule and preventing input corruption and
  NCBI BLAST+ 2.17.0's endless read-back loop (`PD-LOSAT-NCBI-DEFECTS` 1.117).

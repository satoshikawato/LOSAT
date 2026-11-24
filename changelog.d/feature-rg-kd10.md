- `losat blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `makembindex`, `makeprofiledb` and
  `convert2blastmask` now read every option value as its bytes, as NCBI BLAST+ 2.17.0 does: a file or
  database name with a byte that is not UTF-8 (a Latin-1 name, say) opens, a missing one gives NCBI's
  message with that byte, a text value such as `-infmt`, `-range` or `-outfmt` gives NCBI's own error or
  output, and titles and options written into the output keep their bytes. Before, these programs
  refused such a value. Database names, `BLASTDB`, the `[BLAST] BLASTDB` registry entry, the names in
  alias files (`DBLIST`, `MASKLIST`) and the target of a linked alias file are read as bytes too, so the
  search programs, `makeblastdb` and `rpsbproc` now also open a database found through a Latin-1
  `BLASTDB` or alias file, and `makeblastdb` builds in a directory whose name is not UTF-8. Their own
  option values (`-db`, `-query`, `-out` and the like) and `blastn -use_index` with such a path follow
  in a later change.

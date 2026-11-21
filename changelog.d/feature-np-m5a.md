- New program `losat windowmasker` (and `pylosat.windowmasker`): masks the repeats of nucleotide
  sequences with a unit counts file (`-ustat`, the four formats of NCBI's), from FASTA or a BLAST
  database, in the eight output formats, with `-dust`, `-ids`/`-exclude_ids` and the `-t_*_pct`
  settings, writing the bytes of NCBI BLAST+ 2.17.0's `windowmasker`. `-mk_counts` and `-convert`
  are refused for now (manual: What LOSAT refuses).

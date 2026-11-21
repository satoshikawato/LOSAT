- `losat makembindex` (and `makembindex` through `--install-shims`, `losat::api::run_makembindex`,
  `pylosat.makembindex`): the MegaBLAST index of a nucleotide BLAST database, or with
  `-old_style_index true` of FASTA input, with the bytes of NCBI BLAST+ 2.17.0's `.NN.idx`, `.shd`,
  `.map` and `-stat` files and its messages; every option of its `-help` (`-legacy`, `-nmer`,
  `-volsize`, `-stride`, `-ws_hint`, `-db_mask`, `-show_filters`, `-idmap`, `-stat`, `-iformat`,
  `-old_style_index`). The `-legacy false` strides for which NCBI gives no defined index are refused (above `-ws_hint` −
  `-nmer` + 1, where NCBI never finishes, and 32768 or more).

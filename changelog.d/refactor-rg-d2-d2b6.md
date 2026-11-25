- Internal refactor: consolidate the three duplicate implementations of NCBI's lagged Fibonacci
  generator (`CRandom` / `eGetRand_LFG` in `core/blast_encoding.rs`, `utils/dust.rs`, and
  `writedb/convert.rs`) into a single owner (`utils/random_gen.rs::NcbiRandom`), and move
  `find_matrix_file`'s duplicated path-separator and path-joining logic onto `PathRules` in
  `blastinput/ncbi_environment.rs`. Output bytes and runtime behavior are unchanged.

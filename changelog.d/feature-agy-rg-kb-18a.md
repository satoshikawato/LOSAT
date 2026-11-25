- `blastx` and `rpstblastn` now read FASTA queries, `-subject` inputs, XML pre-scans, and Web ABI
  v1 inputs through the shared `blastinput::fasta_reader` (the same reader used by `blastn`,
  `blastp`, `tblastn`, and `tblastx`), removing BLASTX's duplicate line/FASTA reader and unifying
  CLI and library batch cutting onto one `BatchCutter` (`blastx::input_source`). Output bytes are
  identical to `dev` and NCBI BLAST+ 2.17.0, while FASTA reading is faster and uses less memory on
  large `-subject` files; two LOSAT-only refusal messages (sequence exceeding `Int4` length and
  other-molecule Seq-id under `-parse_deflines`) now use the shared reader's wording and program
  name (`RPSTBLASTN` for `rpstblastn`).

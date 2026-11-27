- FASTA records of 4,294,967,296 letters or more (gaps included): `losat dustmasker`, `segmasker`,
  `windowmasker`, `makembindex` and `convert2blastmask` now read them as NCBI BLAST+ 2.17.0 does,
  with their length modulo 2^32 and NCBI's errors where NCBI's object manager or reader stops,
  instead of refusing them. A nucleotide record that NCBI packs into a piece of 4,294,967,293
  bases or more, and in `makembindex -old_style_index true` a record of 4,294,914,101 letters or
  more, are refused: NCBI crashes on them.
- `losat segmasker` prints a record of 2,147,483,648 letters or more without intervals, as NCBI
  BLAST+ does (it passes the length to SEG as a signed 32-bit number), instead of running out of
  memory on it.
- Reading a nucleotide record whose letters change between unambiguous and ambiguous often takes
  5 bytes of memory for each change instead of 21 (`makeblastdb`, the maskers, `makembindex`).
- `losat makeblastdb` builds a nucleotide record with gaps or ambiguous letters without a copy
  of 1 byte a letter: 38 % less peak memory on the longest chromosome of GRCz11 (123 MB, was
  200 MB), the same database.

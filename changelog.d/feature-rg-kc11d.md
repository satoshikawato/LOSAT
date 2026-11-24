- FASTA records of 4,294,967,296 letters or more (gaps included): `losat dustmasker`, `segmasker`,
  `windowmasker`, `makembindex` and `convert2blastmask` now read them as NCBI BLAST+ 2.17.0 does,
  with their length modulo 2^32 and NCBI's errors where NCBI's object manager or reader stops,
  instead of refusing them. A nucleotide record that NCBI packs into a piece of 4,294,967,293
  bases or more, and in `makembindex -old_style_index true` a record of 4,294,914,101 letters or
  more, are refused: NCBI crashes on them.

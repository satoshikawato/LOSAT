- `losat windowmasker -mk_counts` makes the unit counts file of a genome, from FASTA files or BLAST
  databases, in the four formats of NCBI's (`ascii`, `binary`, `oascii`, `obinary`), with `-checkdup`
  (the check for duplicated sequences), `-fa_list`, `-mem`, `-unit`, `-genome_size`, `-smem`, `-meta`,
  `-t_low`/`-t_high` and the `-t_*_pct` percentages; `-convert` writes a counts file in another format.
  Both write the bytes and messages of NCBI BLAST+ 2.17.0's `windowmasker`, so counts for
  `blastn -window_masker_db` and masks for `makeblastdb -mask_data` need no NCBI install. What NCBI
  gives no one result for is refused (manual: What LOSAT refuses).

- `blastp` and `rpsblast` now give NCBI BLAST+'s error of a query batch that fails, such as
  subjects too short for `-searchsp` or `OLD_FSC` or a database whose sequence file is missing,
  after the results of first batches whose queries have no letter to search, in every output
  format. Before, `blastp` refused this case in formats other than `-outfmt` 0 to 4, 6, 7 and 10,
  and `rpsblast` refused it in every format.

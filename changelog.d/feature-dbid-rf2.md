- `tblastn` with `-outfmt` 0 to 4 or `-html`, searching a database whose sequences have parsed ids and tax ids when the
  taxonomy database's `taxdb.btd` is empty (where NCBI BLAST+ crashes), now writes the results of the query batches
  before and then stops, as the other programs and formats do (it stopped before writing them).
- `blastx` and `tblastn` stop such a search with the same error as the other programs, the error of the tax id's lookup
  in `taxdb.btd`, in every output format; where NCBI makes a title from the tax id they named the title instead.

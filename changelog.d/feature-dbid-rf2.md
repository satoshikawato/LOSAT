- `tblastn` with `-outfmt` 0 to 4 or `-html`, searching a database whose sequences have parsed ids and tax ids when the
  taxonomy database's `taxdb.btd` is empty (where NCBI BLAST+ crashes), now writes the results of the query batches
  before and then stops, as the other programs and formats do (it stopped before writing them).
- `blastx` and `tblastn` stop such a search with the same error as the other programs, the error of the tax id's lookup
  in `taxdb.btd`, in every output format; where NCBI makes a title from the tax id they named the title instead.
- `tblastn` no longer refuses a database sequence that NCBI BLAST+ finds as a query read with `-parse_deflines` (whose
  load NCBI would fail on) in the tabular formats and `-outfmt` 18 when a sequence before it in the database is not
  reported: it compared the reported sequences' ranks with their numbers in the database.

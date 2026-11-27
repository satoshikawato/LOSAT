- `-db` searches of `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` with `-outfmt` 6, 7, 10 or 20, `blastn`
  with `-outfmt 17`, and `-outfmt` 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647` now run on a database whose
  stored titles hold bytes outside the printable ASCII range: the output and NCBI BLAST+'s error messages for each
  such byte are NCBI's, where LOSAT refused them before. NCBI writes them where each query batch's prefetch loads
  the sequence, and again at each row whose fields read the title, ids or taxonomy of the definition line. These
  runs are still refused where the prefetch can be off (`PRE_FETCH_SEQS_LIMIT` set to another value), on a
  database with parsed ids, where NCBI's object manager settings (`OBJMGR_*`) are set, where a query batch asks
  for a sequence that NCBI may still keep from an earlier batch or not, depending on the order of objects in its
  memory, and where a batch reports more than 200 sequences.

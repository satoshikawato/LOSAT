- Database search: a database sequence whose id is the id NCBI BLAST+ gives a query of the same
  batch (`lcl|Query_1` for the first query read without `-parse_deflines`, in any letter case,
  or the id of a `-in_pssm` query) is now refused in every format where NCBI shows that query in
  its place: `-outfmt 0` to `5` and `13` to `16`; the tables whose fields show what NCBI finds by
  the id (its title and ids, its length where the query's differs, its letters where the query's
  differ inside an alignment, the default fields outside `blastn`, `staxid` of a sequence with a
  taxid, and the taxonomy names when a taxonomy database is found); the organism report of a
  sequence with a taxid when a taxonomy database is found; SAM (`-outfmt 17`) where the query's
  length differs, or with `SQ` without `SR` its letters inside an alignment; and the `psiblast`
  PSSMs made of such a sequence (the next round's, and the one `-save_pssm_after_last_round`
  saves) where the query's letters differ inside an alignment; also where the taxonomy database
  is damaged, since NCBI loads nothing of such a sequence. With `-query_loc` the query's whole
  record is compared, as NCBI finds it. This holds in all eight search
  programs. LOSAT wrote the sequence's own output there, which differed from NCBI's (before: only
  `-outfmt 1` to `4` of five programs and some `rpsblast` and `rpstblastn` formats were refused).
  The other formats and fields still give NCBI's output, among them a copy of the query in the
  database, or the query with more letters after it, and `rpsblast`'s and `rpstblastn`'s organism
  report of a profile with a taxid when no taxonomy database is found, which LOSAT refused
  before. In `blastp`, `psiblast` and `rpsblast` the refusal now comes after the reports of the
  query batches before, as NCBI writes them.

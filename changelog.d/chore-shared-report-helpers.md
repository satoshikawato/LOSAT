- Report code that BLASTN, BLASTP, BLASTX and TBLASTN repeated is shared in `report`: the per-query
  Seq-align write step (`write_query_to`), the hit-list loop of the BLASTN and BLASTP Seq-align writers
  (`write_hit_list_seqaligns`), the XML `Hit` built from a subject's names (`XmlHit::of`), the archive's
  Karlin block (`ka_block`) and, for BLASTN, BLASTP and TBLASTN, its defline (`defline_of`). Outputs do
  not change.

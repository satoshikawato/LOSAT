- BLASTX `-culling_limit` now has frozen NCBI BLAST+ 2.17.0 reports among the tests (tests only; no output
  or behavior changes): 28 rows in `outfmt0_manifest.tsv` (`-subject`) and 14 in `db_search_manifest/blastx.tsv`
  (`-db`). They cover limits 1, 2 and 5; `-strand minus`; a 20,000-nt query that NCBI splits into two chunks
  and a `-db` search split with `CHUNK_SIZE`; `-max_target_seqs 1`, `-max_hsps 1` and `-subject_besthit`
  with culling; `-task blastx-fast`, `-seg no` and `-ungapped -comp_based_stats F`; ties by subject OID and
  by subject offset; a 200,000-nt query whose frame trees fork; the values 0x2, 2147483645 (no culling) and
  2147483646/7 (the wrap of PD-LOSAT-NCBI-DEFECTS 1.4); and `-outfmt` 0, 5, 6, 7, 8 and 11 (BLAST archive).

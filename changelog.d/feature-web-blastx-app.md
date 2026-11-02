- LOSAT Web offers BLASTX: the BLASTX tab has its search form (Task Selection with blastx and blastx-fast, the
  query's genetic code, NCBI's parameter sections and Other Parameters) with the defaults that the search uses, and
  its runs have the results screen of the other modes. BLASTX shows no query or subject subrange until the engine
  takes `-query_loc` / `-subject_loc` for it. The verification badge says that BLASTX was compared with the frozen
  NCBI rows of the outfmt 0 manifest, not certified with a release.
- `losat blastx -help` gives the defaults of `-word_size`, `-gapopen`, `-gapextend`, `-matrix`, `-threshold`,
  `-window_size`, `-culling_limit`, `-max_target_seqs` and `-sum_stats`, as the other programs do.

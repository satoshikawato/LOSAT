- LOSAT Web adapter: the index scan's parser kind 0 (the `bio::io::fasta` reader of the programs before NCBI's
  reader was ported) is gone, now that the app scans with NCBI's reader kinds; the adapter no longer depends on
  `bio` directly.

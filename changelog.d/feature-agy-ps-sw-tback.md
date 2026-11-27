- `psiblast -use_sw_tback`: supported across round 1 (standard matrix), multi-round PSSM
  iterations (`-num_iterations 0/2/3/...`), `-in_pssm`, `-in_msa`, `-comp_based_stats 0` to `3`,
  `-sum_stats true`/`false`, `-db` and `-subject`, multi-threaded runs (`-num_threads`), all
  supported output formats, and PSSM file exports (`-out_pssm`, `-out_ascii_pssm`), reusing the
  existing square-matrix and position-based Smith-Waterman redo kernels and refusing the non-ending
  and heap-cleanup crash cases as in `blastp` (line `PS-SW`).

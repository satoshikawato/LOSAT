- A run that sets none of LOSAT's diagnostic variables (`LOSAT_X_STATS`, `LOSAT_TIMING`, `LOSAT_X_REACH` and the
  others that `LOSAT/src/utils/xexit.rs` lists) no longer runs the code of their reports at its end. That code was the
  last a run touched, so small runs of `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `rpsbproc` and
  `makeblastdb` now peak 0.5 to 1.2 MiB (3.6 to 14 %) lower than on dev of 2026-11-28: `rpsbproc` of a one-protein
  archive 8.7 -> 7.6 MiB, `makeblastdb` of a 6 kb sequence 8.5 -> 7.3 MiB, `blastp` of one protein against 20
  14.1 -> 13.6 MiB, with 0.02 to 1.3 % fewer instructions. Runs whose peak comes earlier keep it. Outputs do not
  change, and each variable still writes its report.

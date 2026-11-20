- `PD-NCBI-PLATFORM-VARIANCE` 1.4 records the C `qsort` of the oracle hosts: glibc's is a stable merge
  sort in every release, the TBLASTX linking tie class of `s_RevCompareHSPsTbx` is source-underdetermined
  (LOSAT sorts it stably), and macOS / Windows `qsort` is not a merge sort. Documentation only.

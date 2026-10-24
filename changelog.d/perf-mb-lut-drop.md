- BLASTN drops a query batch's lookup tables when its preliminary stage ends, before its traceback,
  where NCBI's traceback stage drops its table (`LOSAT_X_BNLUTDROP`, in the strict set; `=0` holds them to
  the end of the batch as before). CL_PLACEHOLDER Outputs do not change.

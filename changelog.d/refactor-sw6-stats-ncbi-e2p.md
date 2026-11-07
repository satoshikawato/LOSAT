- Engine: BLASTP's composition-based E-value adjustment (`adjust_evalues_for_composition`) converts with the `stats`
  `e_to_p`/`p_to_e` (NCBI's `BLAST_Expm1`/`BLAST_Log1p`) instead of libm `exp_m1`/`ln_1p`, and
  `core/blast_stat/sum_statistics.rs` is gone. The function runs only with unified P-values, which LOSAT refuses, so no
  output changes.

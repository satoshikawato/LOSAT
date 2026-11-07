- Engine: the old statistics layer `core::blast_stat` is gone. `compute_blosum62_ideal_karlin_params` (with its
  once-per-process memo), `compute_lambda_from_score_probs` and `calc_lambda_of_score_probs` moved unchanged to
  `stats/karlin_calc.rs`, where they use `stats`' own score-frequency profile, Newton lambda and Karlin K (the same
  code); the rest of `composition.rs` and `karlin_params.rs` had no users. No output changes.

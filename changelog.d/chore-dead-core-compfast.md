- The composition adjustment's size fallbacks are an assert: `optimize_target_frequencies` asserts
  that the alphabet size is 20 (`COMPO_NUM_TRUE_AA`, as NCBI's only caller passes), so the line-by-line
  ports for a larger vector or matrix (`step_bound_reference`, `scaled_symmetric_product_a_reference`,
  `factor_ltriang_pos_def`, `solve_ltriang_pos_def`) and the flag that chose between the forms were
  deleted. Unused code in `core` (the old `blast_filter` module, one BLASTNA encoder), `utils` (DUST and
  SEG helpers), `report` (four pairwise shims), BLASTP (the test-only twins of the redo and alignment-view
  builders, a query-frame wrapper) and BLASTX (a redo wrapper, the no-observer shims of the HSP stream, an
  unused list copy) was removed or moved under `#[cfg(test)]`. Outputs do not change.

//! TBLASTN CLI, genetic-code, and internal seed-search boundary.
mod args;
mod scoring;
// NCBI c++/src/algo/blast/core/blast_engine.c:804-844:
// for (context=first_context; context<=last_context; context++) {
//     status = s_BlastSearchEngineOneContext(...);
// }
mod search_seed;
// NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
// status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
// Blast_InitHitListSortByScore(init_hitlist);
mod search_init;
// NCBI c++/src/algo/blast/core/blast_engine.c:522-552:
// aux_struct->GetGappedScore(..., init_hitlist, &hsp_list, ...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(...);
// Blast_HSPListSortByScore(hsp_list);
mod search_gapped;
// NCBI c++/src/algo/blast/core/blast_setup.c:729-735,770-847:
// translated db_length /= 3 before per-context length adjustment and search space.
mod stage_d_stats;
// NCBI c++/src/algo/blast/core/link_hsps.c:1602-1810;
// blast_engine.c:870-906: link the preliminary list and reap by prelim E-value.
mod stage_d_linking;
// NCBI core/blast_kappa.c:1896-1957 redoes translated-subject alignments.
mod kappa;
// NCBI composition_adjustment/compo_heap.c:252-275,330-391,439-466:
// Kappa result records enter the NCBI comparator heap and pop in its order.
mod kappa_heap;
// NCBI core/blast_kappa.c:2494-2515 and blast_hits.c:3243-3297,3420-3437:
// popped lists enter Blast_HitListUpdate, then results are reversed.
mod stage_d_results;
pub use args::TblastnArgs;

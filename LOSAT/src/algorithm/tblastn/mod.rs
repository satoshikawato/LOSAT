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
pub use args::TblastnArgs;

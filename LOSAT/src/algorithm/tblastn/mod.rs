//! TBLASTN CLI, genetic-code, and internal seed-search boundary.
mod args;
mod scoring;
// NCBI c++/src/algo/blast/core/blast_engine.c:804-844:
// for (context=first_context; context<=last_context; context++) {
//     status = s_BlastSearchEngineOneContext(...);
// }
mod search_seed;
pub use args::TblastnArgs;

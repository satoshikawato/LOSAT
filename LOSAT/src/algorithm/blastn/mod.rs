pub mod alignment;
pub mod args;
pub mod blast_engine;
pub mod blast_extend;
pub mod constants;
pub mod coordination;
pub mod extension;
pub mod filtering;
pub mod interval_tree;
pub mod lookup;
pub mod ncbi_cutoffs;
pub mod sequence_compare;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
// Blast_HSPReevaluateWithAmbiguitiesGapped(...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
// Blast_IntervalTreeReset(tree);
// ```
pub mod tracing;
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-166
// ```c
// typedef struct BlastHSP {
//    Int4 score;
//    double evalue;
//    BlastSeg query;
//    BlastSeg subject;
//    Int4 context;
// } BlastHSP;
// typedef struct BlastHSPList {
//    Int4 oid;
//    Int4 query_index;
//    BlastHSP** hsp_array;
// } BlastHSPList;
// ```
pub mod hsp;

pub use args::BlastnArgs;
pub use blast_engine::run;
#[cfg(target_arch = "wasm32")]
pub use blast_engine::run_web_pair;

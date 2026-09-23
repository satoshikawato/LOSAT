//! NCBI genetic-code singleton adapter.

// NCBI reference: c++/src/algo/blast/core/gencode_singleton.c:65-69
// ```c
// Uint1* GenCodeSingletonFind(Uint4 gen_code_id) {
//     return DynamicSGenCodeNodeArray_Find(g_theInstance, gen_code_id);
// }
// ```
// One Rust table owner serves both the core and translated-search callers.
pub use crate::utils::genetic_code::GeneticCode;

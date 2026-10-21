pub mod alignment;
pub mod args;
pub mod blast_engine;
pub mod blast_extend;
pub(crate) mod blast_xml;
pub mod constants;
pub mod coordination;
pub mod disc_lookup;
pub mod extension;
pub mod filtering;
pub mod interval_tree;
pub mod lookup;
pub mod ncbi_cutoffs;
pub(crate) mod pairwise;
pub mod query_split;
pub mod scoring;
pub mod sequence_compare;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
// Blast_HSPReevaluateWithAmbiguitiesGapped(...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
// Blast_IntervalTreeReset(tree);
// ```
pub mod tracing;
// NCBI c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520: the subjects' search forms made
// from the scope. EXPERIMENT (LOSAT_X_BNSTORE): the subjects held packed, their letters made when
// read.
pub(crate) mod x_subject_store;
// NCBI c++/src/algo/blast/api/dust_filter.cpp:92-102,171-180: the DUST masks of a query batch, one
// fresh masker per query. EXPERIMENT (LOSAT_X_QDUSTPAR): made on the search pool, in query order.
pub(crate) mod x_query_dust_par;
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
pub(crate) mod seqalign;

pub use args::BlastnArgs;
pub use blast_engine::run;
// EXPERIMENT (LOSAT_X_BNSTORE): the summary line `main` prints at exit. No NCBI counterpart: it
// prints counters and changes no value NCBI computes.
pub use x_subject_store::print_summary as x_subject_store_print_summary;
// EXPERIMENT (LOSAT_X_QDUSTPARSHADOW): the summary line `main` prints at exit. No NCBI counterpart:
// it prints counters and changes no value NCBI computes.
pub use x_query_dust_par::print_summary as x_query_dust_par_print_summary;

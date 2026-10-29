pub mod alignment;
pub(crate) mod archive;
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
// NCBI c++/src/algo/blast/api/winmask_filter.cpp:280-313: the query masks of `-window_masker_db`.
pub(crate) mod window_masker;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
// Blast_HSPReevaluateWithAmbiguitiesGapped(...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
// Blast_IntervalTreeReset(tree);
// ```
pub mod tracing;
// NCBI c++/src/algo/blast/core/blast_hits.c:676-735, blast_parameters.c:346-374: the steps of an
// ungapped (`-ungapped`) search on a subject's HSP list and the word cutoff of its batch.
pub(crate) mod ungapped;
// NCBI c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520: the subjects' search forms made
// from the scope. EXPERIMENT (LOSAT_X_BNSTORE): the subjects held packed, their letters made when
// read.
pub(crate) mod x_subject_store;
// NCBI c++/src/algo/blast/api/dust_filter.cpp:92-102,171-180: the DUST masks of a query batch, one
// fresh masker per query. EXPERIMENT (LOSAT_X_QDUSTPAR): made on the search pool, in query order.
pub(crate) mod x_query_dust_par;
// NCBI c++/src/algo/blast/core/na_ungapped.c:149-242 (`s_NuclUngappedExtendExact`). EXPERIMENT
// (LOSAT_X_BNUNGAPCOUNT): the same loops with their step counts worked out first.
pub(crate) mod x_ungap_count;
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
pub(crate) mod sam;
pub(crate) mod seqalign;
// NCBI c++/src/algo/blast/format/blast_format.cpp:828-832: the tabular rows of each alignment, one
// query after another. EXPERIMENT (LOSAT_X_BNTABPAR): a report's queries formatted on the search pool.
pub(crate) mod x_tab_par;

pub use args::BlastnArgs;
pub use blast_engine::run;
// EXPERIMENT (LOSAT_X_BNSTORE): the summary line `main` prints at exit. No NCBI counterpart: it
// prints counters and changes no value NCBI computes.
pub use x_subject_store::print_summary as x_subject_store_print_summary;
// EXPERIMENT (LOSAT_X_QDUSTPAR): the summary line `main` prints at exit. No NCBI counterpart:
// it prints counters and changes no value NCBI computes.
pub use x_query_dust_par::print_summary as x_query_dust_par_print_summary;
// EXPERIMENT (LOSAT_X_GREEDYROWFUSE): the summary line `main` prints at exit with LOSAT_X_STATS. No
// NCBI counterpart: it prints which build of the greedy row code ran and changes no value NCBI
// computes.
pub use alignment::greedy::x_greedy_row_print_summary;

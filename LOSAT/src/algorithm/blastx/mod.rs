//! BLASTX native serial local-subject preparation, search and reporting.
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
pub mod args;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-326
pub mod input;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:577-587
pub mod query_setup;
pub(crate) mod reader_batch;
pub use args::BlastxArgs;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:960-985
pub mod parameters;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234
pub mod seed;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-545
pub mod preliminary;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:50-62
pub mod split;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:235-271
pub mod search;

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1901-1931
#[cfg(test)]
mod stage_c_tests;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
pub mod linking;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
pub mod statistics;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
pub mod kappa;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
// The redo-stage subject SEG memoised per subject.
mod x_seg_memo;

// `-lcase_masking`: the word scan of a `-subject` with lower-case letters over its scan ranges,
// and its lower case in the outfmt 0 rows.
mod x_soft_subject;

// Diagnostic (LOSAT_TIMING): the stage times of a BLASTX run. No NCBI counterpart.
mod x_timing;

// `-db` read in blocks where it is searched, instead of a copy of
// the whole database.
pub(crate) mod x_db_stream;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1508
pub mod traceback;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-114
pub mod results;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
pub mod runtime;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:386-394
mod culling;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329,377-385
// `-best_hit_overhang` / `-best_hit_score_edge` (series NP, O6a-3): the best hit writer and pipe.
mod besthit;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1777-1784

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:279-297
pub mod report;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1549-1649 (BlastHitList2SeqAlign_OMF): the
// Seq-aligns of each query's final hit list (outfmt 8, 12).
mod anchored;
mod archive;
// The streams of a library call's hold, which the report streams of `cli::NcbiOstream` name.
pub(crate) use blast_xml::HeldFormatWriter;
mod blast_xml;
// The reports of the ungapped HSPs of length 0 that `-window_size 0` finds (NCBI's formatters'
// exceptions, the Seq-loc mapper's rows of `-subject`), whose exceptions TBLASTN's reports share.
pub(crate) mod empty_rows;
mod seqalign;
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp (CBlastTabularInfo): the tabular reports with a
// field list or a delimiter, and outfmt 10 and 20, through the shared writer (`report::tabular`).
mod tabular;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-291
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
pub mod native;
// LOSAT Web's entry (web ABI v2): one search, every requested output format.
mod local;
pub use local::{check_options, run_local, run_local_database};
/// The default size of a query batch, in letters (`GetQueryBatchSize`).
pub(crate) use query_setup::BATCH_SIZE as QUERY_BATCH_SIZE;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod web;
// The subjects of a search by oid: `-subject` records or the stored sequences of `-db`.
pub(crate) mod subjects;
// The database path (`-db`): opening, the stored sequences, titles and ranks per batch.
pub(crate) mod database;

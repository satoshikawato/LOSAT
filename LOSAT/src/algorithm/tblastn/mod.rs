//! TBLASTN CLI, genetic-code, and internal seed-search boundary.
mod args;
// NCBI c++/src/algo/blast/core/blast_engine.c:804-844:
// for (context=first_context; context<=last_context; context++) {
//     status = s_BlastSearchEngineOneContext(...);
// }
mod search_seed;
// EXPERIMENT (LOSAT_X_TBNSCANBF): the lookup scan without a branch on the presence vector.
pub(crate) mod x_scan_bf;
// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the lookup scan on presence bytes, a block of indexes at once.
mod x_scan_byte;
// NCBI c++/src/algo/blast/core/blast_util.c:1072,1084-1086:
// translation_table = s_BlastGetTranslationTable(genetic_code, FALSE);
// length = BLAST_TranslateCompressedSequence(translation_table, nucl_length, nucl_seq, frame, ...);
// The subject frames read from tables.
mod x_subject_translate;
// EXPERIMENT (LOSAT_X_TBNSTORE): the subjects held packed, their letters made when read.
mod x_subject_store;
// NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
// status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
// Blast_InitHitListSortByScore(init_hitlist);
mod search_init;
// NCBI c++/src/algo/blast/core/aa_ungapped.c:516-552: the two-hit tests of every hit.
// EXPERIMENT (LOSAT_X_TBNHITBF): too far and overlap without a branch.
mod x_tbn_hitbf;
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
// NCBI c++/src/algo/blast/core/blast_kappa.c:2352-2390,2418-2479 constructs redo inputs
// from the initial score block, extension options and hit-saving parameters.
mod stage_d_kappa_params;
// NCBI core/blast_kappa.c:1896-1957 redoes translated-subject alignments.
mod kappa;
// NCBI composition_adjustment/compo_heap.c:252-275,330-391,439-466:
// Kappa result records enter the NCBI comparator heap and pop in its order.
mod kappa_heap;
// NCBI core/blast_kappa.c:2494-2515 and blast_hits.c:3243-3297,3420-3437:
// popped lists enter Blast_HitListUpdate, then results are reversed.
mod stage_d_results;
// NCBI c++/src/algo/blast/core/blast_engine.c:870-905;
// blast_traceback.c:1481-1499,1763-1776: preliminary lists enter Kappa,
// then local results pass through the post-traceback pipes.
mod stage_d_pipeline;
// NCBI c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:333-344: each search thread takes the
// next chunk of subjects when it is free.
// Graduated 2026-10-22 (was LOSAT_X_TBNOIDNEXT): the chunks of subjects are handed to threads
// when they are free.
#[cfg(feature = "parallel")]
mod x_oid_flow;
// NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
// PrintOneResultSet dispatches the retained search result to its output format.
mod stage_e_report;
// TBLASTN's BLAST XML (outfmt 5), once per run, from the final hit lists.
mod blast_xml;
// NCBI c++/src/algo/blast/api/blast_seqalign.cpp:1549-1649 (BlastHitList2SeqAlign_OMF): the
// Seq-aligns of each query's final hit list (outfmt 8, 12).
mod seqalign;
// NCBI c++/src/app/blast/tblastn_app.cpp:188-194 and seqsrc_seqdb.cpp: the search of a BLAST
// nucleotide database (`-db`): its sequences, numbers and titles.
mod database;
pub use args::{check_options, run_local, ResolvedTblastnArgs, TblastnArgs};
// The summary lines `main` prints at exit (`LOSAT_X_TBNKEEPFEWSHADOW`, `LOSAT_X_TBNXLATESIMD`).
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub use x_subject_translate::print_summary as x_subject_translate_print_summary;
// EXPERIMENT (LOSAT_X_TBNSTORE): the summary line `main` prints at exit. No NCBI counterpart: it
// prints counters and changes no value NCBI computes.
pub use x_subject_store::print_summary as x_subject_store_print_summary;
// EXPERIMENT (LOSAT_X_TBNSCANBF): the summary line `main` prints at exit. No NCBI counterpart: it
// prints counters and changes no value NCBI computes.
pub use x_scan_bf::print_summary as x_scan_bf_print_summary;
// EXPERIMENT (LOSAT_X_TBNSCANBYTE / LOSAT_X_TBNSCANBYTESHADOW): the summary lines `main` prints at
// exit. No NCBI counterpart: they print counters and change no value NCBI computes.
pub use x_scan_byte::print_summary as x_scan_byte_print_summary;
// EXPERIMENT (LOSAT_X_TBNLINKWIN): the summary line `main` prints at exit. No NCBI counterpart: it
// prints counters and changes no value NCBI computes.
pub use stage_d_linking::x_link_window_print_summary;
// EXPERIMENT (LOSAT_X_TBNHITBF / LOSAT_X_TBNHITBFSHADOW): the summary line `main` prints at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub use x_tbn_hitbf::print_summary as x_tbn_hitbf_print_summary;

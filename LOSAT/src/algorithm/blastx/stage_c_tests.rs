//! Focused Stage C boundary checks, source-derived expectations.
use super::{
    args::ResolvedOptions,
    input::{read_fasta, FastaRecord},
    parameters::{effective_lengths, score_block},
    query_setup::{prepare_queries, PreparedQueryBatch},
    split::{calculate_num_chunks, split_queries, validate_chunk_size},
};
use crate::cli::Commands;
use std::path::Path;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
// ```c++
//                                   "Translated Query-Protein Subject BLAST"));
//     const bool kQueryIsProtein = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastx";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("blastx-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
//
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// ```
pub(super) fn options() -> ResolvedOptions {
    let cli: crate::cli::Cli = crate::cli::try_parse_from([
        "LOSAT",
        "blastx",
        "-query",
        "unused.fna",
        "-subject",
        "unused.faa",
        "-seg",
        "no",
        "-comp_based_stats",
        "0",
    ])
    .unwrap();
    let Commands::Blastx(args) = cli.command else {
        unreachable!()
    };
    args.resolve().unwrap()
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-264
// ```c++
// 	    BLAST_PROF_START( APP.LOOP.PRE );
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
// ```
pub(super) fn fixture(name: &str) -> (Vec<FastaRecord>, PreparedQueryBatch, ResolvedOptions) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../docs/evidence/losatx_stage_a/inputs/generated")
        .join(name);
    let records = read_fasta(&path, false, false).unwrap();
    let o = options();
    let batch = prepare_queries(&records, &o).unwrap();
    (records, batch, o)
}
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1901-1925 CalculateNumberChunks
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1927-1931 InvalidChunkSizeBlastx
// ```c++
// BOOST_AUTO_TEST_CASE(CalculateNumberChunks)
// {
//     EBlastProgramType program = eBlastTypeBlastx;
//     size_t chunk_size = 10002;
//     Uint4 retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size, 10240000, 1);
//     BOOST_REQUIRE_EQUAL(1055, retval);
//
//     retval = SplitQuery_CalculateNumChunks(eBlastTypeBlastx,
//                        &chunk_size, chunk_size/2, 1);
//
//     BOOST_REQUIRE_EQUAL(1, retval);
//
//     retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size,
//                        3*chunk_size-2*SplitQuery_GetOverlapChunkSize(program), 1);
//
//     BOOST_REQUIRE_EQUAL(3, retval);
//
//     retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size,
//                        1+2*chunk_size+SplitQuery_GetOverlapChunkSize(program), 1);
//
//     BOOST_REQUIRE_EQUAL(2, retval);
// }
//
// BOOST_AUTO_TEST_CASE(InvalidChunkSizeBlastx)
// {
//     CAutoEnvironmentVariable tmp_env("CHUNK_SIZE", "40000");
//     BOOST_REQUIRE_THROW(SplitQuery_GetChunkSize(blast::eBlastx), CBlastException);
// }
// ```
#[test]
fn active_n04_mutable_chunk_size_and_invalid_divisibility() {
    let mut size = 10002;
    assert_eq!(calculate_num_chunks(&mut size, 10240000, 1), 1055);
    let length = size / 2;
    assert_eq!(calculate_num_chunks(&mut size, length, 1), 1);
    let length = 3 * size - 2 * 297;
    assert_eq!(calculate_num_chunks(&mut size, length, 1), 3);
    let length = 1 + 2 * size + 297;
    assert_eq!(calculate_num_chunks(&mut size, length, 1), 2);
    assert!(validate_chunk_size(40000).is_err());
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:748-791
// ```c++
//             // The corrections for the contexts corresponding to the plus
//             // strand are always the same, so only calculate the first one
//             if (s_IsPlusStrand(chunk_qinfo[chunk_num], ctx) &&
//                 (ctx % NUM_FRAMES == 1 || ctx % NUM_FRAMES == 2)) {
//                 correction = m_SplitBlk->GetContextOffsets(chunk_num).back();
//                 goto error_check;
//             }
//
//             // If the query length is divisible by CODON_LENGTH, the
//             // corrections for all contexts corresponding to a given strand are
//             // the same, so only calculate the first one
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 0) &&
//                 (ctx % NUM_FRAMES != 0) && (ctx % NUM_FRAMES != 3)) {
//                 correction = m_SplitBlk->GetContextOffsets(chunk_num).back();
//                 goto error_check;
//             }
//
//             // If the query length % CODON_LENGTH == 1, the corrections for the
//             // first two contexts of the negative strand are the same, and the
//             // correction for the last context is one more than that.
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 1) &&
//                 !s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
//
//                 if (ctx % NUM_FRAMES == 4) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back();
//                     goto error_check;
//                 } else if (ctx % NUM_FRAMES == 5) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back() + 1;
//                     goto error_check;
//                 }
//             }
//
//             // If the query length % CODON_LENGTH == 2, the corrections for the
//             // last two contexts of the negative strand are the same, which is
//             // one more that the first context on the negative strand.
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 2) &&
//                 !s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
//
//                 if (ctx % NUM_FRAMES == 4) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back() + 1;
//                     goto error_check;
// ```
#[test]
fn split_modulo_frame_maps_and_shared_positive_corrections() {
    for (len, map, negative) in [
        (19409, [0, 1, 2, 3, 4, 5], [0, 0, 0]),
        (19410, [0, 1, 2, 3, 4, 5], [3136, 3136, 3136]),
        (19411, [0, 1, 2, 4, 5, 3], [3136, 3136, 3137]),
    ] {
        let (records, full, o) = fixture(&format!("S12_split{len}.fna"));
        let chunks = split_queries(&records, &full, &o).unwrap();
        assert_eq!(chunks.len(), if len < 19410 { 1 } else { 2 });
        assert_eq!(chunks[0].absolute_contexts, map);
        assert_eq!(&chunks[0].corrections[3..], negative);
        if chunks.len() > 1 {
            assert_eq!(chunks[1].range, 9705..len);
            assert_eq!(chunks[1].corrections, [3235, 3235, 3235, 0, 0, 0]);
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:76-100
// ```c++
//                        size_t num_queries)
// {
//     // TODO: need to model mem usage and when it's advantageous to split
//     bool retval = true;
//
//     if (program == eBlastTypeMapping) {
//         return false;
//     }
//
//     // if ((concatenated_query_length <= chunk_size+SplitQuery_GetOverlapChunkSize(program)) ||
//    //  if ((concatenated_query_length <= chunk_size) ||
//         // do not split RPS-BLAST
//     if (Blast_SubjectIsPssm(program) ||
//         // the current implementation does NOT support splitting for multiple
//         // blastx queries, loop over queries individually here...
//         (program == eBlastTypeBlastx && num_queries > 1) ||
//         Blast_ProgramIsPhiBlast(program)) {
//         retval = false;
//     }
//
//     return retval;
// }
//
// Uint4
// SplitQuery_CalculateNumChunks(EBlastProgramType program,
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:56-61
// ```c++
//     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
//     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
//     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
//         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
//     /* No split for ungapped mode JIRA SB-1082 */
//     if (!options->GetGappedMode()) m_NumChunks = 1;
// ```
#[test]
fn multiple_queries_and_ungapped_do_not_split() {
    let (mut records, full, mut o) = fixture("S12_split19411.fna");
    o.gapped = false;
    assert_eq!(split_queries(&records, &full, &o).unwrap().len(), 1);
    o.gapped = true;
    records.push(records[0].clone());
    let full = prepare_queries(&records, &o).unwrap();
    assert_eq!(split_queries(&records, &full, &o).unwrap().len(), 1);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:776-785
// ```c++
//
//       /* Effective search space for a given sequence/strand/frame */
//       Int8 effective_search_space =
//           s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
//                                               blast_message);
//
//       kbp = kbp_ptr[index];
//
//       if (query_info->contexts[index].is_valid &&
//           ((query_length = query_info->contexts[index].query_length) > 0) ) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:846-847
// ```c++
//       query_info->contexts[index].eff_searchsp = effective_search_space;
//       query_info->contexts[index].length_adjustment = length_adjustment;
// ```
#[test]
fn invalid_statistics_keep_fixed_space_with_zero_adjustment() {
    let (_, full, o) = fixture("S06_seed.fna");
    let mut p = score_block(&full);
    p[0].valid = false;
    let spaces = vec![37; full.contexts.len()];
    effective_lengths(&full, &mut p, &o, 83, 1, Some(&spaces)).unwrap();
    assert_eq!(p[0].search_space, 37);
    assert_eq!(p[0].length_adjustment, 0);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2781-2787
// ```c++
//       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
//       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
//       if (loop_status) {
//           contexts[context].is_valid = FALSE;
//           sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
//           sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
//           if (!Blast_QueryIsTranslated(program) ) {
// ```
#[test]
fn statistical_validity_is_distinct_from_preparation_validity() {
    let (_, mut full, _) = fixture("S06_seed.fna");
    assert!(full.contexts.iter().all(|c| c.is_valid));
    full.sequence_start.fill(0);
    let p = score_block(&full);
    assert!(p.iter().all(|p| !p.valid));
    assert!(full.contexts.iter().all(|c| c.is_valid));
}
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:282-292 BackboneSequenceTest
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:294-304 SmallboneSequenceTest
// ```c++
// BOOST_AUTO_TEST_CASE(BackboneSequenceTest) {
//   // create a trivial sequence
//   Int4 len = 65534; // 65535 is the maximum possible unsigned short
//   GetSeqBlk(len);
//   FillLookupTable();
//   BOOST_REQUIRE_EQUAL( lookup->bone_type, eBackbone );
//   Int4 num_used = ((AaLookupBackboneCell *)(lookup->thick_backbone))[0].num_used;
//   BOOST_REQUIRE_EQUAL(num_used, len-2);
//   Int4 offset = ((Int4 *)(lookup->overflow))[num_used-1];
//   BOOST_REQUIRE_EQUAL(offset, len-3);
// }
//
// BOOST_AUTO_TEST_CASE(SmallboneSequenceTest) {
//   // create a trivial sequence
//   Int4 len = 65533;
//   GetSeqBlk(len);
//   FillLookupTable();
//   BOOST_REQUIRE_EQUAL( lookup->bone_type, eSmallbone );
//   Int4 num_used = ((AaLookupSmallboneCell *)(lookup->thick_backbone))[0].num_used;
//   BOOST_REQUIRE_EQUAL(num_used, len-2);
//   Int4 offset = ((Uint2 *)(lookup->overflow))[num_used-1];
// ```
#[test]
fn active_ncbi_lookup_overflow_u16_boundary_preserves_all_offsets() {
    use crate::algorithm::tblastx::lookup::{build_ncbi_lookup_from_prepared, QueryContext};
    for len in [65533usize, 65534] {
        let c = QueryContext {
            q_idx: 0,
            f_idx: 0,
            frame: 1,
            aa_seq: vec![0; len + 2],
            aa_seq_nomask: None,
            aa_len: len,
            orig_len: len * 3,
            frame_base: 0,
            is_valid: true,
            karlin_params: Default::default(),
        };
        let lookup = build_ncbi_lookup_from_prepared(
            vec![0; len + 1],
            vec![(0, len as i32 - 1)],
            vec![c],
            0,
        );
        let hits = lookup.get_hits(0);
        assert_eq!(hits.len(), len - 2);
        assert_eq!(hits.last(), Some(&(len as i32 - 3)));
        assert!(hits.iter().copied().eq(0..len as i32 - 2));
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
// ```c++
// /** Split subject sequences if longer than this */
// #define MAX_DBSEQ_LEN 5000000
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
// ```c++
//     if (backup->offset + MAX_DBSEQ_LEN <
//         backup->hard_ranges[backup->hm_index].right) {
//
//         subject->length = MAX_DBSEQ_LEN;
//         backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
// ```
#[test]
fn unported_long_subject_chunking_is_rejected_before_search() {
    let (records, _, o) = fixture("S06_seed.fna");
    let subject = FastaRecord {
        internal_id: "Subject_1".into(),
        title: "long".into(),
        sequence: vec![b'W'; 5_000_001],
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
    };
    let error = super::search::search_preliminary(&records, &[subject], &o)
        .err()
        .expect("unported subject splitting");
    assert_eq!(
        error.to_string(),
        "BLASTX preliminary protein-subject splitting above 5000000 residues is not implemented"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:233-239
// ```c++
//     for (index=0; index<(*seq_info)->GetNumSeqs(); ++index)
//         retval = MIN(retval, (*seq_info)->GetSeqBlk(index)->length);
//
//     if(retval < BLAST_SEQSRC_MINLENGTH)
// 	retval = BLAST_SEQSRC_MINLENGTH;
//
//     return retval;
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205-205
// ```c++
// #define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-987
// ```c++
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        if (Blast_SubjectIsTranslated(program_number)) {
//            min_subject_length/=3;
//        }
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
//
//    if ((status = BlastHitSavingParametersNew(program_number, hit_options, sbp,
// 		                                     query_info, min_subject_length,
// 		                                     (*ext_params)->options->compositionBasedStats,
// 		                                     hit_params)) != 0){
//        return status;
//    }
// ```
#[test]
fn short_subject_cutoff_uses_adapter_minimum_ten() {
    let o = options();
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../docs/evidence/losatx_stage_c/inputs/compressed_overflow.fna");
    let records = read_fasta(&path, false, false).unwrap();
    let full = prepare_queries(&records, &o).unwrap();
    let mut p = score_block(&full);
    super::parameters::subject_parameters(&full, &mut p, &o, 8, 1, 8, None).unwrap();
    // Pinned NCBI CALL0 context0: qlen6468, subject8, min length10 -> max cutoff18.
    assert_eq!(p[0].hit_cutoff_max, 18);
    assert_eq!(p[0].word_cutoff, 18);
}

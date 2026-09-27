//! Comparison-only diagnostic: internal BLASTX preliminary state, no report output.
use anyhow::{bail, Result};
use LOSAT::algorithm::blastx::{input::read_fasta, search::search_preliminary};
use LOSAT::cli::{Cli, Commands};
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
// ```c++
//         // Restore the full query sequence for the traceback stage!
//         if (m_InternalData->m_Queries == NULL) {
//             CRef<ILocalQueryData> query_data
//                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
//             // Query masking info is calculated as a side-effect
//             CBlastScoreBlk sbp
//                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
//                                         NULL, m_Messages, NULL, NULL));
//             m_InternalData->m_Queries = query_data->GetSequenceBlk();
//         }
// ```
fn hex(data: &[u8]) -> String {
    data.iter().map(|b| format!("{b:02x}")).collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
// ```c++
//         // Restore the full query sequence for the traceback stage!
//         if (m_InternalData->m_Queries == NULL) {
//             CRef<ILocalQueryData> query_data
//                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
//             // Query masking info is calculated as a side-effect
//             CBlastScoreBlk sbp
//                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
//                                         NULL, m_Messages, NULL, NULL));
//             m_InternalData->m_Queries = query_data->GetSequenceBlk();
//         }
// ```
fn masks(data: &[(i32, i32)]) -> String {
    if data.is_empty() {
        "-".to_string()
    } else {
        data.iter()
            .map(|(a, b)| format!("{a}:{b}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4071-4077
// ```c++
//              status = Blast_HSPInit(gap_align->query_start,
//                            gap_align->query_stop, gap_align->subject_start,
//                            gap_align->subject_stop,
//                            init_hsp->offsets.qs_offsets.q_off,
//                            init_hsp->offsets.qs_offsets.s_off, context,
//                            query_frame, subject->frame, gap_align->score,
//                            &(gap_align->edit_script), &new_hsp);
// ```
fn hsp(
    stage: &str,
    call: usize,
    i: usize,
    h: &LOSAT::algorithm::blastx::preliminary::PreliminaryHsp,
) {
    println!(
        "{stage}\t{call}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        h.context,
        h.score,
        h.frame,
        h.q_start,
        h.q_end,
        h.q_gapped_start,
        h.s_start,
        h.s_end,
        h.s_gapped_start
    );
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-526
// ```c++
//         BlastInitHitListReset(init_hitlist);
//
//         if (aux_struct->WordFinder) {
//             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                                    word_params, aux_struct->ewp,
//                                    aux_struct->offset_pairs,
//                                    kScanSubjectOffsetArraySize,
//                                    init_hitlist, ungapped_stats);
//
//             if (init_hitlist->total == 0) continue;
//         }
//
//         if (score_options->gapped_calculation) {
//             Int4 prot_length = 0;
//             if (score_options->is_ooframe) {
//                 /* Convert query offsets in all HSPs into the mixed-frame
//                    coordinates */
//                 s_TranslateHSPsToDNAPCoord(program_number, init_hitlist,
//                        query_info, subject->frame, orig_length, backup.offset);
//                 if (kTranslatedSubject) {
//                     prot_length = subject->length;
//                     subject->length = orig_length;
//                 }
//             }
//         /** NB: If queries are concatenated, HSP offsets must be adjusted
//           * inside the following function call, so coordinates are
//           * relative to the individual contexts (i.e. queries, strands or
//           * frames). Contexts should also be filled in HSPs when they
//           * are saved.
//           */
//         /* fence_hit is null, since this is only for prelim stage. */
//         if (aux_struct->GetGappedScore) {
//             status = aux_struct->GetGappedScore(program_number, query,
//                     query_info,
//                     subject, gap_align, score_params, ext_params, hit_params,
//                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
// ```
fn main() -> Result<()> {
    let cli: Cli = LOSAT::cli::try_parse_from(std::env::args_os())?;
    let Commands::Blastx(args) = cli.command else {
        bail!("blastx required");
    };
    let options = args.resolve()?;
    let subjects = read_fasta(&args.subject, true, args.lcase_masking)?;
    let queries = read_fasta(&args.query, false, args.lcase_masking)?;
    let batches = search_preliminary(&queries, &subjects, &options)?;
    for batch in batches {
        println!("SPLIT\t{}\t{}", batch.query_ordinal, batch.chunks.len());
        for (k, stage) in batch.chunks.iter().enumerate() {
            if batch.chunks.len() > 1 {
                println!(
                    "CHUNK\t{k}\t{}\t{}",
                    stage.chunk.range.start, stage.chunk.range.end
                );
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:277-288
                // ```c++
                //             } catch (const CBlastException& e) {
                //                 // This error message is safe to ignore for a given chunk,
                //                 // because the chunks might end up producing a region of
                //                 // the query for which ungapped Karlin-Altschul blocks
                //                 // cannot be calculated
                //                 const string err_msg1("search cannot proceed due to errors "
                //                                      "in all contexts/frames of query "
                //                                      "sequences");
                //                 const string err_msg2(kBlastErrMsg_CantCalculateUngappedKAParams);
                //                 if (e.GetMsg().find(err_msg1) == NPOS && e.GetMsg().find(err_msg2) == NPOS) {
                //                     throw;
                //                 }
                // ```
                println!("CHUNK_STATE\t{k}\t{}", i32::from(stage.statistically_valid));
                for c in 0..stage.chunk.absolute_contexts.len() {
                    println!(
                        "MAP\t{k}\t{c}\t{}\t{}",
                        stage.chunk.absolute_contexts[c], stage.chunk.corrections[c]
                    );
                }
            }
            for s in &stage.subjects {
                println!(
                    "CALL\t{}\t{}\t{}",
                    s.call,
                    s.oid,
                    subjects[s.oid].sequence.len()
                );
                for c in stage.chunk.prepared.first_context..s.parameters.len() {
                    let x = &stage.chunk.prepared.contexts[c];
                    let p = &s.parameters[c];
                    println!(
                        "PARAM\t{}\t{c}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        s.call,
                        x.offset,
                        x.length,
                        x.frame,
                        i32::from(p.valid),
                        p.search_space,
                        p.length_adjustment,
                        p.word_cutoff,
                        p.word_xdrop_init,
                        p.word_xdrop
                    );
                }
                for &(q, sub) in &s.seeds {
                    println!("SEED\t{}\t{q}\t{sub}", s.call);
                }
                for (i, h) in s.initial.iter().enumerate() {
                    println!(
                        "INIT\t{}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}",
                        s.call, h.q_seed, h.s_seed, h.q_start, h.s_start, h.length, h.score
                    );
                }
                println!("WORD_EXIT\t{}\t0\t{}", s.call, s.initial.len());
                if !s.initial.is_empty() {
                    if options.gapped {
                        for c in stage.chunk.prepared.first_context..s.parameters.len() {
                            let p = &s.parameters[c];
                            println!(
                                "HIT_PARAM\t{}\t{c}\t{}\t{}\t{}",
                                s.call, p.hit_cutoff, p.hit_cutoff_max, s.gap_xdrop
                            );
                        }
                    }
                    let raw_stage = if options.gapped { "GAPPED" } else { "UNGAPPED" };
                    println!("{raw_stage}_COUNT\t{}\t{}", s.call, s.raw.len());
                    for (i, h) in s.raw.iter().enumerate() {
                        hsp(raw_stage, s.call, i, h);
                    }
                    if options.gapped {
                        println!("PURGED_COUNT\t{}\t{}", s.call, s.purged.len());
                        for (i, h) in s.purged.iter().enumerate() {
                            hsp("PURGED", s.call, i, h);
                        }
                    }
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
        // ```c++
        //         // Restore the full query sequence for the traceback stage!
        //         if (m_InternalData->m_Queries == NULL) {
        //             CRef<ILocalQueryData> query_data
        //                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
        //             // Query masking info is calculated as a side-effect
        //             CBlastScoreBlk sbp
        //                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
        //                                         NULL, m_Messages, NULL, NULL));
        //             m_InternalData->m_Queries = query_data->GetSequenceBlk();
        //         }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:836-847
        // ```c++
        //         	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
        //
        //         	 // Just in case effective_db_length < 0
        //         	 if (effective_db_length <= 0)
        //         		 effective_db_length = 1;
        //
        //              effective_search_space = effective_db_length *
        //                              (query_length - length_adjustment);
        //          }
        //       }
        //       query_info->contexts[index].eff_searchsp = effective_search_space;
        //       query_info->contexts[index].length_adjustment = length_adjustment;
        // ```
        for (c, x) in batch.restored.contexts.iter().enumerate() {
            let p = &batch.full_parameters[c];
            println!(
                "FULL_CONTEXT\t{}\t{c}\t{}\t{}\t{}\t{}\t{}\t{}",
                batch.query_ordinal,
                x.offset,
                x.length,
                x.frame,
                i32::from(x.is_valid),
                p.search_space,
                p.length_adjustment
            );
        }
        for (q, &len) in batch.restored.original_lengths.iter().enumerate() {
            println!("RESTORED_DNA\t{}\t{len}", batch.query_ordinal + q);
        }
        println!(
            "RESTORED_BUFFER\t{}\t{}\t{}",
            batch.query_ordinal,
            hex(&batch.restored.sequence_start),
            hex(&batch.restored.sequence_start_nomask)
        );
        for (c, x) in batch.restored.contexts.iter().enumerate() {
            println!(
                "RESTORED_MASK\t{}\t{c}\t{}",
                batch.query_ordinal,
                masks(&x.dna_masks)
            );
        }
        for (q, by_oid) in batch.merged_by_query.iter().enumerate() {
            for (oid, hsps) in by_oid.iter().enumerate() {
                println!(
                    "REJOIN_COUNT\t{}\t{oid}\t{}",
                    batch.query_ordinal + q,
                    hsps.len()
                );
                for (i, h) in hsps.iter().enumerate() {
                    hsp(&format!("REJOIN\t{}", batch.query_ordinal + q), oid, i, h);
                }
            }
        }
    }
    Ok(())
}

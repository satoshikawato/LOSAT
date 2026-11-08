//! psiblast's rounds of FASTA queries (`-num_iterations` other than 1; series PS, item PS5a): the
//! round driver (`search_rounds`), which `search` calls for such a run, and the hook of one
//! round's search (`round_outcome`), which `run_resolved_in_pool` calls beside `psiblast_pssms`.
//! Both are cold: blastp's default path only tests `ResolvedBlastpArgs::psi_rounds` and
//! `BlastpReports::psi_round`.
//!
//! The driver runs `run_resolved_in_pool` once per round per query, as NCBI's `CPsiBlast` searches
//! one query per batch and round (psiblast_app.cpp:468, `x_RunLocalPsiBlastIterations`; quoted at
//! `algorithm::psiblast::rounds`), on one search pool and the database opened once: round 1 is
//! blastp's search of the FASTA query, a later round the PSSM search (`algorithm::pssm_search`) of
//! the PSSM that the engine made from the round before (`PsiSearch::from_engine`). The outfmt 0
//! prolog comes with the first call; every call writes no epilog, which the driver writes after
//! the last query (outfmt 7's count is that of the rounds). The PSSM files are written between
//! the calls, where NCBI writes them, after the report of the round before.

use std::io::Write as _;
use std::sync::Arc;

use anyhow::{bail, Result};

use super::{
    blastp_database_name, run_resolved_in_pool, write_pairwise_prologs,
    write_reports_without_queries, BlastpRanges, BlastpReports, BlastpSubjectSource, ShownLetters,
    ShownSubjects, SubjectSet, NCBI_BLASTP_VERSION,
};
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::common::database::HeaderReadError;
use crate::algorithm::psiblast::pssm::{
    ncbistdaa_letters, query_pssm, PssmAlignment, PssmAncillary, PssmQueryRecord, PssmScoring,
    QueryPssm,
};
use crate::algorithm::psiblast::rounds::{
    InclusionSet, IterationState, PsiRoundIo, PsiRoundPlan, RoundFiles, RoundOutcome, SetKey,
};
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::app::CompositionMode;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::query_batch::QueryBatches;
use crate::blastinput::seq_range;
use crate::core::psi::KarlinBlk;
use crate::report::db_subject::ParsedSubject;
use crate::stats::KarlinParams;

/// The run of a psiblast search of FASTA queries with `-num_iterations` other than 1, after the
/// query batches are read (`search`): LOSAT's refusal of a query that NCBI splits, the rounds of
/// each searched batch's query, the epilog, then the stop of the batch after them.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn search_rounds(
    args: &ResolvedBlastpArgs,
    plan: &PsiRoundPlan,
    query_records: &[FastaRecord],
    searched_queries: &[FastaRecord],
    batches: QueryBatches,
    subjects: &SubjectSet<'_>,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196 (the
    // chunks of a split query are set up without the PSSM; `in_pssm::check_in_pssm_options`)
    // A query long enough to be split would have its later rounds' chunks searched without
    // their PSSM: refused before the search, as a PSSM query of that length (DESIGN_ps5 §6).
    let sizes = crate::algorithm::blastp::protein_query_split::protein_query_split_sizes(
        args.app.name(),
        10000,
    )?;
    for &(_, query) in &batches.starts {
        let length = searched_queries[query].seq().len();
        if crate::algorithm::blastp::protein_query_split::split_protein_batch(&[length], sizes)
            .is_some()
        {
            bail!(
                "a query of {length} letters, which NCBI BLAST+ splits into chunks that its later rounds search without the PSSM, with -num_iterations {} is not supported by LOSAT's PSIBLAST",
                plan.num_iterations
            );
        }
    }
    let epilog = batches.stop.is_none();
    if batches.searched_done == 0 {
        // As `search`: the first batch's stop after the prolog, or the reports of no query.
        if batches.stop.is_some() {
            write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
        } else {
            write_reports_without_queries(args, subject_source, subject_label, outputs)?;
        }
    } else {
        crate::utils::threading::with_search_pool(args.num_threads, "blastp", |pool| {
            run_queries(
                args,
                plan,
                query_records,
                searched_queries,
                &batches,
                subjects,
                subject_source,
                subject_label,
                outputs,
                pool,
                epilog,
            )
        })?;
    }
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        return Err(error);
    }
    Ok(())
}

/// The composition-based mode of the options' reader (`in_pssm::setup_messages`).
fn app_mode(mode: BlastpCompositionMode) -> CompositionMode {
    match mode {
        BlastpCompositionMode::NoCompositionBasedStats => CompositionMode::NoCompositionBasedStats,
        BlastpCompositionMode::CompositionBasedStats => CompositionMode::CompositionBasedStats,
        BlastpCompositionMode::CompositionMatrixAdjust => CompositionMode::CompositionMatrixAdjust,
        BlastpCompositionMode::ForceFullMatrixAdjust => CompositionMode::CompoForceFullMatrixAdjust,
    }
}

/// The rounds of every searched batch's query on `pool`, then the epilog (`epilog`).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:508-532
// ```c++
//             ECompoAdjustModes comp_stats_original = opts_hndl->GetOptions().GetCompositionBasedStats();
//
//             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
// 		CRef<CBlastQueryVector> query_batch(input->GetNextSeqBatch(*scope));
//
//                 const bool converged = DoIterations(opts_hndl,
//                            query_batch,
//                            pssm,
//                            db_args,
//                            db_adapter,
//                            scope,
//                            formatter);
//
//                 if (converged && !fmt_args->HasStructuredOutputFormat() &&
//                 	fmt_args->GetFormattedOutputChoice() != CFormattingArgs::eArchiveFormat) {
//                 	out_stream << NcbiEndl << "Search has CONVERGED!" << NcbiEndl;
//             	}
//             	// Reset for next query sequence.
//             	m_AncillaryData.Reset();
//             	pssm.Reset();
//             	opts_hndl->SetOptions().SetCompositionBasedStats(comp_stats_original);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:196-205 (quoted at
// `in_pssm::setup_messages`: a PSSM search of mode 2 or 3 sets the options to 1, so a later
// round of the query neither runs 2 or 3 nor warns again)
#[allow(clippy::too_many_arguments)]
fn run_queries(
    args: &ResolvedBlastpArgs,
    plan: &PsiRoundPlan,
    query_records: &[FastaRecord],
    searched_queries: &[FastaRecord],
    batches: &QueryBatches,
    subjects: &SubjectSet<'_>,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    pool: &crate::utils::threading::SearchPool<'_>,
    epilog: bool,
) -> Result<()> {
    let mut files = RoundFiles::new(
        &plan.request,
        args.out.as_deref().filter(|path| path.as_os_str() != "-"),
    );
    let original_mode = args.comp_based_stats.mode;
    let pssm_mode =
        crate::algorithm::psiblast::in_pssm::position_based_mode(app_mode(original_mode));
    let pssm_round_mode = match pssm_mode {
        CompositionMode::NoCompositionBasedStats => BlastpCompositionMode::NoCompositionBasedStats,
        _ => BlastpCompositionMode::CompositionBasedStats,
    };
    // The rounds formatted so far (`m_QueriesFormatted`: each round's `PrintOneResultSet`).
    let mut calls = 0usize;
    for &(input_start, query) in &batches.starts {
        let records = std::slice::from_ref(&searched_queries[query]);
        let query_input = seq_range::QueryInput::whole(records);
        let ranges = BlastpRanges {
            input_records: std::slice::from_ref(&query_records[input_start]),
            query_input: &query_input,
            subject_placements: &subjects.placements,
            x_ordinal_base: input_start,
        };
        let mut state = IterationState::new(plan.num_iterations);
        // The PSSM that the next round searches (`pssm`, made after the round before).
        let mut pssm: Option<QueryPssm> = None;
        while state.go() {
            if let Some(current) = &pssm {
                files.write(current)?;
            }
            let round = state.number();
            let mut round_args = args.clone();
            if let Some(current) = &pssm {
                // The options' mode when this round's score block is set up: the given one in
                // the query's first PSSM round, 1 (or 0) after it.
                let setup_mode = if round == 2 {
                    app_mode(original_mode)
                } else {
                    pssm_mode
                };
                let mut search =
                    PsiSearch::from_engine(&current.engine, records[0].clone(), Vec::new());
                search.setup_messages = crate::algorithm::psiblast::in_pssm::setup_messages(
                    &search,
                    setup_mode,
                    &args.matrix_name,
                );
                round_args.psi_search = Some(Arc::new(search));
                round_args.comp_based_stats.mode = pssm_round_mode;
            }
            let io = PsiRoundIo {
                state: state.clone(),
                save_last: plan.save_last,
                request: plan.request.clone(),
                outcome: std::sync::Mutex::new(None),
            };
            // The reader's messages of the batch come with its first round.
            let reader_warnings: &[Vec<u8>] = if round == 1 {
                &batches.warnings[query..query + 1]
            } else {
                &[]
            };
            let result = run_resolved_in_pool(
                round_args,
                records,
                subject_source,
                "",
                subject_label,
                outputs,
                pool,
                None,
                &ranges,
                &BlastpReports {
                    reader_warnings,
                    identity_warned: 0,
                    epilog: false,
                    x_prolog: calls == 0,
                    x_queries_before: 0,
                    x_xml2_carry: None,
                    psi_round: Some(&io),
                },
                super::x_lookup_ahead::XAheadBlock::none(),
            );
            // A header file that NCBI cannot read fails the round whose results reach it first,
            // at its prefetch, after the reports of the rounds before (`search`'s rule for the
            // batches): the outfmt 0 prolog if nothing was written, and the reader's messages
            // of a batch whose first round fails.
            if let Err(error) = result {
                if error.is::<HeaderReadError>() {
                    let failure = error
                        .downcast::<HeaderReadError>()
                        .expect("the error was checked to be a HeaderReadError");
                    if calls == 0 {
                        write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
                    }
                    if round == 1 {
                        outputs.diagnostics.write_all(&batches.warnings[query])?;
                    }
                    return Err(failure.error);
                }
                return Err(error);
            }
            calls += 1;
            let outcome = io
                .outcome
                .lock()
                .expect("the round's outcome")
                .take()
                .expect("every round's search reaches its outcome hook");
            if !outcome.has_alignments {
                break;
            }
            state.advance(outcome.set);
            if state.go() || plan.save_last {
                // The PSSM engine's error comes here, after the round's report (no epilog).
                pssm = Some(
                    outcome
                        .next
                        .expect("the hook makes the PSSM where NCBI makes it")?,
                );
            }
        }
        if plan.save_last {
            if let Some(current) = &pssm {
                files.write(current)?;
            }
        }
        if state.has_converged() {
            write_converged(outputs)?;
        }
    }
    if epilog {
        write_round_epilogs(args, subject_source, subject_label, outputs, calls)?;
    }
    Ok(())
}

/// `out_stream << NcbiEndl << "Search has CONVERGED!" << NcbiEndl` for every report that is not
/// structured (psiblast's outfmt 0, 6, 7, 10 and 20).
fn write_converged(outputs: &mut ReportOutputs<'_>) -> Result<()> {
    use crate::blastinput::app::ReportFormat;
    for format in outputs.formats.iter_mut() {
        let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
        if matches!(
            ReportFormat::from_number(number),
            Some(
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::CommaSeparatedValues
                    | ReportFormat::CommaSeparatedValuesWithHeader
            )
        ) {
            let mut writer = format.sink.open()?;
            writer.write_all(b"\nSearch has CONVERGED!\n")?;
            writer.flush()?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2238 (quoted at
// `tabular::write_tabular_without_queries`: outfmt 7's count of the result sets formatted)
/// The epilog of each report after the last round: outfmt 0's, outfmt 7's `# BLAST processed R
/// queries` (R the rounds formatted, `m_QueriesFormatted`), nothing for outfmt 6, 10 and 20.
fn write_round_epilogs(
    args: &ResolvedBlastpArgs,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    rounds: usize,
) -> Result<()> {
    let (num_sequences, total_letters) = subject_source.shown_totals();
    let pairwise_database = subject_source
        .database()
        .map(super::database::BlastpDatabase::pairwise_database);
    let taxdb_found = match subject_source.database() {
        Some(database) => database.taxdb_found(),
        None => args.subject_taxdb_found,
    };
    for format in outputs.formats.iter_mut() {
        if let Some(spec) = super::super::tabular::custom_tabular_spec(format.outfmt) {
            let mut writer = format.sink.open()?;
            super::super::tabular::write_tabular_round_epilog(
                spec,
                taxdb_found,
                args.app.ncbi_name(),
                &mut writer,
                &mut *outputs.diagnostics,
                rounds,
            )?;
            writer.flush()?;
            continue;
        }
        match crate::blastinput::app::parse_formatting_string(format.outfmt)?.number {
            0 => {}
            7 => {
                // The default outfmt 7 (`write_blastp_tabular_hit_lists`' epilog).
                let mut writer = format.sink.open()?;
                writeln!(writer, "# BLAST processed {rounds} queries")?;
                writer.flush()?;
                continue;
            }
            _ => continue,
        }
        let report = crate::report::BlastpPairwiseReport {
            version: NCBI_BLASTP_VERSION.to_string(),
            database_name: blastp_database_name(args, subject_label),
            database_num_sequences: num_sequences,
            database_total_letters: total_letters,
            matrix_name: args.matrix_name.clone(),
            gap_open: args.scoring.gap_open,
            gap_extend: args.scoring.gap_extend,
            word_threshold: args.threshold,
            window_size: args.window_size as i32,
            gapped_karlin: crate::stats::lookup_protein_params(&args.scoring),
            gumbel: None,
            num_descriptions: args.num_descriptions,
            num_alignments: args.num_alignments,
        };
        let mut writer = format.sink.open()?;
        crate::report::pairwise::psiblast_round::write_psiblast_epilog(
            &mut *writer,
            &report,
            pairwise_database.as_ref(),
        )?;
    }
    Ok(())
}

/// A Karlin-Altschul block of the PSSM engine (`core::psi`).
fn engine_block(params: &KarlinParams) -> KarlinBlk {
    KarlinBlk {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
    }
}

/// The outcome of one round's search (`PsiRoundIo::outcome`), from the round's hit list of its
/// one query before the formatter's cut (`psiblast_pssms`' Seq-aligns): whether it has
/// alignments, its inclusion set, and where NCBI makes one the next round's PSSM (an error of
/// the engine is kept for after the round's report). Returns, when the previous round's set is
/// not empty, whether each subject of the round (by subject index, a database's by rank) is in
/// it, for the two-part description table (`psiblast_round::ReportRound`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:529-566
// ```c++
// double GetLowestEvalue(const objects::CDense_seg::TScores& scores,
//                        double* bit_score /* = NULL */)
// {
//     double retval = BLAST_EXPECT_VALUE;
//     ...
//     ITERATE(CDense_seg::TScores, i, scores) {
//         if ( (tmp = s_GetEvalue(**i)) < retval) {
//             retval = tmp;
//         }
//     ...
// void
// CPsiBlastAlignmentProcessor::operator()
//     (const objects::CSeq_align_set& alignments,
//      double evalue_inclusion_threshold,
//      THitIdentifiers& output)
// {
//     output.clear();
//
//     ITERATE(CSeq_align_set::Tdata, hsp, alignments.Get()) {
//         // Look for HSP with score less than inclusion_ethresh
//         double e = GetLowestEvalue((*hsp)->GetScore());
//         if (e < evalue_inclusion_threshold) {
//             CSeq_id_Handle sid =
//                 CSeq_id_Handle::GetHandle((*hsp)->GetSeq_id(1));
//             output.insert(sid);
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:158
// ```c
// #define BLAST_EXPECT_VALUE 10.0 /**< by default, alignments whose expect
// ```
// The Seq-id of a subject is its key (`rounds::SetKey`): a `-subject` record's index
// (`lcl|Subject_N`), a database sequence's ordinal id (`gnl|BL_ORD_ID|oid`, the hit lists name it
// by its rank among the listed ones), or its parsed Seq-align id. The PSSM's multiple alignment
// keeps its rows by subject (`pssm::input_msa`): two Seq-aligns of one id there are copies of one
// record that the engine purges (line DBID admits shared ids for such copies only).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:483-509 (quoted at
// `core::psi::PsiScoreBlk::set_ungapped_statistical_params`: the round's PSI blocks)
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn round_outcome(
    args: &ResolvedBlastpArgs,
    io: &PsiRoundIo,
    hit_lists: &[Option<BlastpHitList>],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
    listed_ordinals: Option<&[u32]>,
    parsed: &[Option<Arc<ParsedSubject>>],
) -> Option<Vec<bool>> {
    const BLAST_EXPECT_VALUE: f64 = 10.0;
    let key = |s_idx: u32| -> SetKey {
        match listed_ordinals {
            None => SetKey::Index(s_idx),
            Some(ordinals) => match parsed.get(s_idx as usize).and_then(Option::as_deref) {
                Some(subject) => {
                    let mut id = Vec::new();
                    crate::seqid::label::fasta(subject.alignment_id(), &mut id);
                    SetKey::Id(id)
                }
                None => SetKey::Index(ordinals[s_idx as usize]),
            },
        }
    };
    let hsps: Vec<&crate::algorithm::blastp::hsp::BlastpHsp> = hit_lists
        .first()
        .and_then(Option::as_ref)
        .map(|hit_list| {
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            hit_list.hsplist_array[..count]
                .iter()
                .flat_map(|hsp_list| &hsp_list.hsps)
                .collect()
        })
        .unwrap_or_default();
    let has_alignments = !hsps.is_empty();
    let lowest_evalue = |e_value: f64| e_value.min(BLAST_EXPECT_VALUE);
    let set: InclusionSet = hsps
        .iter()
        .filter(|hsp| {
            lowest_evalue(crate::report::serial::seqalign::seqalign_evalue(
                hsp.e_value,
            )) < args.inclusion_ethresh
        })
        .map(|hsp| key(hsp.s_idx))
        .collect();
    let previous = io.state.previously_found();
    let found_before = (!previous.is_empty()).then(|| {
        let subjects = hsps
            .iter()
            .map(|hsp| hsp.s_idx as usize + 1)
            .max()
            .unwrap_or(0);
        let mut found = vec![false; subjects];
        for hsp in &hsps {
            found[hsp.s_idx as usize] = previous.contains(&key(hsp.s_idx));
        }
        found
    });
    // NCBI makes the next PSSM after `Advance` when another round follows or the last PSSM is
    // saved (`itr || GetSaveLastPssm()`), and after a round with alignments only.
    let makes_next = if has_alignments {
        let mut after = io.state.clone();
        after.advance(set.clone());
        after.go() || io.save_last
    } else {
        false
    };
    let next =
        makes_next.then(|| next_pssm(args, io, &hsps, ranges, shown_subjects, contexts, gapped));
    *io.outcome.lock().expect("the round's outcome") = Some(RoundOutcome {
        has_alignments,
        set,
        next,
    });
    found_before
}

/// `ComputePssmForNextIteration` of the round's Seq-aligns `hsps` (`psiblast_pssms`' PSSM of the
/// FASTA query, its Seq-align id and title, which every round keeps) with the round's ancillary
/// data: round 1's standard blocks, or the PSSM's `kbp_psi` and `kbp_gap_psi` (and `kbp_gap`,
/// which is `kbp_gap_psi`) after a PSSM round.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:167-178
// ```c++
// CPsiBlastApp::ComputePssmForNextIteration(const CBioseq& bioseq,
//                               CConstRef<CSeq_align_set> sset,
//                               CConstRef<CPSIBlastOptionsHandle> opts_handle,
//                               CRef<CScope> scope,
//                               CRef<CBlastAncillaryData> ancillary_data)
// {
//     CPSIDiagnosticsRequest
//         diags(PSIDiagnosticsRequestNewEx(m_CmdLineArgs->SaveAsciiPssm()));
//     m_AncillaryData = ancillary_data;
//     return PsiBlastComputePssmFromAlignment(bioseq, sset, scope, *opts_handle,
//                                             m_AncillaryData, diags);
// }
// ```
fn next_pssm(
    args: &ResolvedBlastpArgs,
    io: &PsiRoundIo,
    hsps: &[&crate::algorithm::blastp::hsp::BlastpHsp],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
) -> Result<QueryPssm> {
    let subject = |s_idx: u32| -> Vec<u8> {
        match shown_subjects.letters(s_idx as usize) {
            Some(ShownLetters::Ascii(letters)) => ncbistdaa_letters(letters),
            Some(ShownLetters::Stored(residues)) => residues.to_vec(),
            None => Vec::new(),
        }
    };
    let scoring = PssmScoring {
        matrix: args.scoring.matrix,
        matrix_name: &args.matrix_name,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
    };
    let q_shift = ranges.query_input.placements.offset(0);
    let alignments: Vec<PssmAlignment<'_>> = hsps
        .iter()
        .map(|hsp| PssmAlignment {
            subject: hsp.s_idx,
            e_value: crate::report::serial::seqalign::seqalign_evalue(hsp.e_value),
            q_start: hsp.q_start - 1 + q_shift,
            s_start: hsp.s_start - 1,
            q_length: hsp.q_end + 1 - hsp.q_start,
            ops: hsp.gap_info.as_deref(),
        })
        .collect();
    let ordinal = ranges.query_input.ordinal(0);
    let record = &ranges.input_records[ordinal];
    let id = format!(
        "Query_{}",
        args.first_query_ordinal + ranges.x_ordinal_base + ordinal + 1
    );
    let standard_ungapped = engine_block(&contexts[0].karlin_params);
    let table_gapped = engine_block(gapped);
    let ancillary = match args.psi_search.as_deref() {
        None => PssmAncillary {
            ungapped: standard_ungapped,
            gapped: table_gapped,
            psi_ungapped: standard_ungapped,
            psi_gapped: table_gapped,
        },
        Some(psi) => {
            let karlin = super::pssm_hooks::psi_karlin(psi, contexts, gapped);
            let psi_gapped = engine_block(&karlin.gapped);
            PssmAncillary {
                ungapped: standard_ungapped,
                gapped: psi_gapped,
                psi_ungapped: engine_block(&karlin.ungapped),
                psi_gapped,
            }
        }
    };
    query_pssm(
        &io.request,
        scoring,
        table_gapped,
        args.inclusion_ethresh,
        PssmQueryRecord {
            letters: record.seq(),
            id: &id,
            title: record.title(),
        },
        &alignments,
        &subject,
        ancillary,
    )
}

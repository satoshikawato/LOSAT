//! LOSAT Web's BLASTX (web ABI v2): one search over the records that the host registered,
//! written in every requested output format (`run_local`). The query batches are searched and
//! reported by the command line's loop (`native::run_query_batches`).

use super::{
    args::{BlastxArgs, ResolvedOptions, SemanticOptionsError},
    input::FastaRecord,
    native::{
        engine, input_error, options_error, output_error, report_error, run_query_batches,
        write_empty_subject_warnings,
    },
    query_setup::BATCH_SIZE,
    report::{self, LocalReport, RecordsParsed, ShownSubjects},
    runtime::BatchResults,
    subjects::Subjects,
};
use crate::api::local_blast::{FormatProbe, ReportOutputs};
use crate::blastinput::app::{self, Program};
use crate::blastinput::fasta_reader;
use anyhow::{bail, ensure, Context, Result};
use std::io::{self, Write};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// ```c++
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats are several
// CBlastFormat printers over the same result set.
/// Runs one BLASTX search over already parsed records and writes every requested output format
/// from the same result (the entry of web ABI v2; the command line is `native::run`), and gives
/// the caller the final hit list (`ReportOutputs::hits`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the options and
/// of the subjects come in the command line's order, before the search. The `-query` and
/// `-subject` values of `args` are used only as display names. LOSAT Web writes outfmt 0, 6 and 7
/// (without a custom field list or delimiter); the other formats, `-db` and `-parse_deflines` are
/// refused.
///
/// The records are those of NCBI's reader (`fasta_reader`): the queries as its nucleotide reader
/// reads them, the subjects as its protein reader reads them (`blastx_record`). Their messages are
/// written where BLASTX's reader writes them: a subject's as the subjects are read, a query's
/// before the report of the batch that holds it. Queries without a record are NCBI's
/// `Query is Empty!`.
pub fn run_local(
    args: BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    subject_records: &[fasta_reader::FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let requested: Vec<&str> = outputs.formats.iter().map(|format| format.outfmt).collect();
    let (choices, resolutions) = read_options(&args, &requested)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let Some(subject_path) = args.subject.as_deref() else {
        return Err(app::missing_subject_error());
    };
    let ReportOutputs {
        formats,
        diagnostics,
        hits: hits_sink,
        observer,
    } = outputs;
    // NCBI reads the subjects before the queries, each record's messages as it is read.
    let mut subjects = Vec::with_capacity(subject_records.len());
    for record in subject_records {
        subjects.push(blastx_record(record, args.lcase_masking).map_err(input_error)?);
        diagnostics.write_all(&record.warnings)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c++
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    if subjects.is_empty() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    let resolved = check_resolutions(&args, &choices, resolutions, &mut **diagnostics)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:213-217
    // ```c++
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    if query_records.is_empty() {
        diagnostics.write_all(b"Warning: [blastx] Query is Empty!\n")?;
        return Ok(());
    }
    // The letters of a reader batch (`native::run_inner`).
    let batch_size = app::query_batch_size("BLASTX", BATCH_SIZE as u32)? as usize;
    write_empty_subject_warnings(&subjects, &mut **diagnostics)?;
    // The environment of a search of subject records, after the subjects' warnings
    // (`native::run_inner`); none of the formats written here is the archive.
    app::check_search_environment("BLASTX", true, false)?;
    let records_parsed = RecordsParsed::default();
    let shown = ShownSubjects::Fasta {
        path: subject_path,
        parsed: &records_parsed,
    };
    let mut writers = formats
        .iter_mut()
        .map(|format| format.sink.open())
        .collect::<io::Result<Vec<_>>>()?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
    // ```c++
    //         formatter.PrintProlog();
    // ```
    for (index, (writer, options)) in writers.iter_mut().zip(&resolved).enumerate() {
        let mut quiet = io::sink();
        let warnings: &mut dyn Write = if index == 0 {
            &mut **diagnostics
        } else {
            &mut quiet
        };
        report::write_prolog(writer, warnings, options, &subjects, &shown).map_err(report_error)?;
        writer.flush().map_err(output_error)?;
    }
    // NCBI: a batch size of 0 makes the first batch empty, and its query factory refuses it.
    if batch_size == 0 {
        return Err(engine("Empty CBlastQueryVector"));
    }
    // The search's options are those of every format; the formats differ only in the counts of
    // the pairwise report (`run_report_counts`).
    let search_options = &resolved[0];
    let lowercase = args.lcase_masking;
    // The reader batches as BLASTX's reader makes them (`input::emit_fasta_batch`): a record's
    // messages, then the record, and the batch once it holds `batch_size` letters or more; the
    // rest at the end. Only the batch being read is kept.
    let mut batch = Vec::new();
    let mut read = |deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
                    warn: &mut dyn FnMut(&[u8]) -> Result<()>|
     -> Result<Vec<FastaRecord>> {
        let mut letters = 0;
        for record in query_records {
            let query = blastx_record(record, lowercase)?;
            if !record.warnings.is_empty() {
                warn(&record.warnings)?;
            }
            letters += query.sequence.len();
            batch.push(query);
            if letters >= batch_size {
                deliver(&batch)?;
                batch.clear();
                letters = 0;
            }
        }
        if !batch.is_empty() {
            deliver(&batch)?;
        }
        Ok(Vec::new())
    };
    // The errors of a `-subject` search (`native::run_inner`).
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        // LOSAT's refusal, not an error of NCBI's engine.
        if error.is::<super::parameters::SearchSpaceOverflow>() {
            error
        } else {
            engine(error)
        }
    };
    let made_hits = hits_sink.is_some() || choices.iter().any(|choice| choice.number == 0);
    let mut run_hits = Vec::new();
    let mut first_hsp = 0;
    let mut render = |diagnostics: &mut dyn Write,
                      queries: &[FastaRecord],
                      results: &mut Vec<BatchResults>,
                      queries_before: usize|
     -> Result<()> {
        // Every format and the hit list read the lists in the order of the reports.
        report::sort_lists(results);
        let hits = match made_hits {
            true => Some(
                report::pairwise_hits(search_options, queries, &subjects, &shown, results)
                    .map_err(report_error)?,
            ),
            false => None,
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
        // ```c++
        //             	ITERATE(CSearchResultSet, result, *results) {
        //                	    formatter.PrintOneResultSet(**result, query_batch);
        // ```
        for (index, (writer, options)) in writers.iter_mut().zip(&resolved).enumerate() {
            // The warnings of the batch's queries are written once, with the first format.
            let mut quiet = io::sink();
            let warnings: &mut dyn Write = if index == 0 {
                &mut *diagnostics
            } else {
                &mut quiet
            };
            let mut probe = observer
                .as_deref_mut()
                .map(|observer| FormatProbe::new(observer, index));
            report::render_queries_with(
                writer,
                warnings,
                options,
                queries,
                &subjects,
                &shown,
                results,
                queries_before,
                None,
                Some(LocalReport {
                    probe: probe.as_mut(),
                    first_hsp,
                    hits: hits.as_deref(),
                }),
            )
            .map_err(report_error)?;
            writer.flush().map_err(output_error)?;
        }
        first_hsp += match &hits {
            Some(hits) => hits.len(),
            None => report::reported_hsps(results, search_options.gapped),
        };
        // The run's hit list numbers the queries of the run.
        if hits_sink.is_some() {
            for mut hit in hits.into_iter().flatten() {
                hit.hit.q_idx = u32::try_from(queries_before + hit.hit.q_idx as usize)?;
                run_hits.push(hit);
            }
        }
        Ok(())
    };
    let processed = run_query_batches(
        search_options,
        Subjects::Fasta(&subjects),
        false,
        &mut read,
        &mut **diagnostics,
        &search_error,
        &mut render,
    )
    .map_err(input_error)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    for (index, (writer, options)) in writers.iter_mut().zip(&resolved).enumerate() {
        let mut quiet = io::sink();
        let warnings: &mut dyn Write = if index == 0 {
            &mut **diagnostics
        } else {
            &mut quiet
        };
        report::write_epilog(
            writer, warnings, options, &subjects, &shown, processed, None,
        )
        .map_err(report_error)?;
        writer.flush().map_err(output_error)?;
    }
    if let Some(sink) = hits_sink {
        sink(std::sync::Arc::new(run_hits));
    }
    Ok(())
}

/// The checks of the options alone, without the records (the `validate` of web ABI v2), for the
/// output formats `outfmts` of one `run_local` call: `run_local`'s refusals and LOSAT's, the
/// thread count, then NCBI's errors of the options, with the messages of a run. A run reports
/// NCBI's errors after the subjects are read, so for an argv and records that both fail the run
/// reports the records' error.
pub fn check_options(args: &BlastxArgs, outfmts: &[&str]) -> Result<()> {
    let (choices, resolutions) = read_options(args, outfmts)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    check_resolutions(args, &choices, resolutions, &mut io::sink()).map(|_| ())
}

/// The output formats `outfmts` as `run_local` reads them and the options of each, before the
/// records are read: LOSAT Web's refusals and LOSAT's fail here; NCBI's errors of the options
/// (`SemanticOptionsError`) stay in the resolutions, which the command line reports after the
/// subjects are read (`native::run_inner`).
fn read_options(
    args: &BlastxArgs,
    outfmts: &[&str],
) -> Result<(Vec<app::FormatChoice>, Vec<Result<ResolvedOptions>>)> {
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // NCBI parses -outfmt before its option handlers read the subjects (`native::run_inner`).
    let mut choices = Vec::with_capacity(outfmts.len());
    for outfmt in outfmts {
        let choice = app::parse_formatting_string(outfmt)?;
        app::report_format(&choice, Program::Blastx, None)?;
        // No NCBI counterpart: the formats of LOSAT Web (`docs/web/abi_v2.md`).
        if !matches!(choice.number, 0 | 6 | 7) {
            bail!(
                "output format {} is not supported by LOSAT Web's BLASTX",
                choice.number
            );
        }
        if choice.custom() {
            bail!(
                "a custom output format specification (fields or a delimiter) is not supported by LOSAT Web's BLASTX"
            );
        }
        choices.push(choice);
    }
    ensure!(!choices.is_empty(), "no output format was requested");
    // The records keep no ids read from their deflines (the host registers them before it knows
    // the options).
    if args.parse_deflines {
        bail!("-parse_deflines is not supported by LOSAT Web yet");
    }
    // The command line's ranges come to LOSAT Web with line SX's session SXd3.
    if args.query_loc.is_some() || args.subject_loc.is_some() {
        bail!("-query_loc and -subject_loc are not supported by LOSAT Web's BLASTX yet");
    }
    // The options of each format: LOSAT's refusals now, NCBI's errors of the options after the
    // subjects are read (`native::run_inner`).
    let mut resolutions = Vec::with_capacity(choices.len());
    for choice in &choices {
        match args.resolve_with_format(choice) {
            Err(error) if error.downcast_ref::<SemanticOptionsError>().is_none() => {
                return Err(error)
            }
            resolution => resolutions.push(resolution),
        }
    }
    Ok((choices, resolutions))
}

/// The checks of the options that the command line makes after the subjects are read, in its
/// order: the formatting handler's, the hit list size (its warnings and those of an ignored
/// -num_descriptions, the first format's, to `diagnostics`; formats that need different hit lists
/// are refused, `run_report_counts`), then NCBI's errors of each format's options.
fn check_resolutions(
    args: &BlastxArgs,
    choices: &[app::FormatChoice],
    resolutions: Vec<Result<ResolvedOptions>>,
    diagnostics: &mut dyn Write,
) -> Result<Vec<ResolvedOptions>> {
    for choice in choices {
        app::formatting_handler_check(choice, Program::Blastx)?;
    }
    let outfmts: Vec<i32> = choices.iter().map(|choice| choice.number).collect();
    let (_, warnings) = crate::blastinput::report_counts::run_report_counts(
        "blastx",
        &outfmts,
        args.max_target_seqs.map(|n| n as usize),
        args.num_descriptions,
        args.num_alignments,
    )?;
    diagnostics.write_all(&warnings)?;
    resolutions
        .into_iter()
        .map(|resolution| resolution.map_err(options_error))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
// ```c++
//     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
// }
//
// void CFastaReader::x_CloseMask(void)
// {
//     _ASSERT(m_MaskRangeStart != kInvalidSeqPos);
//     m_CurrentMask->SetPacked_int().AddInterval
//         (GetBestID(), m_MaskRangeStart, GetCurrentPos(ePosWithGapsAndSegs) - 1,
//          eNa_strand_plus);
// ```
/// A record of NCBI's reader (`fasta_reader`, the records that LOSAT Web registers) as BLASTX's
/// reader (`input.rs`) reads the same input: the residues in upper case and, with
/// `-lcase_masking` (`lowercase`), the lowercase mask's intervals, which `fasta_reader` keeps as
/// the case of the residues (a `>?` gap inside a masked interval, too); a title that is not UTF-8
/// is BLASTX's reader's error. The reader's messages stay with the caller's record.
fn blastx_record(record: &fasta_reader::FastaRecord, lowercase: bool) -> Result<FastaRecord> {
    let title = String::from_utf8(record.title.clone()).context("invalid UTF-8 FASTA title")?;
    if record.sequence.len() > i32::MAX as usize {
        bail!("BLASTX sequence exceeds Int4 length");
    }
    let mut sequence = record.sequence.clone();
    let mut lowercase_masks: Vec<(i32, i32)> = Vec::new();
    for (offset, residue) in sequence.iter_mut().enumerate() {
        if !residue.is_ascii_lowercase() {
            continue;
        }
        residue.make_ascii_uppercase();
        if lowercase {
            // The length is at most `i32::MAX` (above).
            let offset = offset as i32;
            match lowercase_masks.last_mut().filter(|r| r.1 + 1 == offset) {
                Some(last) => last.1 = offset,
                None => lowercase_masks.push((offset, offset)),
            }
        }
    }
    Ok(FastaRecord {
        internal_id: record.local_id.clone(),
        title,
        sequence,
        lowercase_masks,
        warnings: Vec::new(),
        source: None,
    })
}

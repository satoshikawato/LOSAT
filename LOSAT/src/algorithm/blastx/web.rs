//! Existing reactor ABI adapter; the search, batching and reports use BLASTX owners.
use super::{args::BlastxArgs, input, native, query_setup::BATCH_SIZE, report, runtime};
use anyhow::{bail, Result};
use std::{cell::RefCell, io::Write, path::Path};

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
// The existing CLI v2 parser owns all defaults, lexical/scope validation and
// code32 rejection. Synthetic paths identify memory input; no file is opened.
pub(crate) fn parse_args(extra: &[&str]) -> Result<BlastxArgs> {
    // The CLI's `-db` (a database search) waits for stage DB-7 in LOSAT Web, as BLASTP's.
    if extra
        .iter()
        .any(|key| key.split('=').next().unwrap_or(key) == "-db")
    {
        bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    if extra.iter().any(|key| {
        matches!(
            key.split('=').next().unwrap_or(key),
            "-query" | "-subject" | "-out" | "-help" | "--help" | "--version"
        )
    }) {
        bail!("file/help arguments are unsupported for the BLASTX web API");
    }
    let mut argv = vec!["losat", "blastx", "-query", "query", "-subject", "subject"];
    argv.extend_from_slice(extra);
    // ABI v1 is frozen (plan TD-1): its `-outfmt` keeps the values and the parser errors of
    // the CLI's BLASTX parser before it read `-outfmt` as NCBI (`v1_outfmt`).
    use clap::CommandFactory;
    let command = crate::cli::Cli::command().mut_subcommand("blastx", |blastx| {
        blastx.mut_arg("outfmt", |outfmt| outfmt.value_parser(v1_outfmt))
    });
    match crate::cli::try_parse_from_command::<crate::cli::Cli, _, _>(command, argv)?.command {
        crate::cli::Commands::Blastx(args) => Ok(args),
        _ => unreachable!("fixed BLASTX subcommand"),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2827-2845
// ```c++
//         }
//         int val = 0;
//         try { val = NStr::StringToInt(fmt_choice); }
//         catch (const CStringException&) {   // probably a conversion error
//             CNcbiOstrstream os;
//             os << "'" << fmt_choice << "' is not a valid output format";
//             string msg = CNcbiOstrstreamToString(os);
//             NCBI_THROW(CInputException, eInvalidInput, msg);
//         }
//         if (val < 0 || val >= static_cast<int>(eEndValue)) {
//             string msg("Formatting choice is out of range");
//             throw std::out_of_range(msg);
//         }
//         if (m_IsIgBlast && (val != 3 && val != 4 && val != 7 && val != eAirrRearrangement)) {
//             string msg("Formatting choice is not valid");
//             throw std::out_of_range(msg);
//         }
//         fmt_type = static_cast<EOutputFormat>(val);
//         if ( !(fmt_type == eTabular ||
// ```
/// ABI v1's reading of a BLASTX `-outfmt` value, which is frozen (plan TD-1): the values
/// and the errors of the CLI's BLASTX parser before it read `-outfmt` as NCBI.
fn v1_outfmt(value: &str) -> std::result::Result<String, String> {
    let mut tokens = value.split_whitespace();
    let f = tokens.next().ok_or("empty outfmt")?;
    if !matches!(f, "0" | "6" | "7") {
        return Err("unsupported BLASTX outfmt: use 0, 6 or 7".into());
    }
    let fields: Vec<_> = tokens.collect();
    if f == "0" && !fields.is_empty() {
        return Err("custom fields require outfmt 6/7".into());
    }
    for field in fields {
        if field != "std" && !super::args::FIELDS.contains(&field) {
            return Err(format!("unsupported BLASTX outfmt field '{field}'"));
        }
    }
    Ok(value.to_string())
}

/// The options of ABI v1's BLASTX: the format number and the fields as ABI v1 read them
/// (`v1_outfmt`: the words between white space), before BLASTX read `-outfmt` as NCBI.
fn v1_resolve(args: &BlastxArgs) -> Result<super::args::ResolvedOptions> {
    let mut words = args.outfmt.split_whitespace();
    let outfmt = words.next().unwrap_or_default().parse::<u8>()?;
    args.resolve_with_format(outfmt, words)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:302-314
// ```c++
// CBlastFastaInputSource::CBlastFastaInputSource(const string& user_input,
//                                        const CBlastInputSourceConfig& iconfig)
//     : m_Config(iconfig),
//       m_ReadProteins(iconfig.IsProteinInput())
// {
//     if (user_input.empty()) {
//         NCBI_THROW(CInputException, eEmptyUserInput,
//                    "No sequence input was provided");
//     }
//     m_LineReader.Reset(new CMemoryLineReader(user_input.c_str(),
//                                              user_input.size()));
//     x_InitInputReader();
// }
// ```
// Original bytes are reparsed on each call with this role's gencode/mask
// options; generic FASTA handles never cache translated or masked search state.
pub(crate) fn run_pair(
    args: BlastxArgs,
    query: &[u8],
    subject: &[u8],
    subject_label: &str,
) -> Result<Vec<u8>> {
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let diagnostics = RefCell::new(std::io::stderr());
    let subjects = input::parse_fasta_batches_with_warnings(
        subject,
        true,
        args.lcase_masking,
        usize::MAX,
        &mut |_| Ok(()),
        &mut |warning| native::reader_warning(warning, &mut *diagnostics.borrow_mut()),
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
    // ```c++
    //     if(hitlist_size < 5){
    //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
    //     }
    // ```
    if args.max_target_seqs.unwrap_or(500) < 5 {
        diagnostics
            .borrow_mut()
            .write_all(b"Warning: [blastx] Examining 5 or more matches is recommended\n")?;
    }
    let options = v1_resolve(&args)?;
    if query.is_empty() {
        bail!("No sequence input was provided");
    }
    if subjects.is_empty() {
        bail!("BLASTX requires a protein subject");
    }
    let mut output = Vec::new();
    let label = Path::new(subject_label);
    let shown = report::ShownSubjects::Fasta { path: label };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:784-788
    // ```c++
    //         		if(warning != kEmptyStr){
    //         			warning += ": ";
    //         		}
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    // ```
    for s in subjects.iter().filter(|s| s.sequence.is_empty()) {
        writeln!(
            diagnostics.borrow_mut(),
            "Warning: [blastx] {} {}: Subject sequence contains no data",
            s.internal_id,
            s.title
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
    // ```c++
    //         formatter.PrintProlog();
    // ```
    report::write_prolog(&mut output, &options, &subjects, &shown)?;
    let mut processed = 0;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    //             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
    // ```
    // EXPERIMENT (LOSAT_X_BXPOOL): see native.rs.
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    // ...
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ...
    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
    //                 results = lcl_blast.Run();
    // ...
    //             	ITERATE(CSearchResultSet, result, *results) {
    //                	    formatter.PrintOneResultSet(**result, query_batch);
    //             	}
    // ```
    // Same batch loop as in native.rs, as a closure over the pool its searches run on (see
    // native.rs).
    let mut x_run_batches = |x_pool: Option<&crate::utils::threading::SearchPool<'_>>| {
        // EXPERIMENT (LOSAT_X_BXBATCH): see native.rs.
        #[cfg(feature = "parallel")]
        if let Some(pool) = x_pool.filter(|pool| runtime::x_bx_batch() && pool.enabled()) {
            use rayon::prelude::*;
            enum XEvent {
                Warn(Vec<u8>),
                Batch(Vec<input::FastaRecord>),
            }
            let x_events = RefCell::new(Vec::<XEvent>::new());
            let x_failed = std::cell::Cell::new(false);
            let x_wave = pool.threads() * 4;
            let x_large = |queries: &[input::FastaRecord]| {
                queries.iter().map(|q| q.sequence.len()).sum::<usize>() > 2 * BATCH_SIZE
            };
            // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
            // ```c
            //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
            // ...
            //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
            // ...
            //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
            //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
            //                 results = lcl_blast.Run();
            // ...
            //             	ITERATE(CSearchResultSet, result, *results) {
            //                	    formatter.PrintOneResultSet(**result, query_batch);
            //             	}
            // ```
            // As in native.rs: search the queued batches, then write warnings and results in reader
            // order.
            let mut x_flush = |events: Vec<XEvent>| -> Result<()> {
                let batches: Vec<&[input::FastaRecord]> = events
                    .iter()
                    .filter_map(|event| match event {
                        XEvent::Batch(queries) => Some(queries.as_slice()),
                        XEvent::Warn(_) => None,
                    })
                    .collect();
                let searched = if batches.len() < pool.threads() {
                    batches
                        .iter()
                        .map(|queries| {
                            runtime::x_search_internal_in_pool(queries, &subjects, &options, pool)
                        })
                        .collect::<Vec<_>>()
                } else {
                    pool.install(|| {
                        batches
                            .par_iter()
                            .with_max_len(1)
                            .map(|queries| {
                                if x_large(queries) {
                                    runtime::x_search_internal_in_pool(
                                        queries, &subjects, &options, pool,
                                    )
                                } else {
                                    runtime::x_search_internal_serial(queries, &subjects, &options)
                                }
                            })
                            .collect::<Vec<_>>()
                    })
                };
                let mut searched = searched.into_iter();
                for event in &events {
                    match event {
                        XEvent::Warn(bytes) => diagnostics.borrow_mut().write_all(bytes)?,
                        XEvent::Batch(queries) => {
                            let mut results =
                                searched.next().expect("one search per queued batch")?;
                            report::render_queries(
                                &mut output,
                                &mut *diagnostics.borrow_mut(),
                                &options,
                                queries,
                                &subjects,
                                &shown,
                                &mut results,
                            )?;
                            processed += queries.len();
                        }
                    }
                }
                Ok(())
            };
            let parsed = input::parse_fasta_batches_with_warnings(
                query,
                false,
                args.lcase_masking,
                BATCH_SIZE,
                &mut |queries| {
                    if queries.is_empty() {
                        bail!("Empty CBlastQueryVector");
                    }
                    let large = x_large(queries);
                    let events = {
                        let mut events = x_events.borrow_mut();
                        events.push(XEvent::Batch(queries.to_vec()));
                        let queued = events
                            .iter()
                            .filter(|event| matches!(event, XEvent::Batch(_)))
                            .count();
                        if !large && queued < x_wave {
                            return Ok(());
                        }
                        std::mem::take(&mut *events)
                    };
                    let flushed = x_flush(events);
                    if flushed.is_err() {
                        x_failed.set(true);
                    }
                    flushed
                },
                &mut |warning| {
                    let mut bytes = Vec::new();
                    native::reader_warning(warning, &mut bytes)?;
                    x_events.borrow_mut().push(XEvent::Warn(bytes));
                    Ok(())
                },
            );
            if !x_failed.get() {
                x_flush(x_events.take())?;
            }
            return parsed;
        }
        input::parse_fasta_batches_with_warnings(
            query,
            false,
            args.lcase_masking,
            BATCH_SIZE,
            &mut |queries| {
                if queries.is_empty() {
                    bail!("Empty CBlastQueryVector");
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
                // ```c
                //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
                //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                //                 results = lcl_blast.Run();
                // ```
                // Dispatch point: a kept pool or search_internal; same result either way.
                let mut results = match x_pool {
                    Some(pool) => {
                        runtime::x_search_internal_in_pool(queries, &subjects, &options, pool)
                    }
                    None => runtime::search_internal(queries, &subjects, &options),
                }?;
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
                // ```c++
                //             	ITERATE(CSearchResultSet, result, *results) {
                //                	    formatter.PrintOneResultSet(**result, query_batch);
                // ```
                report::render_queries(
                    &mut output,
                    &mut *diagnostics.borrow_mut(),
                    &options,
                    queries,
                    &subjects,
                    &shown,
                    &mut results,
                )?;
                processed += queries.len();
                Ok(())
            },
            &mut |warning| native::reader_warning(warning, &mut *diagnostics.borrow_mut()),
        )
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
    // ```c
    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
    //                 results = lcl_blast.Run();
    // ```
    // Dispatch point: one pool for the whole run (LOSAT_X_BXPOOL) instead of one per batch, as
    // NCBI's CLocalBlast has per batch. Reuse of threads only.
    if runtime::x_bx_pool() && options.num_threads > 1 {
        crate::utils::threading::x_with_search_pool_local(
            options.num_threads as usize,
            "BLASTX",
            |pool| x_run_batches(Some(pool)),
        )?;
    } else {
        x_run_batches(None)?;
    }
    drop(x_run_batches);
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    report::write_epilog(&mut output, &options, &subjects, &shown, processed)?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ABI v1 is frozen (plan TD-1): its BLASTX `-outfmt` keeps the values and the parser
    // errors of the CLI parser before BLASTX read `-outfmt` as NCBI (dev 0797cc40, whose
    // `LOSAT blastx -outfmt 5` wrote these lines with exit status 2).
    #[test]
    fn web_v1_keeps_its_outfmt_values_and_errors() {
        let error = |outfmt: &str| {
            let error = parse_args(&["-outfmt", outfmt]).unwrap_err();
            crate::cli::render_message(error.downcast_ref::<clap::Error>().unwrap())
        };
        for (outfmt, reason) in [
            ("5", "unsupported BLASTX outfmt: use 0, 6 or 7"),
            ("abc", "unsupported BLASTX outfmt: use 0, 6 or 7"),
            ("13", "unsupported BLASTX outfmt: use 0, 6 or 7"),
            ("0 std", "custom fields require outfmt 6/7"),
            ("6 delim=", "unsupported BLASTX outfmt field 'delim='"),
            ("7 qaccver foo", "unsupported BLASTX outfmt field 'foo'"),
        ] {
            assert_eq!(
                error(outfmt),
                format!(
                    "error: invalid value '{outfmt}' for '-outfmt <OUTFMT>': {reason}\n\nFor more information, try '--help'.\n"
                )
            );
        }
        // White space between the words, read as ABI v1 read it.
        for (outfmt, number, fields) in [
            ("6\tqseqid", 6, &["qseqid"][..]),
            (" 7 ", 7, &[]),
            ("0", 0, &[]),
        ] {
            let args = parse_args(&["-outfmt", outfmt]).unwrap();
            let options = v1_resolve(&args).unwrap();
            assert_eq!(options.outfmt, number, "{outfmt:?}");
            if !fields.is_empty() {
                assert_eq!(options.fields, fields, "{outfmt:?}");
            }
        }
    }
}

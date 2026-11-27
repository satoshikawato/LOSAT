//! Existing reactor ABI adapter; the search, batching and reports use BLASTX owners.
#[cfg(feature = "parallel")]
use super::input;
use super::{args::BlastxArgs, input_source, native, query_setup::BATCH_SIZE, report, runtime};
use crate::api::local_blast::OutputSink;
use anyhow::{bail, Result};
use std::{cell::RefCell, io::Write, path::Path};

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// The existing CLI v2 parser owns all defaults, lexical/scope validation and
// code32 rejection. Synthetic paths identify memory input; no file is opened.
pub(crate) fn parse_args(extra: &[&str]) -> Result<BlastxArgs> {
    // The frozen ABI v1 (plan TD-1) has no database search; LOSAT Web searches a database
    // through ABI v2.
    if extra
        .iter()
        .any(|key| key.split('=').next().unwrap_or(key) == "-db")
    {
        bail!("-db (a database search) is not supported in the browser yet");
    }
    if extra.iter().any(|key| {
        matches!(
            key.split('=').next().unwrap_or(key),
            "-query"
                | "-subject"
                | "-out"
                | "-h"
                | "-help"
                | "-help-full"
                | "-xmlhelp"
                | "--help"
                | "--version"
        )
    }) {
        bail!("file/help arguments are unsupported for the BLASTX web API");
    }
    let mut argv = vec!["losat", "blastx", "-query", "query", "-subject", "subject"];
    argv.extend_from_slice(extra);
    // ABI v1 is frozen (plan TD-1): its `-outfmt` keeps the values and the parser errors of
    // the CLI's BLASTX parser before it read `-outfmt` as NCBI (`v1_outfmt`).
    let command = crate::web_v1_argv::command("blastx").mut_subcommand("blastx", |blastx| {
        blastx.mut_arg("outfmt", |outfmt| outfmt.value_parser(v1_outfmt))
    });
    let argv = argv.into_iter().map(Into::into).collect();
    match crate::web_v1_argv::parse(command, argv).map_err(anyhow::Error::msg)? {
        crate::cli::Commands::Blastx(args) => Ok(args),
        _ => unreachable!("fixed BLASTX subcommand"),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2827-2845
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
/// (`v1_outfmt`: the words between white space, each `std` or one of `args::FIELDS`), before
/// BLASTX read `-outfmt` as NCBI; a list other than the `std` fields is written by the shared
/// tabular writer, as the CLI's (`ResolvedOptions::tabular`), without a taxonomy database
/// (DECISIONS 7 of line CLI C2; ABI v1's fields show no taxonomy).
fn v1_resolve(args: &BlastxArgs) -> Result<super::args::ResolvedOptions> {
    let outfmt = args.outfmt.to_string_lossy();
    let mut words = outfmt.split_whitespace();
    let number = words.next().unwrap_or_default().parse::<u8>()?;
    args.resolve_with_format(&crate::blastinput::app::FormatChoice {
        number: i32::from(number),
        spec: words.collect::<Vec<_>>().join(" ").into_bytes(),
        delimiter: Vec::new(),
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:302-314
// Original bytes are reparsed on each call with this role's gencode/mask
// options; generic FASTA handles never cache translated or masked search state.
#[cfg_attr(
    not(losat_web_abi_v1),
    expect(dead_code, reason = "called only by the web ABI (`web_api`, wasm32)")
)]
pub(crate) fn run_pair(
    args: BlastxArgs,
    query: &[u8],
    subject: &[u8],
    subject_label: &str,
) -> Result<Vec<u8>> {
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let diagnostics = RefCell::new(std::io::stderr());
    let subjects = input_source::read_web_subjects(subject, args.lcase_masking, &mut |message| {
        Ok(diagnostics.borrow_mut().write_all(message)?)
    })?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
    // (The same bytes as before for the options ABI v1 had: a hit list size below 5 from
    // -max_target_seqs; `native::warnings_of` also reads -num_descriptions / -num_alignments,
    // which the CLI parser now accepts.)
    let v1_outfmt_number = args
        .outfmt
        .to_string_lossy()
        .split_whitespace()
        .next()
        .and_then(|word| word.parse::<i32>().ok())
        .unwrap_or(0);
    diagnostics
        .borrow_mut()
        .write_all(&native::warnings_of(&args, v1_outfmt_number))?;
    let options = v1_resolve(&args)?;
    if query.is_empty() {
        bail!("No sequence input was provided");
    }
    if subjects.is_empty() {
        bail!("BLASTX requires a protein subject");
    }
    let mut output = Vec::new();
    let label = Path::new(subject_label);
    let no_parsed = report::RecordsParsed::default();
    let shown = report::ShownSubjects::Fasta {
        path: label,
        parsed: &no_parsed,
    };
    // The run's decision for a matrix name that NCBI's XML does not know
    // (`blast_xml::xml_matrix_failure`): none, as ABI v1 writes outfmt 0, 6 and 7 (`v1_outfmt`);
    // its query would be the run's without scopes of `-parse_deflines`, which ABI v1 does not read
    // (`QueryScopes::xml_failing_query`).
    let xml_matrix = super::blast_xml::xml_matrix_failure(std::slice::from_ref(&options), true)
        .map(|failure| {
            crate::blastinput::query_batch::QueryScopes::xml_failing_query(None, failure)
        });
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:784-788
    super::native::write_empty_subject_warnings(
        &subjects,
        options.app_name,
        &mut *diagnostics.borrow_mut(),
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
    report::write_prolog(
        &mut output,
        &mut *diagnostics.borrow_mut(),
        None,
        &options,
        &subjects,
        &shown,
    )?;
    let mut processed = 0;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
    // One pool for the whole run at 2+ threads: see native.rs.
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
    // Same batch loop as in native.rs, as a closure over the pool its searches run on (see
    // native.rs).
    let mut x_run_batches = |x_pool: Option<&crate::utils::threading::SearchPool<'_>>| {
        // Several query batches searched at a time with a pool of 2+ threads: see native.rs.
        #[cfg(feature = "parallel")]
        if let Some(pool) = x_pool.filter(|pool| pool.enabled()) {
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
                                &mut crate::cli::NcbiOstream::sink(&mut OutputSink::Writer(
                                    &mut output,
                                ))?,
                                &mut *diagnostics.borrow_mut(),
                                &options,
                                queries,
                                &subjects,
                                &shown,
                                &mut results,
                                processed,
                                // No document here holds a report back (`stops`).
                                false,
                                None,
                                xml_matrix,
                                // ABI v1 reads no `-parse_deflines` (`parse_args`).
                                None,
                                // The queued batches are searched: the pool is free.
                                Some(pool),
                            )?;
                            processed += queries.len();
                        }
                    }
                }
                Ok(())
            };
            let parsed = input_source::read_web_queries(
                query,
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
                &mut |message| {
                    x_events.borrow_mut().push(XEvent::Warn(message.to_vec()));
                    Ok(())
                },
            );
            if !x_failed.get() {
                x_flush(x_events.take())?;
            }
            return parsed;
        }
        input_source::read_web_queries(
            query,
            args.lcase_masking,
            BATCH_SIZE,
            &mut |queries| {
                if queries.is_empty() {
                    bail!("Empty CBlastQueryVector");
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
                // Dispatch point: a kept pool or search_internal; same result either way.
                let mut results = match x_pool {
                    Some(pool) => {
                        runtime::x_search_internal_in_pool(queries, &subjects, &options, pool)
                    }
                    None => runtime::search_internal(queries, &subjects, &options),
                }?;
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
                report::render_queries(
                    &mut crate::cli::NcbiOstream::sink(&mut OutputSink::Writer(&mut output))?,
                    &mut *diagnostics.borrow_mut(),
                    &options,
                    queries,
                    &subjects,
                    &shown,
                    &mut results,
                    processed,
                    // No document here holds a report back (`stops`).
                    false,
                    None,
                    xml_matrix,
                    // ABI v1 reads no `-parse_deflines` (`parse_args`).
                    None,
                    x_pool,
                )?;
                processed += queries.len();
                Ok(())
            },
            &mut |message| Ok(diagnostics.borrow_mut().write_all(message)?),
        )
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
    // Dispatch point: with 2+ threads one pool for the whole run instead of one per batch, as
    // NCBI's CLocalBlast has per batch. Reuse of threads only.
    if options.num_threads > 1 {
        crate::utils::threading::x_with_search_pool_local(
            options.num_threads as usize,
            "BLASTX",
            |pool| x_run_batches(Some(pool)),
        )?;
    } else {
        x_run_batches(None)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
    report::write_epilog(
        &mut output,
        &mut *diagnostics.borrow_mut(),
        &options,
        &subjects,
        &shown,
        processed,
        None,
        xml_matrix,
    )?;
    Ok(output)
}

#[cfg(test)]
#[path = "web_tests.rs"]
mod tests;

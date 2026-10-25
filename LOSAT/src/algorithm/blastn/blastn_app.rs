//! BLASTN's application: the command line, the input, and the loop over the query batches.
//!
//! This file is a child module of `blast_engine/run.rs` (declared there with `#[path]`, as
//! `database.rs` is), so that it reaches the run's types and the search of a query batch
//! without changing them.
//!
//! NCBI reference: `c++/src/app/blast/blastn_app.cpp:128-318` (`CBlastnApp::Run`,
//! `x_RunMTBySplitDB`), the option checks of
//! `c++/src/algo/blast/blastinput/blast_args.cpp:1995-3636` and
//! `c++/src/app/blast/blast_app_util.cpp:204-860`, and
//! `c++/src/app/blast/blast_formatter.cpp:429-467`.

use super::blast_format::{
    write_batch_reports, write_before_failed_batch, write_pairwise_prologs, write_reports,
    BatchOutputs, BatchResults, ReportSubjects, FIRST_BATCH_OUTPUTS,
};
use super::{
    database, drop_hsps_of_score_zero, finish_query_hit_list, search_query_batch, x_pair_seed,
    x_qbatch_par, x_query_drop, x_query_stream, x_rpt_stream, BatchStage, PreparedSubjects,
    QueryBatch, QueryBatching, ShownSubjects, SubjectSource, Subjects, MAX_DBSEQ_LEN,
};
use crate::algorithm::blastn::archive;
use crate::algorithm::blastn::args::{BlastnArgs, QueryStrand};
use crate::algorithm::blastn::coordination::{
    configure_task, finalize_task_config, subject_metadata_from_records,
};
use crate::algorithm::blastn::hsp::{BlastnHitList, BlastnOutputFormat};
use crate::algorithm::blastn::lookup::x_mb_lut_keep;
use crate::algorithm::blastn::query_split::SplitSizes;
use crate::algorithm::blastn::scoring::{check_losat_limits, check_scoring_options};
use crate::algorithm::blastn::tracing as blastn_trace;
use crate::algorithm::blastn::x_subject_store::{x_bnstore, x_read_subject_stores, SubjectLetters};
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::app::{self as app, Program};
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords, ReadError,
    ReaderConfig,
};
use crate::blastinput::input_files::{
    check_utf8_file_name, is_blank, open_input, write_empty_subject_warnings_by_length,
};
use crate::blastinput::query_batch::BatchSizeMixer;
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, next_ranged_batch, parse_optional_range, subjects_read, Placements,
    QueryInput, RangeRole, SequenceRange,
};
use crate::core::blast_encoding::encode_subject_ncbi2na_packed;
use crate::report::query_warnings::{invalid_query_warning, query_warning};
use crate::sequence::packed_subject::PackedNucleotide;
use anyhow::{Context, Result};

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::check_ncbi_application_settings`): the data loaders decide whether
/// the readers of the subjects and the queries try a first line as a Seq-id
/// (`ReaderConfig`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
/// ```c++
///     CNcbiApplication* app = CNcbiApplication::Instance();
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         x_LoadDataLoadersConfig(registry);
///         x_LoadBlastDbDataLoaderConfig(registry);
///     }
/// ```
pub fn run(
    mut args: BlastnArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:128-133
    // ```c
    // if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
    // 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
    // }
    // else {
    // 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptions(args));
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // ```c
    //     if (GetExportSearchStrategyStream(args) ||
    //            m_FormattingArgs->ArchiveFormatRequested(args)) {
    //         locality = CBlastOptions::eBoth;
    //     }
    // ```
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too (`cli_output_formats`).
    let output_formats = cli_output_formats(&args.outfmt, args.out.as_deref())?;
    // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot
    // write, before anything is read or written.
    let archive_requested = output_formats.contains(&BlastnOutputFormat::Archive);
    args.archive = archive::handle_for(&args, args.dust_filtering.as_deref(), archive_requested)?;
    // LOSAT's thread capability is checked before any input or output, as before.
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3631-3636
    // ```c
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
    // ```c
    //     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
    //     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
    //     arg.Reset(m_BlastDbArgs);
    //     m_Args.push_back(arg);
    //
    //     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
    //     arg.Reset(m_StdCmdLineArgs);
    //     m_Args.push_back(arg);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // ```c
    //         CRef<blast::CBlastQueryVector> subjects;
    //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
    //                                        subj_range, parse_deflines,
    //                                        use_lcase_masks, subjects, m_IsMapper);
    //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
    // ```
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // before the query and the output are opened and the options are checked. Files are
    // opened when a handler asks for them, and each once: a named pipe gives its bytes to
    // one reader only.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // ```c
    //     } else if (!m_IsIgBlast){
    //         // IgBlast permits use of germline database
    //         NCBI_THROW(CInputException, eInvalidInput,
    //            "Either a BLAST database or subject sequence(s) must be specified");
    //     }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:172-176
    // ```c
    //     catch (const blast::CInputException& e) {                               \
    //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
    //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // The handler of the database arguments raises it, before the query and the output
    // are opened (`CStdCmdLineArgs` comes after it).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    // ```c++
    //     if (args.Exist(kArgDb) && args[kArgDb]) {
    // 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
    //
    //         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
    //                                              mol_type));
    // ```
    // With `-db` the handler of the database arguments only keeps the name: the database
    // is opened after every option handler has run (`search_cli`). With `-subject` it
    // reads the subjects here, before the query and the output are opened.
    let archive = args.archive.clone();
    let mut stderr = std::io::stderr();
    let mut archive_tee;
    let diagnostics: &mut (dyn std::io::Write + Send) = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
            &mut archive_tee
        }
        None => &mut stderr,
    };
    let cli_subjects = if args.db.is_some() {
        None
    } else {
        match read_cli_subjects(&args, &settings, &mut *diagnostics, archive.is_some()) {
            Ok(read) => Some(read),
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        }
    };
    if archive.is_some() {
        if let Some(read) = &cli_subjects {
            archive::check_titles(&read.subjects, "subject")?;
            if let Some(full) = &read.records {
                archive::check_titles(full, "subject")?;
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // ```c
    //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
    //         m_InputStream == NULL) {
    // ...
    //         else {
    //             m_InputStream = &args[kArgQuery].AsInputFile();
    //         }
    //     }
    // ...
    //     else {
    //         m_OutputStream = &args[kArgOutput].AsOutputFile();
    //     }
    // ```
    // NCBI reference: ncbi-blast/c++/src/corelib/ncbiargs.cpp:777-781
    // ```c
    //     if (AsString() == "-") {
    // #if defined(NCBI_OS_MSWIN)
    //         NcbiSys_setmode(NcbiSys_fileno(stdout), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
    // #endif
    //         m_Ios = &cout;
    // ```
    // The output file is created here, before the options are checked, so a run that
    // stops at a check leaves an empty file; `-` is standard output.
    check_utf8_file_name(&args.query, "query", "BLASTN")?;
    // NCBI opens the query file inside `Run` too (`CStdCmdLineArgs` in its handler of the
    // arguments does not): its error is saved and an error archive follows.
    let query_file = match open_input(&args.query, "query", "BLASTN") {
        Ok(file) => file,
        Err(error) => {
            write_run_error_archive(&args, archive.as_ref(), &error);
            return Err(error);
        }
    };
    // The output file is opened once, as NCBI's stream (a named pipe gives one reader one
    // end of file).
    let out_file = open_out_file(&args)?;
    let outfmt = args.outfmt.clone();
    let pairwise = output_formats.contains(&BlastnOutputFormat::Pairwise);
    let mut report = cli_report_stream(out_file);
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut report), &mut *diagnostics);
        search_cli(
            args,
            query_file,
            cli_subjects.as_ref().map(|read| {
                (
                    read.subjects(),
                    SubjectRanges {
                        placements: &read.placements,
                        records: read.records.as_deref(),
                    },
                )
            }),
            &mut outputs,
            output_formats,
            settings,
        )
    };
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status:
    // blastn_app.cpp:331-334).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                if let Some(lines) = archive::error_lines(error, false) {
                    let _ =
                        crate::report::archive::write_error_archive(handle, &lines, &mut report);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut report);
                }
            }
        }
    }
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = std::io::Write::flush(&mut report);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
    // ```c
    // {
    //     m_Outfile.exceptions(NcbiBadbit);
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
    // ```c
    //     catch (const std::ios::failure&) {                                      \
    //         LOG_POST(Error << "BLAST failed to write output");                  \
    //         exit_code = BLAST_OUTPUT_ERROR;                                     \
    //     }                                                                       \
    // ```
    // The outfmt 0 formatter's stream throws when a write fails (oracle: `-out /dev/full`,
    // exit 6). For outfmt 6/7 NCBI aborts instead, and LOSAT reports the error
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).
    if report.failed && pairwise {
        return Err(crate::cli::NativeError {
            exit: 6,
            message: "BLAST failed to write output\n".to_string(),
        }
        .into());
    }
    result?;
    flushed.context("failed to write the output")
}

/// The subjects of `-subject` as the handler of the database arguments reads them (`run`):
/// the searched records (their intervals with `-subject_loc`), where their letters lie in
/// the records, and the records as read (with `-subject_loc`).
/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastnArgs,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    let Some(lines) = archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = cli_report_stream(out_file);
        let _ = crate::report::archive::write_error_archive(handle, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &BlastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if let Some(path) = args.out.as_deref() {
        check_utf8_file_name(path, "out", "BLASTN")?;
    }
    Ok(
        match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
            Some(path) => Some(std::io::BufWriter::new(
                std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
            )),
            None => None,
        },
    )
}

/// The report stream of the command line: the `-out` file, or standard output.
fn cli_report_stream(
    out_file: Option<std::io::BufWriter<std::fs::File>>,
) -> crate::cli::ReportStream {
    crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    }
}

struct CliSubjects {
    subjects: Vec<FastaRecord>,
    placements: Placements,
    records: Option<Vec<FastaRecord>>,
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of `subjects`, packed (`x_read_subject_stores`;
    /// the records keep none).
    stores: Option<Vec<PackedNucleotide>>,
}

impl CliSubjects {
    /// The searched subjects.
    fn subjects(&self) -> Subjects<'_> {
        Subjects::Records(SubjectLetters::new(&self.subjects, self.stores.as_deref()))
    }
}

/// Reads the subjects of `-subject` where NCBI's handler of the database arguments reads
/// them (`run`; NCBI references there), with NCBI's error when `-subject` is missing (a
/// search without `-db`).
fn read_cli_subjects(
    args: &BlastnArgs,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut (dyn std::io::Write + Send),
    keep_deflines: bool,
) -> Result<CliSubjects> {
    let Some(subject_path) = args.subject.as_deref() else {
        return Err(crate::cli::NativeError {
            exit: 1,
            message: "BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified\nPlease refer to the BLAST+ user manual.\n".to_string(),
        }
        .into());
    };
    check_utf8_file_name(subject_path, "subject", "BLASTN")?;
    let subject_file = open_input(subject_path, "subject", "BLASTN")?;
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // ```c++
    //             subj_input_stream = &args[kArgSubject].AsInputFile();
    //         }
    //
    //         TSeqRange subj_range;
    //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
    //             subj_range =
    //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
    //                             "Invalid specification of subject location");
    //         }
    // ```
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "BLASTN")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // ```c++
    //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
    //     CRef<CBlastInput> input(new CBlastInput(fasta));
    //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //     sequences = input->GetAllSeqs(*scope);
    //     return scope;
    // ```
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`).
    let mut subject_source = FastaInputSource::from_argument(
        subject_path,
        subject_file,
        ReaderConfig::subject("BLASTN", false, settings.data_loaders).with_defline(keep_deflines),
    );
    // EXPERIMENT (LOSAT_X_BNSTORE): without a subject range, each record's letters are packed as
    // soon as it has been read, and the records keep none (`x_read_subject_stores`).
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs); with
    // LOSAT_X_BNSTORE the packing of each record is in this interval too. No NCBI counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let (read_subjects, stores) = {
        let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
        // The archive writes the subjects' letters from the records, not from a packed store.
        if x_bnstore() && subject_range.is_none() && !keep_deflines {
            x_read_subject_stores(&mut subject_source, &mut warn)?
        } else {
            (
                read_subjects(&mut subject_source, subject_range.as_ref(), &mut warn)?,
                None,
            )
        }
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    if x_bnstore() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNSTORE", stores.is_some());
    }
    drop(subject_source);
    // NCBI reads records without residues (and intervals without letters) without a message
    // and reports them when it sets up the search (`search`).
    let ranged_subjects = cut_subjects(&read_subjects, subject_range.as_ref())?;
    let (subjects, subject_placements, full_subjects) = match ranged_subjects {
        Some((cut, placements)) => (cut, placements, Some(read_subjects)),
        None => (read_subjects, Placements::default(), None),
    };
    check_subjects_not_empty(&subjects)?;
    Ok(CliSubjects {
        subjects,
        placements: subject_placements,
        records: full_subjects,
        stores,
    })
}

/// The part of `run` after the output is opened: NCBI's filtering handler (`-dust`), the
/// processing of the options, the opening of the database of `-db`, `Query is Empty!`, the
/// reading of the queries and the search. `subjects` are those of `-subject` (`None` with
/// `-db`).
fn search_cli(
    mut args: BlastnArgs,
    query_file: std::fs::File,
    subjects: Option<(Subjects<'_>, SubjectRanges<'_>)>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.resolve_dust()?;
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-219
    // ```c++
    //         try {
    //             // Try to open the BLAST database even for remote searches, as if
    //             // it is available locally, it will be better to fetch the
    //             // sequence data for formatting from this (local) source
    //             CRef<CSeqDB> seqdb = search_db->GetSeqDb();
    //             db_adapter.Reset(new CLocalDbAdapter(*search_db));
    //             dl = RegisterOMDataLoader(seqdb);
    // ```
    // The database of `-db` is opened here: after the option handlers (the query file is
    // open and the output file exists) and before the test for an empty query.
    let database = match args.database_name() {
        Some(name) => Some(database::BlastnDatabase::open(name, &settings)?),
        None => None,
    };
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for serial in output_formats.iter().filter_map(|format| format.serial()) {
            crate::report::serial::seqalign::check_database(
                serial,
                database.name(),
                &database.pairwise_database().definition,
                Program::Blastn,
            )?;
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:201-214
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())) {
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query's data loaders are those of the subjects (`SDataLoaderConfig`; the BLAST
    // database of a database search is added only to the loaders' databases). The position
    // is taken on the opened file, before it is read (`stream_is_empty`). When the subjects
    // were read from standard input too, `cin` has reached its end (a failed stream), so
    // NCBI gets no position for the query either, and its reader is at its end: no query
    // batch is read.
    let shared_standard_input =
        args.query.as_os_str() == "-" && args.subject_path().as_os_str() == "-";
    // EXPERIMENT (LOSAT_X_QSTREAM): what its rule reads (`x_query_stream::Rule`): the task, the
    // largest query batch, and whether LOSAT_X_QBATCHPAR's slots may run: the switch on, the
    // conditions of `x_qbatch_par_off` (`-subject`, two threads or more, the build, the output, fewer
    // subject units than threads), and a budget (LOSAT_X_QBATCHPARMEM when it is set, else
    // `x_qbatch_par::Profile::of_search`'s default, none for a query without a file size, such as
    // standard input) or a slot count set by LOSAT_X_QBATCHPARSLOTS. EXPERIMENT (LOSAT_X_BNQDROP)
    // reads it too (`x_query_drop::Rule`).
    let x_slots_may_run = match (subjects.as_ref(), x_query_file_bytes(&args.query)) {
        (Some((Subjects::Records(letters), _)), bytes) => {
            let switches = x_qbatch_par::Switches::from_env();
            let lengths = (0..letters.records().len()).map(|index| letters.len(index));
            let subject_letters = lengths.clone().map(|length| length as u64).sum();
            let budget = switches.slots.is_some()
                || match switches.mem_mib {
                    Some(mib) => mib > 0,
                    None => x_qbatch_par::ncbi_splits_by_queries(
                        crate::algorithm::blastn::coordination::determine_effective_word_size(
                            &args,
                        ),
                        bytes,
                        subject_letters,
                    ),
                };
            switches.on
                && budget
                && x_qbatch_par_off_of(&args, Some(lengths), args.num_threads).is_none()
        }
        _ => false,
    };
    let x_rule = x_query_stream::Rule {
        megablast: args.task == "megablast",
        batch_letters: x_query_batch_letters(&args.task),
        threads: args.num_threads,
        slots_may_run: x_slots_may_run,
    };
    // EXPERIMENT (LOSAT_X_QSTREAM): whether the query file is read a second time, batch by batch
    // (a second handle on it, at the position its reader starts from), where the reports may be
    // written batch by batch (LOSAT_X_RPTSTREAM).
    let x_opening = x_query_stream::Opening::new(
        &args.query,
        &query_file,
        shared_standard_input,
        x_rpt_stream::switch_on(),
        x_rule,
    );
    let query_config = ReaderConfig::query("BLASTN", false, settings.data_loaders)
        .with_defline(args.archive.is_some());
    let mut query_source = FastaInputSource::from_argument(&args.query, query_file, query_config);
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`run_in_pool`). What ends the reading (the end of the input, an `eEOF`, a reader
    // error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads it in,
    // after the reports of the batches before (`run_in_pool`).
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    // EXPERIMENT (LOSAT_X_QSTREAM): with the switch acting this reading keeps only each record's
    // length and how the reading ends; the batches read their records again (`x_query_stream`).
    let (records, end, x_store) = if x_opening.acts() {
        let (store, end) = x_query_stream::Store::read(
            x_opening,
            &mut query_source,
            query_config,
            query_range.as_ref(),
        )?;
        (Vec::new(), end, store)
    } else {
        let QueryRecords { records, end } = read_queries(&mut query_source);
        (records, end, x_query_stream::Store::off(x_opening))
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    if args.archive.is_some() {
        archive::check_titles(&records, "query")?;
    }
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    // A database has no subject ranges (`-subject_loc` excludes `-db`).
    let no_placements = Placements::default();
    let (subjects, subject_ranges) = match (subjects, database.as_ref()) {
        (Some((subjects, ranges)), _) => (subjects, ranges),
        (None, Some(database)) => (
            Subjects::Database(database),
            SubjectRanges {
                placements: &no_placements,
                records: None,
            },
        ),
        (None, None) => unreachable!("`run` reads the subjects of a search without -db"),
    };
    // EXPERIMENT (LOSAT_X_BNQDROP): where no report reads the queries' letters, the records read
    // before the search move into a store that the batch loop takes them out of, so a batch can drop
    // its letters once it has encoded them (`x_query_drop.rs`).
    let x_taken = !x_store.acts() && x_query_drop::on() && {
        let acts = x_query_drop::Rule {
            formats: &output_formats,
            query_range: query_range.is_some(),
            slots_may_run: x_slots_may_run,
            hit_records: outputs.hits.is_some(),
        }
        .acts();
        crate::utils::xreach::x_reach!("LOSAT_X_BNQDROP", acts);
        acts
    };
    let x_taken_store;
    let queries = if x_store.acts() {
        SearchQueries::Streamed(&x_store)
    } else if x_taken {
        x_taken_store = x_query_drop::Store::new(records);
        SearchQueries::Taken(&x_taken_store)
    } else {
        SearchQueries::Records(&records)
    };
    search(
        args,
        queries,
        reading,
        Ranges {
            query: query_range,
            subjects: subject_ranges,
        },
        subjects,
        outputs,
        output_formats,
    )
}

/// How the reading of the queries of a search ended, after its records.
#[derive(Debug)]
pub(crate) enum QueryReading {
    /// The records of `run_local`: an input without records is NCBI's `Query is Empty!`
    /// (the inputs of the ABI are read as regular files).
    Records,
    /// The CLI's reader reached the end of its input (`CBlastInput::End()`); without
    /// records, no query batch is read.
    End,
    /// The CLI's reader found only blank and comment lines where it expected a record
    /// (`eEOF`), which ends the batch being read.
    BlankLines,
    /// The CLI's reader failed on the record after the last one read (a reader error, or
    /// LOSAT's explicit rejection of a line that NCBI would fetch as a Seq-id), after
    /// writing `warnings`.
    Error { error: ReadError, warnings: Vec<u8> },
}

impl QueryReading {
    /// Whether NCBI's reader reads on after the last record (`End()` is false): an `eEOF`
    /// or an error then comes in the batch being read.
    pub(crate) fn reads_past_records(&self) -> bool {
        matches!(self, QueryReading::BlankLines | QueryReading::Error { .. })
    }
}

/// The query range, read by NCBI's query options handler: after the filtering handler
/// (`-dust`), before the formatting handler and the check of the options.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blastn_args.cpp:100-102
/// ```c++
///     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
///     arg.Reset(m_QueryOptsArgs);
///     m_Args.push_back(arg);
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
/// ```c++
///     // set the sequence range
///     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
///         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
///                                      "Invalid specification of query location");
///     }
/// ```
fn parse_query_range(args: &BlastnArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.query_loc.as_deref(), RangeRole::Query, "BLASTN")
}

/// The ranges of a search: the query range as given, and the subjects' ranges.
struct Ranges<'a> {
    query: Option<SequenceRange>,
    subjects: SubjectRanges<'a>,
}

/// Where each subject's searched letters lie in its record, and the records as read
/// (with `-subject_loc`; the report shows their lowercase letters).
#[derive(Clone, Copy)]
pub(super) struct SubjectRanges<'a> {
    pub(super) placements: &'a Placements,
    pub(super) records: Option<&'a [FastaRecord]>,
}

/// The ranges as the report uses them: the query input and the subjects' ranges.
pub(super) struct ReportRanges<'a> {
    pub(super) queries: &'a QueryInput,
    pub(super) subjects: SubjectRanges<'a>,
}

/// NCBI's error for a subject file without records, raised when it reads the subjects
/// (the first option handler, `run`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// ```c
/// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
///     : m_QueryVector(& queries)
/// {
///     if (queries.Empty()) {
///         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
pub(crate) fn check_subjects_not_empty<R>(subjects: &[R]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// The report of the CLI's `-outfmt`, read as NCBI reads it before its option handlers run
/// (`run`; `blastinput/app.rs` `parse_formatting_string`, `report_format`). A format that
/// NCBI itself rejects later has no report: the run stops at NCBI's error, the formatting
/// handler's in `process_options` (AIRR rearrangement, FASTA) or the XInclude formats' in
/// `search` (13 and 14 to standard output).
fn cli_output_formats(
    outfmt: &str,
    out: Option<&std::path::Path>,
) -> Result<Vec<BlastnOutputFormat>> {
    let choice = app::parse_formatting_string(outfmt)?;
    app::report_format(&choice, Program::Blastn, out)?
        .map(|format| BlastnOutputFormat::of(format, &choice.spec))
        .into_iter()
        .collect()
}

/// The requested output formats of `run_local`, each read as the CLI reads its single
/// `-outfmt` (`cli_output_formats`); a format that NCBI rejects later fails here with NCBI's
/// error, as in the other programs' `run_local` (they have no output file).
fn parse_output_formats<'a>(
    formats: impl IntoIterator<Item = &'a str>,
) -> Result<Vec<BlastnOutputFormat>> {
    formats
        .into_iter()
        .map(|outfmt| {
            let choice = app::parse_formatting_string(outfmt)?;
            match app::report_format(&choice, Program::Blastn, None)? {
                Some(format) => BlastnOutputFormat::of(format, &choice.spec),
                None => {
                    app::formatting_handler_check(&choice, Program::Blastn)?;
                    app::xinclude_check(&choice)?;
                    anyhow::bail!(
                        "output format {} is not supported by LOSAT's BLASTN",
                        choice.number
                    )
                }
            }
        })
        .collect()
}

/// NCBI's processing of the options after the query and the output are opened
/// (`SetOptions`): the formatting handler (its errors for the formats of other programs,
/// then its few-matches warning), and the check of the options.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// ```c
///     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
///     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
/// ...
///     if((m_OutputFormat == eAirrRearrangement) && !(m_FormatFlags & eIsAirrRearrangement) ){
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "AIRR rearrangement format is only applicable to igblastn" );
///     }
///     if (m_OutputFormat == eFasta) {
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "FASTA output format is only applicable to magicblast");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// ```c
///     if(hitlist_size < 5){
///    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3636
/// ```c
///     try { retval->Validate(); }
/// ```
/// The hit list size is the -max_target_seqs value, or 500 when it is omitted.
fn process_options(args: &mut BlastnArgs, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    // The formatting handler comes after the query options handler (blastn_args.cpp:100-105).
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        app::formatting_handler_check(&choice, Program::Blastn)?;
        outfmts.push(choice.number);
    }
    // The hit list size and the counts of the pairwise report, with the formatting
    // handler's warnings (`ExtractAlgorithmOptions`), once for the run.
    let (counts, warnings) = crate::algorithm::blastn::scoring::report_counts_of(args, &outfmts)?;
    args.report_counts = counts;
    outputs.diagnostics.write_all(&warnings)?;
    check_scoring_options(args)
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one BLASTN search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI, web ABI v1 and
/// v2).
///
/// Each requested `-outfmt` is parsed as on the command line, before the search starts.
/// No BLASTN search option depends on the output formats it supports (6 and 7). The
/// `-query` and `-subject` values of `args` are used only as display names. The hit
/// records (`ReportOutputs::hits`) are not produced yet.
///
/// The BLAST archive (outfmt 11) writes each record's DefLine as the reader kept it
/// (`ReaderConfig::with_defline(true)`). Records without a kept defline get `>` and their
/// title, which differs from NCBI for a defline with extra spaces or control bytes.
pub fn run_local(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The BLAST archive (outfmt 11) among the formats: its request, the refusals of what
    // LOSAT cannot write, and a diagnostics stream that saves the messages.
    let archive_requested = outputs.formats.iter().any(|format| {
        app::parse_formatting_string(format.outfmt)
            .is_ok_and(|choice| choice.number == app::output_format::ARCHIVE_FORMAT)
    });
    if archive_requested && args.db.is_none() {
        QueryStrand::parse(&args.strand)?;
        args.archive = archive::handle_for(&args, args.dust_filtering.as_deref(), true)?;
        archive::check_titles(query_records, "query")?;
        archive::check_titles(subject_records, "subject")?;
        if let Some(handle) = args.archive.clone() {
            return outputs.with_archive_diagnostics(handle, |view| {
                run_local_formats(args, query_records, subject_records, view)
            });
        }
    }
    run_local_formats(args, query_records, subject_records, outputs)
}

fn run_local_formats(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The callers of `run_local` (LOSAT Web's adapter, ABI v1) give subject records; a
    // database is read from files by the command line only.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // The command line has checked `-strand` at parsing (clap); a caller of `run_local` has not.
    QueryStrand::parse(&args.strand)?;
    let output_formats = parse_output_formats(outputs.formats.iter().map(|format| format.outfmt))?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "BLASTN")?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((records, placements)) => (records.as_slice(), placements.clone()),
        None => (subject_records, Placements::default()),
    };
    check_subjects_not_empty(searched_subjects)?;
    args.resolve_dust()?;
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    search(
        args,
        SearchQueries::Records(query_records),
        QueryReading::Records,
        Ranges {
            query: query_range,
            subjects: SubjectRanges {
                placements: &subject_placements,
                records: ranged_subjects.is_some().then_some(subject_records),
            },
        },
        Subjects::Records(SubjectLetters::new(searched_subjects, None)),
        outputs,
        output_formats,
    )
}

/// The queries of a search: records read before it, or EXPERIMENT (LOSAT_X_QSTREAM) a query file
/// whose batches read their records when the batch loop needs them (`x_query_stream`).
#[derive(Clone, Copy)]
enum SearchQueries<'a> {
    Records(&'a [FastaRecord]),
    Streamed(&'a x_query_stream::Store),
    /// EXPERIMENT (LOSAT_X_BNQDROP): records read before the search that the batch loop takes out
    /// of their store (`x_query_drop.rs`).
    Taken(&'a x_query_drop::Store),
}

fn search(
    args: BlastnArgs,
    queries: SearchQueries<'_>,
    reading: QueryReading,
    ranges: Ranges<'_>,
    subjects: Subjects<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
) -> Result<()> {
    // EXPERIMENT (LOSAT_X_QSTREAM): a streamed search has no records here (its reading is never
    // `Records`), and its query input comes from the lengths of the first reading; EXPERIMENT
    // (LOSAT_X_BNQDROP) a taken one has its records in its store, with its query input.
    let query_records: &[FastaRecord] = match queries {
        SearchQueries::Records(records) => records,
        SearchQueries::Streamed(_) | SearchQueries::Taken(_) => &[],
    };
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
    // ```c
    //     if(m_IsDbScan) {
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // ```c
    //         catch(CBlastException & e ) {
    //         	// Skip bad subject sequence
    //         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
    //         		seqblk_vec->push_back(subj);
    //         ...
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    //         		continue;
    // ```
    // The formatter, made after `Query is Empty!` and before the outfmt 0 prolog, sets up
    // the subjects: a subject without letters (a record without residues, or an interval
    // that starts just past its record's end) gets its warning, stays in the database
    // statistics with no letters, and never hits. (A database is not set up so.)
    if let Subjects::Records(subject_set) = subjects {
        write_empty_subject_warnings_by_length(
            subject_set.records(),
            |index| subject_set.len(index),
            "blastn",
            outputs.diagnostics,
        )?;
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:255-258
    // ```c
    //         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made (a
    // format 13 or 14 to a file takes the file name as the base of its side files).
    for format in &outputs.formats {
        app::xinclude_check_for(
            &app::parse_formatting_string(format.outfmt)?,
            args.out.as_deref(),
        )?;
    }
    // NCBI reads the queries after `Query is Empty!`, one batch at a time, with its
    // reader's warnings (`run_in_pool`).
    // NCBI's x_CleanAndCompress reads past the end of some titles of punctuation (NCBI
    // crashes when it writes the title of such a subject with hits); LOSAT stops at the end
    // of the title (`report/defline.rs`, approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
    // LOSAT's limits come where NCBI starts the search, after its checks and its
    // `Query is Empty!` success.
    check_losat_limits(&args)?;
    check_unsupported_environment(matches!(subjects, Subjects::Database(_)))?;
    let batching = query_batching_from_environment(
        crate::algorithm::blastn::coordination::task_uses_megablast_chunks(&args.task),
    )?;
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`cut_queries`): a record is searched cut, and one whose interval starts just past
    // its end is a query without data (`run_in_pool`).
    let ranged_queries = ranges
        .query
        .as_ref()
        .filter(|_| matches!(queries, SearchQueries::Records(_)))
        .map(|range| cut_queries(query_records, range));
    let input_records = query_records;
    let (query_records, query_input) = match (ranged_queries, queries) {
        (Some(ranged), _) => (
            std::borrow::Cow::Owned(ranged.records),
            std::borrow::Cow::Owned(ranged.input),
        ),
        (None, SearchQueries::Records(_)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Owned(QueryInput::whole(query_records)),
        ),
        (None, SearchQueries::Streamed(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
        (None, SearchQueries::Taken(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
    };
    let queries = match queries {
        SearchQueries::Records(_) => x_query_stream::Queries::Records {
            input: input_records,
            searched: &query_records,
        },
        SearchQueries::Streamed(store) => x_query_stream::Queries::Streamed(store),
        SearchQueries::Taken(store) => x_query_stream::Queries::Taken(store),
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-138
    // ```c
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    //         }
    // ```
    // The report and the query batches (`GetDbTotalLength`) take the total through an
    // `int`, which LOSAT does not reproduce beyond its range. (A database's total is the
    // `Int8` of `GetBlastDbInfo`, blast_format.cpp:140-154.)
    if let Subjects::Records(subject_set) = subjects {
        let subject_letters: usize = (0..subject_set.records().len())
            .map(|index| subject_set.len(index))
            .sum();
        if subject_letters > i32::MAX as usize {
            anyhow::bail!(
                "the subjects have {subject_letters} letters; NCBI BLAST+ reports a total of 2^31 letters or more through a 32-bit int, which is not supported by LOSAT's BLASTN"
            );
        }
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let x_result = crate::utils::threading::with_search_pool(args.num_threads, "blastn", |pool| {
        run_in_pool(
            args,
            queries,
            reading,
            &query_input,
            subjects,
            ranges.subjects,
            outputs,
            output_formats,
            pool,
            batching,
        )
    });
    // No NCBI counterpart: the PAIRPAR counters with LOSAT_TIMING; nothing without it.
    x_pair_seed::summary();
    x_result
}

/// Rejects the environment variables that put NCBI BLAST+ into a search mode or report that
/// LOSAT does not reproduce.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:204-210
/// ```c
///     if ( (subjects = db_args->GetSubjects(scope)) ) {
///         _ASSERT(search_db.Empty());
/// 	char* bl2seq_legacy = getenv("BL2SEQ_LEGACY");
/// 	if (bl2seq_legacy)
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, false));
/// 	else
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
/// ```
/// `BL2SEQ_LEGACY` (any value, also an empty one) turns off the scan of the subjects as a
/// database: `CLocalBlast` then searches each subject on its own (local_blast.cpp:189,289)
/// and `CBlastFormat` writes the legacy bl2seq report (`m_IsBl2Seq && !m_IsDbScan`). NCBI
/// reads it in that branch only: a database search (`database`) does not read it.
fn check_unsupported_environment(database: bool) -> Result<()> {
    if !database && std::env::var_os("BL2SEQ_LEGACY").is_some() {
        anyhow::bail!(
            "the environment variable BL2SEQ_LEGACY, which makes NCBI BLAST+ search each subject on its own and write its legacy bl2seq report, is not supported by LOSAT's BLASTN"
        );
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:732-737
    // ```c
    // 		char * pre_fetch_limit_str = getenv("PRE_FETCH_SEQS_LIMIT");
    // 		if (pre_fetch_limit_str) {
    // 			int pre_fetch_limit = NStr::StringToInt(pre_fetch_limit_str);
    // 			if(pre_fetch_limit == 0) {
    // 				return false;
    // 			}
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:761-766
    // ```c
    //     if (results.size() == 0) {
    //         return;
    //     }
    //     if(!s_PreFetchSeqs(results, format_type)){
    //     	return;
    //     }
    // ```
    // An integer only decides whether the report's sequences are fetched ahead, which
    // changes no output. Before the report of each query batch NCBI converts the value, and
    // one that `NStr::StringToInt` cannot convert (an empty one too) stops it with a
    // `CStringException` whose text names its build's source files (exit 255).
    if let Some(value) = std::env::var_os("PRE_FETCH_SEQS_LIMIT") {
        if crate::blastinput::query_batch::ncbi_string_to_int(&value).is_none() {
            anyhow::bail!(
                "the environment variable PRE_FETCH_SEQS_LIMIT has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN",
                value.to_string_lossy()
            );
        }
    }
    Ok(())
}

/// An error of a query batch after the first that does not depend on the batch's queries
/// (the Karlin-Altschul error of the scoring options, or the average subject length;
/// `search_query_batch`): the first batch with a valid query meets it, and NCBI raises it
/// after the reports of the batches before, whose queries were all invalid.
#[derive(Debug)]
pub(super) struct ErrorAfterInvalidBatches(pub(super) anyhow::Error);

impl std::fmt::Display for ErrorAfterInvalidBatches {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ErrorAfterInvalidBatches {}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:295-312
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
// 		        lcl_blast.SetBatchNumber( batch_num );
//                 results = lcl_blast.Run();
//     ...
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
/// A database error of a query batch's search (`database::search_blocks`,
/// `database::traceback_listed`) as NCBI reports it: a sequence file that the preliminary
/// search cannot read is NCBI's `BLAST engine error` (`database::PreliminaryReadError`,
/// local_blast.cpp:240-251), one that the traceback cannot read its `BLAST Database
/// error`. NCBI raises it in the batch's search, after the outfmt 0 prolog and the reports
/// of the batches before (`ErrorAfterInvalidBatches`; for the first batch, the prolog of a
/// file sink is written here, the others were written before the search). LOSAT's
/// refusals are errors of their own.
pub(super) fn database_search_error(
    error: anyhow::Error,
    first_batch: bool,
    outputs: BatchOutputs<'_, '_>,
    megablast: bool,
    shown: &ShownSubjects,
) -> anyhow::Error {
    let error = match error.downcast::<database::PreliminaryReadError>() {
        Ok(read_error) => crate::blastinput::app::engine_error(read_error.0),
        Err(error) => error,
    };
    if !error.is::<crate::cli::NativeError>() {
        return error;
    }
    if !first_batch {
        return ErrorAfterInvalidBatches(error).into();
    }
    let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
    match write_pairwise_prologs(outputs, output_formats, megablast, shown, false) {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// Reads the environment variables that set NCBI's query batches and chunks, in NCBI's
/// order: `CHUNK_SIZE` (the `CBatchSizeMixer` maximum), `BATCH_SIZE`, then
/// `OVERLAP_CHUNK_SIZE` (the query splitter of the first batch).
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-262
/// ```c
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
///         int batch_size = m_CmdLineArgs->GetQueryBatchSize();
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:55-59
/// ```c
///     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
///     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
///     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
///     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
///         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
/// ```
/// A value that `NStr::StringToInt` cannot convert makes NCBI stop with its
/// `CStringException`, whose text names the source path of the oracle's build (exit 255);
/// LOSAT rejects it. A negative `CHUNK_SIZE` is a `size_t` near 2^64; LOSAT rejects it where
/// it splits a query batch (`search_query_chunks`), and otherwise searches as NCBI.
fn query_batching_from_environment(megablast: bool) -> Result<QueryBatching> {
    let read = |variable: &str| -> Result<Option<i32>> {
        match std::env::var_os(variable) {
            Some(value) if !is_blank(value.as_encoded_bytes()) => {
                match crate::blastinput::query_batch::ncbi_string_to_int(&value) {
                    Some(number) => Ok(Some(number)),
                    None => anyhow::bail!(
                        "the environment variable {variable} has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN",
                        value.to_string_lossy()
                    ),
                }
            }
            _ => Ok(None),
        }
    };
    let chunk_size = read("CHUNK_SIZE")?;
    let fixed_batch_size = crate::blastinput::query_batch::get_query_batch_size(
        std::env::var_os("BATCH_SIZE").as_deref(),
    )
    .map_err(|value| {
        anyhow::anyhow!(
            "the environment variable BATCH_SIZE has the value {value:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN"
        )
    })?;
    let overlap = read("OVERLAP_CHUNK_SIZE")?;
    let split = SplitSizes::new(megablast, chunk_size, overlap);
    Ok(QueryBatching {
        fixed_batch_size,
        batch_size: fixed_batch_size as u32,
        split,
    })
}

/// Searches the queries in NCBI's query batches and writes the results of all of them.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-300
/// ```c++
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
///         int batch_size = m_CmdLineArgs->GetQueryBatchSize();
///         if (batch_size) {
///             input.SetBatchSize(batch_size);
///         } else {
///             Int8 total_len = formatter.GetDbTotalLength();
///             if (total_len > 0) {
///                 /* the optimal hits per batch scales with total db size */
///                 mixer.SetTargetHits(total_len / 3000);
///             }
///             input.SetBatchSize(mixer.GetBatchSize());
///         }
///         ...
///         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
///             ...
///                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
///                 ...
///                 results = lcl_blast.Run();
///                 if (!batch_size)
///                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
/// ```
/// `GetQueryBatchSize` is the `BATCH_SIZE` value (0 when it is not set) and the chunk size
/// that of `CHUNK_SIZE` or the task (`query_batching_from_environment`).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: BlastnArgs,
    queries: x_query_stream::Queries<'_>,
    reading: QueryReading,
    query_input: &QueryInput,
    subjects: Subjects<'_>,
    subject_ranges: SubjectRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    mut batching: QueryBatching,
) -> Result<()> {
    // EXPERIMENT (LOSAT_X_QBATCHPAR): its summary line comes at the end of every search with the
    // switch and LOSAT_X_STATS set, also one that ends before its batch loop.
    let qbp_shared =
        x_qbatch_par::Shared::<Result<QueryBatch>>::new(x_qbatch_par::Switches::from_env());
    let prepared: PreparedSubjects<'_>;
    let source = match subjects {
        Subjects::Records(subject_set) => {
            let subject_records = subject_set.records();
            let metadata =
                subject_metadata_from_records(subject_records, |index| subject_set.len(index));
            if metadata.db_num_seqs == 0 {
                // EXPERIMENT (LOSAT_X_QSTREAM): the records are read again for their warnings.
                queries.write_input_warnings(outputs.diagnostics)?;
                return Ok(());
            }
            // Diagnostic (LOSAT_TIMING): the subjects' packing for the common line
            // (utils/xtiming.rs). No NCBI counterpart.
            let t_encode = crate::utils::xtiming::clock();
            // The packed subjects hold the packed letters (LOSAT_X_BNSTORE).
            let packed = if subject_set.stores().is_some() {
                Vec::new()
            } else {
                subject_records
                    .iter()
                    .map(|record| encode_subject_ncbi2na_packed(record.seq()))
                    .collect()
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
            prepared = PreparedSubjects {
                records: subject_records,
                stores: subject_set.stores(),
                packed,
                metadata,
            };
            SubjectSource::Prepared(&prepared)
        }
        // A database without sequences fails in the first searched batch
        // (`search_query_batch`).
        Subjects::Database(database) => SubjectSource::Database(database),
    };
    // The subjects as the reports show them, and the stored titles of the database
    // sequences that the final hit lists of the batches name (`database::read_titles`).
    let shown = ShownSubjects::new(&args, source);
    let mut titles = database::DatabaseTitles::new();
    let megablast = args.task == "megablast";
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:258-261
    // ```c
    //         formatter.PrintProlog();
    //
    //         /*** Process the input ***/
    //         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
    // ```
    // NCBI writes the outfmt 0 prolog before it reads the first query batch, and its stream
    // throws at the prolog's first flush when the write fails (blast_format.cpp:118-119),
    // before any query is read or searched.
    write_pairwise_prologs(outputs, &output_formats, megablast, &shown, true)?;
    let mut mixer = BatchSizeMixer::new(batching.split.mixer_max_batch_size());
    if batching.fixed_batch_size == 0 {
        // A database's total is the `Int8` of its index files (`GetDbTotalLength`).
        let total_length = source.metadata().db_len_total as i64;
        if total_length > 0 {
            mixer.set_target_hits((total_length / 3000) as i32);
        }
        // `input.SetBatchSize(mixer.GetBatchSize())` converts the `Int4` to `TSeqPos`.
        batching.batch_size = mixer.batch_size(None) as u32;
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-280
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    //             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    // A reader at the end of its input reads no batch (`!input.End()`): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    if queries.input_count() == 0 && matches!(reading, QueryReading::End) {
        return write_batch_reports(
            &args,
            &[],
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &ReportRanges {
                queries: query_input,
                subjects: subject_ranges,
            },
            outputs,
            &output_formats,
            BatchResults {
                hit_lists: Vec::new(),
                warnings: &[],
                lengths: &[],
                unsearched: Vec::new(),
                query_masks: &[],
                context_karlin: &[],
                query_eff_searchsp: &[],
                batches: &[],
            },
            true,
            None,
        );
    }
    // A batch size of 0 (`CHUNK_SIZE` 1000 without `BATCH_SIZE`: the mixer's maximum is
    // 0) reads no query, so the first batch fails after the prolog (outfmt 0).
    if batching.batch_size == 0 {
        write_pairwise_prologs(outputs, &output_formats, megablast, &shown, false)?;
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }

    // The lengths of the searched queries (their intervals with `-query_loc`; with
    // LOSAT_X_QSTREAM from the first reading's lengths).
    let lengths: Vec<usize> = queries.searched_lengths();
    let report_ranges = ReportRanges {
        queries: query_input,
        subjects: subject_ranges,
    };
    let query_offsets: Vec<usize> = (0..lengths.len())
        .map(|index| query_input.placements.offset(index))
        .collect();
    // EXPERIMENT (LOSAT_X_QBATCHPAR): the batches may be searched ahead of this loop, in slots on
    // other threads (`x_qbatch_par.rs`); the loop computes every batch's key itself and takes a
    // slot's result only under that key. A slot's search is the loop's call for that batch
    // without the outputs (it never searches the first batch).
    let qbp_off = if qbp_shared.switches().on {
        x_qbatch_par_off(&args, source, parallel_pool.threads())
    } else {
        Some(x_qbatch_par::Off::NotSet)
    };
    let qbp_batching = batching;
    let qbp_search = |job: &x_qbatch_par::Job, pool: &crate::utils::threading::SearchPool<'_>| {
        // EXPERIMENT (LOSAT_X_QSTREAM): the batch's records read again (a chunk shared with the
        // window, or a copy); none for a batch the loop has passed, which the loop never takes
        // (cancelled).
        let records = match queries.slot_records(job.searched.clone()) {
            Ok(Some(records)) => records,
            Ok(None) => return Err(x_qbatch_par::Cancelled.into()),
            Err(error) => return Err(error),
        };
        search_query_batch(
            &args,
            x_query_drop::BatchQueries::Shared(&records),
            source,
            None,
            pool,
            QueryBatching {
                batch_size: job.key.batch_size,
                ..qbp_batching
            },
            false,
            BatchStage::Search {
                query_offsets: &query_offsets[job.searched.clone()],
            },
            None,
        )
    };
    let qbp_input = x_qbatch_par::PlanInput {
        query_input,
        lengths: &lengths,
        reading_error: matches!(reading, QueryReading::Error { .. }),
        fixed_batch_size: batching.fixed_batch_size != 0,
        split: batching.split,
    };
    // EXPERIMENT (LOSAT_X_RPTSTREAM): whether the loop writes each batch's reports when it takes
    // the batch (`x_rpt_stream.rs`); the summary lines come when it is dropped, at the end.
    let mut stream = x_rpt_stream::Stream::start(
        &args,
        outputs,
        &output_formats,
        source,
        subject_ranges,
        query_input,
        &lengths,
        batching,
        parallel_pool,
    );
    std::thread::scope(|scope| {
        let mut qbp = x_qbatch_par::Driver::new(
            &qbp_shared,
            &qbp_search,
            scope,
            qbp_input,
            parallel_pool.threads(),
            x_qbatch_par::Profile::of_search(
                &args.task,
                crate::algorithm::blastn::coordination::determine_effective_word_size(&args),
                x_query_file_bytes(&args.query),
                x_subject_letters(source),
            ),
            Box::new(|letters: u64, max_len: usize| {
                let mut config = configure_task(&args);
                let discontig_template = config.mb_template_length > 0;
                finalize_task_config(
                    &mut config,
                    2 * letters as usize,
                    max_len.max(1),
                    discontig_template,
                );
                config.lut_word_length
            }),
            qbp_off,
        );
        // EXPERIMENT (LOSAT_X_BNQDROP): no slot searches the records that the loop takes out of
        // their store (`search_cli`'s rule moved them there where it found that no slot may run).
        if matches!(queries, x_query_stream::Queries::Taken(_)) {
            qbp.x_refuse_taken();
        }
        // EXPERIMENT (LOSAT_X_MBLUTKEEP): the megablast table this loop keeps between its batches
        // (`lookup::x_mb_lut_keep`), freed when the loop ends.
        let x_lut_keep = x_mb_lut_keep::Keep::for_loop();
        let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(lengths.len());
        let mut batch_spans: Vec<std::ops::Range<usize>> = Vec::new();
        let mut unsearched = Vec::with_capacity(lengths.len());
        let mut query_masks = Vec::with_capacity(lengths.len());
        let mut context_karlin = Vec::with_capacity(2 * lengths.len());
        // The warnings written before each query's report (`QueryWarnings`).
        let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); lengths.len()];
        let mut query_eff_searchsp = Vec::with_capacity(2 * lengths.len());
        // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQDROP): the searched records of the batches taken,
        // kept where a writer reads them after their batch: the reports after the loop and the
        // error paths that write the reports of the batches before (LOSAT_X_RPTSTREAM not acting).
        // Without the switches the records read before the search are those.
        let keep = queries.per_batch() && !stream.acting();
        let mut kept: Vec<FastaRecord> = Vec::new();
        // A batch is read from the input records (`input_start..input_end`); its searched
        // queries are `start..end` (all of them without `-query_loc`). A batch size counts the
        // whole records, and a record whose range starts past its end is skipped
        // (`next_ranged_batch`).
        let mut reading = reading;
        // Whether what ended the reading after the records (`reading`) has come.
        let mut reading_ended = false;
        let mut input_start = 0;
        // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-279
        // ```c
        //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
        // 	    BLAST_PROF_START( APP.LOOP.PRE );
        //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
        // ```
        // The reader is at the end of its input after the last record, unless an `eEOF` or an
        // error ends the reading there; that comes in the batch being read when the records
        // run out before the batch reaches its size, and otherwise in a batch of its own.
        while input_start < query_input.input_lengths.len()
            || (reading.reads_past_records() && !reading_ended)
        {
            let (input_end, reached_size) = next_ranged_batch(
                &query_input.input_lengths,
                &query_input.skipped,
                input_start,
                batching.batch_size,
            );
            let std::ops::Range { start, end } = query_input.searched_in(input_start..input_end);
            let reading_ends_here = input_end == query_input.input_lengths.len()
                && !reached_size
                && reading.reads_past_records();
            // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-279
            // ```c
            //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
            // 	    BLAST_PROF_START( APP.LOOP.PRE );
            //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
            // ```
            // `CFastaReader` writes its messages about a batch's queries (their lines and
            // titles) when it reads them, after the report of the batch before: before the
            // report of the batch's first query (`QueryWarnings`). Skipped records are read too.
            // EXPERIMENT (LOSAT_X_QSTREAM): the batch's records, read now from the query file when
            // the switch acts (an error here is the file changed between the two readings).
            let mut batch_records = match queries.batch(input_start..input_end, start..end) {
                Ok(records) => records,
                Err(error) => {
                    qbp.finish();
                    return Err(error);
                }
            };
            let batch_title_warnings: Vec<u8> = batch_records
                .input()
                .iter()
                .flat_map(|record| record.warnings.iter().copied())
                .collect();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:145-152
            // ```c
            //         try { q.Reset(m_Source->GetNextSequence(scope)); }
            //         catch (const CObjReaderParseException& e) {
            //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
            //                 break;
            //             }
            //             throw;
            //         }
            // ```
            // An `eEOF` ends the batch with the records read so far. A reader error (or LOSAT's
            // rejection of a Seq-id line) stops the run while the batch is read: after the
            // reports of the batches before (without the epilog) and the messages of the
            // batch's records, which are neither searched nor reported.
            if reading_ends_here {
                reading_ended = true;
                if let QueryReading::Error {
                    error,
                    warnings: error_warnings,
                } = std::mem::replace(&mut reading, QueryReading::End)
                {
                    // EXPERIMENT (LOSAT_X_QBATCHPAR): no slot runs while the loop writes.
                    qbp.finish();
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                    // written when the loop took them.
                    let before = |outputs: &mut ReportOutputs<'_>| {
                        write_before_failed_batch(
                            &args,
                            queries.before(&kept, start),
                            ReportSubjects {
                                source,
                                shown: &shown,
                                titles: &titles,
                            },
                            &report_ranges,
                            outputs,
                            &output_formats,
                            BatchResults {
                                hit_lists,
                                warnings: &warnings[..start],
                                lengths: &lengths[..start],
                                unsearched,
                                query_masks: &query_masks,
                                context_karlin: &context_karlin,
                                query_eff_searchsp: &query_eff_searchsp,
                                batches: &batch_spans,
                            },
                        )
                    };
                    if stream.acting() {
                        stream.finish_sam(outputs)?;
                    } else {
                        before(outputs)?;
                    }
                    outputs.diagnostics.write_all(&batch_title_warnings)?;
                    outputs.diagnostics.write_all(&error_warnings)?;
                    return Err(error.into_app_error());
                }
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
            // ```c
            // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
            //     : m_QueryVector(& queries)
            // {
            //     if (queries.Empty()) {
            //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
            //     }
            // ```
            // A batch of skipped records only (with `-query_loc`), or a batch that an `eEOF`
            // ends before its first record, fails after the reports of the batches before,
            // without the epilog.
            if start == end {
                // EXPERIMENT (LOSAT_X_QBATCHPAR): no slot runs while the loop writes.
                qbp.finish();
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were written
                // when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                if stream.acting() {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
                outputs.diagnostics.write_all(&batch_title_warnings)?;
                return Err(crate::cli::NativeError {
                    exit: 3,
                    message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
                }
                .into());
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
            // ```c
            //         } catch (const CException& e) {
            //             ...
            //             CRef<CSearchMessage> m
            //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
            //             messages[index].push_back(m);
            //             s_InvalidateQueryContexts(qinfo, index);
            //         }
            //     ...
            //     // Validate that at least one query context is valid
            //     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
            //         NCBI_THROW(CBlastException, eSetup, messages.ToString());
            //     }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
            // ```c
            // TSearchMessages::ToString() const
            // {
            //     string retval;
            //     ITERATE(vector<TQueryMessages>, qm, *this) {
            //         if (qm->empty()) {
            //             continue;
            //         }
            //         ITERATE(TQueryMessages, msg, *qm) {
            //             retval += (*msg)->GetMessage() + " ";
            //         }
            //     }
            //     return retval;
            // }
            // ```
            // A query without letters (`Sequence contains no data`, blast_setup.hpp:190-198) has
            // no valid context, and the strands of every other BLASTN query are valid when the
            // queries are set up: a batch of such queries only stops NCBI with one message per
            // query (`GetMessage` puts the severity first), after the reports of the batches
            // before (`CATCH_ALL`: `BLAST engine error: `, exit 3).
            if batch_records
                .searched()
                .iter()
                .all(|record| record.seq().is_empty())
            {
                // EXPERIMENT (LOSAT_X_QBATCHPAR): no slot runs while the loop writes.
                qbp.finish();
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were written
                // when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                if stream.acting() {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
                outputs.diagnostics.write_all(&batch_title_warnings)?;
                let mut message = String::from("BLAST engine error: ");
                for _ in start..end {
                    message.push_str("Warning: ");
                    message.push_str(NO_DATA_MESSAGE);
                    message.push(' ');
                }
                message.push('\n');
                return Err(crate::cli::NativeError { exit: 3, message }.into());
            }
            warnings[start].extend(batch_title_warnings);
            // EXPERIMENT (LOSAT_X_QBATCHPAR): the batch's result from a slot that searched the same
            // key ahead of the loop, or the loop's own search (on the whole pool, no slot running).
            let key = x_qbatch_par::BatchKey {
                input_start,
                input_end,
                batch_size: batching.batch_size,
            };
            let searched = match qbp.take(key) {
                Some(precomputed) => precomputed,
                None => {
                    // EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop keeps its megablast table from batch
                    // to batch once no slot can search a batch any more.
                    let x_keep = x_lut_keep.as_ref().filter(|_| !qbp.may_speculate());
                    if x_lut_keep.is_some() && x_keep.is_none() {
                        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", false);
                    }
                    search_query_batch(
                        &args,
                        batch_records.batch_queries(),
                        source,
                        Some((&mut *outputs, &output_formats)),
                        parallel_pool,
                        batching,
                        input_start == 0,
                        BatchStage::Search {
                            query_offsets: &query_offsets[start..end],
                        },
                        x_keep,
                    )
                }
            };
            // EXPERIMENT (LOSAT_X_QBATCHPAR): an error ends the loop here, and what it writes
            // comes after every slot has stopped.
            if searched.is_err() {
                qbp.finish();
            }
            let batch = match searched {
                Ok(batch) => batch,
                Err(error) => match error.downcast::<ErrorAfterInvalidBatches>() {
                    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:295-318
                    // ```c
                    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
                    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                    // 		        lcl_blast.SetBatchNumber( batch_num );
                    //                 results = lcl_blast.Run();
                    //     ...
                    //                 ITERATE(CSearchResultSet, result, *results) {
                    //                     formatter.PrintOneResultSet(**result, query_batch);
                    //                 }
                    //     ...
                    //         }
                    //         BLAST_PROF_START( APP.POST );
                    //         formatter.PrintEpilog(opt);
                    // ```
                    // The reports of the batches before are written, then the error, which
                    // skips `PrintEpilog` (the outfmt 0 footer and outfmt 7's last line).
                    Ok(ErrorAfterInvalidBatches(error)) => {
                        // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                        // written when the loop took them.
                        let before = |outputs: &mut ReportOutputs<'_>| {
                            write_batch_reports(
                                &args,
                                queries.before(&kept, start),
                                ReportSubjects {
                                    source,
                                    shown: &shown,
                                    titles: &titles,
                                },
                                &report_ranges,
                                outputs,
                                &output_formats,
                                BatchResults {
                                    hit_lists,
                                    warnings: &warnings[..start],
                                    lengths: &lengths[..start],
                                    unsearched,
                                    query_masks: &query_masks,
                                    context_karlin: &context_karlin,
                                    query_eff_searchsp: &query_eff_searchsp,
                                    batches: &batch_spans,
                                },
                                false,
                                None,
                            )
                        };
                        if stream.acting() {
                            stream.finish_sam(outputs)?;
                        } else {
                            before(outputs)?;
                        }
                        // The failing batch was read (its title warnings) before its search.
                        outputs.diagnostics.write_all(&warnings[start])?;
                        return Err(error);
                    }
                    // The warnings of the batches read so far, in NCBI's order (the outfmt 0
                    // prolog of a Karlin-Altschul error of the first batch is written).
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were written
                    // with their reports; the off path writes no report here (the switch acts only
                    // where LOSAT's refusals cannot come after the first batch, `x_rpt_stream`).
                    Err(error) => {
                        let written = if stream.acting() { start } else { 0 };
                        for query_warnings in &warnings[written..=start] {
                            outputs.diagnostics.write_all(query_warnings)?;
                        }
                        return Err(error);
                    }
                },
            };
            // The batch's LOSAT_TIMING block, where its search wrote it before.
            if !batch.timing.is_empty() {
                eprint!("{}", batch.timing);
            }
            let mut batch_hit_lists = batch.hit_lists;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
            // ```c++
            //     CRef<ILocalQueryData> qdata = m_QueryFactory->MakeLocalQueryData(m_Options);
            //
            //     vector<TSeqLocInfoVector> subj_masks;
            //     TSeqAlignVector aligns =
            //         LocalBlastResults2SeqAlign(hsp_results,
            //                                    *qdata,
            //                                    *m_SeqInfoSrc,
            //                                    m_OptsMemento->m_ProgramType,
            //                                    m_Options->GetGappedMode(),
            //                                    m_Options->GetOutOfFrameMode(),
            //                                    subj_masks,
            //                                    m_ResultType);
            // ```
            // A database: the batch's hit lists are finished now (`finish_query_hit_list`, the
            // end of NCBI's traceback), and the definition lines of the sequences they name are
            // read in NCBI's order where NCBI makes the batch's Seq-aligns
            // (`database::read_titles`), before the HSPs of score 0 are dropped. A header file
            // that NCBI cannot read fails this batch, after the reports of the batches before
            // and the messages of this batch's reading.
            if let SubjectSource::Database(database) = source {
                let hitlist_size = args.report_counts.hitlist_size;
                for (q_idx, hit_list) in batch_hit_lists.iter_mut().enumerate() {
                    if let Some(hit_list) = hit_list {
                        finish_query_hit_list(
                            hit_list,
                            lengths[start + q_idx],
                            hitlist_size,
                            args.max_hsps_per_subject.unwrap_or(0),
                            args.qcov_hsp_perc.unwrap_or(0.0),
                            args.subject_besthit,
                            args.culling_limit_value(),
                        );
                    }
                }
                if let Err(error) = database::read_titles(database, &batch_hit_lists, &mut titles) {
                    // EXPERIMENT (LOSAT_X_QBATCHPAR): no slot runs while the loop writes (the
                    // path is off with a database; for uniformity).
                    qbp.finish();
                    let error = match error.downcast::<database::HeaderReadError>() {
                        Ok(failure) => failure.error,
                        // LOSAT's refusal: the warnings of the batches read so far, as for the
                        // other errors of a batch (below).
                        // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were
                        // written with their reports; the off path writes no report here (the
                        // switch acts only where every title is read before the search,
                        // `x_rpt_stream`).
                        Err(error) => {
                            let written = if stream.acting() { start } else { 0 };
                            for query_warnings in &warnings[written..=start] {
                                outputs.diagnostics.write_all(query_warnings)?;
                            }
                            return Err(error);
                        }
                    };
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                    // written when the loop took them.
                    let before = |outputs: &mut ReportOutputs<'_>| {
                        write_before_failed_batch(
                            &args,
                            queries.before(&kept, start),
                            ReportSubjects {
                                source,
                                shown: &shown,
                                titles: &titles,
                            },
                            &report_ranges,
                            outputs,
                            &output_formats,
                            BatchResults {
                                hit_lists,
                                warnings: &warnings[..start],
                                lengths: &lengths[..start],
                                unsearched,
                                query_masks: &query_masks,
                                context_karlin: &context_karlin,
                                query_eff_searchsp: &query_eff_searchsp,
                                batches: &batch_spans,
                            },
                        )
                    };
                    if stream.acting() {
                        stream.finish_sam(outputs)?;
                    } else {
                        before(outputs)?;
                    }
                    outputs.diagnostics.write_all(&warnings[start])?;
                    return Err(error);
                }
                for hit_list in batch_hit_lists.iter_mut().flatten() {
                    drop_hsps_of_score_zero(hit_list);
                }
            }
            // The batch numbers its queries from 0.
            let first = start as u32;
            hit_lists.extend(batch_hit_lists.into_iter().map(|list| {
                list.map(|mut list| {
                    for hsp_list in &mut list.hsplist_array {
                        hsp_list.query_index += first;
                        for hsp in &mut hsp_list.hsps {
                            hsp.q_idx += first;
                        }
                    }
                    list
                })
            }));
            batch_spans.push(start..end);
            unsearched.extend(std::iter::repeat_n(!batch.searched, end - start));
            query_masks.extend(batch.query_masks);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
            // ```c
            //     if (results.HasWarnings()) {
            //         ERR_POST(Warning << results.GetWarningStrings());
            //     }
            // ```
            // The warning of an invalid query comes with its report (`PrintOneResultSet`),
            // before its preamble (`QueryWarnings`): a query without letters has the message of
            // its set-up (`SetupQueries_OMF`, above), any other invalid query that of its
            // Karlin-Altschul parameters.
            // A structured format (outfmt 8, 12) writes no result warnings
            // (`ReportFormat::writes_query_warnings`); the warnings go with the first format.
            let result_warnings = output_formats
                .first()
                .is_none_or(|format| format.writes_query_warnings());
            // The query's first valid context (`-strand minus`: its minus context).
            let first_context = QueryStrand::from_option(&args.strand).first_context();
            for (index, karlin) in batch
                .context_karlin
                .iter()
                .skip(first_context)
                .step_by(2)
                .enumerate()
            {
                if karlin.is_none() && result_warnings {
                    let query = &batch_records.searched()[index];
                    let ordinal = query_input.ordinal(start + index);
                    warnings[start + index].extend(if lengths[start + index] == 0 {
                        query_warning("blastn", ordinal, query, &[NO_DATA_MESSAGE.to_string()])
                    } else {
                        invalid_query_warning("blastn", ordinal, query)
                    });
                }
            }
            context_karlin.extend(batch.context_karlin);
            query_eff_searchsp.extend(batch.query_eff_searchsp);
            // NCBI's `Int4` count of the batch's successful initial extensions.
            if batching.fixed_batch_size == 0 {
                batching.batch_size =
                    mixer.batch_size(Some(batch.good_init_extends as u32 as i32)) as u32;
            }
            // EXPERIMENT (LOSAT_X_QBATCHPAR): the true gie and next size, from which the batches
            // ahead are planned again.
            qbp.commit(
                key,
                batch.good_init_extends,
                (batching.fixed_batch_size == 0).then_some(&mixer),
                batching.batch_size,
            );
            // EXPERIMENT (LOSAT_X_RPTSTREAM): the batch's reports, written now, after the mixer
            // has its gie, as NCBI's loop prints a batch's result sets (`x_rpt_stream.rs`). The
            // batch's lists and values are taken back from the end of what the loop keeps, numbered
            // from 0 again, and written with the batch's part of the query ranges.
            if stream.acting() {
                let count = end - start;
                let at = hit_lists.len() - count;
                fn take<T>(kept: &mut Vec<T>, at: usize) -> Vec<T> {
                    kept.drain(at..).collect()
                }
                let mut batch_lists = take(&mut hit_lists, at);
                for hsp_list in batch_lists
                    .iter_mut()
                    .flatten()
                    .flat_map(|list| list.hsplist_array.iter_mut())
                {
                    hsp_list.query_index -= first;
                    for hsp in &mut hsp_list.hsps {
                        hsp.q_idx -= first;
                    }
                }
                let batch_unsearched = take(&mut unsearched, at);
                let batch_masks = take(&mut query_masks, at);
                let batch_karlin = take(&mut context_karlin, 2 * at);
                let batch_searchsp = take(&mut query_eff_searchsp, 2 * at);
                let batch_input = x_rpt_stream::batch_input(query_input, start..end);
                let batch_ranges = ReportRanges {
                    queries: &batch_input,
                    subjects: subject_ranges,
                };
                // EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool formats the batch's outfmt 0 report while
                // it is idle (no thread of a LOSAT_X_QBATCHPAR slot exists).
                let x_pool = qbp.x_slots_idle().then_some(parallel_pool);
                stream.write(outputs, count, false, |outputs, written| {
                    write_reports(
                        &args,
                        batch_records.searched(),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &batch_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists: batch_lists,
                            warnings: &warnings[start..end],
                            lengths: &lengths[start..end],
                            unsearched: batch_unsearched,
                            query_masks: &batch_masks,
                            context_karlin: &batch_karlin,
                            query_eff_searchsp: &batch_searchsp,
                            batches: &[],
                        },
                        false,
                        written,
                        x_pool,
                    )
                })?;
                // The written queries' warnings are not read again (the error paths write those
                // of the current batch).
                for query_warnings in &mut warnings[start..end] {
                    *query_warnings = Vec::new();
                }
            }
            // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQDROP): the batch's records are kept for a later
            // writer (moved when the loop took them out of their store), or dropped (the batch is
            // committed: the slots now search batches after it).
            if keep {
                let searched = batch_records.into_searched();
                queries.note_kept(searched.len());
                kept.extend(searched);
            } else {
                drop(batch_records);
            }
            queries.done(input_end);
            input_start = input_end;
        }
        // EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop's kept table is freed with the loop.
        drop(x_lut_keep);

        // EXPERIMENT (LOSAT_X_QBATCHPAR): every slot has stopped before the reports are written.
        qbp.finish();
        let reports = |outputs: &mut ReportOutputs<'_>| {
            write_batch_reports(
                &args,
                queries.all(&kept),
                ReportSubjects {
                    source,
                    shown: &shown,
                    titles: &titles,
                },
                &report_ranges,
                outputs,
                &output_formats,
                BatchResults {
                    hit_lists,
                    warnings: &warnings,
                    lengths: &lengths,
                    unsearched,
                    query_masks: &query_masks,
                    context_karlin: &context_karlin,
                    query_eff_searchsp: &query_eff_searchsp,
                    batches: &batch_spans,
                },
                true,
                Some(parallel_pool),
            )
        };
        if !stream.acting() {
            return reports(outputs);
        }
        // EXPERIMENT (LOSAT_X_RPTSTREAM): every batch's reports are written; the epilogs remain
        // (NCBI's `PrintEpilog` after its loop), written by the report code for no query, with
        // the queries of every batch counted (outfmt 7's last line).
        stream.write(outputs, 0, true, |outputs, written| {
            write_reports(
                &args,
                &[],
                ReportSubjects {
                    source,
                    shown: &shown,
                    titles: &titles,
                },
                &report_ranges,
                outputs,
                &output_formats,
                BatchResults {
                    hit_lists: Vec::new(),
                    warnings: &[],
                    lengths: &[],
                    unsearched: Vec::new(),
                    query_masks: &[],
                    context_karlin: &[],
                    query_eff_searchsp: &[],
                    batches: &[],
                },
                true,
                written,
                None,
            )
        })?;
        Ok(())
    })
}

// No NCBI counterpart: LOSAT_X_QBATCHPAR schedules NCBI's query batches; it changes no value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_QBATCHPAR): why the batch loop of a search searches every batch itself
/// (`None`: slots may search batches ahead of it once its first batch is done;
/// `x_qbatch_par.rs`). A database keeps the loop (each batch reads the database's blocks); so
/// do verbose, debug and trace output, which a batch writes while it is searched, and subjects
/// of at least as many 5 Mb work units as threads, whose batches already keep the pool busy.
/// The subjects' lengths come from their letters (`SubjectLetters::len`), which
/// LOSAT_X_BNSTORE holds packed while the records keep none.
fn x_qbatch_par_off(
    args: &BlastnArgs,
    source: SubjectSource<'_>,
    threads: usize,
) -> Option<x_qbatch_par::Off> {
    let subject_lengths = source.prepared().map(|prepared| {
        let letters = prepared.letters();
        (0..prepared.records.len()).map(move |index| letters.len(index))
    });
    x_qbatch_par_off_of(args, subject_lengths, threads)
}

// No NCBI counterpart (scheduling only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QBATCHPAR): `x_qbatch_par_off` from the letters of each subject
/// (`subject_lengths`; `None` for a database), so that LOSAT_X_QSTREAM's rule reads the same
/// conditions before the subjects are prepared.
fn x_qbatch_par_off_of(
    args: &BlastnArgs,
    subject_lengths: Option<impl Iterator<Item = usize>>,
    threads: usize,
) -> Option<x_qbatch_par::Off> {
    use x_qbatch_par::Off;
    if !cfg!(all(feature = "parallel", not(target_arch = "wasm32"))) {
        return Some(Off::Build);
    }
    if threads < 2 {
        return Some(Off::Threads);
    }
    let Some(subject_lengths) = subject_lengths else {
        return Some(Off::Database);
    };
    if x_debug_output(args) {
        return Some(Off::Debug);
    }
    let units: usize = subject_lengths
        .map(|length| length.div_ceil(MAX_DBSEQ_LEN))
        .sum();
    if units >= threads {
        return Some(Off::Subjects);
    }
    None
}

// No NCBI counterpart (scheduling only): the size NCBI's blastn reads for its split by queries
// (`x_qbatch_par::ncbi_splits_by_queries`); it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QBATCHPAR): the length of the query file in bytes, as NCBI's `CFile`
/// gives it (`stat`, through symbolic links); `None` for standard input (`-`) or a path that
/// cannot be read.
fn x_query_file_bytes(query: &std::path::Path) -> Option<u64> {
    if query == std::path::Path::new("-") {
        return None;
    }
    std::fs::metadata(query).ok().map(|metadata| metadata.len())
}

// No NCBI counterpart (when LOSAT reads the records only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QSTREAM): the letters of the largest query batch of a search of `task`, as
/// `run_in_pool` reads its batches (`query_batching_from_environment`): `BATCH_SIZE` when it is
/// set (a `TSeqPos`), else the `CBatchSizeMixer` maximum (`CHUNK_SIZE`, or the task's chunk size,
/// less 1000: 999,000 letters for blastn, 4,999,000 for megablast and dc-megablast; as the loop's
/// `TSeqPos`, so a maximum below 0 is a batch of about 4.29e9 letters); 0 when the variables cannot
/// be read (the search stops at them before its first batch).
fn x_query_batch_letters(task: &str) -> u64 {
    match query_batching_from_environment(
        crate::algorithm::blastn::coordination::task_uses_megablast_chunks(task),
    ) {
        Ok(batching) if batching.fixed_batch_size != 0 => u64::from(batching.batch_size),
        Ok(batching) => u64::from(batching.split.mixer_max_batch_size() as u32),
        Err(_) => 0,
    }
}

// No NCBI counterpart (scheduling only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QBATCHPAR): the letters of the searched subjects (`-subject`), the
/// database size in `x_qbatch_par::ncbi_splits_by_queries`; 0 for a database, where the slots
/// do not run.
fn x_subject_letters(source: SubjectSource<'_>) -> u64 {
    let Some(prepared) = source.prepared() else {
        return 0;
    };
    let letters = prepared.letters();
    (0..prepared.records.len())
        .map(|index| letters.len(index) as u64)
        .sum()
}

// No NCBI counterpart: LOSAT's verbose, debug and trace output (it changes no value NCBI computes).
/// EXPERIMENT (LOSAT_X_QBATCHPAR, LOSAT_X_RPTSTREAM): whether a search writes verbose, debug or
/// trace output while its batches are searched (`-verbose`, `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`,
/// `LOSAT_DEBUG_COORDS[_START]`, the blastn trace, the thread diagnostics). Both switches then keep
/// the loop as it was: the batches searched one after the other on the loop, and the reports
/// written after the loop, where that output comes before them.
pub(super) fn x_debug_output(args: &BlastnArgs) -> bool {
    args.verbose
        || [
            "BLEMIR_DEBUG",
            "LOSAT_DEBUG_BLASTN",
            "LOSAT_DEBUG_COORDS",
            "LOSAT_DEBUG_COORDS_START",
        ]
        .iter()
        .any(|name| std::env::var_os(name).is_some())
        || blastn_trace::enabled()
        || crate::utils::threading::diagnostics_enabled()
}

// No NCBI counterpart: what LOSAT_X_QBATCHPAR's slots need to know of a batch's result
// (scheduling only).
impl x_qbatch_par::SlotResult for Result<QueryBatch> {
    fn failed(&self) -> bool {
        matches!(self, Err(error) if !error.is::<x_qbatch_par::Cancelled>())
    }

    fn cancelled(&self) -> bool {
        matches!(self, Err(error) if error.is::<x_qbatch_par::Cancelled>())
    }

    /// The memory of the hit lists (lists, HSPs and edit scripts) and the per-query values.
    fn held_bytes(&self) -> u64 {
        let Ok(batch) = self else {
            return 0;
        };
        let hits: usize = batch
            .hit_lists
            .iter()
            .flatten()
            .map(|list| {
                std::mem::size_of_val(list.hsplist_array.as_slice())
                    + list
                        .hsplist_array
                        .iter()
                        .map(|hsp_list| {
                            std::mem::size_of_val(hsp_list.hsps.as_slice())
                                + hsp_list
                                    .hsps
                                    .iter()
                                    .map(|hsp| {
                                        hsp.gap_info
                                            .as_ref()
                                            .map_or(0, |ops| std::mem::size_of_val(ops.as_slice()))
                                    })
                                    .sum::<usize>()
                        })
                        .sum::<usize>()
            })
            .sum();
        let per_query = std::mem::size_of_val(batch.hit_lists.as_slice())
            + batch
                .query_masks
                .iter()
                .map(|masks| std::mem::size_of_val(masks.as_slice()))
                .sum::<usize>()
            + std::mem::size_of_val(batch.context_karlin.as_slice())
            + std::mem::size_of_val(batch.query_eff_searchsp.as_slice())
            + batch.timing.len();
        (hits + per_query) as u64
    }
}

/// NCBI's message for a query or subject without letters.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup.hpp:190-198
/// ```c
///     TSeqPos size() const {
///         TSeqPos retval = x_Size();
///         if (retval == 0) {
///             NCBI_THROW(CBlastException, eInvalidArgument,
///                        "Sequence contains no data");
///         }
///         return retval;
///     }
/// ```
pub(crate) const NO_DATA_MESSAGE: &str = "Sequence contains no data";

#[cfg(test)]
mod tests {
    use super::*;

    /// The NCBI errors of BLASTN's `-outfmt` at NCBI's positions (blast_args.cpp:2745-2748,
    /// 2874-2886, 2975-2977; blastn_app.cpp:255-258).
    #[test]
    fn output_formats_raise_ncbis_errors_where_ncbi_does() {
        use clap::Parser;
        #[derive(Parser)]
        struct Options {
            #[command(flatten)]
            blastn: BlastnArgs,
        }
        let args = |words: &[&str]| {
            let mut argv = vec!["test", "-query", "q", "-subject", "s"];
            argv.extend_from_slice(words);
            crate::cli::try_parse_from::<Options, _, _>(argv)
                .unwrap()
                .blastn
        };
        let native = |error: anyhow::Error| {
            let native = error.downcast::<crate::cli::NativeError>().unwrap();
            (native.exit, native.message)
        };
        let options = |text: &str| {
            (
                1,
                format!(
                    "BLAST query/options error: {text}\nPlease refer to the BLAST+ user manual.\n"
                ),
            )
        };
        // The CLI: a format that NCBI rejects later has no report (`cli_output_formats`).
        let file = Some(std::path::Path::new("o.json"));
        for outfmt in ["13", "14", "19", "21"] {
            assert!(
                cli_output_formats(outfmt, None).unwrap().is_empty(),
                "{outfmt}"
            );
        }
        assert!(cli_output_formats("13", Some(std::path::Path::new("-")))
            .unwrap()
            .is_empty());
        // 13 and 14 to a file are written (OF5c part 1): the XML2 / JSON documents as files.
        for (outfmt, kind) in [
            ("13", crate::report::blast_xml2::Xml2Kind::Files13),
            ("14", crate::report::blast_xml2::Xml2Kind::Files14),
        ] {
            assert_eq!(
                cli_output_formats(outfmt, file).unwrap(),
                [BlastnOutputFormat::Xml2(kind)],
                "{outfmt}"
            );
        }
        for outfmt in ["10"] {
            assert_eq!(
                cli_output_formats(outfmt, file).unwrap_err().to_string(),
                format!("output format {outfmt} is not supported by LOSAT's BLASTN")
            );
        }
        // SAM (OF6) reads its words and ignores the rest and a delimiter.
        assert_eq!(
            cli_output_formats("17", file).unwrap(),
            [BlastnOutputFormat::Sam(
                crate::report::sam::SamSpec::default()
            )]
        );
        assert_eq!(
            cli_output_formats(" 17 delim=, SR SQ foo", None).unwrap(),
            [BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
                subject_as_reference: true,
                sequence: true
            })]
        );
        assert_eq!(
            cli_output_formats("0 qaccver", None).unwrap(),
            [BlastnOutputFormat::Pairwise]
        );
        assert_eq!(
            cli_output_formats("7 qaccver", None).unwrap_err().to_string(),
            "a custom output format specification (fields or a delimiter) is not supported by LOSAT's BLASTN"
        );
        assert_eq!(
            native(cli_output_formats("22", None).unwrap_err()),
            (
                255,
                "Error: Formatting choice is out of range\n".to_string()
            )
        );
        // The formatting handler raises its errors before its warning; 13 passes it.
        for (outfmt, error) in [
            (
                "19",
                Some("AIRR rearrangement format is only applicable to igblastn"),
            ),
            (
                "21",
                Some("FASTA output format is only applicable to magicblast"),
            ),
            ("13", None),
        ] {
            let mut diagnostics = Vec::new();
            let mut sink = Vec::new();
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
            let result = process_options(&mut args(&["-max_target_seqs", "1"]), &mut outputs);
            drop(outputs);
            match error {
                Some(error) => {
                    assert_eq!(native(result.unwrap_err()), options(error), "{outfmt}");
                    assert!(diagnostics.is_empty(), "{outfmt}");
                }
                None => {
                    result.unwrap();
                    assert_eq!(
                        String::from_utf8(diagnostics).unwrap(),
                        "Warning: [blastn] Examining 5 or more matches is recommended\n"
                    );
                }
            }
        }
        // `run_local` has no output file; its formats that NCBI rejects later fail at once.
        for (outfmt, error) in [
            ("13", "Please provide a file name for outfmt 13."),
            ("14", "Please provide a file name for outfmt 14."),
            (
                "19",
                "AIRR rearrangement format is only applicable to igblastn",
            ),
            ("21", "FASTA output format is only applicable to magicblast"),
        ] {
            assert_eq!(
                native(parse_output_formats(["6", outfmt]).unwrap_err()),
                options(error),
                "{outfmt}"
            );
        }
        assert_eq!(
            parse_output_formats(["0", "6", "07"]).unwrap(),
            [
                BlastnOutputFormat::Pairwise,
                BlastnOutputFormat::Tabular,
                BlastnOutputFormat::TabularWithComments
            ]
        );
        assert_eq!(
            parse_output_formats(["17", "17 SQ SR"]).unwrap(),
            [
                BlastnOutputFormat::Sam(crate::report::sam::SamSpec::default()),
                BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
                    subject_as_reference: true,
                    sequence: true
                })
            ]
        );
    }
}

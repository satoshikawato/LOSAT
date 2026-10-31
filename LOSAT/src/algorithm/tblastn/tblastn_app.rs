//! TBLASTN's application (NCBI c++/src/app/blast/tblastn_app.cpp:114-361, `Run` and
//! `x_RunMTBySplitDB`): the command line's `TblastnArgs::run`, web ABI v2's `run_local`, and the
//! search over the query batches that both reach. The options are in `args.rs`.
use anyhow::{Context, Result};
use std::io::Write;

use super::args::{search_settings, ResolvedTblastnArgs, SearchSettings, TblastnArgs};
use super::database::TblastnDatabase;
use super::stage_d_pipeline::{
    run_local_for_report_threads, LocalStageDProfile, SubjectSetNumbers,
};
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{render, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use crate::algorithm::blastn::blast_engine::{QueryReading, NO_DATA_MESSAGE};
use crate::algorithm::common::database::{HeaderReadError, PreliminaryReadError};
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the failure of the batch whose search a read of the database
/// failed (`x_db_stream::XReadError`), as the check's error fails the first batch with a valid
/// context: a reader error of NCBI's own is the engine error; LOSAT's refusal (and any other
/// error) is returned as `Err` to stop the run.
fn x_read_failure(error: anyhow::Error) -> std::result::Result<anyhow::Error, anyhow::Error> {
    let error = match error.downcast::<super::x_db_stream::XReadError>() {
        Ok(super::x_db_stream::XReadError(read)) => read,
        Err(error) => return Err(error),
    };
    match error.downcast::<PreliminaryReadError>() {
        Ok(read_error) => Ok(crate::blastinput::app::engine_error(read_error.0)),
        Err(error) => Err(error),
    }
}

impl TblastnArgs {
    // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
    // results = lcl_blast.Run();
    // formatter.PrintOneResultSet(**result, query);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
    // PrintOneResultSet dispatches the complete result to outfmt 0/6/7.
    /// `settings` are the NCBI application settings that LOSAT reproduces
    /// (`ncbi_environment::check_ncbi_application_settings`): the readers of the subjects
    /// and of the queries (`fasta_reader`) use them.
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
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<()> {
        use crate::blastinput::app;
        use crate::blastinput::fasta_reader::{read_subjects, FastaInputSource, ReaderConfig};
        use crate::blastinput::input_files;
        use crate::blastinput::seq_range;
        let mut args = self;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
        // ```c
        //     if (GetExportSearchStrategyStream(args) ||
        //            m_FormattingArgs->ArchiveFormatRequested(args)) {
        //         locality = CBlastOptions::eBoth;
        //     }
        // ```
        // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
        // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
        // there too.
        let choice = app::parse_formatting_string(&args.outfmt)?;
        let format = app::report_format(&choice, app::Program::Tblastn, args.out.as_deref())?;
        // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot
        // write, before anything is read or written.
        // The error archives follow the `-outfmt` text (`writes_error_archive`).
        args.archive =
            super::archive::handle_for(&args, format == Some(app::ReportFormat::ArchiveFormat))?
                .map(|handle| handle.with_outfmt_text(&args.outfmt));
        let archive = args.archive.clone();
        // LOSAT's thread capability is checked before any input or output, as for BLASTN.
        crate::utils::threading::validate_threads(args.num_threads)?;
        // The diagnostics stream: with the archive it also saves what it writes as the
        // messages of the next archive.
        let mut stderr = std::io::stderr();
        let mut archive_tee;
        let diagnostics: &mut (dyn Write + Send) = match &archive {
            Some(handle) => {
                archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
                &mut archive_tee
            }
            None => &mut stderr,
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2562
        // ```c
        //         CRef<blast::CBlastQueryVector> subjects;
        //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
        //                                        subj_range, parse_deflines,
        //                                        use_lcase_masks, subjects, m_IsMapper);
        //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
        //
        //     } else if (!m_IsIgBlast){
        //         // IgBlast permits use of germline database
        //         NCBI_THROW(CInputException, eInvalidInput,
        //            "Either a BLAST database or subject sequence(s) must be specified");
        //     }
        // ```
        // The handler of the database arguments reads the subjects (an empty subject set
        // fails there), before the query and the output are opened.
        // `-db`: the database is opened with the options (`search_cli`), as NCBI's
        // `InitializeSubject` opens it after `SetOptions`.
        let subjects_read = (|| -> Result<(
            Vec<FastaRecord>,
            seq_range::Placements,
            Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
        )> {
            Ok(if args.db.is_some() {
            (Vec::new(), seq_range::Placements::default(), None)
        } else {
            let Some(subject_path) = args.subject.clone() else {
                return Err(app::missing_subject_error());
            };
            input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTN")?;
            let subject_file = input_files::open_input(&subject_path, "subject", "TBLASTN")?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
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
            let subject_range = seq_range::parse_optional_range(
                args.subject_loc.as_deref(),
                seq_range::RangeRole::Subject,
                "TBLASTN",
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-246
            // ```c++
            //     SDataLoaderConfig dlconfig(read_proteins);
            //     dlconfig.OptimizeForWholeLargeSequenceRetrieval();
            //
            //     CBlastInputSourceConfig iconfig(dlconfig);
            //     iconfig.SetRange(range);
            //     iconfig.SetBelieveDeflines(parse_deflines);
            //     iconfig.SetLowercaseMask(use_lcase_masking);
            //     iconfig.SetSubjectLocalIdMode();
            //     if (!read_proteins && gaps_to_Ns) {
            //         iconfig.SetConvertGapsToNs(true);
            //     }
            //
            //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
            //     CRef<CBlastInput> input(new CBlastInput(fasta));
            //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
            //     sequences = input->GetAllSeqs(*scope);
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2466-2470
            // ```c++
            //     EMoleculeType mol_type = Blast_SubjectIsNucleotide(opts.GetProgramType())
            //         ? CSearchDatabase::eBlastDbIsNucleotide
            //         : CSearchDatabase::eBlastDbIsProtein;
            //     m_IsProtein = (mol_type == CSearchDatabase::eBlastDbIsProtein);
            // ```
            // TBLASTN's subjects are nucleotide (`read_proteins` = `IsProtein()` is false, so the
            // reader assumes nucleotides, `fAssumeNuc`). NCBI reads the records one at a time, writes
            // each one's messages as it reads it, and checks each one's range after reading it: a
            // range that starts past the end of a record stops the reading there
            // (`read_subjects`). Records without residues (and intervals without letters) are
            // read without a message and reported when the formatter sets up the subjects
            // (`before_prolog`).
            let mut subject_source = FastaInputSource::from_argument(
                &subject_path,
                subject_file,
                ReaderConfig::subject("TBLASTN", false, settings.data_loaders)
                    .with_defline(archive.is_some())
                    .with_parse_deflines(args.parse_deflines),
            );
            // EXPERIMENT (LOSAT_X_TBNSTORE): each record's letters packed as it is read
            // (`x_read_subject_stores`), with a subject range after the record is cut to its
            // interval (the record `cut_subjects` makes).
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
            // ```c
            //     virtual SBlastSequence GetCompressedPlusStrand() {
            //         SBlastSequence retval(size());
            //         string ncbi4na = kEmptyStr;
            //         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
            //         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
            // ```
            // Dispatch point: NCBI keeps the subjects in its scope and makes each one's search forms
            // from it. With the switch the records keep their deflines and the letters are held packed,
            // made again when a subject is searched or reported; without it the records keep them.
            // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs);
            // with LOSAT_X_TBNSTORE the packing of each record is in this interval too. No NCBI
            // counterpart.
            let t_read_fasta = crate::utils::xtiming::clock();
            // With the switch the records come back cut to the range (`Some` placements).
            let (read_subject_records, x_subject_stores, x_placements) = {
                let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
                // Both sides of the stores are recorded in x_subject_store.rs. The archive of
                // `-subject` writes every subject's whole letters (and its deflines) in every
                // batch, which the packed stores do not keep (they hold the interval of a
                // `-subject_loc` and no letters), so a run that writes it reads the records
                // as the switch's `=0` side does.
                // Under `-parse_deflines` the reader compares each record's letters
                // (`believed::scope_clash`).
                if super::x_subject_store::x_tbn_store()
                    && archive.is_none()
                    && !args.parse_deflines
                {
                    let (records, stores, placements) =
                        super::x_subject_store::x_read_subject_stores(
                            &mut subject_source,
                            subject_range.as_ref(),
                            &mut warn,
                            args.num_threads,
                        )?;
                    (records, stores, Some(placements))
                } else {
                    (
                        read_subjects(&mut subject_source, subject_range.as_ref(), &mut warn)?,
                        None,
                        None,
                    )
                }
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
            drop(subject_source);
            let (subjects, subject_placements) = match x_placements {
                Some(placements) => (read_subject_records, placements),
                None => {
                    match seq_range::cut_subjects(&read_subject_records, subject_range.as_ref())? {
                        Some((cut, placements)) => (cut, placements),
                        None => (read_subject_records, seq_range::Placements::default()),
                    }
                }
            };
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
            // ```c
            // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
            //     : m_QueryVector(& queries)
            // {
            //     if (queries.Empty()) {
            //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
            //     }
            // ```
            // `m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects))` (blast_args.cpp:2557)
            // fails for a subject input without records.
            if subjects.is_empty() {
                return Err(app::empty_subjects_error());
            }
            (subjects, subject_placements, x_subject_stores)
            })
        })();
        let (subjects, subject_placements, x_subject_stores) = match subjects_read {
            Ok(read) => read,
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        };
        if archive.is_some() {
            super::archive::check_titles(&subjects, "subject")?;
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
        // The query is opened and the output file created before the options are checked;
        // `-` is standard input and standard output.
        input_files::check_utf8_file_name(&args.query, "query", "TBLASTN")?;
        // NCBI opens the query file inside `Run` too: its error is saved and an error archive
        // follows.
        let query_file = match input_files::open_input(&args.query, "query", "TBLASTN") {
            Ok(file) => file,
            Err(error) => {
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        };
        let out_file = open_out_file(&args)?;
        let pairwise = format == Some(app::ReportFormat::Pairwise);
        // The format as the reports read it: the number and the custom fields, or the value as
        // given when it has a delimiter, which `normalized` drops and the custom tabular report
        // reads (`tabular::tabular_spec`, which parses it again as NCBI does).
        let outfmt = if choice.delimiter.is_empty() {
            choice.normalized()
        } else {
            args.outfmt.clone()
        };
        let mut stream = cli_report_stream(out_file);
        let result = {
            let mut outputs =
                ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
            search_cli(
                args,
                query_file,
                &subjects,
                &subject_placements,
                x_subject_stores.as_deref(),
                &mut outputs,
                &choice,
                settings,
            )
        };
        // NCBI writes an error archive after the archives of the batches before when the
        // application stops with one of its own errors, and when the run ends with messages that
        // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
        if let Some(handle) = &archive {
            match &result {
                Err(error) => {
                    if let Some(lines) = super::archive::error_lines(error, false) {
                        let _ = crate::report::archive::write_error_archive(
                            handle,
                            &lines,
                            &mut stream,
                        );
                    }
                }
                Ok(()) => {
                    if handle.has_messages() {
                        let _ =
                            crate::report::archive::write_error_archive(handle, &[], &mut stream);
                    }
                }
            }
        }
        let flushed = std::io::Write::flush(&mut stream);
        // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
        // ```c
        //     catch (const std::ios::failure&) {                                      \
        //         LOG_POST(Error << "BLAST failed to write output");                  \
        //         exit_code = BLAST_OUTPUT_ERROR;                                     \
        //     }                                                                       \
        // ```
        // The outfmt 0 formatter's stream throws when a write fails; for outfmt 6/7 NCBI
        // aborts and LOSAT reports the error (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 3).
        if stream.failed && pairwise {
            return Err(crate::cli::NativeError {
                exit: 6,
                message: "BLAST failed to write output\n".to_string(),
            }
            .into());
        }
        result?;
        flushed.context("failed to write the output")
    }
}

/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &TblastnArgs,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    let Some(lines) = super::archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = cli_report_stream(out_file);
        let _ = crate::report::archive::write_error_archive(handle, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &TblastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if let Some(path) = args.out.as_deref() {
        crate::blastinput::input_files::check_utf8_file_name(path, "out", "TBLASTN")?;
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

/// The part of `run` after the files are opened: NCBI's processing of the options,
/// `Query is Empty!`, the query batch size and the subjects without letters, the reading of
/// the queries and the search.
#[allow(clippy::too_many_arguments)]
fn search_cli(
    mut args: TblastnArgs,
    query_file: std::fs::File,
    subjects: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    outputs: &mut ReportOutputs<'_>,
    choice: &crate::blastinput::app::FormatChoice,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    use crate::blastinput::fasta_reader::{
        read_queries, FastaInputSource, QueryEnd, QueryRecords, ReaderConfig,
    };
    let resolved = args.check_options(outputs.diagnostics)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
    // ```c++
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    // ```
    // The database of `-db` is opened after the options are processed, before
    // `Query is Empty!` (as BLASTN's).
    // The taxid or seqid list (line DI-8) is read first: NCBI's handler of the database
    // arguments reads it (`s_GetTaxIDList`, `CSeqDBFileGiList`) before the database opens.
    let database = match args.database_name() {
        Some(name) => {
            let list = args.db_filter.resolve(&settings, outputs.diagnostics)?;
            let database = super::database::TblastnDatabase::open_restricted(
                name,
                &settings,
                list.as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis);
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                &args.db_filter,
                list.as_ref(),
                database.filtered_empty(),
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "TBLASTN",
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    // ```c++
    // void CBlastTabularInfo::x_CheckTaxDB()
    // {
    //    	if( x_IsFieldRequested(eSubjectSciNames) ||
    //     ...
    //    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    // ```
    // A tabular field list that shows a name of the taxonomy database asks whether `taxdb.bti`
    // is on the search path: its formatters warn without it (the plural names), and LOSAT
    // refuses a name of a subject with a taxid with it (`tabular`), as BLASTN.
    if outputs
        .formats
        .iter()
        .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
        .any(|spec| spec.needs().taxonomy_names)
    {
        args.tabular_taxdb = Some(match database.as_ref() {
            Some(database) => database.taxdb_found(),
            None => crate::algorithm::common::database::subject_taxdb_found(&settings, "TBLASTN")?,
        });
    }
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let (Some(database), Some(serial)) = (
        database.as_ref(),
        crate::blastinput::app::ReportFormat::from_number(choice.number)
            .and_then(crate::report::serial::SerialFormat::of),
    ) {
        crate::report::serial::seqalign::check_database(
            serial,
            database.name(),
            &database.pairwise_database().definition,
            crate::blastinput::app::Program::Tblastn,
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:200-221
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    // ...
    //         if (pssm.Empty()) {
    //             CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                          query_opts->UseLowercaseMasks(),
    //                                          query_opts->GetParseDeflines(),
    //                                          query_opts->GetRange());
    //             if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //                	ERR_POST(Warning << "Query is Empty!");
    //                	return BLAST_EXIT_SUCCESS;
    //             }
    //             fasta.Reset(new CBlastFastaInputSource(
    //                                          m_CmdLineArgs->GetInputStream(),
    //                                          iconfig));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:51,106
    // ```c++
    //     const bool kQueryIsProtein = true;
    // ...
    //     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
    // ```c++
    //     flags += (m_ReadProteins
    //               ? CFastaReader::fAssumeProt
    //               : CFastaReader::fAssumeNuc);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query of tblastn is a protein (`QueryIsProtein()`, `fAssumeProt`) and its data
    // loaders are those of the subjects (`SDataLoaderConfig`). The position is taken on the
    // opened file, before it is read (`stream_is_empty`). When the subjects were read from
    // standard input too, `cin` has reached its end (a failed stream), so NCBI gets no
    // position for the query either, and its reader is at its end: no query batch is read.
    let shared_standard_input = args.query.as_os_str() == "-"
        && args
            .subject
            .as_deref()
            .is_some_and(|path| path.as_os_str() == "-");
    let mut query_source = FastaInputSource::from_argument(
        &args.query,
        query_file,
        ReaderConfig::query("TBLASTN", true, settings.data_loaders)
            .with_defline(args.archive.is_some())
            .with_parse_deflines(args.parse_deflines),
    );
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    let batch_size = before_prolog(subjects, x_subject_stores, outputs.diagnostics)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:267-270
    // ```c
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made.
    crate::blastinput::app::xinclude_check_for(choice, args.out.as_deref())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-277
    // ```c++
    //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //                 query =  input->GetNextSeqBatch(*scope);
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`read_query_batches`). What ends the reading (the end of the input, an `eEOF`, a
    // reader error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads
    // it in, after the reports of the batches before.
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = read_queries(&mut query_source);
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(
        &args,
        resolved,
        batch_size,
        &records,
        reading,
        subjects,
        subject_placements,
        x_subject_stores,
        database.as_ref(),
        outputs,
    )
}

/// NCBI's steps of tblastn between `Query is Empty!` and the outfmt 0 prolog that concern
/// the subjects: the query batch size, then the formatter, which sets up the subjects and
/// warns about each one without letters (a record without residues, or an interval that
/// starts just past its record's end). Such a subject stays in the database statistics
/// with no letters and is never searched (`x_scan_plan` skips subjects
/// shorter than a word). Returns the query batch size.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:217-221
/// ```c
///             fasta.Reset(new CBlastFastaInputSource(
///                                          m_CmdLineArgs->GetInputStream(),
///                                          iconfig));
///             input.Reset(new CBlastInput(&*fasta,
///                                         m_CmdLineArgs->GetQueryBatchSize()));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
/// ```c
///     if(m_IsDbScan) {
/// 	int num_seqs=0;
///         int total_length=0;
/// 	if (!is_remote_search)
///         {
///                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
///                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
///                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
/// ```c
///         catch(CBlastException & e ) {
///         	// Skip bad subject sequence
///         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
///         		seqblk_vec->push_back(subj);
///         ...
///         		warning += "Subject sequence contains no data";
///         		ERR_POST(Warning << warning);
///         		continue;
/// ```
/// NCBI's gapped tblastn never stops for the average subject length, also when every
/// subject is without letters: its shortest subject is raised to `BLAST_SEQSRC_MINLENGTH`
/// before the division by 3 (`run_local_search_with_pool` ports blast_setup.c:969-973 and
/// seqsrc_multiseq.cpp:226-240; LOSAT rejects `-ungapped` and `OLD_FSC`, the two ways to a
/// search without `sbp->gbp`).
fn before_prolog(
    subjects: &[FastaRecord],
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    diagnostics: &mut dyn Write,
) -> Result<u32> {
    let batch_size = crate::blastinput::app::query_batch_size("TBLASTN", 20000)?;
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject's length is its store's.
    crate::blastinput::input_files::write_empty_subject_warnings_by_length(
        subjects,
        |index| x_subject_stores.map_or(subjects[index].seq().len(), |stores| stores[index].len()),
        "tblastn",
        diagnostics,
    )?;
    Ok(batch_size)
}

// NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
// results = lcl_blast.Run();
// formatter.PrintOneResultSet(**result, query);
// NCBI c++/src/app/blast/blast_formatter.cpp:429-467:
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTN search over already parsed records and writes every requested
/// output format from the same result (the entry of web ABI v2; the command line reaches
/// `search` through `search_cli`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the
/// options and of the subjects come before the search. The `-query` and `-subject` values
/// of `args` are used only as display names.
///
/// The records are those of NCBI's reader (`fasta_reader`): the subjects as NCBI's
/// nucleotide reader reads them, the queries as its protein reader reads them. Queries
/// without a record are NCBI's `Query is Empty!`.
pub fn run_local(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // The BLAST archive (outfmt 11) among the formats: its request, the refusals of what LOSAT
    // cannot write, and a diagnostics stream that saves the messages. The records' deflines are
    // those the caller's reader kept; records without one get `>` and their title.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| super::stage_e_report::format_number(format.outfmt) == 11);
    if archive_requested {
        args.archive = super::archive::handle_for(&args, true)?;
        super::archive::check_titles(subject_records, "subject")?;
        if let Some(handle) = args.archive.clone() {
            return outputs.with_archive_diagnostics(handle, |view| {
                run_local_formats(args, query_records, subject_records, view)
            });
        }
    }
    run_local_formats(args, query_records, subject_records, outputs)
}

fn run_local_formats(
    args: TblastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, app::Program::Tblastn, None)?.is_none() {
            app::formatting_handler_check(&choice, app::Program::Tblastn)?;
            app::xinclude_check(&choice)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
        // ```c++
        //     s_ValidateCustomDelim(m_CustomOutputFormatSpec,m_CustomDelim);
        // ```
        // The formatting handler checks the custom delimiter of each format (the command line's
        // in `check_options`).
        app::formatting_handler_check(&choice, app::Program::Tblastn)?;
        outfmts.push(choice.number);
    }
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
        "TBLASTN",
    )?;
    for record in
        &subject_records[..seq_range::subjects_read(subject_records, subject_range.as_ref())]
    {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = seq_range::cut_subjects(subject_records, subject_range.as_ref())?;
    if subject_records.is_empty() {
        return Err(app::empty_subjects_error());
    }
    let resolved = args.check_options_for_formats(&outfmts, outputs.diagnostics)?;
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let whole_subjects = seq_range::Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    // The caller owns the records: their letters are read in place.
    let batch_size = before_prolog(searched_subjects, None, outputs.diagnostics)?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        QueryReading::Records,
        searched_subjects,
        subject_placements,
        None,
        None,
        outputs,
    )
}

/// The outfmt 0 prologs of a run that fails in its first query batch: NCBI writes the prolog
/// (`formatter.PrintProlog()`) before it reads the first batch. `subject_path` names a
/// FASTA subject set of `num_subjects` records (`subject_letters`); a database search shows
/// its `Database:` lines (`database`).
fn write_first_batch_prologs(
    outputs: &mut ReportOutputs<'_>,
    subject_path: Option<&std::path::Path>,
    num_subjects: usize,
    subject_letters: SubjectLetters<'_>,
    database: Option<&crate::report::pairwise::PairwiseDatabase>,
    taxdb_found: bool,
) -> Result<()> {
    let name = subject_path.map_or_else(String::new, |path| {
        format!("User specified sequence set (Input: {})", path.display())
    });
    let total = match database {
        Some(_) => 0,
        None => subject_letters.lengths().sum(),
    };
    for format in outputs.formats.iter_mut() {
        if crate::blastinput::app::parse_formatting_string(format.outfmt)?.number == 0 {
            let mut writer = format.sink.open()?;
            crate::report::pairwise::write_tblastn_pairwise_prolog(
                &mut writer,
                "2.17.0+",
                &name,
                num_subjects,
                total,
                database,
            )?;
            writer.flush()?;
        } else if let Some(kind) = super::stage_e_report::xml2_kind(format.outfmt) {
            // The XML2 and JSON documents start with their prolog too.
            let mut writer = format.sink.open()?;
            writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
            writer.flush()?;
        } else if let Some(spec) = super::tabular::tabular_spec(format.outfmt) {
            // Outfmt 20's line of the field names and its formatter's warning (nothing for the
            // other tabular kinds).
            let mut writer = format.sink.open()?;
            super::tabular::write_tabular_prolog(
                &spec,
                taxdb_found,
                &mut writer,
                &mut *outputs.diagnostics,
            )?;
            writer.flush()?;
        }
    }
    Ok(())
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options, the
/// reading of the subjects and the steps before the prolog (`before_prolog`): LOSAT's
/// limits and the environment, then the query batches (`read_query_batches`), the search of
/// the batches that NCBI searches and their reports.
#[allow(clippy::too_many_arguments)]
fn search(
    args: &TblastnArgs,
    resolved: ResolvedTblastnArgs,
    batch_size: u32,
    query_records: &[FastaRecord],
    reading: QueryReading,
    subject_records: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    // `-db`: the database (the subject records are then none).
    database: Option<&TblastnDatabase>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    // The records are NCBI's reader's records: the subjects of its nucleotide reader (a
    // subject without letters stays in the database, its warning written before the prolog,
    // `before_prolog`) and the queries of its protein reader. The query range applies to
    // every query record as NCBI's batch reader reads it (`seq_range::cut_queries`): a record
    // is searched cut, and one whose interval starts just past its end is a query without
    // data (`Sequence contains no data`, as a record without residues).
    let input_records = query_records;
    let (searched_queries, query_input) = match resolved.query_range.as_ref() {
        Some(range) => {
            let ranged = seq_range::cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            seq_range::QueryInput::whole(query_records),
        ),
    };
    let SearchSettings {
        composition_mode2,
        identity_cbs_reset,
        scoring,
    } = search_settings(&resolved)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_memento_priv.hpp:64-74
    // ```c++
    //     CBlastOptionsMemento(CBlastOptionsLocal* local_opts)
    //     {
    //         m_ProgramType = local_opts->GetProgramType();
    //         m_QueryOpts = local_opts->m_QueryOpts;
    // ...
    //         m_ExtnOpts = local_opts->m_ExtnOpts;
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288
    // ```c++
    //                 CLocalBlast lcl_blast(query_factory, opt_handle, db_adapter);
    // ```
    // IDENTITY's reset (`search_settings`) writes through the memento's pointers into the options
    // that every query batch reuses: the first batch resets the mode and warns its queries, and
    // the later batches start in mode 0 without the warning. These are the searched queries of
    // the first batch.
    let identity_warned = if identity_cbs_reset {
        query_input
            .batches(batch_size)
            .first()
            .map_or(0, |batch| batch.searched.end)
    } else {
        0
    };
    app::check_unsupported_environment("TBLASTN")?;
    app::check_old_fsc("TBLASTN")?;
    app::check_query_split_environment("TBLASTN", false)??;
    let subjects = subject_records;
    // The search reads the residues of the subject records in place (NCBI's subjects stay in the
    // reader's records, `CBlastQueryVector`); they are not copied for it. EXPERIMENT
    // (LOSAT_X_TBNSTORE): or from their packed form.
    let subject_seqs: Vec<&[u8]> = match x_subject_stores {
        Some(_) => Vec::new(),
        None => subjects.iter().map(FastaRecord::seq).collect(),
    };
    let subject_letters = match x_subject_stores {
        Some(stores) => SubjectLetters::Packed(stores),
        None => SubjectLetters::Records(&subject_seqs),
    };
    // What the requested formats need beyond outfmt 0, 6 and 7, decided once for the search
    // (`ReportNeeds`): with the Seq-aligns (outfmt 8, 12, 11), the engine positives of the
    // composition-adjusted matches (`LocalStageDProfile::engine_positives`).
    // The archive's Seq-aligns are those of outfmt 8.
    let engine_positives =
        crate::report::record::ReportNeeds::of(outputs.formats.iter().filter_map(|format| {
            app::ReportFormat::from_number(super::stage_e_report::format_number(format.outfmt))
        }))
        .seqalign
            || args.archive.is_some();
    // A database's numbers and `-dbsize` (`SubjectSetNumbers`), and the `Database:` lines of
    // a database search.
    let subject_set = SubjectSetNumbers {
        database: database.map(TblastnDatabase::numbers),
        dbsize: args.dbsize,
        searchsp: args.searchsp,
    };
    let pairwise_database = database.map(TblastnDatabase::pairwise_database);
    // BLAST XML (5, and the XML2 documents of 16 and 15) cannot carry a database name outside
    // US-ASCII (NCBI's stream throws): refused before the search, and before anything is
    // written.
    if let Some(database) = database {
        if outputs
            .formats
            .iter()
            .any(|format| super::stage_e_report::is_structured_xml(format.outfmt))
        {
            crate::report::blast_xml::names::check_database_name("TBLASTN", database.name())?;
        }
    }
    let subject_path = match database {
        Some(_) => None,
        None => Some(
            args.subject
                .as_ref()
                .context("TBLASTN subject file missing")?,
        ),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-278
    // ```c++
    //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //                 query =  input->GetNextSeqBatch(*scope);
    //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
    // ```
    // NCBI reads the query batches one at a time after the outfmt 0 prolog. The run stops at
    // the first batch that NCBI fails (a reader error, a batch without a query, a batch of
    // queries without letters), after the reports of the batches before (without the
    // epilog) and the messages of that batch's records (`read_query_batches`).
    let batches = crate::blastinput::query_batch::read_query_batches(
        input_records,
        &searched_queries,
        &query_input,
        reading,
        batch_size,
    );
    // The refusal of a title the archive cannot write comes before its first byte: the queries
    // of the batches NCBI completes are written (the records of the batch that a reader error
    // interrupts are not).
    if args.archive.is_some() {
        super::archive::check_titles(&searched_queries[..batches.searched_done], "query")?;
    }
    if batches.searched_done == 0 && batches.stop.is_some() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // ```c
        // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
        //     : m_QueryVector(& queries)
        // {
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // The first batch fails after the outfmt 0 prolog (`formatter.PrintProlog()`).
        write_first_batch_prologs(
            outputs,
            subject_path.map(|path| path.as_path()),
            subjects.len(),
            subject_letters,
            pairwise_database.as_ref(),
            args.tabular_taxdb.unwrap_or(false),
        )?;
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
        // ```c
        //     db_length = BlastSeqSrcGetTotLen(seq_src);
        //
        //     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
        //
        //     /* iterate over all subject sequences */
        //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
        //            != BLAST_SEQSRC_EOF) {
        //        ...
        //        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
        // ```
        // `-db`: the database's sequences, read once for every batch (`database::read_stores`)
        // where NCBI reads them first, in the preliminary search of the first batch with a
        // valid context (a batch without one is not searched, `CheckInternalData`). A database
        // without sequences (BLASTN's `Database memory map file error`,
        // blastn/blast_engine/run.rs) or a sequence file that cannot be read fails that batch,
        // after the reports of the batches before, which are searched without sequences; a run
        // without such a batch never reads the file. EXPERIMENT (LOSAT_X_TBNDBSTREAM): the
        // database is checked as `read_stores` reads it, with the same error, and each sequence
        // is read where it is searched (`x_db_stream.rs`).
        let mut x_stream = None;
        let (database_stores, mut stores_error) = match database {
            Some(database) => {
                let read = if database.printed_num_seqs() == 0 {
                    Err(PreliminaryReadError("Database memory map file error").into())
                } else if super::x_db_stream::x_tbn_db_stream() {
                    database.x_stream(resolved.num_threads).map(|stream| {
                        x_stream = Some(stream);
                        Vec::new()
                    })
                } else {
                    database.read_stores()
                };
                match read {
                    Ok(stores) => (stores, None),
                    Err(error) => (Vec::new(), Some(error)),
                }
            }
            None => (Vec::new(), None),
        };
        let search_letters = match (database, &x_stream) {
            (Some(_), Some(stream)) => SubjectLetters::Stream(stream),
            (Some(_), None) => SubjectLetters::Packed(&database_stores),
            (None, _) => subject_letters,
        };
        // The batches before the stop are searched and reported; a reader at the end of its
        // input reads no batch (`!input->End()`), and the reports have no query (an empty
        // pipe, or `-query -` after `-subject -`).
        let searched_done = batches.searched_done;
        let queries = &searched_queries[..searched_done];
        let query_offsets: Vec<usize> = (0..queries.len())
            .map(|index| query_input.placements.offset(index))
            .collect();
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let seg = resolved.seg.params();
        // NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
        // (*thread)->Run(); (*thread)->Join(&result);
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1417:
        // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
        //                         CConstRef<blast::CBlastQueryVector> queries,
        // NCBI c++/src/app/blast/tblastn_app.cpp:251-252,275-301:
        // input.Reset(new CBlastInput(&*fasta, GetQueryBatchSize()));
        // for (; !input->End(); ...) {
        //     query = input->GetNextSeqBatch(*scope);
        //     CLocalBlast lcl_blast(query_factory, ..., db_adapter);
        //     results = lcl_blast.Run();
        // }
        // NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:70-127:
        // case eTblastn: retval = 20000;
        // Search each accumulated batch independently; preliminary linking
        // and sum-statistics depend on its concatenated query contexts.
        let mut results = Vec::with_capacity(queries.len());
        let mut lengths: Option<LocalSubjectParameters> = None;
        let mut ungapped_karlin = Vec::with_capacity(queries.len());
        let mut query_validity = Vec::with_capacity(queries.len());
        let mut query_batch_skipped = Vec::with_capacity(queries.len());
        // `-db`: the stored titles read so far, and the batch that a header file or the
        // sequence file failed with its error.
        let mut titles = super::database::DatabaseTitles::new();
        let mut failed: Option<(std::ops::Range<usize>, anyhow::Error)> = None;
        // Graduated 2026-10-22 (was LOSAT_X_TBNSSIDE, LOSAT_X_TBNBATCH): the batches of the
        // loop below searched ahead (`x_search_batches`), each with its range; `None` when there
        // is none. A batch whose range does not match is searched by the loop.
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-291
        // ```c
        //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
        // ...
        //                 query =  input->GetNextSeqBatch(*scope);
        //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
        // ...
        //                     CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
        //                     lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
        //                     results = lcl_blast.Run();
        // ```
        // Dispatch point: NCBI searches the query batches one after another (one CLocalBlast per
        // batch). x_search_batches searches them ahead (on one shared pool, or side by side) and
        // returns the reports in batch order; the loop below takes them from it.
        // A database whose sequences could not be read searches its batches one by one below.
        let x_batches = match stores_error {
            Some(_) => Vec::new(),
            None => query_input
                .batches(batch_size)
                .into_iter()
                .map(|batch| batch.searched)
                .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
                .collect(),
        };
        // EXPERIMENT (LOSAT_X_TBNDBSTREAM): a read of the database that failed while the batches
        // were searched ahead: they are searched again one at a time below, so that the batch
        // that reads it fails after the reports of the batches before.
        let x_searched_ahead = super::stage_d_pipeline::x_search_batches(
            &query_seqs,
            search_letters,
            x_batches,
            |range| LocalStageDProfile {
                seg: seg.as_ref(),
                soft_masking: resolved.soft_masking,
                mask_lowercase: resolved.lcase_masking,
                genetic_code: resolved.db_gencode,
                expect_value: resolved.evalue,
                max_target_seqs: resolved.hitlist_size,
                query_offsets: &query_offsets[range],
                subject_set,
                engine_positives,
                culling_limit: resolved.culling_limit,
            },
            composition_mode2,
            resolved.sum_stats,
            scoring,
            resolved.num_threads,
        );
        let mut x_searched = match x_searched_ahead {
            Ok(searched) => searched.map(Vec::into_iter),
            Err(error) if error.is::<super::x_db_stream::XReadError>() => None,
            Err(error) => return Err(error),
        };
        for range in query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
        {
            let (mut batch_results, batch_lengths, batch_karlin, batch_validity) = {
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-291
                // ```c
                //                 query =  input->GetNextSeqBatch(*scope);
                //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
                // ...
                //                     CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
                //                     lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                //                     results = lcl_blast.Run();
                // ```
                // Dispatch point: a report searched ahead for this very range is used; otherwise
                // run_local_for_report_threads searches the batch here, which is the port of
                // CLocalBlast::Run for one batch.
                let x_searched_here = match x_searched.as_mut().and_then(Iterator::next) {
                    Some((x_range, report)) if x_range == range => Ok(report),
                    _ => run_local_for_report_threads(
                        &query_seqs[range.clone()],
                        search_letters,
                        LocalStageDProfile {
                            seg: seg.as_ref(),
                            soft_masking: resolved.soft_masking,
                            mask_lowercase: resolved.lcase_masking,
                            genetic_code: resolved.db_gencode,
                            expect_value: resolved.evalue,
                            max_target_seqs: resolved.hitlist_size,
                            query_offsets: &query_offsets[range.clone()],
                            subject_set,
                            engine_positives,
                            culling_limit: resolved.culling_limit,
                        },
                        composition_mode2,
                        resolved.sum_stats,
                        scoring,
                        resolved.num_threads,
                    ),
                };
                // EXPERIMENT (LOSAT_X_TBNDBSTREAM): a read of the database that failed in this
                // batch fails it as the check's error would (below).
                match x_searched_here {
                    Ok(report) => report,
                    Err(error) => match x_read_failure(error) {
                        Ok(failure) => {
                            failed = Some((range, failure));
                            break;
                        }
                        Err(error) => return Err(error),
                    },
                }
            };
            // `-db`: the first batch with a valid context reads the sequences (above); a
            // reader error of NCBI's own fails it, LOSAT's refusal stops the run.
            if batch_validity.iter().any(|&valid| valid) {
                if let Some(error) = stores_error.take() {
                    match error.downcast::<PreliminaryReadError>() {
                        Ok(read_error) => {
                            failed = Some((range, app::engine_error(read_error.0)));
                            break;
                        }
                        Err(error) => return Err(error),
                    }
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
            // ```c++
            //     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
            //                                                   m_InternalData,
            //                                                   m_Opts,
            //                                                   seqinfo_src,
            //                                                   search_msgs));
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
            // ```c++
            //     if(hit_params->options->query_cov_hsp_perc > 0 || hit_params->options->max_hsps_per_subject > 0 ||
            //        (hit_params->options->hsp_filt_opt != NULL && hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
            //     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
            //     }
            //
            //     /* Re-sort the hit lists according to their best e-values, because they
            //        could have changed. Only do this for a database search. */
            //     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
            //         Blast_HSPResultsSortByEvalue(results);
            //     }
            // ```
            // `-max_hsps` and `-qcov_hsp_perc` filter the batch's final hit lists (each query's
            // context is its searched letters), which the search has sorted by e-value; the
            // sort is a total order (best e-value, first score, OID), so sorting again after
            // the filter gives NCBI's order. A search with hits has subject letters.
            // `-culling_limit` runs NCBI's culling pipe first on each query's hit list
            // (`BlastHSPStreamTBackClose`, before the filters of `s_FilterBlastResults`;
            // hspfilter_culling.c:699-733); its lists leave in the order of the tree, and the
            // sort by e-value below is the one NCBI runs unconditionally after the pipes.
            if resolved.max_hsps > 0
                || resolved.query_cov_hsp_perc > 0.0
                || resolved.culling_limit > 0
            {
                for (hitlist, query) in batch_results.iter_mut().zip(&query_seqs[range.clone()]) {
                    if resolved.culling_limit > 0 {
                        hitlist.cull(resolved.culling_limit, query.len())?;
                    }
                    if resolved.max_hsps > 0 || resolved.query_cov_hsp_perc > 0.0 {
                        hitlist.filter_results(
                            resolved.max_hsps,
                            resolved.query_cov_hsp_perc,
                            query.len(),
                        );
                    }
                    hitlist.sort_by_evalue_for_positive_totlen();
                }
            }
            // `-db`: the definition lines of the batch's final hit lists are read where NCBI
            // makes the batch's Seq-aligns (`database::read_titles`); a header file that NCBI
            // cannot read fails this batch, after the reports of the batches before.
            if let Some(database) = database {
                if let Err(error) =
                    super::database::read_titles(database, &batch_results, range.start, &mut titles)
                {
                    match error.downcast::<HeaderReadError>() {
                        Ok(failure) => {
                            failed = Some((range, failure.error));
                            break;
                        }
                        Err(error) => return Err(error),
                    }
                }
            }
            results.append(&mut batch_results);
            if let Some(all_lengths) = &mut lengths {
                // NCBI tblastn_app.cpp:288-301 formats each batch's own
                // query-context lengths; retain those values in input order.
                all_lengths.lengths.extend(batch_lengths.lengths);
                all_lengths.cutoffs.extend(batch_lengths.cutoffs);
            } else {
                lengths = Some(batch_lengths);
            }
            ungapped_karlin.extend(batch_karlin);
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
            // if (CheckInternalData() != 0) each result in this Run() batch
            // receives an unsearched ancillary state and null align set.
            let skipped = batch_validity.iter().all(|&valid| !valid);
            query_batch_skipped.extend(std::iter::repeat(skipped).take(batch_validity.len()));
            query_validity.extend(batch_validity);
        }
        // A run without a searched query has no parameters of a query set.
        let mut lengths = lengths.unwrap_or(LocalSubjectParameters {
            lengths: Vec::new(),
            cutoffs: Vec::new(),
            link: None,
            hit_cutoff_min: 0,
            word_cutoff_min: 0,
            prelim_evalue: 0.0,
        });
        // A batch that a database file failed is not reported, nor are the batches after it.
        let queries = match &failed {
            Some((range, _)) => {
                results.truncate(range.start);
                lengths.lengths.truncate(range.start);
                lengths.cutoffs.truncate(range.start);
                ungapped_karlin.truncate(range.start);
                query_validity.truncate(range.start);
                query_batch_skipped.truncate(range.start);
                &queries[..range.start]
            }
            None => queries,
        };
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1443-1451:
        // if (results.HasWarnings()) ERR_POST(Warning << results.GetWarningStrings());
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
        // ```c
        //         } catch (const CException& e) {
        //             ...
        //             CRef<CSearchMessage> m
        //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
        //             messages[index].push_back(m);
        //             s_InvalidateQueryContexts(qinfo, index);
        //         }
        // ```
        // The warnings of each query are written before its report, once (`QueryWarnings`):
        // the reader's messages of the query's batch (before its first query,
        // `read_query_batches`), then the query's own, on one line in NCBI's order
        // (`query_warning`): the O characters read as X, and the message of its set-up (a
        // query without letters) or of its Karlin-Altschul parameters (an invalid query
        // with letters).
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
        // ```c++
        //     if (m_FormatType == CFormattingArgs::eAsnText
        //       ...
        //       || m_FormatType == CFormattingArgs::eJsonSeqalign
        //       || m_FormatType == CFormattingArgs::eSAM)
        //     {
        //         x_PrintStructuredReport(results, queries);
        //         return;
        //     }
        //     ...
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        // A structured format (outfmt 8, 12) writes none of a query's own warnings
        // (`ReportFormat::writes_query_warnings`); the warnings go with the first format.
        let result_warnings = outputs.formats.first().is_none_or(|format| {
            crate::blastinput::app::ReportFormat::from_number(super::stage_e_report::format_number(
                format.outfmt,
            ))
            .is_none_or(crate::blastinput::app::ReportFormat::writes_query_warnings)
        });
        let filtered_empty = database.is_some_and(|database| database.filtered_empty());
        let mut query_warning_lines: Vec<Vec<u8>> = queries
            .iter()
            .zip(&query_validity)
            .enumerate()
            .map(|(index, (query, valid))| {
                let mut lines = batches.warnings[index].clone();
                if !result_warnings {
                    return lines;
                }
                let mut messages = Vec::new();
                messages.extend(crate::report::query_warnings::replaced_o_message(
                    query.seq(),
                ));
                if index < identity_warned {
                    messages.push(
                        crate::report::query_warnings::IDENTITY_CBS_RESET_MESSAGE.to_string(),
                    );
                }
                if query.seq().is_empty() {
                    messages.push(NO_DATA_MESSAGE.to_string());
                } else if !valid {
                    messages.push(crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string());
                }
                if filtered_empty && !query_batch_skipped[index] {
                    messages
                        .push(crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string());
                }
                lines.extend(crate::report::query_warnings::query_warning(
                    "tblastn",
                    query_input.ordinal(index),
                    query,
                    &messages,
                ));
                lines
            })
            .collect();
        // `-db`: the reported sequences by rank, with their records and packed subjects; the
        // database's other sequences are dropped before the reports.
        let reported = match database {
            Some(database) => Some(
                super::database::reported_subjects(
                    database,
                    &titles,
                    &mut results,
                    search_letters,
                    parsed_shown(&outputs.formats),
                )
                // EXPERIMENT (LOSAT_X_TBNDBSTREAM): a read of a reported sequence that
                // failed: the run's error, as the check's error would be.
                .map_err(|error| x_read_failure(error).unwrap_or_else(|error| error))?,
            ),
            None => None,
        };
        drop(database_stores);
        super::x_db_stream::x_print_shadow_summary();
        // The query batches that were read and searched, by their searched queries (the batches
        // before a database failure or a reader stop): each is reported.
        let batch_spans: Vec<std::ops::Range<usize>> = query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= queries.len())
            .collect();
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:298-302
        // ```c++
        //                     BlastFormatter_PreFetchSequenceData(*results, scope,
        //                     		                            fmt_args->GetFormattedOutputChoice());
        //                     ITERATE(CSearchResultSet, result, *results) {
        //                         formatter.PrintOneResultSet(**result, query);
        //                     }
        // ```
        // Before the reports of a batch NCBI fetches every subject that the batch reports, and
        // with `taxdb.bti` on the search path the taxonomy database writes its messages there
        // (`TblastnDatabase::prefetch_messages`: the texts of opening it at the run's first
        // lookup, a volume's lookup of a taxid once per run): after the reader's messages of the
        // batch and before the warnings of its first query. They go with the warnings, where the
        // first format is outfmt 0, 5, 6 or 7 (`database::prefetches_subjects`); XML2 and JSON
        // write theirs at each `Hit` (`TaxdbPrinter`). They are read here, after the refusals of
        // the reported subjects and before any format is written, so a refusal (a taxonomy
        // database or a taxid that LOSAT cannot read as NCBI does) comes before any byte, as
        // those refusals do; without `taxdb.bti` nothing is read. A failed batch is not reported
        // and prefetches nothing (`failed`, a header file or the sequence file that NCBI cannot
        // read during the batch's search).
        if let (Some(reported), Some(database)) = (&reported, database) {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:716-725
            // ```c++
            // s_IsPrefetchFormat(blast::CFormattingArgs::EOutputFormat format_type)
            // {
            // 	if ((format_type == CFormattingArgs::eAsnText) ||
            // 	    (format_type == CFormattingArgs::eAsnBinary) ||
            // 	    (format_type == CFormattingArgs::eArchiveFormat)||
            // 	    (format_type == CFormattingArgs::eJsonSeqalign)) {
            // 		return false;
            // 	}
            // 	return true;
            // ```
            // Outfmt 10 and 20 prefetch as the tabular 6 and 7 do (a name of the taxonomy
            // database that the prefetch would look up is LOSAT's refusal, `tabular`).
            if outputs.formats.first().is_some_and(|format| {
                crate::algorithm::common::database::prefetches_subjects(Some(
                    match super::stage_e_report::format_number(format.outfmt) {
                        10 | 20 => 6,
                        number => number,
                    },
                ))
            }) && database.taxdb_found()
            {
                // A custom tabular report with a plural name looks every taxid of its rows up,
                // taxid 0 too, which the prefetch does not: the rows of a batch with hits may open
                // the taxonomy database before a later batch's prefetch
                // (`TaxdbRun::open_for_rows`).
                let plural_names = outputs
                    .formats
                    .first()
                    .and_then(|format| super::tabular::tabular_spec(format.outfmt))
                    .is_some_and(|spec| spec.checks_taxdb());
                for span in &batch_spans {
                    let hitlists = &results[span.clone()];
                    let messages = database.prefetch_messages(reported, hitlists)?;
                    if plural_names
                        && hitlists
                            .iter()
                            .flat_map(super::stage_d_results::KappaResultHitList::lists)
                            .any(|list| !list.hsps.hsps.is_empty())
                    {
                        database.open_for_tabular_rows()?;
                    }
                    if messages.is_empty() {
                        continue;
                    }
                    let reader = batches.warnings[span.start].len();
                    query_warning_lines[span.start].splice(reader..reader, messages);
                }
            }
        }
        // `-db`, XML2 and JSON: the taxid of each reported sequence, and what the taxonomy
        // database says about it (the refusals of both come before any byte).
        let xml_taxa = match (&reported, database) {
            (Some(reported), Some(database))
                if outputs
                    .formats
                    .iter()
                    .any(|format| super::stage_e_report::xml2_kind(format.outfmt).is_some()) =>
            {
                database.read_taxids(reported, &results)?
            }
            _ => Default::default(),
        };
        // The custom tabular reports (`tabular`): whether `taxdb.bti` was found, the taxids of the
        // reported database sequences where a field shows them and the taxonomy database's names
        // that the rows look up (read here, before any byte, so that their refusals come before
        // it), and the `-query_loc` range's length (`qlen`, `qcovs`).
        let tabular_specs: Vec<_> = outputs
            .formats
            .iter()
            .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
            .collect();
        let taxonomy = match (&reported, database) {
            (Some(reported), Some(database))
                if tabular_specs.iter().any(|spec| spec.needs().taxonomy) =>
            {
                database.tabular_taxonomy(reported, &results)?
            }
            _ => Vec::new(),
        };
        let names = match (tabular_specs.first(), database) {
            (Some(spec), Some(database)) => database.tabular_names(spec, &taxonomy)?,
            _ => None,
        };
        let tabular_inputs = super::tabular::TabularInputs {
            taxdb_found: args.tabular_taxdb.unwrap_or(false),
            taxonomy,
            names,
            range_length: resolved
                .query_range
                .as_ref()
                .map(crate::blastinput::seq_range::SequenceRange::length),
        };
        // The records of a subject set whose ids `-parse_deflines` read, by subject index.
        let records_parsed = if args.parse_deflines && reported.is_none() {
            crate::report::db_subject::parsed_records(
                subjects
                    .iter()
                    .map(|record| (record.believed(), record.title.as_slice())),
                args.show_gis,
                false,
                crate::blastinput::app::parse_formatting_string(&args.outfmt)?.number == 0,
                "TBLASTN",
            )?
        } else {
            Vec::new()
        };
        let (report_records, report_letters, shown) = match (&reported, database) {
            (Some(reported), Some(database)) => (
                reported.records.as_slice(),
                SubjectLetters::Packed(&reported.stores),
                ShownSubjects::Database {
                    name: database.name(),
                    pairwise: pairwise_database
                        .as_ref()
                        .expect("a database search has its Database: lines"),
                    reported,
                },
            ),
            _ => (
                subjects,
                subject_letters,
                ShownSubjects::Fasta {
                    path: subject_path.expect("a search of -subject has its path"),
                    parsed: &records_parsed,
                },
            ),
        };
        // The BLAST archive (outfmt 11): one document per query batch that was read and searched
        // (`batch_spans`).
        let archive_render = args
            .archive
            .as_ref()
            .map(|handle| super::archive::ArchiveRender {
                handle,
                spans: &batch_spans,
            });
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); it is
        // written after all batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        render(
            outputs,
            queries,
            report_records,
            report_letters,
            &shown,
            subject_set,
            &mut results,
            &lengths,
            &ungapped_karlin,
            &query_validity,
            &query_batch_skipped,
            scoring,
            &resolved.matrix_name,
            crate::blastinput::report_counts::ReportCounts {
                hitlist_size: resolved.hitlist_size,
                num_descriptions: resolved.num_descriptions,
                num_alignments: resolved.num_alignments,
            },
            resolved.db_gencode,
            seg.as_ref(),
            resolved.lcase_masking,
            &query_warning_lines,
            &super::stage_e_report::ReportRanges {
                queries: &query_input.placements,
                subjects: subject_placements,
                // NCBI reference: c++/src/app/blast/tblastn_app.cpp:342
                // ```c++
                //         formatter.PrintEpilog(opt);
                // ```
                // A run that stops at a query batch writes no epilog.
                epilog: batches.stop.is_none() && failed.is_none(),
            },
            super::blast_xml::XmlValues {
                expect: resolved.evalue,
                soft_masking: resolved.soft_masking,
                gapped: resolved.gapped,
                // IDENTITY's reset is NCBI's option of the search (`search_settings`).
                cbs: if identity_cbs_reset {
                    0
                } else {
                    resolved.composition_mode as i32
                },
                db_gencode: i32::from(resolved.db_gencode),
                identity_warned,
                filtered_empty,
            },
            xml_taxa,
            args.out.as_deref(),
            archive_render.as_ref(),
            &tabular_inputs,
        )?;
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
        // The failed batch's reading messages, then its error.
        if let Some((range, error)) = failed {
            for warnings in &batches.warnings[range] {
                outputs.diagnostics.write_all(warnings)?;
            }
            return Err(error);
        }
    }
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        return Err(error);
    }
    Ok(())
}

/// What the reports `formats` show of the reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`common::database::ParsedShown`): the tabular formats (6, 7, 10,
/// 20), outfmt 0, and the formats that join the alignments of following hits whose Seq-align
/// ids match (0, 5, 13 to 16).
fn parsed_shown(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> crate::algorithm::common::database::ParsedShown {
    let numbers: Vec<i32> = formats
        .iter()
        .map(|format| super::stage_e_report::format_number(&format.outfmt))
        .collect();
    let pairwise = numbers.contains(&0);
    crate::algorithm::common::database::ParsedShown {
        tabular: numbers
            .iter()
            .any(|&number| matches!(number, 6 | 7 | 10 | 20)),
        pairwise,
        joined: pairwise
            || formats
                .iter()
                .any(|format| super::stage_e_report::is_structured_xml(&format.outfmt)),
        sam: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // EXPERIMENT (LOSAT_X_TBNDBSTREAM): a read that failed during the search fails its batch with
    // NCBI's engine error (exit status 3) when the reader's error is NCBI's own, and stops the run
    // with LOSAT's refusal or any other error, as the check's error does.
    #[test]
    fn a_read_failed_in_the_search_fails_as_the_check_s_error() {
        use super::super::x_db_stream::XReadError;
        let engine = x_read_failure(
            XReadError(PreliminaryReadError("Database memory map file error").into()).into(),
        )
        .expect("the engine error");
        let native = engine
            .downcast::<crate::cli::NativeError>()
            .expect("a NativeError");
        assert_eq!(native.exit, 3);
        assert_eq!(
            native.message,
            "BLAST engine error: Database memory map file error\n"
        );
        let refusal = x_read_failure(XReadError(anyhow::anyhow!("refused")).into()).unwrap_err();
        assert_eq!(refusal.to_string(), "refused");
        let other = x_read_failure(anyhow::anyhow!("other")).unwrap_err();
        assert_eq!(other.to_string(), "other");
    }
}

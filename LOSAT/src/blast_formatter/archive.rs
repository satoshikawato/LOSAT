//! blast_formatter's run over the objects of a BLAST archive (`-archive`): NCBI's
//! `CBlastFormatterApp::Run` with `CRemoteBlast::LoadFromArchive` and `PrintFormattedOutput`
//! (c++/src/app/blast/blast_formatter.cpp:329-531) over the archive reader
//! (`report::archive::read`), and the formats 8, 9 and 12, which write each query's Seq-aligns as
//! the archive holds them, cut to the hit list size, inside the one Seq-annot wrapper of the
//! search programs (`report::serial::seqalign::write_document_to`, decision F3-2 of line RG).
//!
//! The object loop keeps what NCBI's `CRemoteBlast` keeps from one object to the next: the
//! options of the first search object (`GetSearchOptions` builds them once: the program, the
//! taxids, the database's filtering algorithm), the hit list size that `ExtractAlgorithmOptions`
//! set last, the database that the last database object named (`m_Dbs`, which an object of
//! `-subject` sequences does not clear) and the `-out` stream, opened when the first object is
//! formatted. Where NCBI's result would depend on its servers (a database that cannot be opened
//! here), on memory its input does not decide, or on a form LOSAT does not read, the run stops
//! with LOSAT's refusal (`refusals::Refusal`) at the point where NCBI meets it.

use std::io::{BufWriter, Read, Write};
use std::path::Path;

use super::formatting::{extract_formatting, sam_check};
use super::refusals::{self, Refusal};
use super::search_program::{search_program, SearchProgram};
use super::{Options, PROGRAM};
use crate::algorithm::common::database::{
    database_error, database_error_of, open_search_database, pairwise_database, refusal,
    DatabaseOpenError,
};
use crate::blastinput::app::{
    check_search_environment, engine_error, options_error_bytes, output_format,
    prefetch_seqs_limit, xinclude_check_for,
};
use crate::blastinput::db_mask::name_mask_error;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::{command_line_argument_error, NativeError};
use crate::ncbi_app::files::{self, ends_like_sigpipe, OutputFile};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::NcbiException;
use crate::report::archive::read::options::{archive_request, BuiltOptions};
use crate::report::archive::read::request::{Queries, Subject};
use crate::report::archive::read::{
    seq_align_id, Archive, ArchiveReader, Message, Search, Severity, Stop,
};
use crate::report::serial::asn_types::stream::{guess_encoding, Encoding, GuessStop};
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::write::write_value;
use crate::report::serial::asn_types::ReadError;
use crate::report::serial::encoder::Encoder;
use crate::report::serial::seqalign::{
    write_document_to, AnnotDesc, DatabaseDesc, SeqAnnotDocument, SerialFormat,
};
use crate::seqdb::{DbList, ListWarning, SeqDb, SeqDbError, SeqType, TaxidList};
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::ParsedSeqId;

/// The size of the buffer in front of the output (the reports are written in large blocks; the
/// buffer is flushed before each diagnostic line where the output is the standard output).
const OUTPUT_BUFFER: usize = 1 << 16;

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:506-529
// ```c++
//         if (args[kArgArchive].HasValue()) {
//             CNcbiIstream& istr = args[kArgArchive].AsInputFile();
//             try { m_RmtBlast.Reset(new CRemoteBlast(istr)); }
//             catch (const CBlastException& e) {
//                 if (e.GetErrCode() == CBlastException::eInvalidArgument) {
//                     NCBI_RETHROW(e, CInputException, eInvalidInput,
//                                  "Invalid input format for BLAST Archive.");
//                 }
//             }
//
// 	    m_LoadFromArchive = true;
//             try {
//                 while (m_RmtBlast->LoadFromArchive()) {
//                 	if(!m_RmtBlast->IsErrMsgArchive()) {
//                 		status = PrintFormattedOutput();
//                 	}
//                 }
//             } catch (const CSerialException& e) {
//                 NCBI_RETHROW(e, CInputException, eInvalidInput,
//                              "Invalid input format for BLAST Archive.");
//             }
//             x_AddCmdOptions();
//             m_UsageReport.AddParam(CBlastUsageReport::eExitStatus, status);
//     	    return status;
// ```
/// The formatting of the archive `archive` with the command line's `options` and the
/// application's `settings` (the registry's `[BLAST] BLASTDB`), its reports written to `-out`
/// (`out` is the standard output), its diagnostics on `diag`; the exit status. The errors are
/// `CATCH_ALL`'s (`super::catch_all`), after what was written before them.
pub(crate) fn format_archive(
    options: &Options,
    settings: &ApplicationSettings,
    archive: &mut dyn Read,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let mut run = Run {
        options,
        settings,
        output: Output {
            stdout: Some(out),
            file: None,
        },
        search_options: None,
        hitlist: 0,
        database: None,
        opened: None,
        buf: Vec::new(),
    };
    let outcome = run.objects(archive, diag);
    // `cerr` is tied to `cout`: the output written before the error comes first.
    let flushed = run.output.flush();
    let outcome = match outcome {
        Ok(()) => flushed.map_err(anyhow::Error::from),
        Err(error) => Err(error),
    };
    match outcome {
        Ok(()) => 0,
        Err(error) => super::catch_all(diag, write_failure(error)),
    }
}

/// What the loop keeps from one archive object to the next (`CRemoteBlast`'s members, and the
/// `-out` stream of `CArgs`).
struct Run<'r, 'o> {
    options: &'r Options,
    settings: &'r ApplicationSettings,
    output: Output<'o>,
    /// The options of the first search object (`m_CBOH`, `m_TaxidList`,
    /// `m_DbFilteringAlgorithmId`, `m_Task`), which every later object is formatted with.
    search_options: Option<BuiltOptions>,
    /// The options' hit list size: the first object's, then what `ExtractAlgorithmOptions` set
    /// for the last formatted object (`opts.SetHitlistSize` on the kept options).
    hitlist: i32,
    /// `m_Dbs`: the name and residue type of the database of the last database object.
    database: Option<(Vec<u8>, bool)>,
    /// The database opened for an earlier object (`opened_database`).
    opened: Option<OpenedDatabase>,
    /// The buffer of one Seq-annot.
    buf: Vec<u8>,
}

/// The `-out` stream: the standard output until `-out` is opened (`AsOutputFile`, at the first
/// formatted object), then the file or the standard output behind a buffer.
struct Output<'o> {
    stdout: Option<&'o mut dyn Write>,
    file: Option<BufWriter<OutputFile<'o>>>,
}

impl<'o> Output<'o> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:335
    // ```c++
    //     CNcbiOstream& out = args[kArgOutput].AsOutputFile();
    // ```
    /// `args[kArgOutput].AsOutputFile()` of `-out` `value`: opened the first time (NCBI's
    /// `Argument "out". File is not accessible` there), the same stream after.
    fn open(&mut self, value: &str) -> anyhow::Result<()> {
        if self.file.is_some() {
            return Ok(());
        }
        let stdout = self
            .stdout
            .take()
            .expect("the standard output until -out is opened");
        let file = files::open_output_file("out", value, stdout)
            .map_err(|error| command_line_argument_error(error.message()))?;
        self.file = Some(BufWriter::with_capacity(OUTPUT_BUFFER, file));
        Ok(())
    }

    /// The opened stream.
    fn writer(&mut self) -> &mut BufWriter<OutputFile<'o>> {
        self.file
            .as_mut()
            .expect("-out opened at the object's start")
    }

    /// Writes out what the buffer holds.
    fn flush(&mut self) -> std::io::Result<()> {
        match &mut self.file {
            Some(file) => file.flush(),
            None => Ok(()),
        }
    }

    /// `cerr`'s tie to `cout` before a diagnostic line: the standard output's buffer is written
    /// out first (a file of `-out` is not; a failure comes back at the next write).
    fn before_diagnostic(&mut self) {
        if let Some(file) = &mut self.file {
            if matches!(file.get_ref(), OutputFile::Stdout(_)) {
                let _ = file.flush();
            }
        }
    }
}

/// The diagnostics as a stream, for the posts of the bytes that NCBI's object stream replaces
/// while it writes a Seq-annot (`write_document_to`).
struct DiagWriter<'d, 'a>(&'d mut Diag<'a>);

impl Write for DiagWriter<'_, '_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write_text(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The database that an earlier object opened: its name, residue type and title, and the
/// warnings of its taxid list, which NCBI posts each time it opens it again.
struct OpenedDatabase {
    name: Vec<u8>,
    protein: bool,
    title: Vec<u8>,
    warnings: Vec<u8>,
}

impl<'r, 'o> Run<'r, 'o> {
    /// `CRemoteBlast(istr)` and the loop of `LoadFromArchive`: the encoding guessed (an archive
    /// in ASN.1 binary or XML is refused by name), then each object read, its messages posted
    /// and, unless it is the error archive, formatted.
    fn objects(&mut self, archive: &mut dyn Read, diag: &mut Diag) -> anyhow::Result<()> {
        let guessed = guess_encoding(archive);
        match &guessed.encoding {
            Ok(Encoding::Text) => {}
            Ok(Encoding::Binary) => return Err(Refusal::BinaryArchive.into()),
            Ok(Encoding::Xml) => return Err(Refusal::XmlArchive.into()),
            Ok(Encoding::Other) => return Err(invalid_archive()),
            Err(GuessStop::Exception(exception)) => return Err(exception_error(exception)),
            Err(GuessStop::Undecided(what)) => return Err(Refusal::Undecided(what).into()),
        }
        let mut reader = ArchiveReader::new(guessed.into_text());
        while let Some((read, posts)) = reader.next_archive() {
            for post in &posts {
                post_line(&mut self.output, diag, &post.line(None));
            }
            let archive = read.map_err(read_error)?;
            self.object(archive, diag)?;
        }
        Ok(())
    }

    /// One object as `LoadFromArchive` reads it (its messages, then the request: `CImportStrategy`
    /// and, for the first search object, `GetSearchOptions`) and `PrintFormattedOutput` formats
    /// it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1727-1739
    // ```c++
    //         if (strategy.GetSubject()->IsDatabase())
    //             x_SetDatabase(strategy.GetSubject()->GetDatabase());
    //         else
    //             m_SubjectSequences = strategy.GetSubject()->SetSequences();
    //
    //         if(m_Service == "psi")
    //         {
    //         	// Would have errored out in CImportStrategy if we can't get queue search
    //         	m_FormatOpts.Reset( strategy.GetWebFormatOptions());
    //         }
    //         // Ignore return value, want side effect of setting fields.
    //         GetSearchOptions();
    // ```
    fn object(&mut self, archive: Archive, diag: &mut Diag) -> anyhow::Result<()> {
        self.post_messages(&archive.messages, diag)?;
        let search = match archive.search.map_err(stop_error)? {
            Some(search) => search,
            None => return Ok(()),
        };
        if let Subject::Database(name) = &search.request.subject {
            self.database = Some((name.clone(), search.request.db_is_protein()));
        }
        if self.search_options.is_none() {
            let options = archive_request(&search.request).map_err(|e| exception_error(&e))?;
            self.hitlist = options.hitlist_size();
            self.search_options = Some(options);
        }
        self.print(&search, diag)
    }

    /// The stored messages as `LoadFromArchive` posts them at the diagnostic level `Warning`:
    /// the text alone, in LOSAT's words (`wording::losat_bytes`); `Info` and `Trace` are below the
    /// level. A `Fatal` message ends NCBI's program (`abort`): refused.
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:505
    // ```c++
    //         SetDiagPostLevel(eDiag_Warning);
    // ```
    // (`LoadFromArchive` posts each message with `LOG_POST`, quoted at
    // `report::archive::read::messages`.)
    fn post_messages(&mut self, messages: &[Message], diag: &mut Diag) -> anyhow::Result<()> {
        for message in messages {
            match message.severity {
                Severity::Info | Severity::Trace => {}
                Severity::Fatal => return Err(Refusal::FatalMessage.into()),
                Severity::Warning | Severity::Error | Severity::Critical => {
                    let mut line = crate::wording::losat_bytes(&message.text).into_owned();
                    line.push(b'\n');
                    post_line(&mut self.output, diag, &line);
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:335-359
    // ```c++
    //     CNcbiOstream& out = args[kArgOutput].AsOutputFile();
    //     CFormattingArgs fmt_args(false, CFormattingArgs::eIsSAM) ;
    //
    //     CRef<CBlastOptionsHandle> opts_handle = m_RmtBlast->GetSearchOptions();
    //     CBlastOptions& opts = opts_handle->SetOptions();
    //     fmt_args.ExtractAlgorithmOptions(args, opts);
    //     ...
    //     const EBlastProgramType p = opts.GetProgramType();
    //     if((fmt_args.GetFormattedOutputChoice() == CFormattingArgs::eSAM) &&
    //        (p != eBlastTypeBlastn )) {
    // 		NCBI_THROW(CInputException, eInvalidInput,
    // 		                        "SAM format is only applicable to blastn results" );
    //     }
    //
    //     CRef<CBlastQueryVector> queries =
    // 		x_ExtractQueries(Blast_QueryIsProtein(p)?true:false);
    //     CRef<CScope> scope = queries->GetScope(0);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:426-430
    // ```c++
    //     if(UseXInclude(fmt_args, args[kArgOutput].AsString())) {
    //        	formatter.SetBaseFile(args[kArgOutput].AsString());
    //     }
    //     CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
    //     formatter.PrintProlog();
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:453-467
    // ```c++
    //         while (1)
    // 	{
    //         	BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                        fmt_args.GetFormattedOutputChoice());
    //     		ITERATE(CSearchResultSet, result, *results) {
    //     			if(isPsiBlast)
    //     			{
    //     				formatter.PrintOneResultSet(**result, queries,
    // 					m_RmtBlast->GetPsiNumberOfIterations());
    //     			}
    //     			else
    //     			{
    //     				formatter.PrintOneResultSet(**result, queries);
    //     			}
    //     		}
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2957-2962
    // ```c++
    //     	if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //     		hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    //     	}
    //     	else if (args.Exist(kArgNumAlignments) && args[kArgNumAlignments]) {
    //     		hitlist_size = args[kArgNumAlignments].AsInteger();
    //     	}
    // ```
    // (the formats up to 4: `-max_target_seqs`, or `max(m_NumDescriptions, m_NumAlignments)`
    // where `-num_descriptions` or `-num_alignments` is given, quoted at
    // `blastinput::report_counts::report_counts_from`.)
    /// Whether the command line's options set the hit list size of `ExtractAlgorithmOptions` for
    /// the output format `outfmt`, so that the options' own size is not kept.
    fn options_set_hitlist(&self, outfmt: i32) -> bool {
        let options = self.options;
        options.max_target_seqs.is_some()
            || options.num_alignments.is_some()
            || (outfmt <= output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES
                && options.num_descriptions.is_some())
    }

    /// `PrintFormattedOutput` of one search object, in NCBI's order: `-out` opened, the
    /// formatting options of the object (`extract_formatting` with the kept hit list size, its
    /// warnings), the check of SAM, the queries (`check_queries`), the database or the subject
    /// sequences (`subject_database`), the check of the XInclude formats (13 and 14 without
    /// `-out`); the output formats that LOSAT does not write are refused there, after NCBI's
    /// errors that come before any output. Then the results (`GetResultSet`: the alignments of
    /// each query, the masks), the prefetch's variable, and the Seq-annot of each query that has
    /// alignments (`x_PrintStructuredReport`).
    fn print(&mut self, search: &Search, diag: &mut Diag) -> anyhow::Result<()> {
        self.output.open(&self.options.out)?;
        // A negative hit list size (an archive that NCBI's searches do not write, probes
        // P9.005, P9.009): `ExtractAlgorithmOptions` keeps it where no option sets the size, and
        // posts its warning (`hitlist_size < 5`) as for 0, which the handler reads in its place.
        let formatting = extract_formatting(self.options, self.hitlist.max(0) as usize)?;
        post_line(&mut self.output, diag, &formatting.warnings);
        if self.hitlist >= 0 || self.options_set_hitlist(formatting.choice.number) {
            self.hitlist = i32::try_from(formatting.counts.hitlist_size)
                .expect("a hit list size that the handler read from an int");
        }
        let search_options = self
            .search_options
            .as_ref()
            .expect("the options of the first search object");
        let program = search_program(search_options.program);
        sam_check(
            &formatting.choice,
            matches!(program, SearchProgram::Named(program) if program.blastn),
        )?;
        let request = &search.request;
        check_queries(&request.queries)?;
        let database = match &self.database {
            Some((name, protein)) => Some(opened_database(
                &mut self.opened,
                name,
                *protein,
                search_options,
                self.settings,
            )),
            None => {
                check_subjects(&request.subject)?;
                // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:276-283
                // ```c++
                //         char* bl2seq_legacy = getenv("BL2SEQ_LEGACY");
                //         if (bl2seq_legacy) {
                //         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, false));
                //         }
                //         else {
                //             db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
                //         }
                // ```
                // (One rule for the variable: the search programs' refusal of it with
                // `-subject`, decision F3-6 (5) of line RG.)
                check_search_environment(PROGRAM, true)?;
                None
            }
        };
        let database = match database {
            Some((outcome, warnings)) => {
                post_line(&mut self.output, diag, &warnings);
                Some(outcome?)
            }
            None => None,
        };
        xinclude_check_for(&formatting.choice, Some(Path::new(&self.options.out)))?;
        if let Some(refusal) = refusals::output_format(formatting.choice.number) {
            return Err(refusal.into());
        }
        let format = match formatting.choice.number {
            output_format::ASN_TEXT => SerialFormat::AsnText,
            output_format::ASN_BINARY => SerialFormat::AsnBinary,
            _ => SerialFormat::Json,
        };
        let sets = search
            .results
            .alignment_sets(&request.queries)
            .map_err(stop_error)?;
        search
            .results
            .masks(request.queries.num_queries(), || {
                request.query_is_translated()
            })
            .map_err(stop_error)?;
        // `BlastFormatter_PreFetchSequenceData`: a value of `PRE_FETCH_SEQS_LIMIT` that is no
        // integer stops NCBI here (the search programs' one reader of the variable, decision
        // F3-6 (4)); NCBI's formats 8, 9 and 12 fetch nothing for any other value but `INT_MAX`
        // (`s_IsPrefetchFormat`, quoted at `common::database::prefetch_reads_seq_aligns`), with
        // which NCBI loads each subject through its data loaders: refused (decision C3-1).
        if prefetch_seqs_limit()? == Some(i32::MAX) {
            return Err(Refusal::Prefetch.into());
        }
        let alignments = search.results.alignments();
        // `PruneSeqalign`'s `unsigned int number`: a negative size is its value modulo 2^32.
        let number = self.hitlist as u32;
        for set in sets {
            let set = &alignments[set];
            if set.is_empty() {
                continue;
            }
            let pruned = prune_seqalign(set, number).map_err(stop_error)?;
            let program = match program {
                SearchProgram::Named(program) => program,
                SearchProgram::Other(task) => return Err(Refusal::Program(task).into()),
            };
            let desc = AnnotDesc {
                task: program.task,
                program: program.number,
                database: match &database {
                    Some(database) => DatabaseDesc::of_name(
                        &database.name,
                        &database.title,
                        program.subject_nucleotide,
                    ),
                    None => DatabaseDesc::None,
                },
            };
            let document = ArchiveDocument {
                desc,
                aligns: &pruned,
            };
            write_document_to(
                format,
                &document,
                None,
                &mut self.buf,
                self.output.writer(),
                &mut DiagWriter(&mut *diag),
            )?;
        }
        Ok(())
    }
}

/// Posts `text` (whole diagnostic lines, none when empty) after the standard output written
/// before it (`Output::before_diagnostic`).
fn post_line(output: &mut Output<'_>, diag: &mut Diag, text: &[u8]) {
    if text.is_empty() {
        return;
    }
    output.before_diagnostic();
    diag.write_text(text);
}

/// The Seq-annot of one query of an archive: the program's type and database in its `desc`, and
/// the query's Seq-aligns as read, cut by `prune_seqalign`.
struct ArchiveDocument<'d> {
    desc: AnnotDesc<'d>,
    aligns: &'d [V<'d>],
}

impl SeqAnnotDocument for ArchiveDocument<'_> {
    fn desc(&self) -> AnnotDesc<'_> {
        self.desc
    }

    /// Each Seq-align as the archive holds it (`asn_types::write::write_value`): NCBI's
    /// formatter writes the objects it read.
    fn write_aligns<E: Encoder>(&self, enc: &mut E) {
        for align in self.aligns {
            enc.element();
            write_value(enc, align.ty, align.value);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
// ```c++
// void CAlignFormatUtil::PruneSeqalign(const CSeq_align_set& source_aln,
//                                      CSeq_align_set& new_aln,
//                                      unsigned int number)
// {
//     CConstRef<CSeq_id> previous_id, subid;
//     bool is_first_aln = true;
//     unsigned int num_align = 0;
//     ITERATE(CSeq_align_set::Tdata, iter, source_aln.Get()){
//
//         if ((*iter)->GetSegs().IsDisc()) {
//             ++num_align;
//         } else {
//             subid = &((*iter)->GetSeq_id(1));
//             if(is_first_aln || (!is_first_aln && !subid->Match(*previous_id))){
//                 ++num_align;
//             }
//
//             if(num_align > number) {
//                  break;
//             }
//
//             is_first_aln = false;
//             previous_id = subid;
//         }
//         new_aln.Set().push_back(*iter);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
// ```c++
//         if (results.HasAlignments()) {
//     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
//             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
//             if(m_FormatType == CFormattingArgs::eAsnText)
//             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//             else
//             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//         }
// ```
/// `PruneSeqalign(aln_set, number)` (`number` is the hit list size, `m_HitlistSize`): the
/// Seq-aligns of the first `number` subjects, a new subject where the subject's id (row 1) does
/// not match the one before; a `disc` Seq-align counts as one more and is kept without the
/// check. NCBI writes the Seq-annot of a query with alignments also when nothing is left.
fn prune_seqalign<'a>(set: &[V<'a>], number: u32) -> Result<Vec<V<'a>>, Stop> {
    let mut kept = Vec::with_capacity(set.len());
    let mut previous: Option<&ParsedSeqId> = None;
    let mut count: u32 = 0;
    for &align in set {
        if align.field("segs").choice().0 == "disc" {
            count = count.wrapping_add(1);
        } else {
            let subject = seq_align_id(align, 1)?;
            if previous.is_none_or(|previous| compare(subject, previous) != Sic::Yes) {
                count = count.wrapping_add(1);
            }
            if count > number {
                break;
            }
            previous = Some(subject);
        }
        kept.push(align);
    }
    Ok(kept)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:170-180
// ```c++
// CBlastFormatterApp::x_BuildQueryFromPssm(const CPssmWithParameters& pssm)
// {
//     if ( !pssm.HasQuery() ) {
//         throw runtime_error("PSSM has no query");
//     }
//     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
//     const CSeq_entry& seq_entry = pssm.GetQuery();
//     if ( !seq_entry.IsSeq() ) {
//         throw runtime_error("Cannot have multiple queries in a PSSM");
//     }
//     SSeqLoc ssl = x_QueryBioseqToSSeqLoc(seq_entry.GetSeq(), scope);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:199-221
// ```c++
//     if (b4_queries->IsPssm()) {
//         retval->AddQuery(x_BuildQueryFromPssm(b4_queries->GetPssm()));
//     } else if (b4_queries->IsSeq_loc_list()) {
//         CRef<CScope> scope = m_QueryScopeSource->NewScope();
//         ITERATE(CBlast4_queries::TSeq_loc_list, seqloc,
//                 b4_queries->GetSeq_loc_list()) {
//             _ASSERT( !(*seqloc)->GetId()->IsLocal() );
//             CRef<CBlastSearchQuery> query(new CBlastSearchQuery(**seqloc,
//                                                                 *scope));
//             retval->AddQuery(query);
//         }
//     } else if (b4_queries->IsBioseq_set()) {
//         CTypeConstIterator<CBioseq> itr(ConstBegin(b4_queries->GetBioseq_set(),
//                                                    eDetectLoops));
//         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
//         for (; itr; ++itr) {
//             SSeqLoc ssl = x_QueryBioseqToSSeqLoc(*itr, scope);
//             CRef<CBlastSearchQuery> query(new CBlastSearchQuery(*ssl.seqloc,
//                                                                 *ssl.scope));
//             retval->AddQuery(query);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:164-165
// ```c++
//     CRef<CSeq_loc> seqloc(new CSeq_loc);
//     seqloc->SetWhole().Assign(*bioseq.GetFirstId());
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/sseqloc.hpp:261-266
// ```c++
//     /// Currently we only support whole or int.  Throw exception otherwise
//     void x_Validate() {
//         if (seqloc->IsWhole() || seqloc->IsInt()) return;
//         NCBI_THROW(CBlastException, eInvalidArgument,
//                    "Only whole or int typed seq_loc is supported for CBlastQueryVector");
//     }
// ```
/// `x_ExtractQueries` and the scope of its first query (`GetScope(0)`): NCBI's errors of a PSSM
/// without a query or with a set of them (`std::runtime_error`, exit status 255) and of a
/// location other than `whole` or `int` (`BLAST engine error`, exit status 3); a Bioseq without
/// an id (`*GetFirstId()` of a null pointer) and a request without a query (`GetScope(0)` of an
/// empty vector) are refused.
fn check_queries(queries: &Queries) -> anyhow::Result<()> {
    let count = match queries {
        Queries::Pssm(_) => {
            let query = queries
                .pssm()
                .and_then(|pssm| pssm.field("pssm").get("query"));
            let Some(query) = query else {
                return Err(NativeError::std_exception("PSSM has no query").into());
            };
            match query.choice() {
                ("seq", bioseq) if bioseq.field("id").item(0).is_none() => {
                    return Err(Refusal::Undecided(QUERY_WITHOUT_ID).into())
                }
                ("seq", _) => 1,
                _ => {
                    return Err(NativeError::std_exception(
                        "Cannot have multiple queries in a PSSM",
                    )
                    .into())
                }
            }
        }
        Queries::SeqLocList(_) => {
            let mut count = 0;
            for loc in queries.seq_locs() {
                if !matches!(loc.choice().0, "whole" | "int") {
                    return Err(engine_error(
                        "Only whole or int typed seq_loc is supported for CBlastQueryVector",
                    ));
                }
                count += 1;
            }
            count
        }
        Queries::BioseqSet(_) => {
            let bioseqs = queries.bioseqs();
            if bioseqs.iter().any(|bioseq| bioseq.ids().next().is_none()) {
                return Err(Refusal::Undecided(QUERY_WITHOUT_ID).into());
            }
            bioseqs.len()
        }
    };
    if count == 0 {
        return Err(Refusal::Undecided(
            "the scope of the first query of a request without a query",
        )
        .into());
    }
    Ok(())
}

/// What NCBI dereferences for a query Bioseq without an id.
const QUERY_WITHOUT_ID: &str = "the first id of a query that has no id";

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:234-246
// ```c++
//     if (remote_blast->GetSubjectSeqLocs().empty()) {
//         const list<CRef<CBioseq> > subjects =
//             remote_blast->GetSubjectSequences();
//         ITERATE(list<CRef<CBioseq> >, bioseq, subjects) {
//             subj_scope->AddBioseq(**bioseq);
//             CRef<CSeq_id> seqid = FindBestChoice((*bioseq)->GetId(),
//                                                  CSeq_id::BestRank);
//             const TSeqPos length = (*bioseq)->GetInst().GetLength();
//             CRef<CSeq_loc> sl(new CSeq_loc(*seqid, 0, length-1));
//             retval.push_back(SSeqLoc(sl, subj_scope));
//         }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:390-392
// ```c++
//         TSeqLocVector subjects = s_ConvertSubjects2TSeqLocVector(m_RmtBlast);
//         CRef<IQueryFactory> subject_factory(new CObjMgr_QueryFactory(subjects));
//         CRef<CScope> subj_scope = subjects.front().scope;
// ```
/// `s_ConvertSubjects2TSeqLocVector` of the object's subject sequences (no database object came
/// before): a Bioseq without its length is NCBI's unset member (invalid input); one without an
/// id NCBI ranks (`*seqid` of a null `CRef`) and no subject at all (`front()` of an empty vector;
/// the archive reader makes a `seq-loc-list` subject an empty list) are refused.
fn check_subjects(subject: &Subject) -> anyhow::Result<()> {
    let bioseqs = subject.bioseqs();
    for bioseq in &bioseqs {
        if bioseq.best_id().is_none() {
            return Err(Refusal::Undecided("the best id of a subject that has none").into());
        }
        if bioseq.length().is_none() {
            return Err(invalid_archive());
        }
    }
    if bioseqs.is_empty() {
        return Err(Refusal::Undecided("the first subject of a request without one").into());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:364-387
// ```c++
//     if (m_RmtBlast->IsDbSearch())
//     {
//         CRef<CBlast4_database> db = m_RmtBlast->GetDatabases();
//         _ASSERT(db);
//         _TRACE("Fetching results for " + Blast_ProgramNameFromType(p) + " on "
//                + db->GetName());
//         filtering_algorithm = m_RmtBlast->GetDbFilteringAlgorithmId();
//         CRef<CSearchDatabase> search_db(new CSearchDatabase(db->GetName(), db->IsProtein()
//                               ? CSearchDatabase::eBlastDbIsProtein
//                               : CSearchDatabase::eBlastDbIsNucleotide));
//
//         if(m_RmtBlast->GetTaxidList().size() > 0) {
//         	 CSeqDBGiList *gilist = new CSeqDBGiList();
//         	 gilist->AddTaxIds(m_RmtBlast->GetTaxidList());
//         	 search_db->SetGiList(gilist);
//         }
//
//         if(m_RmtBlast->GetNegativeTaxidList().size() > 0) {
//         	 CSeqDBGiList *gilist = new CSeqDBGiList();
//         	 gilist->AddTaxIds(m_RmtBlast->GetTaxidList());
//         	 search_db->SetNegativeGiList(gilist);
//         }
//         db_args->SetSearchDatabase(search_db);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:285-297
// ```c++
//         try {
//             // Try to open the BLAST database even for remote searches, as if
//             // it is available locally, it will be better to fetch the
//             // sequence data for formatting from this (local) source
//             CRef<CSeqDB> seqdb = search_db->GetSeqDb();
//             db_adapter.Reset(new CLocalDbAdapter(*search_db));
//             scope->AddDataLoader(RegisterOMDataLoader(seqdb), CBlastDatabaseArgs::kSubjectsDataLoaderPriority);
//             LOG_POST(Info <<"Add local loader " << search_db->GetDatabaseName());
//         } catch (const CSeqDBException&) {
//            	SetDiagPostLevel(eDiag_Critical);
//             string remote_loader = kEmptyStr;
//             try {
//             db_adapter.Reset(new CLocalDbAdapter(*search_db));
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:380-398
// ```c++
//     CRef<CSeqDB> seqdb(new CSeqDB(dbname, info.is_protein
//                           ? CSeqDB::eProtein : CSeqDB::eNucleotide));
//     if ( !seqdb ) {
//         return false;
//     }
//     info.name = seqdb->GetDBNameList();
//     info.definition = seqdb->GetTitle();
//     if (info.definition.empty())
//         info.definition = info.name;
//     ...
//     if (dbfilt_algorithm == -1) {
//         return true;
//     }
// ```
/// The database of a database search as `PrintFormattedOutput` opens it, through the search
/// programs' open path (`common::database::open_search_database`, the list of
/// `SeqDb::restrict`, the alias files' counts, the formatter's reopening of each name and its
/// filtering algorithm): the database (its title for the Seq-annot) and the warnings of the
/// taxid list, which come before the outcome. NCBI opens the database with the taxid list of the
/// first object's options (`GetSeqDb`); an error of NCBI's database reader there makes NCBI ask
/// its servers for the database (`CRemoteBlastDbDataLoader`): refused (`Refusal::Database`,
/// decisions F3-5 and F3-6 (1) of line RG). The negative taxid list is built of the positive one
/// (an empty list where the request has only a negative one), so it restricts nothing. The
/// formatter's reopening of each name without the list (`GetBlastDbInfo`) and the filtering
/// algorithm (`Filtering algorithm ID N is not supported.`) end with `BLAST Database error`,
/// exit status 2. An object whose database was opened for an earlier one reuses it: the options
/// are the first object's, so NCBI opens it again the same way and posts the same warnings.
fn opened_database<'d>(
    opened: &'d mut Option<OpenedDatabase>,
    name: &[u8],
    protein: bool,
    options: &BuiltOptions,
    settings: &ApplicationSettings,
) -> (anyhow::Result<&'d OpenedDatabase>, Vec<u8>) {
    let reused = opened
        .as_ref()
        .is_some_and(|database| database.name == name && database.protein == protein);
    if !reused {
        let mut warnings = Vec::new();
        match open_database(name, protein, options, settings, &mut warnings) {
            Ok(title) => {
                *opened = Some(OpenedDatabase {
                    name: name.to_vec(),
                    protein,
                    title,
                    warnings,
                });
            }
            Err(error) => return (Err(error), warnings),
        }
    }
    let database = opened.as_ref().expect("opened above or before");
    (Ok(database), database.warnings.clone())
}

/// `opened_database` of a database not opened before: its title (`pairwise_database`), the
/// warnings of its taxid list appended to `warnings`.
fn open_database(
    name: &[u8],
    protein: bool,
    options: &BuiltOptions,
    settings: &ApplicationSettings,
    warnings: &mut Vec<u8>,
) -> anyhow::Result<Vec<u8>> {
    let text = std::str::from_utf8(name)
        .map_err(|_| Refusal::Database(String::from_utf8_lossy(name).into_owned()))?;
    let seqtype = if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let mut db: SeqDb = match open_search_database(text, seqtype, settings) {
        Ok(db) => db,
        Err(DatabaseOpenError::Registry(reason)) => return Err(refusal(&reason, PROGRAM)),
        Err(DatabaseOpenError::Reader(error)) => return Err(remote_database(&error, text)),
    };
    if let Some(taxids) = options
        .taxid_list
        .as_ref()
        .filter(|taxids| !taxids.is_empty())
    {
        let list = DbList::Taxids(TaxidList::Positive(taxids.iter().copied().collect()));
        let mut posted = Vec::new();
        let restricted = db.restrict(&list, &mut posted);
        for warning in posted {
            list_warning(warning, warnings);
        }
        restricted.map_err(|error| remote_database(&error, text))?;
    }
    db.check_alias_values()
        .map_err(|error| remote_database(&error, text))?;
    // `GetSearchOptions` keeps the filtering algorithm of the options only when it is not -1
    // (remote_blast.cpp:2056-2059).
    let filtering = options.db_filtering_algorithm_id.filter(|&id| id != -1);
    db.read_name_records(text, filtering.is_some())
        .map_err(|error| database_error(error, PROGRAM))?;
    if let Some(id) = filtering {
        db.name_mask_details(id)
            .map_err(|error| name_mask_error(&error, PROGRAM))?;
    }
    Ok(pairwise_database(&db).definition)
}

/// A warning of the taxid list as blast_formatter posts it (no prefix: its `Run` sets none).
fn list_warning(warning: ListWarning, out: &mut Vec<u8>) {
    match warning {
        ListWarning::Log(text) => out.extend_from_slice(&text),
        ListWarning::Posted(text) => {
            out.extend_from_slice(b"Warning: ");
            out.extend_from_slice(text.as_bytes());
        }
        ListWarning::Error(text) => {
            out.extend_from_slice(b"Error: ");
            out.extend_from_slice(&text);
        }
        ListWarning::SeqId(warning) => {
            crate::writedb::fasta_source::defline_ids::post_seq_id_warning_of(
                &PROGRAM.to_ascii_lowercase(),
                warning,
                &mut |text| out.extend_from_slice(text),
            );
        }
    }
    out.push(b'\n');
}

/// An error of NCBI's database reader where `GetSeqDb` opens the database `name`: a
/// `CSeqDBException` makes NCBI ask its servers (refused); another exception of NCBI's leaves
/// for `CATCH_ALL`; a database that LOSAT does not read is LOSAT's refusal.
fn remote_database(error: &SeqDbError, name: &str) -> anyhow::Error {
    match error.ncbi_exception() {
        Some(exception) if !exception.kind().starts_with("CSeqDBException::") => {
            NativeError::exception(exception).into()
        }
        Some(_) => Refusal::Database(name.to_string()).into(),
        None if error.ncbi_message().is_some() => Refusal::Database(name.to_string()).into(),
        None => database_error_of(error, PROGRAM),
    }
}

/// `Invalid input format for BLAST Archive.`: the `CInputException` that `Run` throws for an
/// input of another format and for a `CSerialException` of the loop (`BLAST query/options
/// error`, exit status 1).
fn invalid_archive() -> anyhow::Error {
    options_error_bytes(b"Invalid input format for BLAST Archive.")
}

/// The archive reader's error, as the loop ends with it: a `CSerialException` is
/// `invalid_archive`; the end of the text inside the root object (`CEofException`, a
/// `CIOException` that `CATCH_ALL` ends without a post) keeps the exit status 0.
fn read_error(error: ReadError) -> anyhow::Error {
    match error {
        ReadError::Serial(_) => invalid_archive(),
        ReadError::Eof => NativeError::silent_end().into(),
        ReadError::Exception(exception) => exception_error(&exception),
        ReadError::Std(what) => NativeError::std_exception(what).into(),
        ReadError::Unsupported(text) => Refusal::UnreadText(text).into(),
    }
}

/// Where NCBI's code stops on the archive (`Stop`), as the run ends there: the serial exceptions
/// of a CHOICE or an unset member are `invalid_archive` (the loop's `catch`); an exception is
/// `exception_error`'s; a null pointer, an `abort` and a read of memory the input does not decide
/// are refused.
fn stop_error(stop: Stop) -> anyhow::Error {
    match stop {
        Stop::Exception(exception) => exception_error(&exception),
        Stop::InvalidSelection { .. } | Stop::Unassigned { .. } => invalid_archive(),
        Stop::NullPointer(what) | Stop::Abort(what) | Stop::Undefined(what) => {
            Refusal::Undecided(what).into()
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-228
// ```c++
//     catch (const blast::CBlastException& e) {                               \
//         const string& msg = e.GetMsg();                                     \
//         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
//             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
//             exit_code = BLAST_INPUT_ERROR;                                  \
//         } else if ((NStr::Find(msg, "Out of memory") != NPOS) ||            \
//             (NStr::Find(msg, "Failed to allocate") != NPOS)) {              \
//             LOG_POST(Error << "BLAST ran out of memory: " << e.GetMsg());   \
//             exit_code = BLAST_OUT_OF_MEMORY;                                \
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
//         }                                                                   \
// ```
/// An exception of NCBI's that leaves the loop, by its class: a serial exception is
/// `invalid_archive`; `CInputException` and `CArgException` are their `CATCH_ALL` branches; the
/// others are `NativeError::exception`'s (a `CBlastException`, `CRemoteBlastException` too, by its
/// code and message; a `CSeqDBException`; `Error: ` with `what()`).
fn exception_error(exception: &NcbiException) -> anyhow::Error {
    let kind = exception.kind();
    let class = kind.split("::").next().unwrap_or(kind);
    let msg = exception.get_msg();
    match class {
        "CSerialException" => invalid_archive(),
        class if crate::wording::SERIAL_EXCEPTION_CLASSES.contains(&class) => invalid_archive(),
        "CInputException" => options_error_bytes(msg),
        "CArgException" => command_line_argument_error(msg),
        _ => NativeError::exception(exception).into(),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
// ```c++
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
/// A failed write of the output, as the run ends with it: the formatter's stream throws on a
/// failure (`m_Outfile.exceptions(NcbiBadbit)`, blast_format.cpp:119): `BLAST failed to write
/// output`, exit status 6; a closed pipe ends the run at once with exit status 0 and no message
/// (where NCBI's process ends by SIGPIPE, `files::ends_like_sigpipe`). Any other error is
/// passed on.
fn write_failure(error: anyhow::Error) -> anyhow::Error {
    let Some(io) = error.downcast_ref::<std::io::Error>() else {
        return error;
    };
    if ends_like_sigpipe(io) {
        return NativeError::silent_end().into();
    }
    NativeError {
        exit: 6,
        message: b"BLAST failed to write output\n".to_vec(),
    }
    .into()
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod archive_tests;

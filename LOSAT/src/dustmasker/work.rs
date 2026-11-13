//! `CDustMaskApplication::Run`: the writer of `-outfmt`, the reader of `-infmt`, and the masks of
//! each nucleotide sequence.

use std::io::Write;

use crate::ncbi_app::files::{open_input_file, Cin, Ostream};
use crate::ncbi_exception::{err_post_x, NcbiException};
use crate::seqdb::mask::FilterProgram;
use crate::seqdb::SeqType;
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{prints_title, BdbError, MaskBdbReader};
use crate::seqmasks_io::run::{finish, serial_write_failure, write_failure, RunFailure};
use crate::seqmasks_io::writer::{MaskWriter, MaskedBioseq, Serial, WriteError};

use super::masker::dust_masks_skip_ns;

/// The arguments of the run, as NCBI's `CArgs` gives them (`AsInteger`, `AsString`,
/// `AsBoolean`).
pub(super) struct Options<'a> {
    pub input: &'a str,
    pub infmt: &'a str,
    pub outfmt: &'a str,
    pub window: i32,
    pub level: i32,
    pub linker: i32,
    pub parse_seqids: bool,
    pub hard_masking: bool,
    /// The `BLASTDB` of the configuration (`[BLAST] BLASTDB`), for `CSeqDB`'s search path
    /// (`ApplicationSettings::blastdb`: `Err` when LOSAT does not read the entry as NCBI does).
    pub blastdb_registry: Result<Option<&'a str>, &'a str>,
}

/// One sequence of the input (`CBioseq`): its Seq-ids, its title descriptor, its IUPAC letters.
pub(super) struct Record {
    pub ids: Vec<ParsedSeqId>,
    pub title: Option<Vec<u8>>,
    pub iupac: Vec<u8>,
}

/// How a run stops early.
pub(super) enum Failure {
    /// An NCBI exception that leaves `Run`: `NCBI_REPORT_EXCEPTION_X(16, "Application's
    /// execution failed", e)`, exit status 3.
    Ncbi(NcbiException),
    /// A `std::exception` that leaves `Run`: `ERR_POST_X(7, "Application's execution failed: "
    /// << e.what())`, exit status 3.
    Std(String),
    /// A case LOSAT does not reproduce.
    Refused(String),
    /// A `CArgException` of `AsInputFile` in `Run` (the reader's report of it).
    Arg(crate::ncbi_app::args::ArgException),
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`Ostream::closed`): the
    /// run ends at once, exit status 0, no message.
    ClosedPipe,
}

impl From<WriteError> for Failure {
    fn from(error: WriteError) -> Self {
        write_failure(error, super::PROGRAM)
    }
}

impl RunFailure for Failure {
    fn refused(message: String) -> Self {
        Failure::Refused(message)
    }

    fn closed_pipe() -> Self {
        Failure::ClosedPipe
    }

    fn std_exception(what: String) -> Self {
        Failure::Std(what)
    }

    fn stops_writing(&self) -> bool {
        matches!(self, Failure::Refused(_) | Failure::ClosedPipe)
    }
}

/// LOSAT's refusal of `what` (a noun phrase).
pub(super) fn unsupported(what: &str) -> String {
    crate::seqmasks_io::run::unsupported(what, super::PROGRAM)
}

/// `CMaskReader::GetNextSequence` of a source of sequences.
pub(super) trait SequenceSource {
    fn next_record(&mut self) -> Result<Option<Record>, Failure>;
}

/// What the application posts and returns for a failure: the diagnostic lines and the exit
/// status, or LOSAT's refusal.
pub(super) enum Report {
    Posts { lines: String, status: i32 },
    Refused(String),
    Arg(crate::ncbi_app::args::ArgException),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:668-690 (AppMain's handlers of `Run`)
// ```c++
//         catch (CException& e) {
//             NCBI_REPORT_EXCEPTION_X(16,
//                                     "Application's execution failed", e);
//             *got_exception = true;
//             *exit_code = 3;
//         }
//         catch (exception& e) {
//             ERR_POST_X(7, "Application's execution failed: " << e.what());
//             *got_exception = true;
//             *exit_code = 3;
//         }
// ```
// (`NCBI_USE_ERRCODE_X` of ncbiapp.cpp is `Corelib_App`, 106.)
/// The report of `failure`.
pub(super) fn report(failure: Failure) -> Report {
    match failure {
        Failure::Ncbi(exception) => Report::Posts {
            lines: exception.report_default_ex(106, 16, "Application's execution failed"),
            status: 3,
        },
        Failure::Std(what) => Report::Posts {
            lines: err_post_x(106, 7, &format!("Application's execution failed: {what}")),
            status: 3,
        },
        Failure::Refused(reason) => Report::Refused(reason),
        Failure::Arg(error) => Report::Arg(error),
        Failure::ClosedPipe => Report::Posts {
            lines: String::new(),
            status: 0,
        },
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:180-202
// ```c++
// string BuildAlgorithmParametersString(const CArgs& args)
// {
//     CNcbiOstrstream os;
//     if (args.Exist("locut") &&
//     ...
//     } else if (args.Exist("level") &&
//                args.Exist("linker") &&
//                args.Exist("window")) {
//         // DUST
//         os << "window=" << args["window"].AsInteger() << "; "
//            << "level=" << args["level"].AsInteger() << "; "
//            << "linker=" << args["linker"].AsInteger();
//     }
//     return CNcbiOstrstreamToString(os);
// }
// ```
/// `x_GetWriter`: the writer of `-outfmt`, or `Hard masking can only be applied for fasta
/// output`.
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:116-164
/// ```c++
///     if (args["hard_masking"].AsBoolean() && (format != "fasta")) {
///         throw runtime_error("Hard masking can only be applied for fasta output");
///     }
///     if (format == "interval") {
///         ...
///         retval = new CMaskWriterInt(output);
///     } else if (format == "acclist") {
///         ...
///         retval = new CMaskWriterTabular(output);
///     } else if (format == "fasta") {
///         ...
///         retval = new CMaskWriterFasta(output, hard_masking);
///     } else if (NStr::StartsWith(format, "seqloc_asn1_binary")) {
///         ...
///     } else if (NStr::StartsWith(format, "seqloc_")) {
///         ...
///         retval = new CMaskWriterSeqLoc(output, format);
///     } else if (NStr::StartsWith(format, "maskinfo_asn1_bin")) {
///         ...
///         retval =
///             new CMaskWriterBlastDbMaskInfo(output, format, 2,
///                                eBlast_filter_program_dust,
///                                BuildAlgorithmParametersString(args));
///     } else if (NStr::StartsWith(format, "maskinfo_")) {
///         ... (the same)
/// ```
/// (`-outfmt` takes only the nine values of its constraint, and the `seqloc_asn1_binary` branch
/// opens the stream in binary mode, which changes nothing on POSIX.)
fn writer(options: &Options<'_>) -> Result<MaskWriter, Failure> {
    if options.hard_masking && options.outfmt != "fasta" {
        return Err(Failure::Std(
            "Hard masking can only be applied for fasta output".to_string(),
        ));
    }
    let mask_info = |serial| {
        let parameters = format!(
            "window={}; level={}; linker={}",
            options.window, options.level, options.linker
        );
        MaskWriter::mask_info(
            Some(serial),
            2,
            FilterProgram::Dust,
            parameters.into_bytes(),
        )
    };
    Ok(match options.outfmt {
        "interval" => MaskWriter::Interval,
        "acclist" => MaskWriter::Tabular,
        "fasta" => MaskWriter::Fasta {
            hard_masking: options.hard_masking,
        },
        "seqloc_asn1_bin" => MaskWriter::SeqLoc(Serial::AsnBinary),
        "seqloc_asn1_text" => MaskWriter::SeqLoc(Serial::AsnText),
        "seqloc_xml" => MaskWriter::SeqLoc(Serial::Xml),
        "maskinfo_asn1_bin" => mask_info(Serial::AsnBinary),
        "maskinfo_asn1_text" => mask_info(Serial::AsnText),
        "maskinfo_xml" => mask_info(Serial::Xml),
        other => unreachable!("-outfmt {other} passed its constraint"),
    })
}

/// `CDustMaskApplication::Run` after the arguments: `stdin` is the standard input, which the
/// FASTA reader of `-in -` reads (`CMaskFastaReader`). The `maskinfo_*` writer writes its object even when the run
/// fails after the writer was made (its destructor runs as the exception leaves `Run`).
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:281-317
/// ```c++
///     CNcbiOstream& output_stream = GetArgs()[kOutput].AsOutputFile();
///     ...
///     unique_ptr<CMaskWriter> writer(x_GetWriter());
///     CMaskReader * reader = x_GetReader();
///
///     while( (aSeqEntry = reader->GetNextSequence()).NotEmpty() )
///     {
///         CScope scope( *om );
///         CSeq_entry_Handle seh = scope.AddTopLevelSeqEntry( *aSeqEntry );
///         CBioseq_CI bs_iter(seh, CSeq_inst::eMol_na);
///         for ( ;  bs_iter;  ++bs_iter)
///         {
///             CBioseq_Handle bsh = *bs_iter;
///             if (bsh.GetBioseqLength() == 0)
///                 continue;
///             CSeqVector data
///                 = bsh.GetSeqVector( CBioseq_Handle::eCoding_Iupac );
///             std::unique_ptr< duster_type::TMaskList > res = GetDustMasks_SkipNs(data, level, window, linker);
///             if (res.get()) {
///                 writer->Print(bsh, *res, GetArgs()["parse_seqids"] );
///             }
///         }
///     }
///
///     output_stream << flush;
///     return 0;
/// ```
/// `posts` gets the messages that the serial documents post as they are written
/// (`MaskWriter::print`), which the caller writes before the run's failure.
pub(super) fn run(
    options: &Options<'_>,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    posts: &mut Vec<u8>,
) -> Result<(), Failure> {
    let mut writer = writer(options)?;
    // A failed write leaves NCBI's stream bad, and the text writers go on silently (exit 0):
    // `out` drops the bytes after the first failed write.
    let mut out = Ostream::new(out);
    let result = masks_of_input(options, stdin, &mut writer, &mut out, posts);
    finish(
        writer,
        &mut out,
        result,
        options.outfmt,
        super::PROGRAM,
        posts,
    )
}

fn masks_of_input(
    options: &Options<'_>,
    stdin: &mut Cin<'_>,
    writer: &mut MaskWriter,
    out: &mut Ostream<'_>,
    posts: &mut Vec<u8>,
) -> Result<(), Failure> {
    let mut source = reader(options, stdin)?;
    let level = options.level as u32;
    let window = options.window as u32;
    let linker = options.linker as u32;
    while let Some(record) = source.next_record()? {
        if record.iupac.is_empty() {
            continue;
        }
        let mask = dust_masks_skip_ns(&record.iupac, level, window, linker);
        let seq = MaskedBioseq {
            ids: &record.ids,
            title: record.title.as_deref(),
            iupac: &record.iupac,
            protein: false,
        };
        writer.print(&seq, &mask, options.parse_seqids, out, posts)?;
        serial_write_failure::<Failure>(out, options.outfmt, super::PROGRAM)?;
    }
    Ok(())
}

/// `x_GetReader`: `CMaskFastaReader` of `-in`, `CMaskBDBReader` of the database `-in`, or
/// `Unknown input format`.
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:167-182
/// ```c++
///     if( format == "fasta" ) {
///         CNcbiIstream& input_stream  = GetArgs()[kInput].AsInputFile();
///         return new CMaskFastaReader(
///                 input_stream, true, args["parse_seqids"] );
///     }
///     else if( format == "blastdb" ) {
///         return new CMaskBDBReader( args[kInput].AsString() );
///     }
///     else {
///         throw runtime_error( "Unknown input format" );
///     }
/// ```
///
/// `-in` is opened here, as `x_GetReader` asks for it (`AsInputFile`, fasta input only), after
/// the writer: an `-in` that does not open still lets a `maskinfo_*` writer write its object.
fn reader<'a>(
    options: &Options<'_>,
    stdin: &'a mut Cin<'_>,
) -> Result<Box<dyn SequenceSource + 'a>, Failure> {
    match options.infmt {
        "fasta" => {
            let input = open_input_file("in", options.input, stdin).map_err(Failure::Arg)?;
            super::fasta::source(input, options.parse_seqids)
        }
        "blastdb" => Ok(Box::new(DatabaseSource::open(
            options.input,
            options.blastdb_registry,
            prints_title(options.outfmt),
        )?)),
        _ => Err(Failure::Std("Unknown input format".to_string())),
    }
}

/// `CMaskBDBReader(name)` of a nucleotide database (`seqmasks_io::bdb_reader`).
struct DatabaseSource {
    reader: MaskBdbReader,
}

impl DatabaseSource {
    fn open(
        name: &str,
        registry: Result<Option<&str>, &str>,
        title_printed: bool,
    ) -> Result<Self, Failure> {
        let reader = MaskBdbReader::open(
            name,
            SeqType::Nucleotide,
            registry,
            title_printed,
            super::PROGRAM,
        )
        .map_err(database_failure)?;
        Ok(DatabaseSource { reader })
    }
}

/// An error of the database reader, reported as `AppMain` reports an exception of `Run`
/// (`report`), or LOSAT's refusal.
fn database_failure(error: BdbError) -> Failure {
    match error {
        BdbError::Ncbi(exception) => Failure::Ncbi(exception),
        BdbError::Refused(reason) => Failure::Refused(reason),
        BdbError::Unsupported(what) => Failure::Refused(unsupported(&what)),
    }
}

impl SequenceSource for DatabaseSource {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        let record = self.reader.next_record().map_err(database_failure)?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: Some(record.title),
            iupac: record.residues,
        }))
    }
}

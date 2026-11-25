//! `CDustMaskApplication::Run`: the writer of `-outfmt`, the reader of `-infmt`, and the masks of
//! each nucleotide sequence.

use std::io::Write;

use crate::ncbi_app::files::{open_input_file, Cin, Ostream};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{err_post_x, NcbiException, ThrowSite};
use crate::seqdb::mask::FilterProgram;
use crate::seqdb::SeqType;
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{prints_title, BdbError, MaskBdbReader};
use crate::seqmasks_io::run::{
    finish, post_read_texts, serial_write_failure, write_failure, RunFailure,
};
use crate::seqmasks_io::writer::{MaskWriter, MaskedBioseq, Serial, WriteError};

use super::masker::dust_masks_skip_ns;

/// The arguments of the run, as NCBI's `CArgs` gives them (`AsInteger`, `AsString`,
/// `AsBoolean`).
pub(super) struct Options<'a> {
    /// `-in`: the bytes of a FASTA file's name or of a database's.
    pub input: &'a std::ffi::OsStr,
    /// `-infmt`, which NCBI compares as its bytes.
    pub infmt: std::borrow::Cow<'a, [u8]>,
    pub outfmt: &'a str,
    pub window: i32,
    pub level: i32,
    pub linker: i32,
    pub parse_seqids: bool,
    pub hard_masking: bool,
    /// The settings of the database reader (`ApplicationSettings::database`): the `BLASTDB` of
    /// the configuration (`[BLAST] BLASTDB`) for `CSeqDB`'s search path, and the checks of its
    /// reads.
    pub database: crate::blastinput::ncbi_environment::DatabaseSettings<'a>,
}

/// One sequence of the input (`CBioseq`): its Seq-ids, its title descriptor, its IUPAC letters.
pub(super) struct Record {
    pub ids: Vec<ParsedSeqId>,
    pub title: Option<Vec<u8>>,
    /// The name of a database sequence's source, from which its empty title is made
    /// (`seqmasks_io::bdb_reader::BdbRecord::taxname`); none for a FASTA record.
    pub taxname: Option<std::sync::Arc<[u8]>>,
    /// Whether a database sequence's Bioseq has a source (`BdbRecord::has_source`).
    pub has_source: bool,
    pub iupac: Vec<u8>,
    /// The taxonomy database's texts of reading a database sequence
    /// (`seqmasks_io::bdb_reader::BdbRecord::texts`); none for a FASTA record.
    pub texts: Vec<String>,
}

/// How a run stops early.
pub(super) enum Failure {
    /// An NCBI exception that leaves `Run`: reported by `NCBI_REPORT_EXCEPTION_X` (subcode 16,
    /// title "Application's execution failed"), exit status 3.
    Ncbi(NcbiException),
    /// A `std::exception` that leaves `Run`: posted by `ERR_POST_X` (subcode 7) as "Application's
    /// execution failed: " and `e.what()`, exit status 3.
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
    Posts { lines: Vec<u8>, status: i32 },
    Refused(String),
    Arg(crate::ncbi_app::args::ArgException),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:668-690,725-731 (AppMain's handlers of `Run`)
// (`NCBI_USE_ERRCODE_X` of ncbiapp.cpp is `Corelib_App`, 106.)
/// The report of `failure`.
pub(super) fn report(failure: Failure) -> Report {
    match failure {
        Failure::Ncbi(exception) => Report::Posts {
            lines: exception.report_default_ex(106, 16, "Application's execution failed"),
            status: 3,
        },
        Failure::Std(what) => Report::Posts {
            lines: err_post_x(106, 7, &format!("Application's execution failed: {what}"))
                .into_bytes(),
            status: 3,
        },
        Failure::Refused(reason) => Report::Refused(reason),
        Failure::Arg(error) => Report::Arg(error),
        Failure::ClosedPipe => Report::Posts {
            lines: Vec::new(),
            status: 0,
        },
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:180-202
/// `x_GetWriter`: the writer of `-outfmt`, or `Hard masking can only be applied for fasta
/// output`.
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:116-164
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
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:276-317
pub(super) fn run(
    options: &Options<'_>,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let mut writer = writer(options)?;
    // A failed write leaves NCBI's stream bad, and the text writers go on silently (exit 0):
    // `out` drops the bytes after the first failed write.
    let mut out = Ostream::new(out);
    let result = masks_of_input(options, stdin, &mut writer, &mut out, diag);
    finish(
        writer,
        &mut out,
        diag,
        result,
        options.outfmt,
        super::PROGRAM,
    )
}

fn masks_of_input(
    options: &Options<'_>,
    stdin: &mut Cin<'_>,
    writer: &mut MaskWriter,
    out: &mut Ostream<'_>,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let mut source = reader(options, stdin)?;
    let level = options.level as u32;
    let window = options.window as u32;
    let linker = options.linker as u32;
    while let Some(record) = source.next_record()? {
        post_read_texts::<Failure>(&record.texts, out, diag)?;
        if record.iupac.is_empty() {
            continue;
        }
        let mask = dust_masks_skip_ns(&record.iupac, level, window, linker);
        let seq = MaskedBioseq {
            ids: &record.ids,
            handle_seq_id: None,
            title: record.title.as_deref(),
            taxname: record.taxname.as_ref(),
            has_source: record.has_source,
            iupac: &record.iupac,
            protein: false,
        };
        writer.print(&seq, &mask, options.parse_seqids, out, diag)?;
        serial_write_failure::<Failure>(out, options.outfmt, super::PROGRAM)?;
    }
    Ok(())
}

/// `x_GetReader`: `CMaskFastaReader` of `-in`, `CMaskBDBReader` of the database `-in`, or
/// `Unknown input format`.
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:167-182
///
/// `-in` is opened here, as `x_GetReader` asks for it (`AsInputFile`, fasta input only), after
/// the writer: an `-in` that does not open still lets a `maskinfo_*` writer write its object.
fn reader<'a>(
    options: &Options<'_>,
    stdin: &'a mut Cin<'_>,
) -> Result<Box<dyn SequenceSource + 'a>, Failure> {
    match &*options.infmt {
        b"fasta" => {
            let input = open_input_file("in", options.input, stdin).map_err(Failure::Arg)?;
            super::fasta::source(input, options.parse_seqids)
        }
        b"blastdb" => Ok(Box::new(DatabaseSource::open(
            &crate::utils::os_bytes::os_str_bytes(options.input),
            options.database,
            prints_title(options.outfmt),
        )?)),
        _ => Err(Failure::Std("Unknown input format".to_string())),
    }
}

/// `CMaskBDBReader` of `name`, a nucleotide database (`seqmasks_io::bdb_reader`).
struct DatabaseSource {
    reader: MaskBdbReader,
}

impl DatabaseSource {
    fn open(
        name: &[u8],
        database: crate::blastinput::ncbi_environment::DatabaseSettings<'_>,
        title_printed: bool,
    ) -> Result<Self, Failure> {
        let reader = MaskBdbReader::open(
            name,
            SeqType::Nucleotide,
            database,
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

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:292-296
// NCBI reference (598d8ae6): c++/src/objmgr/seq_entry_info.cpp:112-123
// `CMaskBDBReader` gives an empty Seq-entry for an OID whose `GetBioseq` gives no Bioseq (no
// letters, or every definition line erased by the alias files' `OID_MASK_TYPE`), and the object
// manager's check of its choice throws on it: the run ends there, after the output of the earlier
// sequences (probes dust_omnx_int and dust_ntall_fasta of line RG D1b-2, exit status 3).
// windowmasker's masking skips such an entry (`windowmasker::work`); segmasker's too.
const UNSET_SEQ_ENTRY: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_entry_info.cpp",
    line: 119,
    kind: "CUnassignedMember::eGet",
    function: None,
};

impl SequenceSource for DatabaseSource {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        let record = self.reader.next_record().map_err(database_failure)?;
        if record
            .as_ref()
            .is_some_and(|record| record.residues.is_empty())
        {
            return Err(Failure::Ncbi(NcbiException::new(
                UNSET_SEQ_ENTRY,
                "Seq_entry.set",
            )));
        }
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: Some(record.title),
            taxname: record.taxname,
            has_source: record.has_source,
            iupac: record.residues,
            texts: record.texts,
        }))
    }
}

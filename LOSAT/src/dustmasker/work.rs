//! `CDustMaskApplication::Run`: the writer of `-outfmt`, the reader of `-infmt`, and the masks of
//! each nucleotide sequence.

use std::io::Write;

use crate::ncbi_exception::{err_post_x, NcbiException, ThrowSite};
use crate::report::db_subject::bioseq_ids;
use crate::report::defline::bioseq_title_from_source;
use crate::seqdb::alias::SearchPath;
use crate::seqdb::{SeqDb, SeqType};
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::writer::{FilterProgram, MaskWriter, MaskedBioseq, Serial, WriteError};

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
    /// The `BLASTDB` of the configuration (`[BLAST] BLASTDB`), for `CSeqDB`'s search path.
    pub blastdb_registry: Option<&'a str>,
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
}

impl From<WriteError> for Failure {
    fn from(error: WriteError) -> Self {
        match error {
            WriteError::Io(error) => {
                Failure::Refused(format!("writing the output failed: {error}"))
            }
            WriteError::Refused(reason) => Failure::Refused(reason),
        }
    }
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
        MaskWriter::mask_info(serial, 2, FilterProgram::Dust, parameters.into_bytes())
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

/// `CDustMaskApplication::Run` after the arguments: `fasta_source` makes the FASTA reader of
/// `-in` (`CMaskFastaReader`). The `maskinfo_*` writer writes its object even when the run
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
pub(super) fn run(
    options: &Options<'_>,
    fasta_source: &mut dyn FnMut() -> Result<Box<dyn SequenceSource>, Failure>,
    out: &mut dyn Write,
) -> Result<(), Failure> {
    let mut writer = writer(options)?;
    let result = masks_of_input(options, fasta_source, &mut writer, out);
    let finished = writer.finish(out);
    result?;
    finished?;
    out.flush()
        .map_err(|error| Failure::Refused(format!("writing the output failed: {error}")))
}

fn masks_of_input(
    options: &Options<'_>,
    fasta_source: &mut dyn FnMut() -> Result<Box<dyn SequenceSource>, Failure>,
    writer: &mut MaskWriter,
    out: &mut dyn Write,
) -> Result<(), Failure> {
    let mut source = reader(options, fasta_source)?;
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
        };
        writer.print(&seq, &mask, options.parse_seqids, out)?;
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
fn reader(
    options: &Options<'_>,
    fasta_source: &mut dyn FnMut() -> Result<Box<dyn SequenceSource>, Failure>,
) -> Result<Box<dyn SequenceSource>, Failure> {
    match options.infmt {
        "fasta" => fasta_source(),
        "blastdb" => Ok(Box::new(DatabaseSource::open(
            options.input,
            options.blastdb_registry,
        )?)),
        _ => Err(Failure::Std("Unknown input format".to_string())),
    }
}

/// `CMaskBDBReader`: the sequences of a nucleotide BLAST database in OID order, each as
/// `CSeqDB::GetBioseq` makes it (the Seq-ids of the first definition line, the title joined
/// from all of them).
///
/// NCBI reference (598d8ae6): c++/include/objtools/seqmasks_io/mask_bdb_reader.hpp:55-61
/// ```c++
///     CMaskBDBReader( const string & dbname, bool is_nucl = true)
///         : CMaskReader( cin ),
///           seqdb_( new CSeqDB( dbname, (is_nucl
///                                        ? CSeqDB::eNucleotide
///                                        : CSeqDB::eProtein)) ),
///           oid_( 0 )
///     {}
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp:41-55
/// ```c++
/// CRef< CSeq_entry > CMaskBDBReader::GetNextSequence()
/// {
///     CRef< CSeq_entry > entry( null );
///
///     if( oid_ < seqdb_->GetNumOIDs() ) {
///         CRef< CBioseq > seq = seqdb_->GetBioseq( oid_++ );
///         entry.Reset( new CSeq_entry );
///         if( seq.NotEmpty() ) entry->SetSeq( *seq );
///     }
///
///     return entry;
/// }
/// ```
struct DatabaseSource {
    db: SeqDb,
    oid: u32,
}

impl DatabaseSource {
    fn open(name: &str, registry: Option<&str>) -> Result<Self, Failure> {
        let path = SearchPath::from_environment(registry).map_err(database_failure)?;
        let db = SeqDb::open_in(name, SeqType::Nucleotide, &path).map_err(database_failure)?;
        Ok(DatabaseSource { db, oid: 0 })
    }
}

/// An error of the database reader: NCBI's `CSeqDBException` (`SeqDbError::ncbi_message`), or
/// LOSAT's refusal.
fn database_failure(error: crate::seqdb::SeqDbError) -> Failure {
    match error.ncbi_message() {
        Some(text) => Failure::Ncbi(NcbiException::new(
            ThrowSite::generated("CSeqDBException::eInvalid"),
            text,
        )),
        None => Failure::Refused(error.to_string()),
    }
}

impl SequenceSource for DatabaseSource {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        if self.oid >= self.db.num_seqs() {
            return Ok(None);
        }
        let oid = self.oid;
        self.oid += 1;
        let set = self.db.defline_set(oid).map_err(database_failure)?;
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err(Failure::Refused(format!(
                "sequence {oid} of the database: a definition line with bytes NCBI reads as '#' \
                 (it writes an error each time it reads them), which is not supported by \
                 LOSAT's dustmasker"
            )));
        }
        let ids = set
            .deflines
            .first()
            .map(|defline| defline.ids.clone())
            .unwrap_or_default();
        let title = set.bioseq_title();
        // NCBI makes the title of a sequence without one from its BioSource (the taxonomy
        // names of a definition line with a taxid; seqdbvol.cpp:1055-1160).
        let has_source = set
            .deflines
            .iter()
            .any(|defline| defline.taxid.is_some_and(|taxid| taxid > 0));
        if has_source && bioseq_title_from_source(&title, &bioseq_ids(&ids)) {
            return Err(Failure::Refused(format!(
                "sequence {oid} of the database has no title and a taxid (NCBI makes its title \
                 from the taxonomy), which is not supported by LOSAT's dustmasker"
            )));
        }
        let length = self.db.nucleotide_length(oid).map_err(database_failure)?;
        let iupac = self
            .db
            .nucleotide_iupac_range(oid, 0, length)
            .map_err(database_failure)?;
        Ok(Some(Record {
            ids,
            title: Some(title),
            iupac,
        }))
    }
}

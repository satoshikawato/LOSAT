//! The PSSM file of `-in_pssm` (psiblast, series PS item PS4a; tblastn's PSI-TBLASTN, item PS7): the file read
//! as NCBI's `CPsiBlastArgs` reads it for both programs, NCBI's errors for files it cannot read, NCBI's checks of
//! the PSSM (`CPsiBlastValidate::Pssm`), the query record made from the PSSM's Bioseq, the score rows
//! (`CScorematPssmConverter::GetScores`) and the statistical values of `PsiBlastSetupScoreBlock`, the frequency
//! ratios (item PS4b), the scores of a PSSM of frequency ratios only (the PSSM engine's score block), and LOSAT's
//! refusals of what the PSSM search does not read (DECISIONS 11 of line PS). The reader belongs to the PSSM search
//! (`PsiSearch`), which psiblast and tblastn both run; their own options stay with them
//! (`algorithm::psiblast::in_pssm`, `algorithm::tblastn::in_pssm`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1363

use std::path::Path;

use anyhow::{bail, Result};

use super::rows::QueryPssm;
use super::{PsiSearch, PssmQueryRequest, PssmStats};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::pssm_text::PssmReadError;
use crate::config::ScoringMatrix;
use crate::core::blast_util::{BLASTAA_SIZE, BLAST_SCORE_MIN};
use crate::core::psi::{KarlinBlk, PsiScoreBlk};
use crate::ncbi_exception::{NcbiException, Severity};
use crate::report::serial::asn_types::ReadError;
use crate::report::serial::pssm::{ObjectId, PssmFinalData, PssmWithParameters, SeqData, SeqEntry};
use crate::seqid::ParsedSeqId;
use crate::utils::format_guess::{Format, FormatGuess, Stop};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-118
/// The letters of NCBIstdaa, the one table of the PSSM reader and of psiblast's PSSM files (the ASCII PSSM, the query
/// letters of a PSSM file; `algorithm::psiblast::in_pssm_files`).
pub(crate) const NCBISTDAA_TO_AMINOACID: &[u8; BLASTAA_SIZE] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
/// `BLAST_SCORE_MAX` (`INT2_MAX`, blast_stat.h:122): the largest score of one letter pair.
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

/// No NCBI counterpart: LOSAT's bound on the sum of a PSSM's columns' largest positive scores
/// (`read_in_pssm`'s range rule), half of `Int4`'s range as NCBI's `MININT`.
const PSSM_SCORE_SUM_LIMIT: i64 = 1 << 30;

/// What `-in_pssm` makes of the file: the `switch` on the format that `CFormatGuess` guesses for the stream.
#[derive(Debug)]
enum PssmFileFormat {
    /// `eTextASN`: the text reader reads the stream from `start` (the start of the file, or its end after
    /// `TestFormatNewick` read a NEXUS file through).
    TextAsn { start: usize },
    /// `eBinaryASN`: the binary reader reads the stream from `start` (the stream read as `MSerial_AsnBinary`).
    BinaryAsn { start: usize },
    /// `eXml`: the XML reader reads the stream from `start` (the stream read as `MSerial_Xml`).
    Xml { start: usize },
    /// A format that the `switch` sends to "Unsupported format for PSSM" (`eUnknown` and every format but the
    /// three).
    Unsupported,
    /// The guess threw (`TestFormatGlimmer3`'s `NStr::StringToInt`): `CATCH_ALL`'s report (exit status 255).
    Exception(NcbiException),
    /// LOSAT's refusal (`what`, then `detail` and the path in parentheses): a file whose guess reads memory the
    /// file does not decide.
    Refused {
        what: &'static str,
        detail: &'static str,
    },
}

/// The refusal of a guess that reads memory its file does not decide.
const UNDECIDED_GUESS: &str =
    "a PSSM file whose format NCBI's format guess decides by reading memory past the file";

/// The format that `CFormatGuess` guesses for the stream (`utils::format_guess`, the whole order of `sm_CheckOrder`
/// without hints) and the `switch` of `CPsiBlastArgs` on it. A test that reads past its data decides nothing here when
/// the formats it could make the guess are all unsupported (`TestFormatFasta` before `eTextASN` and `eBinaryASN`).
fn guess_format(bytes: &[u8]) -> PssmFileFormat {
    let guessed = match FormatGuess::over(bytes).guess_format() {
        Ok(guessed) => guessed,
        Err(Stop::Exception(exception)) => return PssmFileFormat::Exception(exception),
        Err(Stop::Undecided(detail)) => {
            return PssmFileFormat::Refused {
                what: UNDECIDED_GUESS,
                detail,
            }
        }
    };
    match guessed.format {
        Format::TextAsn | Format::BinaryAsn if guessed.fasta_undecided => PssmFileFormat::Refused {
            what: UNDECIDED_GUESS,
            detail: "TestFormatFasta searches for '|' past a buffer of 10 bytes or fewer",
        },
        Format::TextAsn => PssmFileFormat::TextAsn {
            start: guessed.stream_pos,
        },
        Format::BinaryAsn => PssmFileFormat::BinaryAsn {
            start: guessed.stream_pos,
        },
        Format::Xml => PssmFileFormat::Xml {
            start: guessed.stream_pos,
        },
        _ => PssmFileFormat::Unsupported,
    }
}

/// Who reads `-in_pssm`: the program that LOSAT's refusals name, the options' matrix, with
/// which NCBI's PSSM engine computes the scores of a PSSM of frequency ratios only
/// (`scores_from_freq_ratios`), and whether a search strategy's request holds the PSSM
/// (`request`). psiblast reads with `PSIBLAST` (BLOSUM62: its `-in_pssm` refuses the other
/// matrices, DECISIONS 11 of line PS); tblastn's PSI-TBLASTN with its own options
/// (`algorithm::tblastn::in_pssm`, DECISIONS 35).
#[derive(Clone, Copy, Debug)]
pub(crate) struct PssmReader<'a> {
    /// The program of LOSAT's refusals (`PSIBLAST`, `TBLASTN`).
    pub(crate) program: &'static str,
    /// The options' matrix name (`opts.GetMatrixName()`).
    pub(crate) matrix_name: &'a str,
    /// Whether the run writes a search strategy (`-export_search_strategy`) and searches a database
    /// (`-db`, not `-subject`): with both its request holds the PSSM (`super::pssm_request_written`),
    /// so the search keeps the PSSM as read (`PsiSearch::request`; only such a run pays for that
    /// copy); with the strategy alone NCBI stops before any report (`PsiSearch::stop_before_search`).
    pub(crate) strategy: bool,
    pub(crate) database: bool,
    /// The run's BLAST archive (outfmt 11), which saves each of the reader's posts as one message
    /// of its severity (`diag_log::Posts::write`); `None` without one.
    pub(crate) archive: Option<&'a crate::report::archive::ArchiveHandle>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
/// Whether `error` is NCBI's error of a `-in_pssm` file that cannot be opened (`cli::inaccessible`),
/// which NCBI raises inside `Run`: the options' handler opens the file there, so its `CATCH_ALL`
/// writes the error archive of outfmt 11 for it (NCBI BLAST+ 2.17.0, tblastn's and psiblast's
/// `-in_pssm` of a missing file with `-subject` and with `-db`), as for the files that
/// `report::archive::error_lines` takes as the subjects' (`subject_file`).
pub(crate) fn file_error_in_run(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<crate::cli::NativeError>()
        .is_some_and(|native| {
            native
                .message
                .starts_with(b"Command line argument error: Argument \"in_pssm\".")
        })
}

impl PssmReader<'static> {
    /// psiblast's reader.
    pub(crate) const PSIBLAST: Self = Self {
        program: "PSIBLAST",
        matrix_name: "BLOSUM62",
        strategy: false,
        database: true,
        archive: None,
    };
}

impl PssmReader<'_> {
    /// The refusal of the PSSM search's limits ("... is not supported by LOSAT's PSIBLAST").
    fn refuse(&self, what: impl std::fmt::Display) -> anyhow::Error {
        anyhow::anyhow!("{what} is not supported by LOSAT's {}", self.program)
    }
}

/// The bytes of an input file argument that is read whole (`-in_pssm`, psiblast's `-in_msa`;
/// `name` is the argument's name in NCBI's message): `-` is the standard input; a file that
/// cannot be opened is NCBI's "File is not accessible"; a stream that opens and cannot be read (a
/// directory, as the file or as the standard input) gives the bytes it gave before it failed,
/// none for a directory, as NCBI's stream does.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
/// `-` is the standard input (`kArgPSIInputChkPntFile` and `kArgMSAInputFile` are `eInputFile`s),
/// as `-query -` is; NCBI's messages name the file `-`. `cin` is always there, so `-` is never
/// "not accessible": a standard input that cannot be read (`< .`, a directory) is an empty
/// stream to the reader (NCBI BLAST+ 2.17.0: `-in_pssm - < .` ends with "Unsupported format for
/// PSSM", `-in_msa - < .` with the alignment reader's failure; review probes `manual/closed` of
/// line PS item PS6, O4).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-732
/// Only a file that cannot be opened is "not accessible": one that opens and cannot be read (a
/// directory) gives the reader no bytes, which NCBI's format guess takes as an unsupported
/// format (probe `pf.dir` of line PS PS7).
pub(crate) fn read_input_argument(name: &str, path: &Path) -> Result<Vec<u8>> {
    if path.as_os_str() == "-" {
        return Ok(bytes_until_failure(std::io::stdin().lock()));
    }
    std::fs::File::open(path)
        .map(bytes_until_failure)
        .map_err(|_| crate::cli::inaccessible(name, path))
}

/// The bytes that `stream` gives before it ends or fails: a failed read is the end of the stream
/// (none from a directory), as for NCBI's readers of an `istream` whose read failed.
fn bytes_until_failure(mut stream: impl std::io::Read) -> Vec<u8> {
    let mut bytes = Vec::new();
    let _ = stream.read_to_end(&mut bytes);
    bytes
}

/// Reads `-in_pssm path`: NCBI's errors for a file it cannot read or guess, the messages of NCBI's
/// reader (to `diagnostics`), LOSAT's refusals of what it does not read, then the PSSM search's
/// query, rows and values, and the PSSM as NCBI's object holds it when the search starts
/// (`InPssmFile`, which psiblast's rounds and PSSM files write back).
pub(crate) fn read_in_pssm(
    path: &Path,
    reader: PssmReader<'_>,
    diagnostics: &mut dyn std::io::Write,
) -> Result<(PsiSearch, InPssmFile)> {
    let bytes = read_input_argument("in_pssm", path)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-263 (the clauses of
    // `CATCH_ALL` that an error of the reader other than its `CSerialException` reaches)
    // NCBI's reader posts its messages as it reads (`Post`), before its error when it stops. Its
    // `CSerialException` is the handler's "Unrecognized format"; a `CUtilException` of a number
    // or a `CStringException` of a tag leaves the handler and `CATCH_ALL` reports it with the
    // lines of the reader's stack (exit status 255), and so does a `std::length_error` of
    // Dense-seg's reserve hook (`Error: vector::reserve`) or of the reserve of a binary string or
    // byte block; a `bad_alloc` of such a reserve is "BLAST ran out of memory" (exit status 4),
    // each after the byte block's post 803.4 where a byte block's reserve threw; a
    // `CEofException` (a `CIOException`, the input ending directly in the root object) is caught
    // without a word and the run ends with exit status 0 (NCBI BLAST+ 2.17.0, `PssmWithParameters
    // ::=` alone; line RG D3-1 probe e01; a binary file cut after its last byte, line RG D3c-3a).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
    // (`CBlastAppDiagHandler::Post`, cited in `report::diag_log`): the BLAST archive saves each post as one message of its
    // severity; the reader's posts are errors (the post 803.7 of a replaced byte, whose text has
    // three lines, the error post of `PostRead` and 803.4 of a byte block; NCBI BLAST+ 2.17.0
    // tblastn -outfmt 11 of a PSSM whose title has a control byte and whose ratios overflow: `code 2`), and
    // the binary reader's warning of a UTF8String tag for a VisibleString.
    // NCBI's reader recurses once per value on the main thread's stack and crashes when a PSSM's
    // values nest deeper than that stack holds (NCBI BLAST+ 2.17.0 with 8 MiB: from a few thousand
    // levels, by the stack limit, the command line, the environment and the stack's random start).
    // LOSAT's reader and the object manager's step keep the values they are in on the heap, so
    // every depth is read, and the output is NCBI's with a stack large enough
    // (PD-LOSAT-NCBI-DEFECTS exception 15).
    let read = match guess_format(&bytes) {
        PssmFileFormat::TextAsn { start } => {
            crate::blastinput::pssm_text::read_pssm(&bytes[start..])
        }
        PssmFileFormat::BinaryAsn { start } => {
            crate::blastinput::pssm_text::read_pssm_binary(&bytes[start..])
        }
        PssmFileFormat::Xml { start } => {
            crate::blastinput::pssm_text::read_pssm_xml(&bytes[start..])
        }
        PssmFileFormat::Unsupported => {
            return Err(crate::blastinput::app::options_error(
                "Unsupported format for PSSM",
            ))
        }
        PssmFileFormat::Exception(exception) => {
            return Err(crate::cli::NativeError::exception(&exception).into())
        }
        PssmFileFormat::Refused { what, detail } => {
            return Err(reader.refuse(format!("{what} ({detail}; {})", path.display())))
        }
    };
    let program = reader.program.to_ascii_lowercase();
    let mut posts = crate::report::diag_log::Posts::new();
    for post in &read.posts {
        let severity = match post.severity() {
            Severity::Warning => crate::report::diag_log::WARNING,
            Severity::Error => crate::report::diag_log::ERROR,
        };
        posts.post(severity, post.line(Some(&program)));
    }
    posts.write(diagnostics, reader.archive)?;
    let pssm = match read.result {
        Ok(pssm) => pssm,
        Err(PssmReadError::Ncbi(ReadError::Serial(_))) => {
            return Err(crate::blastinput::app::options_error(&format!(
                "Unrecognized format for PSSM in {} (must be PssmWithParameters)",
                path.display()
            )))
        }
        Err(PssmReadError::Ncbi(ReadError::Exception(exception))) => {
            return Err(crate::cli::NativeError::exception(&exception).into())
        }
        Err(PssmReadError::Ncbi(ReadError::Std(what))) => {
            return Err(crate::cli::NativeError::std_exception(what).into())
        }
        Err(PssmReadError::Ncbi(ReadError::OutOfMemory)) => {
            return Err(crate::cli::NativeError::out_of_memory().into())
        }
        Err(PssmReadError::Ncbi(ReadError::Eof)) => {
            return Err(crate::cli::NativeError::silent_end().into())
        }
        // The readers' own limits (`read_pssm`): LOSAT's refusal, after NCBI's reader read
        // the file. (What NCBI reads and LOSAT's search does not take is refused after the object
        // manager's step, `psi_search_and_object`.)
        Err(PssmReadError::Refused(what))
        | Err(PssmReadError::Ncbi(ReadError::Unsupported(what))) => {
            return Err(reader.refuse(format!(
                "a PSSM that LOSAT's reader does not read ({}: {what})",
                path.display()
            )))
        }
    };
    psi_search_and_object(pssm, reader)
}

/// `CPsiBlastValidate::Pssm`'s error for `pssm`, which NCBI raises (`CBlastException`, "BLAST
/// engine error") when the search is made, after the outfmt 0 prolog.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:569-625
/// (The query checks between them are never reached from psiblast's application, which reads
/// the query of the PSSM before; LOSAT refuses those PSSMs, `psi_search_of`.)
fn validation_error(pssm: &PssmWithParameters) -> Option<String> {
    let scores = pssm.pssm.final_data.as_ref().map(|data| &data.scores);
    let missing_scores = scores.is_none_or(|scores| scores.is_empty());
    let missing_freq_ratios = pssm
        .pssm
        .intermediate_data
        .as_ref()
        .is_none_or(|data| data.freq_ratios.is_empty());
    if missing_freq_ratios && missing_scores {
        return Some("PSSM data must contain either scores or frequency ratios".to_string());
    }
    if !missing_scores {
        if let Some(factor) = pssm
            .pssm
            .final_data
            .as_ref()
            .and_then(|data| data.scaling_factor)
            .filter(|&factor| factor != 1)
        {
            return Some(format!(
                "PSSM has a scaling factor of {factor}. PSI-BLAST does not accept scaled PSSMs"
            ));
        }
    }
    if pssm.pssm.is_protein == Some(false) {
        return Some("PSSM does not represent protein scoring matrix".to_string());
    }
    None
}

/// The PSSM search of `pssm`: its query record, rows and values, or LOSAT's refusal.
#[cfg(test)]
pub(crate) fn psi_search_of(pssm: PssmWithParameters) -> Result<PsiSearch> {
    psi_search_and_object(pssm, PssmReader::PSIBLAST).map(|(search, _)| search)
}

/// The query's residues as NCBI's query set-up reads them (`query_residues`).
struct QueryResidues {
    /// One letter per residue: the NCBIstdaa letter of its code (`-` for the gap code 0), `X` for
    /// an invalid one (NCBI leaves it out of the query and stops, `invalid`, so no search reads it).
    letters: Vec<u8>,
    /// The 0-based positions of the residues whose code is no NCBIstdaa code.
    invalid: Vec<usize>,
    /// The Bioseq's letters as NCBI's reports read them where they are not `letters`
    /// (`PsiSearch::bioseq_letters`): iupacaa letters with one in lower case, as they are.
    bioseq: Option<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/bioseq_extract_data_priv.cpp:105-121
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1707 (`CIupacaaToStdaa`)
// and 1748-1781 (`CEaaToStdaa`): 0xff for every byte but these (rows 0x28-0x2f and 0x40-0x7f)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:867-902 and 894-904
/// The query's residues (`Seq-inst.seq-data`) as NCBI's PSSM search sets the query up
/// (`CBlastSeqVectorFromCSeq_data`, then `GetSequenceProtein`): each coding converts to NCBIstdaa
/// with its own table, and a code of 28 or more is an invalid residue. `iupacaa` knows the 26
/// letters of either case and neither `*` nor `-`; `ncbieaa` knows the 26 capital letters, `*`
/// (25) and `-` (the gap code 0, a residue the search reads); `ncbistdaa` is its own code. The
/// code of a letter is its capital's place in `NCBISTDAA_TO_AMINOACID` (the tables give every
/// letter that code). No residues when the data is missing or of another coding (refused before,
/// `psi_search_and_object`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:530-571
/// The reports read the Bioseq's letters through its `CSeqVector` in Iupac coding (alnvec.cpp:
/// 116-121, cited at `blastp::blast_engine::restore_displayed_o`): iupacaa data is copied as it
/// is, so a lower-case letter, which the search reads as its capital, is shown and compared in
/// lower case (`bioseq`); ncbieaa's and ncbistdaa's letters are the record's (`letters`: the
/// capitals, `*`, and `-` for the gap residue).
fn query_residues(seq_data: Option<&SeqData>) -> QueryResidues {
    let mut residues = QueryResidues {
        letters: Vec::new(),
        invalid: Vec::new(),
        bioseq: None,
    };
    // Each letter of iupacaa and ncbieaa that converts is its capital's place in
    // `NCBISTDAA_TO_AMINOACID`, which holds all 26 capitals: its letter is the capital itself.
    let push = |residues: &mut QueryResidues, position: usize, letter: Option<u8>| match letter {
        Some(letter) => residues.letters.push(letter),
        None => {
            residues.letters.push(b'X');
            residues.invalid.push(position);
        }
    };
    match seq_data {
        Some(SeqData::Ncbistdaa(codes)) => {
            residues.letters.reserve(codes.len());
            for (position, &code) in codes.iter().enumerate() {
                let letter = NCBISTDAA_TO_AMINOACID.get(usize::from(code)).copied();
                push(&mut residues, position, letter);
            }
        }
        Some(SeqData::Iupacaa(text)) => {
            residues.letters.reserve(text.len());
            for (position, &c) in text.iter().enumerate() {
                let letter = c.is_ascii_alphabetic().then(|| c.to_ascii_uppercase());
                push(&mut residues, position, letter);
            }
            if text.iter().any(u8::is_ascii_lowercase) {
                residues.bioseq = Some(text.clone());
            }
        }
        Some(SeqData::Ncbieaa(text)) => {
            residues.letters.reserve(text.len());
            for (position, &c) in text.iter().enumerate() {
                let letter = (c.is_ascii_uppercase() || c == b'*' || c == b'-').then_some(c);
                push(&mut residues, position, letter);
            }
        }
        Some(
            SeqData::Iupacna(_) | SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Other(_),
        )
        | None => {}
    }
    residues
}

/// The query record of the PSSM's Bioseq (`PsiSearch::query`, the reasons in
/// `psi_search_and_object`): one local id is the record's id (`CObject_id::AsString`), which an
/// integer keeps as the Bioseq's Seq-id for the Seq-aligns (`FastaRecord::keeping_bioseq_ids`:
/// `local id 5`, where the record's local id is a string); other ids are believed; a Bioseq
/// without one (only where the object manager stops the run first,
/// `stopped_as_query_joins_scope`) has no id.
fn query_record_of(ids: &[ParsedSeqId], title: &[u8], letters: &[u8]) -> FastaRecord {
    match ids {
        [] => FastaRecord::new(String::new(), title, letters),
        [ParsedSeqId::Local(ObjectId::Str(id))] => {
            FastaRecord::new(String::from_utf8_lossy(id).into_owned(), title, letters)
        }
        [ParsedSeqId::Local(ObjectId::Id(id))] => {
            FastaRecord::new(id.to_string(), title, letters).keeping_bioseq_ids(ids.to_vec())
        }
        _ => FastaRecord::new(String::new(), title, letters).believing(ids.to_vec()),
    }
}

/// The PSSM search of a PSSM whose query NCBI's object manager stops at as it joins the scope
/// (`AnnotOutcome::abort`, `psi_search_and_object`): the query record, the object manager's
/// messages and exception; no rows, values, request or PSSM files (NCBI reads none of them).
fn stopped_as_query_joins_scope(
    pssm: PssmWithParameters,
    annot: super::annot_index::AnnotOutcome,
    defline_reads_annot: bool,
) -> (PsiSearch, InPssmFile) {
    // A query Seq-entry set has no Bioseq of its own: an empty record.
    let empty = crate::report::serial::pssm::PssmQuery {
        ids: Vec::new(),
        descr: None,
        inst: crate::report::serial::pssm::SeqInst::new(("raw", 2), ("aa", 3), None, None),
        annot: None,
    };
    let query = pssm
        .pssm
        .query
        .as_ref()
        .and_then(SeqEntry::seq)
        .unwrap_or(&empty);
    let title = query.titles().first().copied().unwrap_or(&[]).to_vec();
    let search = PsiSearch {
        pssm: QueryPssm { rows: Vec::new() },
        freq_ratios: Vec::new(),
        stats: PssmStats::default(),
        query: query_record_of(
            &query.ids,
            &title,
            &query_residues(query.inst.seq_data.as_ref()).letters,
        ),
        // No report reads the query of a run that stops as its query joins the scope.
        bioseq_letters: None,
        query_ids: query.ids.clone(),
        validation_error: None,
        scores_error: None,
        setup_error: None,
        setup_messages: Vec::new(),
        warning_title: Some(Vec::new()),
        request: None,
        archive_pssm: None,
        archive_matrix_name_unwritable: false,
        annot,
        defline_reads_annot,
        round_query: None,
        split_setup: None,
        from_msa: false,
        query_other_coding: None,
    };
    let file = InPssmFile {
        object: pssm,
        computed: false,
        valid: false,
    };
    (search, file)
}

// The throw site NCBI BLAST+ 2.17.0 prints for `CSeq_entry::GetSeq()` of another variant: the
// `CInvalidChoiceSelection` of datatool-generated code (`Seq_entry_.cpp`, line 99, as its binary
// reports it; the file is not in the source tree, so this is a reference to the printed text only).
const SEQ_ENTRY_GET_SEQ: crate::ncbi_exception::ThrowSite = crate::ncbi_exception::ThrowSite {
    file: "src/objects/seqset/Seq_entry_.cpp",
    line: 99,
    kind: "CInvalidChoiceSelection::eFail",
    function: None,
};

/// The PSSM search of a PSSM that has no query Seq-entry (an `-in_pssm` PSSM without one, the
/// alignment's `-in_msa` PSSM whose query NCBI's `x_ExtractQueryForPssm` did not find): the
/// application adds the empty entry to its scope and `GetSeq()` of it ends the run after the
/// prolog, before any report.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537-540 (cited at
/// `psi_search_and_object`)
/// (NCBI BLAST+ 2.17.0: `CSeq_entry::GetSeq(): Invalid choice selection: NCBI-Seqset::Seq-entry.not
/// set`, exit status 255; PS6 corpus `e.stop`, `r.ch.q.star`.)
pub(crate) fn psi_search_unset_query(pssm: PssmWithParameters) -> (PsiSearch, InPssmFile) {
    let annot = super::annot_index::AnnotOutcome {
        abort: Some(super::annot_index::Abort::of(
            crate::cli::NativeError::exception(&NcbiException::new(
                SEQ_ENTRY_GET_SEQ,
                "CSeq_entry::GetSeq(): Invalid choice selection: NCBI-Seqset::Seq-entry.not set",
            )),
            false,
        )),
        ..super::annot_index::AnnotOutcome::default()
    };
    stopped_as_query_joins_scope(pssm, annot, false)
}

/// The PSSM search of `pssm` and the PSSM as NCBI's `CPsiBlastImpl` leaves it before the first
/// round: a PSSM of frequency ratios only gets the scores and gapped values that the PSSM engine
/// computes for it (`with_computed_scores`).
pub(crate) fn psi_search_and_object(
    pssm: PssmWithParameters,
    reader: PssmReader<'_>,
) -> Result<(PsiSearch, InPssmFile)> {
    psi_search_and_object_of(pssm, reader, IDS_OF_THE_PSSM)
}

/// What LOSAT's refusals of a query Seq-id (`believed::refused_id`) call the place the id comes
/// from, for a PSSM read from `-in_pssm`: the PSSM's query.
const IDS_OF_THE_PSSM: &str = "the PSSM";

/// The same for the PSSM that psiblast makes of `-in_msa`'s alignment (`psiblast::in_msa`): the
/// user gave no PSSM, and the query's ids are those of the master sequence's row (N1 of the PS6
/// review).
pub(crate) const IDS_OF_THE_MASTER_SEQUENCE: &str = "the alignment's master sequence";

/// [`psi_search_and_object`] of a PSSM whose query's Seq-ids the refusals name as those of
/// `ids_of` ([`IDS_OF_THE_PSSM`], [`IDS_OF_THE_MASTER_SEQUENCE`]).
pub(crate) fn psi_search_and_object_of(
    pssm: PssmWithParameters,
    reader: PssmReader<'_>,
    ids_of: &str,
) -> Result<(PsiSearch, InPssmFile)> {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:535-540
    // NCBI's application adds the query's Seq-entry to its scope and takes its Bioseq here: a
    // PSSM without one (no `query`, or `query NULL`) gets an empty Seq-entry, whose `GetSeq()`
    // ends the run as the alignment's PSSM without a query does (NCBI BLAST+ 2.17.0, psiblast and
    // tblastn, probe D3c-3b q_absent_ps, q_absent_tbn, t_query_null, q_null_tbn, b_query_null).
    let Some(entry) = pssm.pssm.query.as_ref() else {
        return Ok(psi_search_unset_query(pssm));
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537 (c++/src/app/blast/tblastn_app.cpp:226)
    // The query's Seq-entry joins NCBI's scope with its annotations, which NCBI's object manager
    // checks and indexes, and with its Bioseqs' ids (`annot_index::index_entry`: its messages, the
    // diagnostic stream it can leave failed, its exceptions), after the reader and the subjects and
    // before NCBI reads anything else of the PSSM; the run writes them where NCBI does
    // (`PsiSearch::write_annot_posts`). What NCBI does not do the same way on every run (memory
    // addresses, crashes, undefined behaviour) is refused here, before any output. With `-db` no
    // other sequence is in NCBI's scope before the query (`PssmReader::database`).
    let program = reader.program.to_ascii_lowercase();
    let mut annot = super::annot_index::index_entry(&program, entry, reader.database)
        .map_err(|refusal| reader.refuse(refusal.text()))?;
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:539-540 (c++/src/app/blast/tblastn_app.cpp:228)
    // A query Seq-entry set that the object manager took ends the run here: `GetSeq()` of a set is
    // NCBI's `CInvalidChoiceSelection` (NCBI BLAST+ 2.17.0, line RG D3-1 fix 4 `r4/probe` s06, s07,
    // s11: after the object manager's messages, psiblast after its prolog).
    let query = match entry {
        SeqEntry::Seq(query) => Some(query.as_ref()),
        SeqEntry::Set(_) => {
            if annot.abort.is_none() {
                annot.abort = Some(super::annot_index::Abort::of(
                    crate::cli::NativeError::exception(&NcbiException::new(
                        SEQ_ENTRY_GET_SEQ,
                        "CSeq_entry::GetSeq(): Invalid choice selection: NCBI-Seqset::Seq-entry.set",
                    )),
                    annot.muted,
                ));
            }
            None
        }
    };
    let defline_reads_annot = query.is_some_and(defline_reads_annot);
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537, 548 and
    // c++/src/app/blast/tblastn_app.cpp:226, 307, 324
    // The object manager's exception ends the run as the query joins the scope, before NCBI reads
    // any other value of the PSSM (its check, scores, sizes, statistics, letters and ids; its
    // search strategy's request and the PSSM files): NCBI BLAST+ 2.17.0 writes plain t04's stop
    // for every such PSSM (line RG D3-1 review 4 S2 and `r3/s2/`: a local id among other ids, the
    // sizes, a score outside one letter's range, the scores' sum, a negative lambda, no Seq-id, a
    // `BL_ORD_ID` or `giim` id, an id without an accession or a name, an empty local id, a
    // nucleotide Seq-inst, letters of another length, frequency ratios only, by row). LOSAT's
    // refusals of those values are for the runs that read them: also of a query Seq-entry set (it
    // ends at `GetSeq()`, above) and of Seq-data that LOSAT's search does not read (NCBI's object
    // manager does not read it, line RG D3-1 review 4b S1).
    let Some(query) = query.filter(|_| annot.abort.is_none()) else {
        return Ok(stopped_as_query_joins_scope(
            pssm,
            annot,
            defline_reads_annot,
        ));
    };
    if let Some(
        data
        @ (SeqData::Iupacna(_) | SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Other(_)),
    ) = &query.inst.seq_data
    {
        let (name, _, _) = data.variant();
        return Err(reader.refuse(format!(
            "a PSSM query whose Seq-data is {name} (LOSAT's PSSM search reads iupacaa, ncbieaa and ncbistdaa)"
        )));
    }
    let validation_error = validation_error(&pssm);
    // A PSSM of frequency ratios only: NCBI computes its scores with the PSSM engine
    // (`x_CreatePssmScoresFromFrequencyRatios`, below).
    let scores: &[i32] = match pssm.pssm.final_data.as_ref() {
        Some(data) if !data.scores.is_empty() => data.scores.as_slice(),
        _ => &[],
    };
    let ratios_only = validation_error.is_none() && scores.is_empty();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:466-473
    // The query's residues (`query_residues`) and NCBI's error of its invalid ones
    // (`report::query_warnings::invalid_residues_error`), which the query set-up
    // (`SetupQueries_OMF`) throws after NCBI's check of the PSSM: in `CPsiBlastImpl`'s constructor
    // for a PSSM of frequency ratios only, whose PSSM engine reads the query's sequence block (the
    // engine then never runs: `scores_error`), else where `CPsiBlast::Run` sets the search's query
    // up (`setup_error`, `PsiSearch::stop_before_search`). Either way NCBI reads none of the PSSM's
    // values but its check's: the refusals of those values below are for the runs that read them
    // (`reads_values`).
    let residues = query_residues(query.inst.seq_data.as_ref());
    let invalid_residues = validation_error
        .is_none()
        .then(|| crate::report::query_warnings::invalid_residues_error(&residues.invalid))
        .flatten();
    // A nucleotide Bioseq (`mol` dna, rna or na) is the query set-up's first error
    // (`SetupQueries_OMF` checks each query's molecule before its letters), also for a PSSM of
    // frequency ratios only, whose engine sets the query up first too.
    let nucleotide_query = matches!(query.inst.mol.1, 1 | 2 | 4);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:532-543
    // Then the set-up takes the query's title (`GetTitle`, before its letters), which reads every
    // descriptor's variant: an empty descriptor is the query's error
    // (`report::query_warnings::null_descriptor_error`; NCBI BLAST+ 2.17.0, psiblast and tblastn,
    // on every report path and with a PSSM of frequency ratios only, probe D3c-3b `probe4`).
    let query_error = if validation_error.is_some() {
        None
    } else if nucleotide_query {
        Some("Nucleotide Bioseq specified in program which expects protein query".to_string())
    } else if query.has_empty() {
        Some(crate::report::query_warnings::null_descriptor_error())
    } else {
        invalid_residues
    };
    let (mut scores_error, setup_error) = if ratios_only {
        (query_error, None)
    } else {
        (None, query_error)
    };
    let reads_values =
        validation_error.is_none() && scores_error.is_none() && setup_error.is_none();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:648-660
    // The engine's scores are column by column, which NCBI's search then reads by the PSSM's
    // own `byRow` (`PsiBlastComputePssmScores` copies the scores, not `byRow`): a ratio-only
    // PSSM written row by row would be searched with its scores transposed.
    if ratios_only && reads_values && pssm.pssm.by_row == Some(true) {
        return Err(reader.refuse(
            "a PSSM of frequency ratios only written by row (NCBI reads the scores it computes for it in the wrong order)",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537-540 (cited above)
    // The query's Seq-ids: NCBI's scope holds the PSSM's Bioseq with all of them, and the query's
    // Seq-loc (so its Seq-aligns, the PSSMs of the later rounds and its warnings) takes the first.
    // One local id (`local str "Query_N"` in a psiblast PSSM) is shown as a FASTA query's local id
    // is: the reports show the title alone, and with `-parse_deflines` the id
    // (`PsiSearch::query_record`); the record's local id is the Seq-id's `GetSeqIdString()`, the
    // string or the integer's decimal digits (`CObject_id::AsString`). Any other id (`gi 4504347`
    // of a `-parse_deflines` PSSM, `swissprot`, `pdb`, ...), and several ids, are shown with or
    // without `-parse_deflines`, as `-parse_deflines` shows the ids of a FASTA query (NCBI BLAST+
    // 2.17.0: `Query= gi|4504347 <title>`, `qseqid` `gi|4504347`; `GetSeqIdString` with `cbs` and
    // `believe_local_id` asks `believe_local_id` of a local id only, align_format_util.cpp).
    let ids = query.ids.as_slice();
    if ids.is_empty() {
        return Err(reader.refuse(
            "a PSSM query without a Seq-id (NCBI BLAST+ makes the query's Seq-loc of its first Seq-id)",
        ));
    }
    if let Some(refusal) =
        crate::blastinput::fasta_reader::believed::refused_id(ids, "query", ids_of, reader.program)
    {
        return Err(refusal);
    }
    // A local id among other ids: without `-parse_deflines` NCBI's reports show the title alone
    // and its tabular ids put the title's first word in the local id's place (`qseqid`
    // `gi|4504347|hemoglobin`, `s_ReplaceLocalId`), a mix that no FASTA query has (one without
    // `-parse_deflines` has its local id only); NCBI's psiblast writes the one Seq-id of the
    // query's Seq-aligns.
    if ids.len() > 1 && ids.iter().any(|id| matches!(id, ParsedSeqId::Local(_))) {
        return Err(reader.refuse(
            "a PSSM query with a local Seq-id among other Seq-ids (NCBI's reports show it as no FASTA query is shown)",
        ));
    }
    // A Textseq-id with neither an accession nor a name: NCBI's label of it is empty, and what NCBI
    // BLAST+ 2.17.0 does with it depends on the database and the report: `BLAST Database error:
    // Accessions to Oids lookup error in tfp` (exit 2, `-db tfp`; line RG Q1b's review, S1), the
    // pairwise report's `Failed to resolve SeqId: gb||` (exit 255, after its prolog), and `Unknown`
    // for `qseqid` in a tabular report (exit 0; line RG D3-1 review 4c `out_z/`), which LOSAT's
    // reports do not write.
    if ids.iter().any(|id| {
        id.textseq()
            .is_some_and(|text| text.accession.is_none() && text.name.is_none())
    }) {
        return Err(reader.refuse(
            "a PSSM query whose Seq-id has neither an accession nor a name (NCBI's run then depends on the database and the report: a database error, 'Failed to resolve SeqId' or the query 'Unknown')",
        ));
    }
    // A local id that is empty: NCBI's reports show it as no id. (A byte outside the visible range
    // of a string of the PSSM is NCBI's `#` with its message, made by the reader,
    // `report::serial::asn_types`.)
    if ids
        .iter()
        .any(|id| matches!(id, ParsedSeqId::Local(ObjectId::Str(name)) if name.is_empty()))
    {
        return Err(reader.refuse("a PSSM query whose Seq-id holds an empty local id"));
    }
    let titles = query.titles();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:669-683
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
    // NCBI's `Query=` and `# Query:` lines join every title of the Bioseq and trim the spaces at
    // the line's ends, `qseqid` reads the first title's first word (a title that starts with a
    // space has an empty first word, so `qseqid` is the local id), and the warnings read the first
    // title without its trailing `.` and spaces (`GetTitle`, below). A FASTA
    // defline gives one title without white space at its ends (psiblast's own PSSMs hold such a
    // title); the PSSM search reads that one only. These refusals concern the reports: a run that
    // NCBI stops before its search writes none of them (`PsiSearch::stop_before_search`, below).
    let title_refusal = if titles.len() > 1 {
        Some(reader.refuse(format!(
            "a PSSM query with {} titles (NCBI's Query= line joins them, its tabular query id and warnings read the first)",
            titles.len()
        )))
    } else {
        None
    };
    let title: &[u8] = titles.first().copied().unwrap_or(&[]);
    // `NStr::TruncateSpaces` (`eTrunc_Both`) of the `Query=` line would cut the title's ends.
    let title_refusal = title_refusal.or_else(|| {
        (crate::ncbi_str::truncate_spaces(title).len() != title.len()).then(|| {
            reader.refuse(
                "a PSSM query whose title starts or ends with white space (NCBI's Query= line, tabular query id and warnings each read such a title their own way)",
            )
        })
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/bioseq_extract_data_priv.cpp:344-372
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:115-118 (cited at
    // `CPsiBlastImpl::x_ExtractQueryFromPssm`: the PSSM's Bioseq goes through
    // `CObjMgrFree_QueryFactory`, whose query source is `CBlastQuerySourceBioseqSet`)
    // The warnings of a PSSM query (`blast_setup_cxx.cpp:533-543`, cited at
    // `report::query_warnings::warning_query_id`) take this title, in every round (the query
    // factory is fixed for the run), while the reports take the Bioseq's title as it is.
    let warning_title: Vec<u8> = if query.has_molinfo() {
        Vec::new()
    } else {
        let end = title
            .iter()
            .rposition(|&c| c != b'.' && c != b' ')
            .map_or(0, |last| last + 1);
        title[..end].to_vec()
    };
    // `other` (255): psiblast and tblastn search it as a protein (NCBI's query set-up refuses
    // only the nucleotide molecules; item PS15b-2 of line PS for tblastn).
    let mol_other = if query.inst.mol.1 == 255 {
        match query.inst.seq_data.as_ref() {
            // The coding that NCBI's object manager converts (`CSeq_data::E_Choice`).
            Some(SeqData::Iupacaa(_)) => Some(2),
            Some(SeqData::Ncbieaa(_)) => Some(8),
            Some(SeqData::Ncbistdaa(_)) => Some(10),
            _ => None,
        }
    } else {
        None
    };
    if !(matches!(query.inst.mol.1, 1..=4) || mol_other.is_some()) || query.inst.repr.1 != 2 {
        return Err(reader.refuse(format!(
            "a PSSM query whose Seq-inst is not a raw protein (repr {}, mol {})",
            query.inst.repr.0, query.inst.mol.0
        )));
    }
    // A gap letter (`-` of ncbieaa, the code 0 of ncbistdaa) is a residue that NCBI's search
    // reads (`query_residues`); the record keeps it as `-`, which the reports show and count as a
    // residue (`report::row_letters`). Lower-case iupacaa letters are searched as their capitals,
    // which the record holds; the reports read the Bioseq's own (`PsiSearch::bioseq_letters`).
    let letters = residues.letters;
    let length = query.inst.length.map_or(-1, i64::from);
    if letters.is_empty() || letters.len() as i64 != length {
        return Err(reader.refuse(
            "a PSSM query whose letters are missing or of another length than its Seq-inst",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:240-261
    // NCBI only asserts the sizes (its release build reads or writes past the matrix): LOSAT
    // refuses a PSSM whose sizes do not match its data, its query, or the 28 rows.
    let rows_n = pssm.pssm.num_rows;
    let columns_n = pssm.pssm.num_columns;
    let cells = i64::from(rows_n) * i64::from(columns_n);
    let ratios_len = pssm
        .pssm
        .intermediate_data
        .as_ref()
        .map_or(0, |data| data.freq_ratios.len());
    if reads_values
        && (!(1..=BLASTAA_SIZE as i32).contains(&rows_n)
            || i64::from(columns_n) != letters.len() as i64
            || (!ratios_only && scores.len() as i64 != cells)
            || (ratios_len != 0 && ratios_len as i64 != cells))
    {
        return Err(reader.refuse(format!(
            "a PSSM whose dimensions ({rows_n} rows, {columns_n} columns) do not match its scores ({}), frequency ratios ({ratios_len}) or query ({} letters)",
            scores.len(),
            letters.len()
        )));
    }
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:121-122
    // NCBI copies the scores into the matrix without a check and adds them in `Int4` (a score of
    // 40,000 runs; one near `INT4_MAX` would overflow, undefined behaviour). LOSAT keeps the scores
    // to one letter's range and the sums of a search below 2^30 (below): range rules of LOSAT's,
    // not NCBI errors.
    if reads_values {
        if let Some(score) = scores
            .iter()
            .find(|&&score| !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&score))
        {
            return Err(reader.refuse(format!(
                "a PSSM score of {score}, outside BLAST_SCORE_MIN..BLAST_SCORE_MAX ({BLAST_SCORE_MIN}..{BLAST_SCORE_MAX}, one letter's score range in NCBI BLAST+),"
            )));
        }
    }
    let by_row = pssm.pssm.by_row == Some(true);
    let mut rows = if reads_values {
        convert_2_matrix(
            scores,
            rows_n as usize,
            columns_n as usize,
            by_row,
            BLAST_SCORE_MIN,
        )
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:58
    // An alignment adds each column's score at most once, so its score and the ungapped and gapped
    // sums before it stay below the sum of the columns' largest positive scores. LOSAT refuses a
    // PSSM whose sum reaches 2^30 (where NCBI's `Int4` sums can reach `MININT`'s distance from 0 and
    // overflow, undefined behaviour): every PSSM of fewer than 19,800 columns passes, and a PSSM of
    // NCBI's engine (scores of a few tens) up to millions of columns.
    if reads_values {
        let sum: i64 = rows
            .iter()
            .map(|column| i64::from(column.iter().copied().max().unwrap_or(0).max(0)))
            .sum();
        if sum >= PSSM_SCORE_SUM_LIMIT {
            return Err(reader.refuse(format!(
                "a PSSM whose columns' largest scores add up to {sum} (2^30 or more, where NCBI BLAST+'s Int4 sums can overflow),"
            )));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:264-284
    // The frequency ratios (`psi_matrix->freq_ratios`, 28 per position, 0.0 past the file's
    // rows), which the position-based composition-based statistics read; none when the PSSM has
    // none (NCBI's then stay zero).
    let freq_ratios = match pssm.pssm.intermediate_data.as_ref() {
        Some(data) if reads_values && !data.freq_ratios.is_empty() => convert_2_matrix(
            &data.freq_ratios,
            rows_n as usize,
            columns_n as usize,
            by_row,
            0.0,
        ),
        _ => Vec::new(),
    };
    let final_data = pssm.pssm.final_data.as_ref();
    let stat = |value: Option<f64>| value.filter(|&value| value != 0.0);
    let mut stats = PssmStats {
        lambda_ungapped: stat(final_data.and_then(|data| data.lambda_ungapped)),
        kappa_ungapped: stat(final_data.and_then(|data| data.kappa_ungapped)),
        h_ungapped: stat(final_data.and_then(|data| data.h_ungapped)),
        lambda: stat(final_data.map(|data| data.lambda)),
        kappa: stat(final_data.map(|data| data.kappa)),
        h: stat(final_data.map(|data| data.h)),
    };
    let mut computed = None;
    if ratios_only && reads_values {
        match scores_from_freq_ratios(&letters, &freq_ratios, reader)? {
            Ok((scores, gapped)) => {
                rows = scores;
                stats.lambda = stat(Some(gapped.lambda));
                stats.kappa = stat(Some(gapped.k));
                stats.h = stat(Some(gapped.h));
                computed = Some(gapped);
            }
            // The engine's error, raised by `CPsiBlastImpl`'s constructor after the outfmt 0
            // prolog and after the search strategy's request (`PsiSearch::stop_before_search`).
            Err(engine_error) => scores_error = Some(engine_error),
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:110-112
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:43-50
    // NCBI takes a negative value of the file into `kbp_psi` or `kbp_gap_psi` (0 is
    // `kInvalidStat`: the query's BLOSUM62 value). A `kbp_psi` that is not valid stops the search's
    // setup with NCBI's error of a query without Karlin-Altschul parameters
    // (`blastp::blast_engine::psiblast_stop_at_invalid_query`, line PS DECISIONS 30). A
    // `kbp_gap_psi` that is not valid does not stop NCBI: it runs with it (no gapped X-drop, a
    // cutoff of 1, E-values of -1), an output without meaning, refused as a PSSM whose scaling gives
    // no Karlin-Altschul parameters is (line PS DECISIONS 28), unless the setup stops first. NCBI's
    // check of the PSSM's file (`validation_error`) comes first.
    //
    // A value that is not a finite number (ASN.1 `PLUS-INFINITY`, `MINUS-INFINITY`,
    // `NOT-A-NUMBER`, in text or binary) NCBI takes in the same way: an ungapped one that is not
    // above 0 (`NOT-A-NUMBER` too) makes `kbp_psi` invalid, as `PsiSearch::ungapped_karlin_valid`
    // holds it. Where NCBI then goes on, LOSAT's search does not yet compute what NCBI computes
    // from a gapped `lambda` of `PLUS-INFINITY`, a gapped `NOT-A-NUMBER` or psiblast's
    // `lambdaUngapped` of `PLUS-INFINITY` (line RG D3c-3a, `D/d3/c3a/probe/realtext`: NCBI's
    // E-values 0.0 and bit scores `inf` or -1.0, its hits; handed to line PS,
    // `D/d3/c3a/HANDOFF_PS_nonfinite_lambda.md`): refused with the negative values. Gapped `kappa`
    // and `h` of `PLUS-INFINITY` and the other ungapped values of `PLUS-INFINITY` give NCBI's
    // output; tblastn reads `lambdaUngapped` only for its word X-drop, where it gives NCBI's too.
    let ungapped_invalid = [
        stats.lambda_ungapped,
        stats.kappa_ungapped,
        stats.h_ungapped,
    ]
    .into_iter()
    .flatten()
    .any(|value| value.is_nan() || value <= 0.0);
    let unusable = |name: &'static str, value: Option<f64>| match value {
        Some(value) if value < 0.0 => Some((name, "negative")),
        Some(value) if value.is_nan() => Some((name, "not a number")),
        Some(value) if value == f64::INFINITY && name == "lambda" => Some((name, "infinite")),
        _ => None,
    };
    let gapped_unusable = unusable("lambda", stats.lambda)
        .or_else(|| unusable("kappa", stats.kappa))
        .or_else(|| unusable("h", stats.h));
    let ungapped_unusable = (reader.program == PssmReader::PSIBLAST.program
        && stats.lambda_ungapped == Some(f64::INFINITY))
    .then_some(("lambdaUngapped", "infinite"));
    if validation_error.is_none()
        && scores_error.is_none()
        && setup_error.is_none()
        && !ungapped_invalid
    {
        match (gapped_unusable, ungapped_unusable) {
            (Some((name, "negative")), _) => {
                return Err(reader.refuse(format!(
                    "a PSSM whose gapped Karlin-Altschul parameter {name} is negative (NCBI BLAST+ 2.17.0 goes on with a Karlin-Altschul block that it holds invalid)"
                )));
            }
            (Some((name, what)), _) | (None, Some((name, what))) => {
                let kind = if name == "lambdaUngapped" {
                    "ungapped"
                } else {
                    "gapped"
                };
                return Err(reader.refuse(format!(
                    "a PSSM whose {kind} Karlin-Altschul parameter {name} is {what} (LOSAT does not yet compute the hits, E-values and bit scores that NCBI BLAST+ 2.17.0 computes from it)"
                )));
            }
            (None, None) => {}
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:533-557
    // The search strategy's request holds the object that the application read, which
    // `SaveSearchStrategy` serializes before `CPsiBlastImpl` computes the scores of a PSSM of
    // frequency ratios only (psiblast_app.cpp:340, tblastn_app.cpp:307): the object as read, kept
    // for a request that NCBI writes (a PSSM that NCBI's check refuses writes none; one whose
    // scores the engine cannot compute is written first), with every member it read (the PSSM
    // writer writes them all, `report::serial::pssm`).
    let request = (super::pssm_request_written(reader.strategy, reader.database)
        && validation_error.is_none())
    .then(|| PssmQueryRequest {
        pssm: std::sync::Arc::new(pssm.clone()),
    });
    // NCBI reference (598d8ae6): c++/src/serial/delaybuf.cpp:63-84
    // That ASN.1 text write reads each delay buffer of another format (`CDelayBuffer::Update`):
    // one whose read fails stops NCBI with its exception after part of the strategy file (rc
    // 255; line RG-F3, PR 2: probe AR5 SO), refused before any output as `-out_pssm` refuses it
    // (`algorithm::psiblast::in_pssm_files`).
    if let Some(request) = &request {
        let text = crate::report::archive::ArchiveEncoding::Text;
        // (A list element stated as NULL is written: the strategy's writer throws where NCBI's
        // stream does, `report::archive::PssmGap::NullElement`.)
        if crate::report::archive::pssm_unwritten(text, &request.pssm)
            == Some(crate::report::archive::PssmGap::ReadAgain)
        {
            return Err(reader.refuse(
                "-export_search_strategy with a PSSM whose author names (a member read into a delay buffer of another format) NCBI reads again with an error (NCBI BLAST+ writes part of the file and stops with a C++ exception)",
            ));
        }
    }
    let search = PsiSearch {
        pssm: QueryPssm { rows },
        freq_ratios,
        stats,
        query: query_record_of(ids, title, &letters),
        bioseq_letters: residues.bioseq.map(std::sync::Arc::from),
        query_ids: ids.to_vec(),
        validation_error,
        scores_error,
        setup_error,
        setup_messages: Vec::new(),
        warning_title: Some(warning_title),
        request,
        archive_pssm: None,
        archive_matrix_name_unwritable: false,
        annot,
        defline_reads_annot,
        round_query: None,
        split_setup: None,
        from_msa: false,
        query_other_coding: mol_other,
    };
    // The title's refusals (above), for a run whose reports NCBI writes (line RG F1b-2, review 4
    // R5-6: with `-export_search_strategy -subject`, an invalid PSSM or the PSSM engine's error, NCBI
    // writes the prolog and its error whatever the title; so it does when its object manager's
    // exception ends the run as the query joins the scope, `PsiSearch::annot_abort`, line RG D3-1).
    if let Some(refusal) = title_refusal {
        if search
            .stop_before_search(reader.strategy, reader.database)
            .is_none()
            && search.annot.abort.is_none()
        {
            return Err(refusal);
        }
    }
    let mut pssm = pssm;
    if let Some(gapped) = computed {
        with_computed_scores(&mut pssm, &search.pssm.rows, gapped);
    }
    let file = InPssmFile {
        object: pssm,
        computed: computed.is_some(),
        valid: search.validation_error.is_none() && search.scores_error.is_none(),
    };
    Ok((search, file))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:559 and 705,
// blastxml2_format.cpp:590 (the query definition of the XML and JSON reports:
// `sequence::CDeflineGenerator` on the query's Bioseq)
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1977-2077 (`x_SetTitleFromProtein`:
// a protein without a title takes the name of its `prot` feature from the scope's annotations;
// `x_SetBioSrc`, :1410, reads its source features)
/// Whether NCBI's definition line of the PSSM's query (the XML and JSON reports' query
/// definition) reads its annotations: a query without a title whose annotations hold a feature (a
/// feature table, or a feature Seq-table, whose rows the lookup makes into features). NCBI BLAST+
/// 2.17.0 then names the query after its `prot` feature (line RG D3-1 `om/` cases x01-x04, outfmt 5
/// and 16); LOSAT does not port that lookup (`om/SPEC.md` §1.3) and refuses those reports.
fn defline_reads_annot(query: &crate::report::serial::pssm::PssmQuery) -> bool {
    use crate::report::serial::asn_types::{member_index, members, table::named, Value};
    if !query.titles().is_empty() {
        return false;
    }
    let data = member_index(named::SEQ_ANNOT, "data");
    let kinds = members(members(named::SEQ_ANNOT)[data].ty);
    query
        .annot
        .iter()
        .flatten()
        .any(|annot| match annot.member(data) {
            Some(Value::Choice(variant, value)) => match kinds[usize::from(*variant)].name() {
                "ftable" => matches!(value.as_ref(), Value::List(feats) if !feats.is_empty()),
                "seq-table" => {
                    let feat_type = value
                        .member(member_index(named::SEQ_TABLE, "feat-type"))
                        .map_or(0, |value| match value {
                            Value::Int(n) => *n,
                            _ => 0,
                        });
                    (1..23).contains(&feat_type)
                }
                _ => false,
            },
            _ => false,
        })
}

/// The refusal of a report whose query definition NCBI makes from the annotations of a PSSM's
/// query (`PsiSearch::defline_reads_annot`): XML and JSON (`outfmt` 5, 13, 14, 15, 16).
pub(crate) fn refuse_annot_defline(search: &PsiSearch, outfmt: i32, program: &str) -> Result<()> {
    if search.defline_reads_annot && matches!(outfmt, 5 | 13..=16) {
        bail!(
            "output format {outfmt} with a PSSM query without a title whose annot holds features (NCBI BLAST+ names the query in XML and JSON from those features, which LOSAT does not yet) is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The scores of a PSSM of frequency ratios only, and the gapped values of its scaling:
/// `PsiBlastComputePssmScores` runs the PSSM engine on the ratios (padded to 28 rows) of the
/// query's letters as the search holds them (O as X), with the options' matrix (`reader`:
/// BLOSUM62 for psiblast, which runs no other matrix with `-in_pssm`) and the gap costs 11/1 (the
/// defaults of `CPsiBlastInputFreqRatios`), and keeps the scores and `lambda`, `kappa`, `h`; the
/// ungapped values stay absent (`kbp_psi` then takes `kbp_std`'s, `PsiSearch::ungapped_karlin`).
/// A matrix whose table has no 11/1 row stops the engine with NCBI's error (the inner `Err`,
/// `CBlastException` "BLAST engine error" when the search is made); a matrix NCBI does not know
/// is the options' error that follows the reading (`validate_protein_options`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:466-485
fn scores_from_freq_ratios(
    letters: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    reader: PssmReader<'_>,
) -> Result<std::result::Result<(Vec<[i32; BLASTAA_SIZE]>, crate::core::psi::KarlinBlk), String>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:979-985
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3531-3561
    // The engine's score block takes the matrix with the gap costs 11/1.
    let matrix_name = reader.matrix_name.to_ascii_uppercase();
    let Ok(matrix) = matrix_name.parse::<crate::config::ScoringMatrix>() else {
        return Ok(Err(format!("{matrix_name} is not a supported matrix")));
    };
    if crate::stats::protein_options::karlin_blk_gapped_load_from_tables(11, 1, &matrix_name, false)
        .is_err()
    {
        return Ok(Err(format!(
            "Gap existence and extension values of 11 and 1 not supported for {matrix_name}"
        )));
    }
    let query: Vec<u8> = letters
        .iter()
        .map(
            |&letter| match crate::utils::matrix::aa_char_to_ncbistdaa(letter) {
                crate::utils::matrix::ncbistdaa::O => crate::utils::matrix::ncbistdaa::X,
                code => code,
            },
        )
        .collect();
    let spec = crate::config::ProteinScoringSpec {
        matrix,
        gap_open: 11,
        gap_extend: 1,
    };
    let gapped = crate::stats::tables::lookup_protein_params(&spec);
    let sbp = score_block(
        spec.matrix,
        crate::core::psi::KarlinBlk {
            lambda: gapped.lambda,
            k: gapped.k,
            h: gapped.h,
        },
    )?;
    crate::core::psi::pssm_from_freq_ratios(&query, freq_ratios, &sbp)
        .map(Ok)
        .map_err(|_| {
            reader.refuse("a PSSM of frequency ratios from which NCBI's PSSM engine makes no scaled scores (NCBI BLAST+ 2.17.0 does not finish or stops with an error)")
        })
}

/// `GetScores` or `GetFreqRatios` + the copy of `PsiBlastSetupScoreBlock`:
/// `rows[position][residue]` of 28 cells, those past the file's rows `fill` (`BLAST_SCORE_MIN`
/// for the scores, 0.0 for the frequency ratios).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:215-237
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:140-144
pub(crate) fn convert_2_matrix<T: Copy>(
    values: &[T],
    num_rows: usize,
    num_columns: usize,
    by_row: bool,
    fill: T,
) -> Vec<[T; BLASTAA_SIZE]> {
    let mut rows: Vec<[T; BLASTAA_SIZE]> = (0..num_columns).map(|_| [fill; BLASTAA_SIZE]).collect();
    if by_row {
        let mut values = values.iter().copied();
        for r in 0..num_rows {
            for row in rows.iter_mut() {
                row[r] = values.next().unwrap_or(fill);
            }
        }
    } else if num_rows > 0 {
        // Column by column (NCBI's engine writes `byRow` false): each column's values are the next
        // `num_rows` of the file's, copied as one slice; the cells past the values keep `fill`.
        let column_rows = num_rows.min(BLASTAA_SIZE);
        for (row, column) in rows.iter_mut().zip(values.chunks(column_rows)) {
            row[..column.len()].copy_from_slice(column);
        }
    }
    rows
}

/// The PSSM file of `-in_pssm` as `read_in_pssm` reads it: the
/// object after `CPsiBlastImpl`'s constructor, whether the scores
/// were computed from the frequency ratios, and whether NCBI's constructor of the search passes
/// (`CPsiBlastValidate::Pssm` and the PSSM engine's scores of a PSSM of frequency ratios only; it
/// throws before any PSSM file is written otherwise).
pub(crate) struct InPssmFile {
    pub(crate) object: PssmWithParameters,
    pub(crate) computed: bool,
    pub(crate) valid: bool,
}

impl InPssmFile {
    /// Whether NCBI cannot write the object's parameters: a PSSM of frequency ratios only without
    /// `params.rpsdbparams` gets those parameters from `CPsiBlastImpl`'s constructor (the gap
    /// costs, `PsiBlastAddAncillaryPssmData`) without the matrix name, and the object stream
    /// stops with "cannot write unassigned member matrixName" (a C++ exception, rc 255) wherever
    /// it writes the object: `-out_pssm` (`algorithm::psiblast::in_pssm_files`) and the archive's
    /// request (`-outfmt 11`, `algorithm::psiblast::in_pssm`). The one owner of that condition.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:426-433
    pub(crate) fn matrix_name_unwritable(&self) -> bool {
        self.computed
            && self
                .object
                .params
                .as_ref()
                .is_none_or(|params| params.rpsdbparams.is_none())
    }
}

/// The scores and gapped values that `CPsiBlastImpl`'s constructor sets in a PSSM of frequency
/// ratios only (`rows[position][residue]` and `gapped`, the PSSM engine's), with the frequency
/// ratios padded to 28 rows when the file has fewer, and the gap costs of the search (BLOSUM62's
/// 11 and 1: LOSAT's psiblast runs no other, `scores_from_freq_ratios`) in its
/// parameters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:486-503
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:440-463
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:426-433
/// The engine's scores are column by column (`x_PSIMatrix2Asn1`, cited at
/// `psi_search_and_object`). A PSSM of frequency ratios written by row is LOSAT's refusal
/// (`psi_search_and_object`), so the ratios are padded column by column. A file without `rpsdbparams` gets
/// parameters without a matrix name, which NCBI cannot write: refused with `-out_pssm`
/// (`algorithm::psiblast::in_pssm_files::rounds_input`), and the object is left without them.
pub(crate) fn with_computed_scores(
    pssm: &mut PssmWithParameters,
    rows: &[[i32; BLASTAA_SIZE]],
    gapped: KarlinBlk,
) {
    let matrix = &mut pssm.pssm;
    let rows_n = usize::try_from(matrix.num_rows).unwrap_or(0);
    if rows_n < BLASTAA_SIZE {
        if let Some(data) = matrix.intermediate_data.as_mut() {
            let mut padded = Vec::with_capacity(rows.len() * BLASTAA_SIZE);
            for column in data.freq_ratios.chunks(rows_n.max(1)) {
                padded.extend_from_slice(column);
                padded.resize(padded.len() + BLASTAA_SIZE - rows_n, 0.0);
            }
            data.freq_ratios = padded;
        }
        matrix.num_rows = BLASTAA_SIZE as i32;
    }
    let scores: Vec<i32> = rows.iter().flatten().copied().collect();
    match matrix.final_data.as_mut() {
        Some(data) => {
            data.scores = scores;
            data.lambda = gapped.lambda;
            data.kappa = gapped.k;
            data.h = gapped.h;
        }
        None => {
            matrix.final_data = Some(PssmFinalData {
                scores,
                lambda: gapped.lambda,
                kappa: gapped.k,
                h: gapped.h,
                scaling_factor: None,
                lambda_ungapped: None,
                kappa_ungapped: None,
                h_ungapped: None,
                word_score_threshold: None,
            })
        }
    }
    if let Some(rps) = pssm
        .params
        .as_mut()
        .and_then(|params| params.rpsdbparams.as_mut())
    {
        rps.gap_open = Some(11);
        rps.gap_extend = Some(1);
    }
}

/// The matrix as the score block holds it, with its smallest and largest scores
/// (`sbp->loscore` and `sbp->hiscore`, `BlastScoreBlkMatrixFill`).
pub(crate) fn score_table(
    matrix: ScoringMatrix,
) -> ([[i32; BLASTAA_SIZE]; BLASTAA_SIZE], i32, i32) {
    let mut table = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];
    let (mut loscore, mut hiscore) = (i32::MAX, i32::MIN);
    for (i, row) in table.iter_mut().enumerate() {
        for (j, score) in row.iter_mut().enumerate() {
            *score = crate::utils::matrix::protein_score(matrix, i as u8, j as u8);
            if *score <= BLAST_SCORE_MIN || *score >= i16::MAX as i32 {
                continue;
            }
            loscore = loscore.min(*score);
            hiscore = hiscore.max(*score);
        }
    }
    (table, loscore, hiscore)
}

/// The score block of `CPssmEngine::x_InitializeScoreBlock` that the engine reads.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:521-578
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2843-2848
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1506-1528 (the matrix's
// `loscore` and `hiscore`: its smallest and largest scores between BLAST_SCORE_MIN and
// BLAST_SCORE_MAX)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/matrix_freq_ratios.c:1661-1736
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:104-120
// `gapped` is the gapped block of the matrix and gap costs (`Blast_KarlinBlkGappedCalc`), the
// search's: `kbp_gap_std[0]`, copied to `kbp_gap_psi[0]`.
pub(crate) fn score_block(matrix: ScoringMatrix, gapped: KarlinBlk) -> Result<PsiScoreBlk> {
    // psiblast's `Validate` refuses IDENTITY, the one matrix without frequency ratios
    // (`stats::protein_options`, the standard matrices only).
    let data = crate::core::composition_adjustment::adjust_scores::compo_matrix_data(matrix)
        .ok_or_else(|| anyhow::anyhow!("no frequency ratios for {}", matrix.ncbi_name()))?;
    let (table, loscore, hiscore) = score_table(matrix);
    let std_comp = crate::stats::karlin_calc::compute_std_aa_composition();
    let sfp = crate::stats::karlin_calc::compute_score_freq_profile_for_matrix(
        &std_comp, &std_comp, loscore, hiscore, matrix,
    );
    let ideal = crate::stats::karlin_calc::compute_karlin_params_ungapped(&sfp)
        .map_err(|error| crate::blastinput::app::engine_error(&error))?;
    Ok(PsiScoreBlk {
        matrix: table,
        ideal_lambda: ideal.lambda,
        ideal_k: ideal.k,
        gap_std_k: gapped.k,
        gap_psi: gapped,
        freq_ratios: data.freq_ratios,
        bit_scale_factor: if matrix == ScoringMatrix::Blosum45 {
            3
        } else {
            2
        },
        background: data.background,
    })
}

/// The score block of `CPssmEngine::x_InitializeScoreBlock` for the matrix named `matrix_name`
/// (`GetMatrixName()`), the gap costs and the query (NCBIstdaa), psiblast's `-in_msa` and
/// makeprofiledb's PSSMs of frequency ratios only: the inner `Err` is the message of NCBI's
/// `CBlastException` for a matrix it does not know and for gap costs that the matrix's table does
/// not hold; the block holds the gapped block
/// `kbp_gap_psi` copies. `None` for IDENTITY, whose table has no frequency ratios (the engine
/// stops there).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:498-520
/// The ungapped calculation writes its message for a query whose parameters it cannot compute
/// and its status is then replaced by the gapped calculation's:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2783-2790
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3557-3562
/// The engine throws the first message that was written, so with gap costs that the matrix does
/// not have the error is the ungapped calculation's for such a query and the gap costs' for any
/// other (whose allowed values follow as further messages); with gap costs that the matrix has,
/// the status is 0 and the engine goes on without the query's block (`ungapped_karlin_valid`).
/// NCBI BLAST+ 2.17.0, PS6 corpus `b.gap.bad`, `o.ok.matrix`; review probes `k.X20.badgap`,
/// `z.x20.go0.ge5`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:566-570
#[allow(clippy::type_complexity)]
pub(crate) fn engine_score_block(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    query: &[u8],
) -> Result<std::result::Result<(ScoringMatrix, Option<PsiScoreBlk>), String>> {
    let Ok(matrix) = matrix_name.parse::<ScoringMatrix>() else {
        return Ok(Err("Unknown error code -1".to_string()));
    };
    let spec = crate::config::ProteinScoringSpec {
        matrix,
        gap_open,
        gap_extend,
    };
    if !crate::stats::tables::protein_scoring_supported(&spec) {
        if query_ungapped_block(matrix, query).is_none() {
            return Ok(Err(
                crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()
            ));
        }
        return Ok(Err(format!(
            "Gap existence and extension values of {} and {} not supported for {}",
            spec.gap_open,
            spec.gap_extend,
            matrix_name.to_ascii_uppercase()
        )));
    }
    let params = crate::stats::tables::lookup_protein_params(&spec);
    let gapped = KarlinBlk {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
    };
    let sbp = if matrix == ScoringMatrix::Identity {
        None
    } else {
        Some(score_block(matrix, gapped)?)
    };
    Ok(Ok((matrix, sbp)))
}

/// The query's ungapped Karlin-Altschul block in the engine's score block of `matrix`
/// (`kbp_psi[0]`, `Blast_ScoreBlkKbpUngappedCalc` of the query's composition), `None` where it
/// cannot be computed: NCBI's score block leaves the query's context invalid there (a query of
/// `X` alone, or of letters that no positive score pairs with the standard composition), and the
/// engine then reads a null block.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2786
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:514-520 (the ungapped
/// calculation's status is replaced by the gapped one's, so the failure goes unreported and the
/// run goes on with `kbp_psi[0]` null; `s_PSISavePssm` and the scaling read it: NCBI BLAST+ 2.17.0
/// ends with SIGSEGV, PS6 corpus `r.ch.q.allX`, `x.in.allx.60`, `x.in.allO`)
pub(crate) fn query_ungapped_block(matrix: ScoringMatrix, query: &[u8]) -> Option<KarlinBlk> {
    let (_, loscore, hiscore) = score_table(matrix);
    let mut sentinelled = Vec::with_capacity(query.len() + 2);
    sentinelled.push(0);
    sentinelled.extend_from_slice(query);
    sentinelled.push(0);
    let composition = crate::stats::karlin_calc::compute_aa_composition(&sentinelled, query.len());
    let standard = crate::stats::karlin_calc::compute_std_aa_composition();
    let profile = crate::stats::karlin_calc::compute_score_freq_profile_for_matrix(
        &composition,
        &standard,
        loscore,
        hiscore,
        matrix,
    );
    crate::stats::karlin_calc::compute_karlin_params_ungapped(&profile)
        .ok()
        .map(|params| KarlinBlk {
            lambda: params.lambda,
            k: params.k,
            h: params.h,
        })
}

#[cfg(test)]
#[path = "in_pssm_tests.rs"]
mod tests;

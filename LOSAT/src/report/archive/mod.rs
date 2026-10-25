//! The BLAST archive (outfmt 11, `Blast4-archive`), shared by the programs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1051-1064
//! ```c++
//! void CBlastFormat::PrintArchive(CRef<objects::CBlast4_archive> archive, CNcbiOstream& out)
//! {
//!     if (archive.Empty()) {
//!         return;
//!     }
//!     string outfmt = CNcbiEnvironment().Get("ARCHIVE_FORMAT");
//!     if (outfmt.empty()) {
//!         out << MSerial_AsnText << *archive;
//!     } else if (!NStr::CompareNocase(outfmt, "xml")) {
//!         out << MSerial_Xml << *archive;
//!     } else if (NStr::StartsWith(outfmt, "bin", NStr::eNocase)) {
//!         out << MSerial_AsnBinary << *archive;
//!     }
//! }
//! ```
//! One `Blast4-archive` per query batch is written after the batch's search, one error
//! archive when the application stops with a diagnostic (`PrintErrorArchive`). Each archive
//! is written straight to the sink, element by element (`writer`), from data the program
//! holds for its other formats.
//!
//! - `options`: the request's `algorithm-options` (the fold of NCBI's `Set*` calls);
//! - `bioseq`: the Bioseq of a query or subject record;
//! - `masks`: the `Blast4-mask` lists (results and the `LCaseMask` options);
//! - `writer`: the document.
//!
//! What a run keeps from one archive to the next is `ArchiveState`: the diagnostics saved
//! since the last archive (`DiagLog`) and the last `EffectiveSearchSpace`.

pub mod bioseq;
pub mod masks;
pub mod options;
pub mod writer;

use std::borrow::Cow;
use std::ffi::OsStr;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::diag_log::{DiagLog, Message, ERROR};
use crate::stats::tables::KarlinParams;
use options::OptionInputs;

pub use writer::{
    write_archive, write_error_archive, AlignSource, ArchiveOut, BatchArchive, KaBlock,
    NoAlignments,
};

/// The Karlin-Altschul block of the results (`Blast4-ka-block`) that `params` makes.
pub(crate) fn ka_block(params: &KarlinParams, gapped: bool) -> KaBlock {
    KaBlock {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
        gapped,
    }
}

/// The defline line of the Bioseq: the line as read, else `>` and the title.
pub(crate) fn defline_of(record: &FastaRecord) -> Cow<'_, [u8]> {
    match record.archive_defline() {
        Some(defline) => Cow::Borrowed(defline),
        None => {
            let mut line = vec![b'>'];
            line.extend_from_slice(&record.title);
            Cow::Owned(line)
        }
    }
}

/// How `ARCHIVE_FORMAT` makes NCBI write the archive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArchiveEncoding {
    /// ASN.1 text (the variable is unset or empty).
    Text,
    /// ASN.1 binary (a value that starts with `bin`, in any case).
    Binary,
    /// XML (the value `xml` in any case): NCBI's XML object stream (`XmlEncoder::archive`).
    Xml,
    /// Any other value: NCBI writes nothing.
    Nothing,
}

impl ArchiveEncoding {
    /// `PrintArchive`'s reading of `ARCHIVE_FORMAT`.
    pub fn of(value: Option<&OsStr>) -> Self {
        let Some(value) = value else {
            return ArchiveEncoding::Text;
        };
        let value = value.to_string_lossy();
        if value.is_empty() {
            ArchiveEncoding::Text
        } else if value.eq_ignore_ascii_case("xml") {
            ArchiveEncoding::Xml
        } else if value.len() >= 3 && value.as_bytes()[..3].eq_ignore_ascii_case(b"bin") {
            ArchiveEncoding::Binary
        } else {
            ArchiveEncoding::Nothing
        }
    }

    /// The encoding the environment asks for.
    pub fn from_env() -> Self {
        Self::of(std::env::var_os("ARCHIVE_FORMAT").as_deref())
    }

    /// Whether this encoding writes every byte of `text` into a VisibleString of the archive
    /// as NCBI does, without a diagnostic.
    ///
    /// Text and binary: only `0x20..=0x7E` (`GoodVisibleChar`); NCBI writes `#` for any other
    /// byte and posts a `Bad char` diagnostic for each, whose text and place LOSAT does not
    /// reproduce. XML: the archive's stream has no encoding, so every byte is written as it
    /// is (TAB, LF and CR as character references) except `0x01..=0x1F` without those three,
    /// which `x_VerifyChar` (`BAD_CHAR`) replaces by `#` with a diagnostic
    /// (objostrxml.cpp:548-565).
    pub fn writes_text(self, text: &[u8]) -> bool {
        match self {
            ArchiveEncoding::Xml => !text
                .iter()
                .any(|byte| matches!(byte, 0x01..=0x08 | 0x0b | 0x0c | 0x0e..=0x1f)),
            _ => text.iter().all(|&byte| is_visible(byte)),
        }
    }
}

/// Whether NCBI's `PrintErrorArchive` writes the error archive for this `-outfmt` text.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-929
/// ```c++
/// void PrintErrorArchive(const CArgs & a, const list<CRef<CBlast4_error> > & msg)
/// {
///     try {
///         if(NStr::StringToInt(NStr::TruncateSpaces(a[align_format::kArgOutputFormat].AsString())) == CFormattingArgs::eArchiveFormat) {
///             ...
///             CBlastFormat::PrintArchive(archive, a[kArgOutput].AsOutputFile());
///         }
///     } catch (...) {}
/// }
/// ```
/// The conversion throws (and nothing is written) for text after the number, such as
/// `11 qseqid`; the batches' archives accept such a text (`-outfmt "11 <spec>"`), only the
/// error archive does not. `TruncateSpaces` strips `isspace` characters at both ends,
/// `StringToInt` takes an optional sign and decimal digits (leading zeros allowed) that fit an
/// `int`.
pub fn writes_error_archive(outfmt: &str) -> bool {
    let is_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r');
    let text = outfmt.trim_matches(is_space);
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'+') => (false, &text[1..]),
        Some(b'-') => (true, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return false;
    }
    let mut value: i64 = 0;
    for b in digits.bytes() {
        value = value * 10 + i64::from(b - b'0');
        if value > i64::from(i32::MAX) + 1 {
            return false;
        }
    }
    let value = if negative { -value } else { value };
    value == i64::from(crate::blastinput::app::output_format::ARCHIVE_FORMAT)
}

/// The `queue-search` request of a run, as the options handle and the arguments make it.
#[derive(Clone, Debug)]
pub struct ArchiveRequest {
    /// `blastn`, `blastp`, ...
    pub program: &'static str,
    /// `megablast` or `plain` (`options::service`).
    pub service: &'static str,
    /// The `algorithm-options` without `EffectiveSearchSpace`.
    pub algorithm: Vec<options::Param>,
    /// `subject database "<-db as typed>"`; `None` for `-subject` (the sequences are written
    /// from the records).
    pub database: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Default)]
struct RunState {
    log: DiagLog,
    /// The `EffectiveSearchSpace` of the last searched batch (`None` before one).
    last_search_space: Option<i64>,
    /// Whether the diagnostics that are written now are not the archive's (the warnings that
    /// another format of the same search writes, which a run of the archive alone does not).
    muted: bool,
}

#[derive(Default)]
struct Inner {
    real: RunState,
}

impl Inner {
    fn state(&mut self) -> &mut RunState {
        &mut self.real
    }
}

struct Shared {
    encoding: ArchiveEncoding,
    request: ArchiveRequest,
    /// Whether the error archives are written (`writes_error_archive` of the `-outfmt` text;
    /// true when the run has no such text).
    error_archive: AtomicBool,
    state: Mutex<Inner>,
}

/// What an archive-writing run shares between its diagnostics stream and its writer.
#[derive(Clone)]
pub struct ArchiveHandle(Arc<Shared>);

impl std::fmt::Debug for ArchiveHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ArchiveHandle")
            .field("encoding", &self.0.encoding)
            .finish_non_exhaustive()
    }
}

impl ArchiveHandle {
    /// The handle of a run of `inputs.program` that writes the archive (`-outfmt 11`): the
    /// request from the recorded options, `ARCHIVE_FORMAT` read once. `database` is the
    /// `-db` argument as typed (`None` for `-subject`); `label` names the program in a refusal.
    /// Refuses, before anything is written, what LOSAT cannot write: a `-db` name or a
    /// `-window_masker_db` path with a byte that `ArchiveEncoding::writes_text` refuses (NCBI
    /// writes `#` and diagnostics there). These are the only strings the archive copies from
    /// the command line: the `-matrix` name is one of the supported matrices (the parser
    /// refuses any other before the search, in letters, digits and case only), and the other
    /// recorded options are numbers or switches.
    pub fn for_program(
        inputs: &OptionInputs,
        database: Option<&str>,
        label: &str,
    ) -> anyhow::Result<Self> {
        Self::for_encoding(ArchiveEncoding::from_env(), inputs, database, label)
    }

    /// `for_program` for a given encoding (`ARCHIVE_FORMAT` read by the caller).
    pub fn for_encoding(
        encoding: ArchiveEncoding,
        inputs: &OptionInputs,
        database: Option<&str>,
        label: &str,
    ) -> anyhow::Result<Self> {
        for (what, text) in [
            ("-db name", database),
            ("-window_masker_db path", inputs.window_masker_db.as_deref()),
        ] {
            if let Some(text) = text {
                if !encoding.writes_text(text.as_bytes()) {
                    anyhow::bail!(
                        "a {what} with a character the archive cannot write as NCBI does (outside the visible ASCII range in text and binary, a control character in XML) is not supported by LOSAT's {label} with -outfmt 11"
                    );
                }
            }
        }
        let program = match inputs.program {
            Program::Blastn => "blastn",
            Program::Blastp => "blastp",
            Program::Blastx => "blastx",
            Program::Tblastn => "tblastn",
            Program::Tblastx => "tblastx",
        };
        Ok(Self::new(
            encoding,
            ArchiveRequest {
                program,
                service: options::service(inputs.program, &inputs.task),
                algorithm: options::algorithm_options(inputs),
                database: database.map(|name| name.as_bytes().to_vec()),
            },
        ))
    }

    pub fn new(encoding: ArchiveEncoding, request: ArchiveRequest) -> Self {
        Self(Arc::new(Shared {
            encoding,
            request,
            error_archive: AtomicBool::new(true),
            state: Mutex::new(Inner::default()),
        }))
    }

    pub fn encoding(&self) -> ArchiveEncoding {
        self.0.encoding
    }

    /// Ties the error archives to the `-outfmt` text of the command line
    /// (`writes_error_archive`): NCBI's `PrintErrorArchive` writes none for `11 <spec>`.
    /// The archives of the batches are not affected.
    #[must_use]
    pub fn with_outfmt_text(self, outfmt: &str) -> Self {
        self.0
            .error_archive
            .store(writes_error_archive(outfmt), Ordering::Relaxed);
        self
    }

    /// Whether `write_error_archive` writes (see `with_outfmt_text`).
    pub fn error_archive_enabled(&self) -> bool {
        self.0.error_archive.load(Ordering::Relaxed)
    }

    pub fn request(&self) -> &ArchiveRequest {
        &self.0.request
    }

    fn with_state<R>(&self, f: impl FnOnce(&mut RunState) -> R) -> R {
        let mut inner = self
            .0
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        f(inner.state())
    }

    /// Saves the lines that a write to the diagnostics stream completes.
    pub fn feed(&self, bytes: &[u8]) {
        self.with_state(|state| {
            if !state.muted {
                state.log.feed(bytes);
            }
        });
    }

    /// While muted, what is written to the diagnostics stream is not saved. The report of a
    /// format other than the archive is written with the handle muted: a search with several
    /// formats writes the warnings of the first one once, and the archive of the same search
    /// alone has none of them.
    pub fn set_muted(&self, muted: bool) {
        self.with_state(|state| state.muted = muted);
    }

    /// Saves a message (the error lines of an aborted run).
    pub fn post(&self, severity: i32, line: &[u8]) {
        self.with_state(|state| state.log.post(severity, line));
    }

    /// Whether messages are saved since the last archive.
    pub fn has_messages(&self) -> bool {
        self.with_state(|state| !state.log.is_empty())
    }

    /// The messages saved since the last archive.
    pub fn take_messages(&self) -> Vec<Message> {
        self.with_state(|state| state.log.take())
    }

    /// The `EffectiveSearchSpace` option of the next archive: the batch's own value when
    /// the batch was searched (kept for the batches after it), else the one of the last
    /// searched batch, which is absent before one.
    pub fn search_space_option(&self, searched: Option<i64>) -> Option<i64> {
        self.with_state(|state| {
            if searched.is_some() {
                state.last_search_space = searched;
            }
            state.last_search_space
        })
    }
}

/// Whether `byte` is one of NCBI's visible characters (`GoodVisibleChar`: `' '` to `'~'`).
/// NCBI writes `#` for any other byte of a VisibleString and posts a diagnostic for each,
/// whose text and place LOSAT does not reproduce.
fn is_visible(byte: u8) -> bool {
    (0x20..=0x7E).contains(&byte)
}

/// The refusal of a title with a byte outside `0x20..=0x7E`, which the archive's VisibleString
/// cannot hold (see `is_visible`). `role` is `query` or `subject`.
pub fn check_titles(records: &[FastaRecord], role: &str, label: &str) -> anyhow::Result<()> {
    if let Some(record) = records
        .iter()
        .find(|record| !record.title.iter().all(|&byte| is_visible(byte)))
    {
        return Err(anyhow::Error::new(TitleRefusal(format!(
            "the title of the {role} {} has a character outside the visible ASCII range, which is not supported by LOSAT's {label} with -outfmt 11",
            record.local_id
        ))));
    }
    Ok(())
}

/// The error of `check_titles`, so that a caller can tell this refusal from the others (to
/// remove an `-out` file it created for nothing).
#[derive(Debug)]
pub struct TitleRefusal(String);

impl std::fmt::Display for TitleRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TitleRefusal {}

/// The error lines of a run that stopped with `error`, when NCBI writes an error archive for
/// it (`PrintErrorArchive` after its `CATCH_ALL`): the errors of the application's `Run`,
/// i.e. NCBI's own error texts, one message per diagnostic, severity `eDiag_Error`. LOSAT's own
/// refusals (and the argument errors raised before the application runs) write none.
/// `subject_file` is whether the error came from opening the `-query` file or reading the
/// `-subject` file, whose "Command line argument error" NCBI raises inside `Run`.
pub fn error_lines(error: &anyhow::Error, subject_file: bool) -> Option<Vec<Message>> {
    let native = error.downcast_ref::<crate::cli::NativeError>()?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    // `Error: NCBI C++ Exception:` is that report of an exception of NCBI's C++ Toolkit
    // (BLASTN's `-window_masker_db` file that cannot be opened).
    const PREFIXES: [&str; 5] = [
        "BLAST query error:",
        "BLAST Database error:",
        "BLAST engine error:",
        "BLAST query/options error:",
        "Error: NCBI C++ Exception:",
    ];
    let known = PREFIXES
        .iter()
        .any(|prefix| native.message.starts_with(prefix))
        || (subject_file && native.message.starts_with("Command line argument error:"));
    if !known {
        return None;
    }
    // One diagnostic per `ERR_POST`: the text of an error is one message even when it has
    // several lines (an exception's chain), and the macro of an invalid input posts the line
    // "Please refer to the BLAST+ user manual." as a second one.
    const REFER: &str = "Please refer to the BLAST+ user manual.";
    let text = native.message.trim_end_matches('\n');
    let (first, refer) = match text.strip_suffix(REFER) {
        Some(head) if head.ends_with('\n') => (head, Some(REFER)),
        _ => (text, None),
    };
    Some(
        std::iter::once(first)
            .chain(refer)
            .filter_map(|part| {
                let text = crate::report::diag_log::sanitize(part.as_bytes());
                (!text.is_empty()).then_some((ERROR, text))
            })
            .collect(),
    )
}

/// A diagnostics stream that saves what it writes as the messages of the next archive.
pub struct DiagTee<'a> {
    inner: &'a mut (dyn Write + Send),
    handle: ArchiveHandle,
}

impl<'a> DiagTee<'a> {
    pub fn new(inner: &'a mut (dyn Write + Send), handle: ArchiveHandle) -> Self {
        Self { inner, handle }
    }
}

impl Write for DiagTee<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.handle.feed(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn archive_format_is_read_like_print_archive() {
        let of = |value: &str| ArchiveEncoding::of(Some(OsStr::new(value)));
        assert_eq!(ArchiveEncoding::of(None), ArchiveEncoding::Text);
        assert_eq!(of(""), ArchiveEncoding::Text);
        assert_eq!(of("xml"), ArchiveEncoding::Xml);
        assert_eq!(of("XML"), ArchiveEncoding::Xml);
        assert_eq!(of("Xml"), ArchiveEncoding::Xml);
        assert_eq!(of("bin"), ArchiveEncoding::Binary);
        assert_eq!(of("BIN"), ArchiveEncoding::Binary);
        assert_eq!(of("binary"), ArchiveEncoding::Binary);
        assert_eq!(of("bin foo"), ArchiveEncoding::Binary);
        assert_eq!(of(" bin"), ArchiveEncoding::Nothing);
        assert_eq!(of("asn"), ArchiveEncoding::Nothing);
        assert_eq!(of("text"), ArchiveEncoding::Nothing);
        assert_eq!(of("foo"), ArchiveEncoding::Nothing);
        assert_eq!(of("xmlx"), ArchiveEncoding::Nothing);
        assert_eq!(of("bi"), ArchiveEncoding::Nothing);
    }

    /// The strings the archive copies from the command line are refused before the first byte
    /// where the encoding cannot write them as NCBI does: text and binary write `#` for a byte
    /// outside `0x20..=0x7E`, XML writes every byte but a control one as it is.
    #[test]
    fn a_db_name_or_masker_path_the_encoding_cannot_write_is_refused() {
        let inputs = |path: Option<&str>| {
            let mut inputs = OptionInputs::new(Program::Blastn, "megablast");
            inputs.window_masker_db = path.map(str::to_string);
            inputs
        };
        let ok = |encoding, path: Option<&str>, db: Option<&str>| {
            ArchiveHandle::for_encoding(encoding, &inputs(path), db, "BLASTN").is_ok()
        };
        for encoding in [
            ArchiveEncoding::Text,
            ArchiveEncoding::Binary,
            ArchiveEncoding::Xml,
        ] {
            assert!(ok(encoding, None, None));
            assert!(ok(encoding, Some("wm&<\"'q/mg.obinary"), Some("d b&<")));
            // A control byte: NCBI writes `#` and a diagnostic, but TAB, LF and CR in XML.
            for bad in ["a\x01b", "a\x1fb", "a\x0bb", "a\x0eb"] {
                assert!(!ok(encoding, Some(bad), None), "{encoding:?} {bad:?}");
                assert!(!ok(encoding, None, Some(bad)), "{encoding:?} {bad:?}");
            }
        }
        for good in ["a\tb", "a\nb", "a\rb"] {
            assert!(ok(ArchiveEncoding::Xml, Some(good), Some(good)));
            assert!(!ok(ArchiveEncoding::Text, Some(good), None));
        }
        // Bytes from 0x7F up: refused in text and binary, written as they are in XML.
        for high in ["wm\u{e9}/mg.obinary", "d\u{7f}b", "\u{1f600}"] {
            assert!(ok(ArchiveEncoding::Xml, Some(high), Some(high)), "{high:?}");
            for encoding in [ArchiveEncoding::Text, ArchiveEncoding::Binary] {
                assert!(!ok(encoding, Some(high), None), "{encoding:?} {high:?}");
                assert!(!ok(encoding, None, Some(high)), "{encoding:?} {high:?}");
            }
        }
        let refusal = ArchiveHandle::for_encoding(
            ArchiveEncoding::Text,
            &inputs(Some("w\u{e9}")),
            None,
            "BLASTN",
        )
        .unwrap_err()
        .to_string();
        assert!(refusal.contains("-window_masker_db path"), "{refusal}");
        assert!(refusal.contains("-outfmt 11"), "{refusal}");
    }

    #[test]
    fn the_search_space_option_carries_over_unsearched_batches() {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastp",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
            },
        );
        assert_eq!(handle.search_space_option(None), None);
        assert_eq!(handle.search_space_option(Some(8738)), Some(8738));
        assert_eq!(handle.search_space_option(None), Some(8738));
        assert_eq!(handle.search_space_option(Some(0)), Some(0));
        assert_eq!(handle.search_space_option(None), Some(0));
    }

    #[test]
    fn the_tee_saves_lines_for_the_next_archive() {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastn",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
            },
        );
        let mut seen: Vec<u8> = Vec::new();
        {
            let mut tee = DiagTee::new(&mut seen, handle.clone());
            tee.write_all(b"Warning: [blastn] a  b\n").unwrap();
        }
        assert_eq!(seen, b"Warning: [blastn] a  b\n");
        assert_eq!(
            handle.take_messages(),
            vec![(1, b"Warning: [blastn] a b".to_vec())]
        );
        assert!(handle.take_messages().is_empty());
    }

    /// NCBI 2.17.0 (`PrintErrorArchive`): each text below was run with a missing `-query` and
    /// wrote the error archive (`true`) or nothing (`false`).
    #[test]
    fn the_error_archive_follows_the_outfmt_text_like_print_error_archive() {
        for (text, writes) in [
            ("11", true),
            (" 11", true),
            ("11 ", true),
            ("\t11\n", true),
            ("\x0b11\x0c", true),
            ("011", true),
            ("+11", true),
            ("00000000000011", true),
            ("11 qseqid", false),
            (" 11 qseqid", false),
            ("11qseqid", false),
            ("11.0", false),
            ("1 1", false),
            ("-11", false),
            ("0x0B", false),
            ("1_1", false),
            ("11,0", false),
            ("4294967307", false),
            ("2147483648", false),
            ("", false),
            (" ", false),
            ("+", false),
            ("6", false),
        ] {
            assert_eq!(writes_error_archive(text), writes, "{text:?}");
        }
    }

    #[test]
    fn a_handle_tied_to_a_spec_writes_no_error_archive() {
        let make = || {
            ArchiveHandle::new(
                ArchiveEncoding::Text,
                ArchiveRequest {
                    program: "tblastx",
                    service: "plain",
                    algorithm: Vec::new(),
                    database: None,
                },
            )
        };
        let lines = vec![(ERROR, b"Error: x".to_vec())];
        // No text (the API): the error archive is written.
        let mut sink = Vec::new();
        write_error_archive(&make(), &lines, &mut sink).unwrap();
        assert!(!sink.is_empty());
        // `11`: written; `11 qseqid`: nothing, and the messages stay unread.
        let mut sink = Vec::new();
        write_error_archive(&make().with_outfmt_text(" 11 "), &lines, &mut sink).unwrap();
        assert!(!sink.is_empty());
        let mut sink = Vec::new();
        let spec = make().with_outfmt_text("11 qseqid");
        assert!(!spec.error_archive_enabled());
        write_error_archive(&spec, &lines, &mut sink).unwrap();
        assert!(sink.is_empty());
    }
}

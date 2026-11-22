//! NCBI's `CFormatGuess` (c++/src/util/format_guess.cpp, NCBI BLAST+ 2.17.0), the one port of it that LOSAT's callers
//! share: makeblastdb's FASTA check (`writedb::format_guess::appears_fasta`, the preferred formats eBinaryASN,
//! eTextASN and eFasta) and its mask files (`writedb::mask_info::guess`, eBinaryASN and eTextASN), and the
//! `-in_pssm` file of psiblast and tblastn, whose `CFormatGuess` (its `Format` of the input) has no hints
//! (`algorithm::pssm_search::in_pssm`): every test of `sm_CheckOrder` up to eBinaryASN, then ePhrapAce and
//! eUCSCRegion.
//!
//! The port keeps the state of NCBI's object as it is: the stream and what its reads put back
//! (`CStreamUtils::Stepback`), the test buffer and its size (doubled while the sample is all comment), the split lines
//! (which `IsInputRepeatMaskerWithHeader` trims in place, `TestFormatNewick` replaces when the split fails, and an
//! all-comment pass can leave from a shorter buffer), the statistics of `EnsureStats` (counted once, on the buffer
//! of the first call), and the stream position that `TestFormatNewick`'s NEXUS scan moves. C's `isalpha`, `isdigit`,
//! `isspace`, `isprint` and `isgraph` are the "C" locale's (ASCII only; a `char` of 0x80 or more is negative, and
//! glibc's tables give it no class). The C++ works on bytes, so this port works on bytes.
//!
//! Not ported: `TestFormatJson` (eJSON comes after every format a caller reads, so `guess_format` answers
//! [`Format::UnknownOrJson`]), the format hints (each caller of the hints runs the tests it prefers itself), and
//! `x_TestInput` (each caller opens its file first). Three tests read memory that their input does not decide (an
//! element past the end of a vector, a C string past the end of the test buffer): [`Stop::Undecided`], or
//! [`Guessed::fasta_undecided`] for `TestFormatFasta`, whose answer can only be eFasta or the guess without it. One
//! throws: `TestFormatGlimmer3`'s `NStr::StringToInt` of a frame beyond `int` ([`Stop::Exception`]).

use std::io::Read;
use std::ops::Range;

use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::{is_space, truncate_spaces, truncate_spaces_begin};

mod line_formats;

#[cfg(test)]
mod guess_tests;

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:174-183
pub(crate) const F_DNA_MAIN_ALPHABET: u8 = 1 << 0;
pub(crate) const F_DNA_AMBIG_ALPHABET: u8 = 1 << 1;
pub(crate) const F_PROTEIN_ALPHABET: u8 = 1 << 2;
pub(crate) const F_LINE_END: u8 = 1 << 3;
pub(crate) const F_ALPHA: u8 = 1 << 4;
pub(crate) const F_DIGIT: u8 = 1 << 5;
pub(crate) const F_SPACE: u8 = 1 << 6;
pub(crate) const F_INVALID: u8 = 1 << 7;

/// C `isprint` in the "C" locale.
pub(crate) fn c_isprint(c: u8) -> bool {
    (0x20..=0x7e).contains(&c)
}

/// C `isgraph` in the "C" locale.
fn c_isgraph(c: u8) -> bool {
    (0x21..=0x7e).contains(&c)
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:249-286
/// `symbol_type_table` after `init_symbol_type_table()`.
pub(crate) static SYMBOL_TYPE_TABLE: [u8; 256] = init_symbol_type_table();

const fn init_symbol_type_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    let main = b"ACGNTU";
    let mut i = 0;
    while i < main.len() {
        table[main[i] as usize] |= F_DNA_MAIN_ALPHABET;
        table[main[i].to_ascii_lowercase() as usize] |= F_DNA_MAIN_ALPHABET;
        i += 1;
    }
    let ambig = b"BDHKMRSVWY";
    let mut i = 0;
    while i < ambig.len() {
        table[ambig[i] as usize] |= F_DNA_AMBIG_ALPHABET;
        table[ambig[i].to_ascii_lowercase() as usize] |= F_DNA_AMBIG_ALPHABET;
        i += 1;
    }
    let protein = b"ACDEFGHIKLMNPQRSTVWYBZX";
    let mut i = 0;
    while i < protein.len() {
        table[protein[i] as usize] |= F_PROTEIN_ALPHABET;
        table[protein[i].to_ascii_lowercase() as usize] |= F_PROTEIN_ALPHABET;
        i += 1;
    }
    table[b'-' as usize] |= F_DNA_MAIN_ALPHABET | F_PROTEIN_ALPHABET;
    table[b'*' as usize] |= F_PROTEIN_ALPHABET;
    table[b'\r' as usize] |= F_LINE_END;
    table[b'\n' as usize] |= F_LINE_END;
    let mut c = 1;
    while c < 256 {
        let b = c as u8;
        if b.is_ascii_alphabetic() {
            table[c] |= F_ALPHA;
        }
        if b.is_ascii_digit() {
            table[c] |= F_DIGIT;
        }
        if is_space(b) {
            table[c] |= F_SPACE;
        }
        c += 1;
    }
    table[0] |= F_INVALID;
    table
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:2509-2510
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7439-7545
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-321
// `NStr::Split` of `data` at `delims` into lines, with `NStr::fSplit_Tokenize`, with `delims` a set of single characters
// (`fSplit_ByPattern` is not set, so every character of the delimiter string is a delimiter): the non-empty runs of
// bytes between delimiters. Leading, trailing and repeated delimiters give no empty token.

/// The tokens a test reads by position: `IsLinePsl`'s 21 after its first column is the most.
const KEPT_TOKENS: usize = 24;

/// The tokens of `NStr::Split` of `str` at `delims` with `NStr::fSplit_Tokenize` (that is,
/// `fSplit_MergeDelimiters` and `fSplit_Truncate` together, the same flags): how many there are, and the first
/// [`KEPT_TOKENS`] of them, which are all the tests read by position (a test that reads every token goes through the line itself).
#[derive(Clone, Copy)]
pub(crate) struct Tokens<'a> {
    kept: [&'a [u8]; KEPT_TOKENS],
    len: usize,
}

impl<'a> Tokens<'a> {
    /// `tokens.size()`.
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// `tokens.empty()`.
    pub(crate) fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// `tokens[i]` of a token that exists (`i < len()`).
    pub(crate) fn at(&self, i: usize) -> &'a [u8] {
        assert!(i < self.len.min(KEPT_TOKENS), "token {i} of {}", self.len);
        self.kept[i]
    }

    /// `tokens[i]` when it exists.
    pub(crate) fn get(&self, i: usize) -> Option<&'a [u8]> {
        (i < self.len.min(KEPT_TOKENS)).then(|| self.kept[i])
    }
}

/// [`Tokens`] of `str` split at the bytes `delims`. The one port of `NStr::Split` with
/// `fSplit_Tokenize`: the format guess's tests and the Clustal scanner of psiblast's `-in_msa`
/// (`blastinput::aln_reader::clustal`) read their lines through it.
pub(crate) fn tokens<'a>(str: &'a [u8], delims: &[u8]) -> Tokens<'a> {
    let mut tokens = Tokens {
        kept: [&[]; KEPT_TOKENS],
        len: 0,
    };
    for token in str
        .split(|c| delims.contains(c))
        .filter(|token| !token.is_empty())
    {
        if tokens.len < KEPT_TOKENS {
            tokens.kept[tokens.len] = token;
        }
        tokens.len += 1;
    }
    tokens
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3574-3597
/// `NStr::SplitInTwo` of `str` at `delims` with `flags`: `str1` (the bytes before the first delimiter, all of `str`
/// when there is none), `str2` (after it; after the run of delimiters with `fSplit_MergeDelimiters`, `merge`) and
/// whether a delimiter was found. `CStrTokenizeBase::Advance` (cited at [`tokens`]) without
/// `fSplit_Truncate_Begin` keeps an empty `str1` before a leading delimiter.
pub(crate) fn split_in_two<'a>(
    str: &'a [u8],
    delims: &[u8],
    merge: bool,
) -> (&'a [u8], &'a [u8], bool) {
    match str.iter().position(|c| delims.contains(c)) {
        None => (str, &[], false),
        Some(at) => {
            let mut rest = at + 1;
            if merge {
                while rest < str.len() && delims.contains(&str[rest]) {
                    rest += 1;
                }
            }
            (&str[..at], &str[rest..], true)
        }
    }
}

/// `NStr::Find` of `pattern` in `str` (case-sensitive) is not `NPOS`.
fn contains(str: &[u8], pattern: &[u8]) -> bool {
    find(str, pattern).is_some()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2887-2920 (`NStr::Find`, `eCase`, `eForwardSearch`)
/// The position of the first `byte` in `hay`, read eight bytes at a time (a byte equal to `byte` is a zero byte of
/// the word XOR `byte` in every lane; the lowest flagged lane is the first zero, a borrow only flags lanes above
/// one).
fn find_byte(hay: &[u8], byte: u8) -> Option<usize> {
    const LOW: u64 = 0x0101_0101_0101_0101;
    const HIGH: u64 = 0x8080_8080_8080_8080;
    let lanes = LOW * u64::from(byte);
    let mut chunks = hay.chunks_exact(8);
    let mut at = 0;
    for chunk in &mut chunks {
        let word = u64::from_le_bytes(chunk.try_into().expect("eight bytes")) ^ lanes;
        let zero = word.wrapping_sub(LOW) & !word & HIGH;
        if zero != 0 {
            return Some(at + (zero.trailing_zeros() / 8) as usize);
        }
        at += 8;
    }
    chunks
        .remainder()
        .iter()
        .position(|&c| c == byte)
        .map(|rest| at + rest)
}

/// `NStr::Find` of `pattern` in `str`: the position of the first `pattern` (one byte or more) in `str`, found from its
/// first byte.
fn find(str: &[u8], pattern: &[u8]) -> Option<usize> {
    let (&first, rest) = pattern.split_first()?;
    if pattern.len() > str.len() {
        return None;
    }
    let last = str.len() - pattern.len();
    let mut from = 0;
    while from <= last {
        let at = from + find_byte(&str[from..=last], first)?;
        if &str[at + 1..at + pattern.len()] == rest {
            return Some(at);
        }
        from = at + 1;
    }
    None
}

/// The C string that starts `bytes` (up to its first NUL): what the `std::string` and `CTempString`
/// constructors from a `const char*` and `sscanf` read of a buffer that holds a NUL.
fn c_string(bytes: &[u8]) -> &[u8] {
    &bytes[..bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len())]
}

/// Reads from `input` until `data` holds `want` bytes, the input ends, or a read fails. A read error counts as end of data.
/// `ended` remembers that the input is at its end (NCBI's stream stays at EOF as well).
fn fill_to<R: Read>(input: &mut R, data: &mut Vec<u8>, want: usize, ended: &mut bool) {
    if *ended || data.len() >= want {
        return;
    }
    let before = data.len();
    let missing = (want - before) as u64;
    // `read_to_end` retries interrupted reads and keeps the bytes read before an error.
    let _ = input.by_ref().take(missing).read_to_end(data);
    if data.len() < want {
        *ended = true;
    }
}

/// The bytes of the stream `CFormatGuess` reads: all of them for a stream over bytes in memory, those read so far
/// for a reader.
pub(crate) trait Source {
    /// The bytes from the start of the stream.
    fn bytes(&self) -> &[u8];
    /// Reads on until [`Source::bytes`] holds `want` bytes or the input ends.
    fn fill_to(&mut self, want: usize);
}

impl Source for &[u8] {
    fn bytes(&self) -> &[u8] {
        self
    }

    fn fill_to(&mut self, _want: usize) {}
}

/// A [`Source`] that reads a reader as far as the guess reads it.
pub(crate) struct ReadSource<R> {
    reader: R,
    bytes: Vec<u8>,
    ended: bool,
}

impl<R: Read> Source for ReadSource<R> {
    fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    fn fill_to(&mut self, want: usize) {
        fill_to(&mut self.reader, &mut self.bytes, want, &mut self.ended);
    }
}

/// The stream `CFormatGuess` reads (`m_Stream`): its bytes (a read that is stepped back leaves them for the next one,
/// as `CStreamUtils::Stepback` puts them back into the stream), the read position and the stream state.
struct Stream<S> {
    source: S,
    pos: usize,
    /// `good()`: neither `eofbit` nor `failbit` is set.
    good: bool,
}

impl<S: Source> Stream<S> {
    /// `read` of `n` bytes from `m_Stream` and its `gcount()`: the bytes read, as a range of `bytes`. A read of fewer than
    /// `n` bytes sets `eofbit` and `failbit`; a stream that is not good reads nothing (its `sentry` fails).
    fn read(&mut self, n: usize) -> Range<usize> {
        let start = self.pos;
        if !self.good {
            return start..start;
        }
        self.source.fill_to(start.saturating_add(n));
        let end = self.source.bytes().len().min(start.saturating_add(n));
        if end - start < n {
            self.good = false;
        }
        self.pos = end;
        start..end
    }

    /// `m_Stream.clear()`.
    fn clear(&mut self) {
        self.good = true;
    }

    /// `CStreamUtils::Stepback` of `m_Stream` with the bytes just read: the next read reads them again.
    fn stepback(&mut self, read: &Range<usize>) {
        self.pos = read.start;
    }
}

/// What `CFormatGuess::EnsureStats` counts (`m_iStatsCount*`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Stats {
    /// `m_iStatsCountAlNumChars`
    pub(crate) al_num_chars: usize,
    /// `m_iStatsCountBraces`
    pub(crate) braces: usize,
    /// `m_iStatsCountData`
    pub(crate) data: usize,
    /// `m_iStatsCountDnaChars`
    pub(crate) dna_chars: usize,
    /// `m_iStatsCountAaChars`
    pub(crate) aa_chars: usize,
}

/// The non-empty runs of `data` between the bytes for which `ends` holds (`NStr::Split` with `fSplit_Tokenize`, the
/// lines of `NcbiGetline`): the bytes that are no line end and follow one or the start.
fn count_runs(data: &[u8], ends: impl Fn(u8) -> bool + Copy) -> usize {
    let mut runs = 0;
    let mut after_end = true;
    for &c in data {
        let end = ends(c);
        runs += usize::from(after_end & !end);
        after_end = end;
    }
    runs
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:687-749
// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:55-134

/// `EnsureStats` on the sample. `NcbiGetline` with `\r\n` as the delimiters ends a line at `\r` or
/// `\n`; a `\r` followed by `\n` is one delimiter and every other pair (`\n\r`, `\r\r`, `\n\n`) is two, so an empty line only adds an empty `strLine`, which
/// `EnsureStats` does not count. The lines are therefore the non-empty runs between `\r`/`\n` bytes, each with one `\n`
/// appended; the loop ends at the first read that finds nothing left (`fail()`), after the last line.
fn stats_of(data: &[u8]) -> Stats {
    let line_end = |c: u8| c == b'\n' || c == b'\r';
    // Every byte counted by value, then those of the lines that start with `>` (found from their `>`) moved to
    // their own counts; the line ends themselves are no line's bytes, and each line counts a `\n` (a space).
    let mut all = [0u32; 256];
    let mut lines = 0;
    let mut after_end = true;
    for &c in data {
        all[usize::from(c)] += 1;
        let end = line_end(c);
        lines += usize::from(after_end & !end);
        after_end = end;
    }
    let mut header = [0u32; 256];
    let mut from = 0;
    while let Some(at) = find_byte(&data[from..], b'>').map(|at| from + at) {
        if at == 0 || line_end(data[at - 1]) {
            let end = data[at..]
                .iter()
                .position(|&c| line_end(c))
                .map_or(data.len(), |len| at + len);
            for &c in &data[at..end] {
                header[usize::from(c)] += 1;
            }
            from = end;
        } else {
            from = at + 1;
        }
        if from >= data.len() {
            break;
        }
    }
    let mut stats = Stats {
        al_num_chars: lines,
        ..Stats::default()
    };
    for (c, &type_) in SYMBOL_TYPE_TABLE.iter().enumerate() {
        if c == usize::from(b'\n') || c == usize::from(b'\r') {
            continue;
        }
        let (in_header, in_data) = (header[c] as usize, (all[c] - header[c]) as usize);
        if type_ & (F_ALPHA | F_DIGIT | F_SPACE) != 0 {
            stats.al_num_chars += in_header + in_data;
        } else if c == usize::from(b'{') || c == usize::from(b'}') {
            stats.braces += in_header + in_data;
        }
        if type_ & F_SPACE == 0 {
            stats.data += in_data;
        }
        if type_ & F_DNA_MAIN_ALPHABET != 0 {
            stats.dna_chars += in_data;
        }
        if type_ & F_PROTEIN_ALPHABET != 0 {
            stats.aa_chars += in_data;
        }
    }
    stats
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:46-87
/// `CFormatGuess::EFormat`, of the formats `sm_CheckOrder` lists.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    Bam,
    Zip,
    Zstd,
    GZip,
    BZip2,
    Lzo,
    Sra,
    Psl,
    Rmo,
    Vcf,
    Gvf,
    Gff3,
    Gtf,
    GffAugustus,
    Gff2,
    Glimmer3,
    Agp,
    Xml,
    Wiggle,
    Newick,
    Bed,
    Bed15,
    Hgvs,
    DistanceMatrix,
    FlatFileSequence,
    FlatFileUniProt,
    FlatFileEna,
    FlatFileGenbank,
    FiveColFeatureTable,
    SnpMarkers,
    Fasta,
    TextAsn,
    Alignment,
    Taxplot,
    Table,
    BinaryAsn,
    PhrapAce,
    UcscRegion,
    /// `eJSON` or `eUnknown`: no test up to `eUCSCRegion` takes the input (`TestFormatJson` is not ported), or
    /// `eUnknown` when there is no test buffer and `TestFormatNewick` does not take the input.
    UnknownOrJson,
}

/// `sm_CheckOrder` without `eJSON`.
const CHECK_ORDER: [Format; 38] = [
    Format::Bam,
    Format::Zip,
    Format::Zstd,
    Format::GZip,
    Format::BZip2,
    Format::Lzo,
    Format::Sra,
    Format::Psl,
    Format::Rmo,
    Format::Vcf,
    Format::Gvf,
    Format::Gff3,
    Format::Gtf,
    Format::GffAugustus,
    Format::Gff2,
    Format::Glimmer3,
    Format::Agp,
    Format::Xml,
    Format::Wiggle,
    Format::Newick,
    Format::Bed,
    Format::Bed15,
    Format::Hgvs,
    Format::DistanceMatrix,
    Format::FlatFileSequence,
    Format::FlatFileUniProt,
    Format::FlatFileEna,
    Format::FlatFileGenbank,
    Format::FiveColFeatureTable,
    Format::SnpMarkers,
    Format::Fasta,
    Format::TextAsn,
    Format::Alignment,
    Format::Taxplot,
    Format::Table,
    Format::BinaryAsn,
    Format::PhrapAce,
    Format::UcscRegion,
];

/// A test of the guess that does not give an answer.
#[derive(Debug)]
pub(crate) enum Stop {
    /// The test throws (`TestFormatGlimmer3`'s `NStr::StringToInt`): `GuessFormat` passes the exception on.
    Exception(NcbiException),
    /// The test reads memory its input does not decide (named here): NCBI's answer, or whether it crashes, is not
    /// known.
    Undecided(&'static str),
}

/// What `CFormatGuess::GuessFormat` answers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Guessed {
    pub(crate) format: Format,
    /// Where the guess leaves the stream: its start, or past what `TestFormatNewick` read of a NEXUS file.
    pub(crate) stream_pos: usize,
    /// `TestFormatFasta` searched past its test buffer (`NStr::Find` of a C string that the buffer does not end), so
    /// it may have answered eFasta; `format` is the guess had it not.
    pub(crate) fasta_undecided: bool,
}

/// The line ends `EnsureSplitLines` split its buffer at.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum LineEnds {
    /// No split: the buffer is the one line (it holds the whole file).
    #[default]
    Whole,
    /// `\r` and `\n` (the buffer holds a `\r\n`).
    CrLf,
    /// `\n`.
    Lf,
    /// `\r`.
    Cr,
}

impl LineEnds {
    /// Whether `c` ends a line.
    fn ends(self, c: u8) -> bool {
        match self {
            LineEnds::Whole => false,
            LineEnds::CrLf => c == b'\r' || c == b'\n',
            LineEnds::Lf => c == b'\n',
            LineEnds::Cr => c == b'\r',
        }
    }
}

/// `m_TestLines`, without a copy of a line: the non-empty runs of the stream's bytes `range` between `ends`, the
/// first `count` of them (`EnsureSplitLines` drops a last line that the buffer may have cut), the first `trimmed` of
/// those without the spaces at their ends (`IsInputRepeatMaskerWithHeader` trims them in place), then `extra`, the C
/// string that `TestFormatNewick` (or `TestFormatHgvs`) adds when the split fails. The tests read the lines again
/// from the bytes each time.
#[derive(Clone, Debug, Default)]
struct Lines {
    range: Range<usize>,
    ends: LineEnds,
    count: usize,
    trimmed: usize,
    extra: Option<Range<usize>>,
}

impl Lines {
    /// `m_TestLines.size()`.
    fn len(&self) -> usize {
        self.count + usize::from(self.extra.is_some())
    }

    /// `m_TestLines.empty()`.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `m_TestLines.clear()`.
    fn clear(&mut self) {
        *self = Lines::default();
    }
}

/// `CFormatGuess` on a stream (the `CFormatGuess` constructor from a `CNcbiIstream&`).
pub(crate) struct FormatGuess<S> {
    stream: Stream<S>,
    /// `m_pTestBuffer` (`m_iTestDataSize` bytes), as the bytes of the stream it holds.
    buffer: Option<Range<usize>>,
    /// A NUL follows the buffer (`TestFormatNewick`'s and `TestFormatHgvs`'s buffers; `EnsureTestBuffer`'s has none).
    terminated: bool,
    /// `m_iTestBufferSize`
    buffer_size: usize,
    /// `m_TestLines`.
    lines: Lines,
    /// `m_bSplitDone`
    split_done: bool,
    /// `m_iStatsCount*`, once `m_bStatsAreValid`.
    stats: Option<Stats>,
    fasta_undecided: bool,
}

impl<R: Read> FormatGuess<ReadSource<R>> {
    /// `CFormatGuess` on the input (`Initialize`), here a reader.
    pub(crate) fn new(input: R) -> Self {
        FormatGuess::on(ReadSource {
            reader: input,
            bytes: Vec::new(),
            ended: false,
        })
    }

    /// The bytes the guess took from its input (makeblastdb reads the mask file on from them).
    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.stream.source.bytes
    }

    /// The bytes the guess took from its input and the input after them: a reader of the stream
    /// reads on from the bytes ([`Guessed::stream_pos`] on), as NCBI's reader reads the bytes
    /// that the guess put back into the stream.
    pub(crate) fn into_parts(self) -> (Vec<u8>, R) {
        (self.stream.source.bytes, self.stream.source.reader)
    }
}

impl<'a> FormatGuess<&'a [u8]> {
    /// `CFormatGuess` on the input, a stream over `bytes` (the whole input, which the guess does not copy).
    pub(crate) fn over(bytes: &'a [u8]) -> Self {
        FormatGuess::on(bytes)
    }
}

impl<S: Source> FormatGuess<S> {
    /// `CFormatGuess` on the input (`Initialize`).
    fn on(source: S) -> Self {
        FormatGuess {
            stream: Stream {
                source,
                pos: 0,
                good: true,
            },
            buffer: None,
            terminated: false,
            buffer_size: 0,
            lines: Lines::default(),
            split_done: false,
            stats: None,
            fasta_undecided: false,
        }
    }

    /// The test buffer (`m_pTestBuffer`, `m_iTestDataSize` bytes; empty when there is none).
    pub(crate) fn data(&self) -> &[u8] {
        match &self.buffer {
            Some(range) => &self.stream.source.bytes()[range.clone()],
            None => &[],
        }
    }

    /// `m_TestLines`, read from the stream's bytes.
    fn lines(&self) -> impl Iterator<Item = &[u8]> + Clone + '_ {
        let bytes = self.stream.source.bytes();
        let Lines {
            range,
            ends,
            count,
            trimmed,
            extra,
        } = &self.lines;
        let (ends, trimmed) = (*ends, *trimmed);
        bytes[range.clone()]
            .split(move |&c| ends.ends(c))
            .filter(|line| !line.is_empty())
            .take(*count)
            .enumerate()
            .map(move |(at, line)| {
                if at < trimmed {
                    truncate_spaces(line)
                } else {
                    line
                }
            })
            .chain(extra.clone().map(|extra| &bytes[extra]))
    }

    /// Whether `pattern` (of bytes that end no line) is in the bytes the lines are cut from: a line can hold it,
    /// or begin with it, only then (the test reads the lines after that).
    fn lines_may_hold(&self, pattern: &[u8]) -> bool {
        let bytes = self.stream.source.bytes();
        find(&bytes[self.lines.range.clone()], pattern).is_some()
            || self
                .lines
                .extra
                .clone()
                .is_some_and(|extra| find(&bytes[extra], pattern).is_some())
    }

    /// `m_TestLines[at]`.
    fn line_at(&self, at: usize) -> Option<&[u8]> {
        self.lines().nth(at)
    }

    /// The `push_back` of one line to `m_TestLines` in `TestFormatNewick` and `TestFormatHgvs`: their
    /// split failed, so this is the one line they add (the split lines left are those of a buffer whose split failed too, or of a pass that
    /// doubled it).
    fn push_line(&mut self, line: Range<usize>) {
        debug_assert!(self.lines.extra.is_none(), "one line added");
        self.lines.extra = Some(line);
    }

    /// The line `m_pTestBuffer` holds as a C string (added to `m_TestLines` by `push_back` in
    /// `TestFormatNewick` and `TestFormatHgvs`).
    fn c_string_line(&self) -> Range<usize> {
        let buffer = self.buffer.clone().unwrap_or_default();
        buffer.start..buffer.start + c_string(self.data()).len()
    }

    /// The statistics of `EnsureStats` (after a call that returned true).
    pub(crate) fn stats(&self) -> Stats {
        self.stats.unwrap_or_default()
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:452-491
    /// `GuessFormat()` without hints.
    pub(crate) fn guess_format(&mut self) -> Result<Guessed, Stop> {
        let format = if !self.ensure_test_buffer() {
            if self.test_format_newick() {
                Format::Newick
            } else {
                Format::UnknownOrJson
            }
        } else {
            let mut found = Format::UnknownOrJson;
            for format in CHECK_ORDER {
                if self.test_format(format)? {
                    found = format;
                    break;
                }
            }
            found
        };
        Ok(Guessed {
            format,
            stream_pos: self.stream.pos,
            fasta_undecided: self.fasta_undecided,
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:513-601 (`x_TestFormat`: each format
    // calls its own test, and the UCSC region format always answers false)
    /// `x_TestFormat` of `format`.
    fn test_format(&mut self, format: Format) -> Result<bool, Stop> {
        Ok(match format {
            Format::Bam => false,
            Format::Zip => self.test_format_zip(),
            Format::Zstd => self.test_magic(&[0x28, 0xb5, 0x2f, 0xfd]),
            Format::GZip => self.test_magic(&[0x1f, 0x8b]),
            Format::BZip2 => self.test_format_bzip2(),
            Format::Lzo => self.test_format_lzo(),
            Format::Sra => self.test_format_sra(),
            Format::Psl => self.test_format_psl(),
            Format::Rmo => self.test_format_repeat_masker(),
            Format::Vcf => self.test_format_vcf(),
            Format::Gvf => self.test_format_gvf()?,
            Format::Gff3 => self.test_format_gff3(),
            Format::Gtf => self.test_format_gtf(),
            Format::GffAugustus => self.test_format_augustus(),
            Format::Gff2 => self.test_format_gff2(),
            Format::Glimmer3 => self.test_format_glimmer3()?,
            Format::Agp => self.test_format_agp(),
            Format::Xml => self.test_format_xml(),
            Format::Wiggle => self.test_format_wiggle(),
            Format::Newick => self.test_format_newick(),
            Format::Bed => self.test_format_bed(),
            Format::Bed15 => self.test_format_bed15(),
            Format::Hgvs => self.test_format_hgvs(),
            Format::DistanceMatrix => self.test_format_distance_matrix(),
            Format::FlatFileSequence => self.test_format_flat_file_sequence(),
            Format::FlatFileUniProt => self.test_format_flat_file_uniprot(),
            Format::FlatFileEna => self.test_format_flat_file_ena(),
            Format::FlatFileGenbank => self.test_format_flat_file_genbank(),
            Format::FiveColFeatureTable => self.test_format_five_col_feature_table(),
            Format::SnpMarkers => self.test_format_snp_markers(),
            Format::Fasta => match self.test_format_fasta() {
                Some(answer) => answer,
                None => {
                    self.fasta_undecided = true;
                    false
                }
            },
            Format::TextAsn => self.test_format_text_asn(),
            Format::Alignment => self.test_format_alignment(),
            Format::Taxplot => false,
            Format::Table => self.test_format_table(),
            Format::BinaryAsn => self.test_format_binary_asn(),
            Format::PhrapAce => self.test_format_phrap_ace()?,
            Format::UcscRegion | Format::UnknownOrJson => false,
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:625-684
    /// `EnsureTestBuffer()`: false for an empty stream (which stays at its end) and for an all-comment sample that
    /// ended before the buffer was full. Each pass reads the stream again from where the last one put its bytes
    /// back (`CStreamUtils::Stepback`), so the sample is the next 8096 bytes times `Multiplier`.
    pub(crate) fn ensure_test_buffer(&mut self) -> bool {
        if self.buffer.is_some() {
            return true;
        }
        if !self.stream.good {
            return false;
        }
        const TEST_BUFFER_GRANULARITY: usize = 8096;
        let mut multiplier: usize = 1;
        loop {
            self.buffer_size = multiplier * TEST_BUFFER_GRANULARITY;
            let read = self.stream.read(self.buffer_size);
            if read.is_empty() {
                self.buffer = None;
                self.buffer_size = 0;
                return false; // empty file
            }
            self.stream.clear(); // in case we reached eof
            self.stream.stepback(&read);
            let data_size = read.len();
            self.buffer = Some(read);
            self.terminated = false;
            if self.is_all_comment() {
                if multiplier >= 1024 {
                    return true;
                }
                multiplier *= 2;
                self.buffer = None;
                if data_size < self.buffer_size {
                    return false;
                }
                continue;
            }
            return true;
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3747-3764
    fn is_ascii_text(&self) -> bool {
        const REQUIRED_ASCII_RATIO: f64 = 0.9;
        let data = self.data();
        let count = data.len();
        let count_print = data.iter().filter(|&&c| c_isprint(c)).count();
        (count_print as f64) >= (count as f64) * REQUIRED_ASCII_RATIO
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3767-3792
    /// Lines are not trimmed: " # x" is not a comment. When `EnsureSplitLines` gives no line (it returned false
    /// early), the loop does not run and the sample counts as all comment.
    fn is_all_comment(&mut self) -> bool {
        if !self.is_ascii_text() {
            return false;
        }
        self.split_done = false;
        self.lines.clear();
        self.ensure_split_lines();
        self.lines()
            .all(|line| line.is_empty() || line.starts_with(b"#") || line.starts_with(b"--"))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3692-3744
    /// `EnsureSplitLines()`: an early `return false` leaves `m_TestLines` as they were (empty after
    /// `IsAllComment`'s `clear()`, and at the start).
    pub(crate) fn ensure_split_lines(&mut self) -> bool {
        if self.split_done {
            return !self.lines.is_empty();
        }
        self.split_done = true;
        const MIN_HIGH_RATIO: usize = 20;
        let Some(range) = self.buffer.clone() else {
            self.lines.clear();
            return false;
        };
        let base = range.start;
        let data = &self.stream.source.bytes()[range];
        let high_count = data.iter().filter(|&&c| 0x80 & c != 0).count();
        if 0 < high_count && data.len() / high_count < MIN_HIGH_RATIO {
            return false;
        }
        self.lines.clear();
        // `NStr::Split` of `data` at `delims` into `m_TestLines` with `NStr::fSplit_Tokenize`: the non-empty
        // runs between delimiters.
        let ends = if data.contains(&b'\r') && find(data, b"\r\n").is_some() {
            LineEnds::CrLf
        } else if data.contains(&b'\n') {
            LineEnds::Lf
        } else if data.contains(&b'\r') {
            LineEnds::Cr
        } else if data.len() == self.buffer_size {
            // most likely single truncated line
            return false;
        } else {
            // test buffer contains the entire file
            LineEnds::Whole
        };
        let mut count = match ends {
            LineEnds::Whole => 1,
            LineEnds::CrLf => count_runs(data, |c| c == b'\r' || c == b'\n'),
            LineEnds::Lf => count_runs(data, |c| c == b'\n'),
            LineEnds::Cr => count_runs(data, |c| c == b'\r'),
        };
        if data.len() == self.buffer_size && count > 1 {
            // multiple lines, last likely truncated
            count -= 1;
        }
        self.lines = Lines {
            range: base..base + data.len(),
            ends,
            count,
            trimmed: 0,
            extra: None,
        };
        !self.lines.is_empty()
    }

    /// `EnsureStats()`: the statistics are counted once, on the test buffer of the first call.
    pub(crate) fn ensure_stats(&mut self) -> bool {
        if self.stats.is_some() {
            return true;
        }
        if !self.ensure_test_buffer() {
            return false;
        }
        self.stats = Some(stats_of(self.data()));
        true
    }

    /// `EnsureTestBuffer()` and then `EnsureSplitLines()`, the first check of most tests.
    fn ensure_lines(&mut self) -> bool {
        self.ensure_test_buffer() && self.ensure_split_lines()
    }

    /// `EnsureStats()` and then `EnsureSplitLines()`.
    fn ensure_stats_and_lines(&mut self) -> bool {
        self.ensure_stats() && self.ensure_split_lines()
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1928-1950
    /// `TestFormatZip`.
    fn test_format_zip(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        data.len() >= 4
            && data.starts_with(b"PK")
            && matches!((data[2], data[3]), (1, 2) | (3, 4) | (5, 6) | (7, 8))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1953-1991 (`TestFormatGZip`, `TestFormatZstd`)
    /// `TestFormatGZip`, `TestFormatZstd`: the buffer starts with `magic`.
    fn test_magic(&mut self, magic: &[u8]) -> bool {
        self.ensure_test_buffer() && self.data().starts_with(magic)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1994-2016
    /// `TestFormatBZip2`.
    fn test_format_bzip2(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        data.len() >= 4 && data.starts_with(b"BZh") && (b'1'..=b'9').contains(&data[3])
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2019-2044
    /// `TestFormatLzo`.
    fn test_format_lzo(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        (data.starts_with(b"LZO") && data.get(3).is_none_or(|&c| c == 0))
            || (data.len() >= 4 && &data[1..4] == b"LZO" && data.get(4).is_none_or(|&c| c == 0))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2047-2063
    /// `TestFormatSra`.
    fn test_format_sra(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        data.len() >= 16
            && data.starts_with(b"NCBI.sra")
            && matches!(&data[8..12], b"\x05\x03\x19\x88" | b"\x88\x19\x03\x05")
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1277-1313
    /// `TestFormatXml`. `eTrunc_Begin` drops C's `isspace` bytes, the vertical tab included.
    fn test_format_xml(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let input = truncate_spaces_begin(self.data());
        let starts_nocase = |prefix: &[u8]| {
            input.len() >= prefix.len() && input[..prefix.len()].eq_ignore_ascii_case(prefix)
        };
        starts_nocase(b"<?XML")
            || starts_nocase(b"<!DOCTYPE")
            || input.starts_with(b"<Blast4-request>")
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1161-1186
    /// `TestFormatBinaryAsn`: a byte that is neither `isgraph` nor `isspace`, other than `\1` (`conf` never becomes
    /// `eYes`).
    pub(crate) fn test_format_binary_asn(&mut self) -> bool {
        self.ensure_test_buffer()
            && self
                .data()
                .iter()
                .any(|&c| !c_isgraph(c) && !is_space(c) && c != 1)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1678-1716
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3680-3689
    /// `TestFormatTextAsn`. The statistics are those of `EnsureStats`'s buffer, the size and the lines those of the
    /// buffer now (another one after `TestFormatNewick` replaced it). `NcbiGetline` with `\n\r` as the delimiters (cited at
    /// [`stats_of`]) ends a line at either byte, and the empty lines of a pair have no field: comments.
    pub(crate) fn test_format_text_asn(&mut self) -> bool {
        if !self.ensure_stats() {
            return false;
        }
        let data = self.data();
        if data.is_empty() || data[0] == b'>' {
            return false;
        }
        let stats = self.stats();
        if ((stats.al_num_chars + stats.braces) as f64) / (data.len() as f64) < 0.80 {
            return false;
        }
        for line in data.split(|&c| c == b'\n' || c == b'\r') {
            let fields = tokens(line, b" \t");
            if fields.is_empty() || fields.at(0).starts_with(b"--") {
                continue;
            }
            return fields.len() >= 2
                && fields.at(1) == b"::="
                && fields.at(0)[0].is_ascii_alphabetic();
        }
        false
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1640-1675
    // NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:340-347 (`CTempString::operator[]`, cited at
    // `skip_comment_and_blank`: NUL past the end)
    /// `TestFormatFasta`: `Some(answer)`, or `None` when its `NStr::Find` of `|` in `m_pTestBuffer` decides and the input
    /// does not. `m_pTestBuffer` is a `char*` there: the C string ends at its first NUL, which `EnsureTestBuffer`'s
    /// buffer does not hold after its bytes. The first `|` of the bytes decides when it comes at position 10 or before
    /// (and before a NUL), and so does a NUL; a buffer of 11 bytes or more without either gives a position past 10.
    /// Only a shorter buffer without a NUL or `|` leaves the answer to the memory after it (makeblastdb's
    /// `appears_fasta` takes that for false).
    pub(crate) fn test_format_fasta(&mut self) -> Option<bool> {
        if !self.ensure_stats() {
            return Some(false);
        }
        let data = self.data();
        let stats = self.stats();
        // reject obvious misfits:
        let header = skip_comment_and_blank(data);
        if data.is_empty() || header.is_empty() || header[0] != b'>' {
            return Some(false);
        }
        if stats.data == 0 {
            if 0.75 > (stats.al_num_chars as f64) / (data.len() as f64) {
                return Some(false);
            }
            let searched = c_string(data);
            if let Some(bar) = searched.iter().position(|&c| c == b'|') {
                if bar <= 10 {
                    return Some(true);
                }
            }
            let decided = self.terminated || searched.len() < data.len() || data.len() > 10;
            return decided.then_some(false);
        }

        // remaining decision based on text stats:
        let al_num_fraction = stats.al_num_chars as f64 / data.len() as f64;
        let dna_fraction = stats.dna_chars as f64 / stats.data as f64;
        let aa_fraction = stats.aa_chars as f64 / stats.data as f64;

        // want at least 80% text-ish overall:
        if al_num_fraction < 0.8 {
            return Some(false);
        }

        // want more than 91 percent of either DNA content or AA content in what we
        // presume is data:
        Some(dna_fraction > 0.91 || aa_fraction > 0.91)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1061-1158
    /// `TestFormatNewick`. When the split fails it reads a buffer of its own (8096 bytes from where the stream is,
    /// a NUL after them) and adds it to the lines as a C string, which the later tests read. The NEXUS scan reads
    /// the stream on without stepping back: the text reader then starts where it stopped. The `CTempString` of `test_buf` is
    /// the C string of the chunk with the 12 bytes `strncpy` carried before it (a NUL among them hides the chunk);
    /// a short chunk is the last (`eofbit`), so its overlapping `strncpy` is never read.
    fn test_format_newick(&mut self) -> bool {
        if !self.ensure_test_buffer() || !self.ensure_split_lines() {
            const BUFFSIZE: usize = 8096;
            let read = self.stream.read(BUFFSIZE);
            self.stream.clear(); // in case we reached eof
            self.stream.stepback(&read);
            self.buffer = Some(read);
            self.terminated = true;
            let line = self.c_string_line();
            self.push_line(line);
        }

        let is_nexus =
            self.lines_may_hold(b"#NEXUS") && self.lines().any(|line| contains(line, b"#NEXUS"));
        if is_nexus {
            const READ_SIZE: usize = 16384;
            const CHECK_SIZE: usize = 12;
            let mut test_buf = vec![b' '; CHECK_SIZE + READ_SIZE + 1];
            const MAX_READS: usize = 32768; // max read to locate tree: 512 MB
            for _ in 0..MAX_READS {
                let read = self.stream.read(READ_SIZE);
                let num_read = read.len();
                if num_read > 0 {
                    test_buf[CHECK_SIZE..CHECK_SIZE + num_read]
                        .copy_from_slice(&self.stream.source.bytes()[read]);
                    test_buf[CHECK_SIZE + num_read] = 0; // null terminator
                    if find_nocase(c_string(&test_buf), b"begin trees;") {
                        self.stream.clear(); // in case we reached eof
                        return true;
                    }
                    // `strncpy` of `check_size` bytes from offset `num_read` of `test_buf` to its
                    // start: up to a NUL, the rest NUL.
                    let carried = c_string(&test_buf[num_read..num_read + CHECK_SIZE]).len();
                    test_buf.copy_within(num_read..num_read + carried, 0);
                    test_buf[carried..CHECK_SIZE].fill(0);
                }
                if !self.stream.good {
                    self.stream.clear(); // clear eof
                    break;
                }
            }
            return false;
        }

        const MAX_SAMPLE_SIZE: usize = 8 * 1024 - 1;
        let read = self.stream.read(MAX_SAMPLE_SIZE);
        self.stream.clear(); // in case we reached eof
        self.stream.stepback(&read);
        if read.is_empty() {
            return false;
        }
        is_sample_newick(&self.stream.source.bytes()[read])
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1624-1637
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:340-347
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:537-549
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:362-373
/// `SkipCommentAndBlank`. `text[0]` of an empty text is NUL, which is not in `;#!`. A comment line without a following
/// line end (`find_first_of` gives `npos`) makes `substr` from `npos` the empty string.
fn skip_comment_and_blank(mut text: &[u8]) -> &[u8] {
    loop {
        text = truncate_spaces_begin(text);
        let first = text.first().copied().unwrap_or(0);
        if matches!(first, b';' | b'#' | b'!') {
            // `find_first_of` of `NEW_LINE_SYMBOLS` from index 1
            let pos = text
                .iter()
                .skip(1)
                .position(|&c| c == b'\r' || c == b'\n')
                .map(|p| p + 1);
            text = match pos {
                Some(p) => &text[p..],
                None => &[],
            };
        } else {
            break;
        }
    }
    text
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2989-2999 (`NStr::FindNoCase`: the first letter in both
// cases, then the pattern compared without case)
/// `NStr::FindNoCase` of `pattern` in `str` is not `NPOS` (ASCII case).
fn find_nocase(str: &[u8], pattern: &[u8]) -> bool {
    str.len() >= pattern.len()
        && str
            .windows(pattern.len())
            .any(|window| window.eq_ignore_ascii_case(pattern))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2946-3071
/// `IsSampleNewick` of `pSample`, the sample read (`cline` is its C string, up to its first NUL; a sample whose first
/// byte after the spaces is not `(` is no tree before the C string is cut).
fn is_sample_newick(sample: &[u8]) -> bool {
    if truncate_spaces_begin(sample).first() != Some(&b'(') {
        return false;
    }
    let line = truncate_spaces(c_string(sample));
    if line.first() != Some(&b'(') {
        return false;
    }
    // Strip out comments:
    let mut trimmed = Vec::with_capacity(line.len());
    let mut in_comment = false;
    for &c in line {
        if !in_comment {
            if c != b'[' {
                trimmed.push(c);
            } else {
                in_comment = true;
            }
        } else if c == b']' {
            in_comment = false;
        }
    }
    // Compress quoted labels:
    let mut quoted = Vec::with_capacity(trimmed.len());
    let mut in_quote = false;
    for &c in &trimmed {
        if !in_quote {
            if c != b'\'' {
                quoted.push(c);
            } else {
                in_quote = true;
                quoted.push(b'A');
            }
        } else if c == b'\'' {
            in_quote = false;
        }
    }
    // Strip distance markers (`c_str()[ii]` is NUL at the end):
    let at = |ii: usize| quoted.get(ii).copied().unwrap_or(0);
    let mut stripped = Vec::with_capacity(quoted.len());
    let mut ii = 0;
    while at(ii) != 0 {
        if at(ii) != b':' {
            stripped.push(at(ii));
            ii += 1;
        } else {
            ii += 1;
            if at(ii) == b'-' || at(ii) == b'+' {
                ii += 1;
            }
            while at(ii).is_ascii_digit() {
                ii += 1;
            }
            if at(ii) == b'.' {
                ii += 1;
                while at(ii).is_ascii_digit() {
                    ii += 1;
                }
            }
        }
    }
    // Rough lexical analysis of what's left. Bail immediately on fault:
    if stripped.first() != Some(&b'(') {
        return false;
    }
    let mut paren_count: usize = 1;
    for &c in &stripped[1..] {
        match c {
            b'(' => paren_count += 1,
            b')' => {
                if paren_count == 0 {
                    return false;
                }
                paren_count -= 1;
            }
            b',' if paren_count == 0 => return false,
            _ => {}
        }
    }
    true
}

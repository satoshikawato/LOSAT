//! The objects of a serial stream read one after another (`CObjectIStreamAsn` over an `istream`),
//! and the choice of the stream's encoding that NCBI's readers of a stream make before they read
//! (`CFormatGuess`'s `Format` of the stream): blast_formatter's `-archive` (Blast4-archive objects until the
//! end of the data, `CRemoteBlast::LoadFromArchive`) and the search strategy of
//! `-import_search_strategy` (one Blast4-get-search-strategy-reply, or else one Blast4-request,
//! `ExtractBlast4Request`).
//!
//! The guess is the one port of `CFormatGuess` (`utils::format_guess`); the text reader is the
//! generic reader of the type table (`read`) on a stream read in blocks (`StreamSource`), so one
//! object at a time is in memory. Binary ASN.1 and XML are recognised and left to their callers:
//! a front end that reads them into the table's values would take the input from
//! [`Guessed::into_parts`].

use std::io::{Read as IoRead, Seek, SeekFrom};

use super::read::{read_object, Read};
use super::Value;
use crate::ncbi_exception::NcbiException;
use crate::report::serial::asn_text_reader::{AsnSource, Reader, StreamSource};
use crate::utils::format_guess::{Format, FormatGuess, Stop};

/// The encoding of a serial stream, as the readers of NCBI's archives and strategies take
/// `CFormatGuess`'s answer.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:796-818
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Encoding {
    /// `eTextASN`: read by [`ObjectStream`].
    Text,
    /// `eBinaryASN`: recognised; LOSAT has no binary reader of the table's types yet.
    Binary,
    /// `eXml`: recognised; LOSAT has no XML reader of the table's types yet.
    Xml,
    /// Any other answer (`eUnknown`, `eJSON`, another format): the caller's invalid input.
    Other,
}

/// A guess that gives no answer.
#[derive(Debug)]
pub(crate) enum GuessStop {
    /// A test of the guess throws (`TestFormatGlimmer3`'s `NStr::StringToInt`): NCBI's caller
    /// gets the exception.
    Exception(NcbiException),
    /// A test of the guess reads memory its input does not decide, so NCBI's answer is not
    /// known: named here.
    Undecided(&'static str),
}

/// The encoding of a stream and the stream itself, its first bytes taken by the guess.
pub(crate) struct Guessed<R> {
    pub(crate) encoding: Result<Encoding, GuessStop>,
    /// The bytes the guess took from the input.
    taken: Vec<u8>,
    /// Where the guess left the stream in `taken` (its start, or past what `TestFormatNewick`
    /// read of a NEXUS text).
    start: usize,
    input: R,
}

/// `CFormatGuess`'s `Format` of the input (the whole order of `sm_CheckOrder`, no hints), its answer
/// taken as [`Encoding`]. `TestFormatFasta`, which comes before `eTextASN` and `eBinaryASN` in
/// the order, may search past its buffer (`Guessed::fasta_undecided`): its answer would be
/// `eFasta`, an invalid input, so a text or binary answer after it is not decided.
pub(crate) fn guess_encoding<R: IoRead>(input: R) -> Guessed<R> {
    let mut guess = FormatGuess::new(input);
    let answer = guess.guess_format();
    let (taken, input) = guess.into_parts();
    let (encoding, start) = match answer {
        Ok(guessed) => {
            let encoding = match guessed.format {
                Format::TextAsn | Format::BinaryAsn if guessed.fasta_undecided => {
                    Err(GuessStop::Undecided(
                        "TestFormatFasta searches for '|' past a buffer of 10 bytes or fewer",
                    ))
                }
                Format::TextAsn => Ok(Encoding::Text),
                Format::BinaryAsn => Ok(Encoding::Binary),
                Format::Xml => Ok(Encoding::Xml),
                _ => Ok(Encoding::Other),
            };
            (encoding, guessed.stream_pos)
        }
        Err(Stop::Exception(exception)) => (Err(GuessStop::Exception(exception)), 0),
        Err(Stop::Undecided(what)) => (Err(GuessStop::Undecided(what)), 0),
    };
    Guessed {
        encoding,
        start: start.min(taken.len()),
        taken,
        input,
    }
}

impl<R: IoRead> Guessed<R> {
    /// The text reader of the stream (a new `CObjectIStreamAsn` of the file), from where the guess left it.
    pub(crate) fn into_text(self) -> ObjectStream<R> {
        ObjectStream::new(self.input, &self.taken[self.start..])
    }

    /// The stream from where the guess left it: the bytes the guess took from there, and the
    /// input after them (the front end of another encoding reads from here).
    pub(crate) fn into_parts(mut self) -> (Vec<u8>, R) {
        self.taken.drain(..self.start);
        (self.taken, self.input)
    }
}

/// The objects of an ASN.1 text stream, read one at a time: one reader for the whole stream, as
/// NCBI keeps one object stream (its line count and its buffer go on from one object to the
/// next).
pub(crate) struct ObjectStream<R: IoRead> {
    r: Reader<StreamSource<R>>,
    /// A read failed: NCBI's stream has its fail flags set, and its `EndOfData` is true from
    /// then on.
    failed: bool,
}

impl<R: IoRead> ObjectStream<R> {
    /// The objects of `input`, whose first bytes `first` were already taken from it.
    pub(crate) fn new(input: R, first: &[u8]) -> Self {
        ObjectStream {
            r: Reader::new(
                StreamSource::new(input, first),
                "ASN.1 text",
                "ASN.1 text reader",
            ),
            failed: false,
        }
    }

    /// `EndOfData()`: whether no object follows. True after a failed read and at the end of the
    /// input; otherwise the white space and comments are skipped, and the end of the input in
    /// them (or an unclosed comment) is the end too. Anything else, garbage included, is the start
    /// of an object for the next read.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:98-109
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:588-597
    // (`m_Input.EndOfData()` and `HasMore()` of `CIStreamBuffer` ask whether a byte follows,
    // reading the stream for one: the end of the input here.)
    pub(crate) fn end_of_data(&mut self) -> bool {
        if self.failed || self.r.src.byte_at(self.r.pos).is_none() {
            return true;
        }
        self.r.skip_white_space().is_err()
    }

    /// `in >> object` of the root type `root`: one object, from its file header (`Read`), with
    /// the messages posted. The stream stands after the object's last byte; a failed read leaves
    /// it failed (`end_of_data`).
    pub(crate) fn read(&mut self, root: u16) -> Read<Value> {
        let read = read_object(&mut self.r, root);
        if read.0.is_err() {
            self.failed = true;
        }
        read
    }

    /// The text after what was read: the bytes the stream holds from there, then the input
    /// after them.
    pub(crate) fn into_rest(self) -> (Vec<u8>, R) {
        let pos = self.r.pos;
        self.r.src.into_rest(pos)
    }

    /// The bytes the stream keeps.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.r.src.buffered()
    }
}

impl<R: IoRead + Seek> ObjectStream<R> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2302
    /// The seek to the start that follows a failed read of a strategy reply: the input from its
    /// start, for a new guess and read; the error of an input that cannot seek (a pipe), whose
    /// stream NCBI's `seekg` leaves failed.
    pub(crate) fn rewind(self) -> std::io::Result<R> {
        let (_, mut input) = self.into_rest();
        input.seek(SeekFrom::Start(0))?;
        Ok(input)
    }
}

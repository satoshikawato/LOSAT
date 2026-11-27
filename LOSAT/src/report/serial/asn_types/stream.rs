//! The objects of a serial stream read one after another (`CObjectIStreamAsn` over an `istream`),
//! and the choice of the stream's encoding that NCBI's readers of a stream make before they read
//! (`CFormatGuess`'s `Format` of the stream): blast_formatter's `-archive` (Blast4-archive objects until the
//! end of the data, `CRemoteBlast::LoadFromArchive`) and the search strategy of
//! `-import_search_strategy` (one Blast4-get-search-strategy-reply, or else one Blast4-request,
//! `ExtractBlast4Request`).
//!
//! The guess is the one port of `CFormatGuess` (`utils::format_guess`); each encoding is read by the
//! generic reader of the type table (`read::ObjectStream`) over its stream primitives (text,
//! binary, NCBI's ASN.1-style XML) on a stream read in blocks (`StreamSource`), so one object at a
//! time is in memory. A caller that does not read an encoding takes the input from
//! [`Guessed::into_parts`].

use std::io::{Read as IoRead, Seek, SeekFrom};

use super::read::ObjectStream;
use super::read_xml::{XmlInput, XmlStyle};
use crate::ncbi_exception::NcbiException;
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{Names, TextInput};
use crate::report::serial::asn_text_reader::StreamSource;
use crate::utils::format_guess::{Format, FormatGuess, Stop};

/// The encoding of a serial stream, as the readers of NCBI's archives and strategies take
/// `CFormatGuess`'s answer.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:796-818
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Encoding {
    /// `eTextASN`: read by [`ObjectStream`].
    Text,
    /// `eBinaryASN`: read by [`Guessed::into_binary`], or refused by the caller.
    Binary,
    /// `eXml`: read by [`Guessed::into_xml`], or refused by the caller.
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
    /// The text reader of the stream (a new `CObjectIStreamAsn` of the file), from where the
    /// guess left it, whose refusals call the text by `names`.
    pub(crate) fn into_text(self, names: Names) -> TextStream<R> {
        text_stream(self.input, &self.taken[self.start..], names)
    }

    /// The binary reader of the stream (a new `CObjectIStreamAsnBinary` of the file), from where
    /// the guess left it, whose refusals call the input by `names`.
    pub(crate) fn into_binary(self, names: Names) -> BinaryStream<R> {
        let (first, input) = self.into_parts();
        ObjectStream::binary(StreamSource::new(input, &first), names)
    }

    /// The XML reader of the stream (`CObjectIStream::Open` with `eSerial_Xml`: ASN.1-style XML,
    /// `XmlStyle::Asn`), from where the guess left it, whose refusals call the input by `names`.
    pub(crate) fn into_xml(self, names: Names) -> XmlStream<R> {
        let (first, input) = self.into_parts();
        ObjectStream::xml(StreamSource::new(input, &first), XmlStyle::Asn, names)
    }

    /// The stream from where the guess left it: the bytes the guess took from there, and the
    /// input after them (the front end of another encoding reads from here).
    pub(crate) fn into_parts(mut self) -> (Vec<u8>, R) {
        self.taken.drain(..self.start);
        (self.taken, self.input)
    }
}

/// The objects of an ASN.1 binary stream read in blocks, one at a time ([`TextStream`]'s binary
/// twin).
pub(crate) type BinaryStream<R> = ObjectStream<BerReader<StreamSource<R>>>;

/// The objects of an ASN.1-style XML stream read in blocks, one at a time ([`TextStream`]'s XML
/// twin).
pub(crate) type XmlStream<R> = ObjectStream<XmlInput<StreamSource<R>>>;

/// The objects of an ASN.1 text stream read in blocks, one at a time: one reader for the whole
/// stream (`read::ObjectStream`), as NCBI keeps one object stream (its buffer goes on from one
/// object to the next).
pub(crate) type TextStream<R> = ObjectStream<TextInput<StreamSource<R>>>;

/// The objects of `input`, whose first bytes `first` were already taken from it, which LOSAT's
/// refusals call by `names`.
pub(crate) fn text_stream<R: IoRead>(input: R, first: &[u8], names: Names) -> TextStream<R> {
    ObjectStream::text(StreamSource::new(input, first), names)
}

impl<R: IoRead> TextStream<R> {
    /// The text after what was read: the bytes the stream holds from there, then the input
    /// after them.
    pub(crate) fn into_rest(self) -> (Vec<u8>, R) {
        let r = self.into_input().r;
        let pos = r.pos;
        r.src.into_rest(pos)
    }

    /// The bytes the stream keeps.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.input().r.src.buffered()
    }
}

impl<R: IoRead + Seek> TextStream<R> {
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

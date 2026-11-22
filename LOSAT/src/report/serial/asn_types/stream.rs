//! The objects of a serial stream read one after another (`CObjectIStreamAsn`,
//! `CObjectIStreamAsnBinary` or `CObjectIStreamXml` over an `istream`), and the choice of the
//! stream's encoding that NCBI's readers of a stream make before they read (`CFormatGuess`'s
//! `Format` of the stream): blast_formatter's `-archive` (Blast4-archive objects until the
//! end of the data, `CRemoteBlast::LoadFromArchive`) and the search strategy of
//! `-import_search_strategy` (one Blast4-get-search-strategy-reply, or else one Blast4-request,
//! `ExtractBlast4Request`).
//!
//! The guess is the one port of `CFormatGuess` (`utils::format_guess`); each encoding is read by
//! the generic reader of the type table (`read`) over that encoding's front end (`read_text`,
//! `read_ber`, `read_xml`) on a stream read in blocks (`StreamSource`), so one object at a time
//! is in memory.

use std::io::{Read as IoRead, Seek, SeekFrom};

use super::read::{read_object, Read};
use super::read_ber::Ber;
use super::read_text::Text;
use super::read_xml::Xml;
pub(crate) use super::read_xml::XmlStyle;
use super::Value;
use crate::ncbi_exception::NcbiException;
use crate::report::serial::asn_text_reader::{Reader, StreamSource};
use crate::utils::format_guess::{Format, FormatGuess, Stop};

/// The encoding of a serial stream, as the readers of NCBI's archives and strategies take
/// `CFormatGuess`'s answer.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:796-818
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Encoding {
    /// `eTextASN`: `CObjectIStreamAsn`.
    Text,
    /// `eBinaryASN`: `CObjectIStreamAsnBinary`.
    Binary,
    /// `eXml`: the stream `CObjectIStream::Open` makes for the format `eSerial_Xml`,
    /// `CObjectIStreamXml`.
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

    /// The object stream of the encoding guessed (a new `CObjectIStreamAsn` or
    /// `CObjectIStreamAsnBinary` over the input, or the stream `CObjectIStream::Open` makes for
    /// `eSerial_Xml`: the ASN.1-style XML, `SetEnforcedStdXml` not called), from where the guess
    /// left the input; `None` for an encoding that is not one of them (the caller's invalid
    /// input).
    pub(crate) fn into_stream(self) -> Option<ObjectStream<R>> {
        let first = &self.taken[self.start..];
        match self.encoding {
            Ok(Encoding::Text) => Some(ObjectStream::new(self.input, first)),
            Ok(Encoding::Binary) => Some(ObjectStream::binary(self.input, first)),
            Ok(Encoding::Xml) => Some(ObjectStream::xml(self.input, first, XmlStyle::Asn)),
            _ => None,
        }
    }

    /// The stream from where the guess left it: the bytes the guess took from there, and the
    /// input after them.
    pub(crate) fn into_parts(mut self) -> (Vec<u8>, R) {
        self.taken.drain(..self.start);
        (self.taken, self.input)
    }
}

/// The objects of a serial stream, read one at a time: one reader for the whole stream, as NCBI
/// keeps one object stream (its position, line count and buffer go on from one object to the
/// next).
pub(crate) struct ObjectStream<R: IoRead> {
    front: Fronts<R>,
    /// A read failed: NCBI's stream has its fail flags set, and its `EndOfData` is true from
    /// then on.
    failed: bool,
}

/// The front end of the stream's encoding.
enum Fronts<R: IoRead> {
    Text(Text<StreamSource<R>>),
    Binary(Ber<StreamSource<R>>),
    Xml(Xml<StreamSource<R>>),
}

impl<R: IoRead> ObjectStream<R> {
    /// The ASN.1 text objects of `input`, whose first bytes `first` were already taken from it.
    pub(crate) fn new(input: R, first: &[u8]) -> Self {
        ObjectStream {
            front: Fronts::Text(Text::new(Reader::new(
                StreamSource::new(input, first),
                "ASN.1 text",
                "ASN.1 text reader",
            ))),
            failed: false,
        }
    }

    /// The binary ASN.1 objects of `input`, whose first bytes `first` were already taken from it.
    pub(crate) fn binary(input: R, first: &[u8]) -> Self {
        ObjectStream {
            front: Fronts::Binary(Ber::new(StreamSource::new(input, first))),
            failed: false,
        }
    }

    /// The XML objects of `input` in the `style`, whose first bytes `first` were already taken
    /// from it.
    pub(crate) fn xml(input: R, first: &[u8], style: XmlStyle) -> Self {
        ObjectStream {
            front: Fronts::Xml(Xml::new(StreamSource::new(input, first), style)),
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
    // The binary stream keeps `EndOfData` of its base: any byte after an object starts the next.
    // The XML stream's (`CObjectIStreamXml::EndOfData`, cited at `Xml::end_of_data`) skips white
    // space, comments and processing instructions as the text stream does, and anything that
    // fails there is the end too.
    #[cfg(test)]
    pub(crate) fn end_of_data(&mut self) -> bool {
        use crate::report::serial::asn_text_reader::AsnSource;
        if self.failed {
            return true;
        }
        match &mut self.front {
            Fronts::Text(text) => {
                let r = &mut text.r;
                r.src.byte_at(r.pos).is_none() || r.skip_white_space().is_err()
            }
            Fronts::Binary(ber) => ber.src.byte_at(ber.pos).is_none(),
            Fronts::Xml(xml) => xml.end_of_data(),
        }
    }

    /// One object of the root type `root` extracted from the stream (`Read`: its file header,
    /// then the value), with the messages posted. The stream stands after the object's last
    /// byte; a failed read leaves it failed (`end_of_data`).
    pub(crate) fn read(&mut self, root: u16) -> Read<Value> {
        let read = match &mut self.front {
            Fronts::Text(text) => read_object(text, root),
            Fronts::Binary(ber) => read_object(ber, root),
            Fronts::Xml(xml) => read_object(xml, root),
        };
        if read.0.is_err() {
            self.failed = true;
        }
        read
    }

    /// The input after what was read: the bytes the stream holds from there, then the input
    /// after them.
    pub(crate) fn into_rest(self) -> (Vec<u8>, R) {
        match self.front {
            Fronts::Text(text) => {
                let pos = text.r.pos;
                text.r.src.into_rest(pos)
            }
            Fronts::Binary(ber) => ber.src.into_rest(ber.pos),
            Fronts::Xml(xml) => xml.x.src.into_rest(xml.x.pos),
        }
    }

    /// The bytes the stream keeps.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        match &self.front {
            Fronts::Text(text) => text.r.src.buffered(),
            Fronts::Binary(ber) => ber.src.buffered(),
            Fronts::Xml(xml) => xml.x.src.buffered(),
        }
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

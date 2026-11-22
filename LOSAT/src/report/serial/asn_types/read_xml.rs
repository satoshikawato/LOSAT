//! The XML front end of the generic reader (`read`): NCBI's `CObjectIStreamXml` as
//! `CObjectIStream::Open` with the format `eSerial_Xml` makes it for blast_formatter's `-archive`
//! (`SetEnforcedStdXml` is not called): the ASN.1-style XML that `ARCHIVE_FORMAT=xml` writes,
//! where each type of an ASN.1 module is a tag (`<Seq-id>`), each member of a SEQUENCE or
//! variant of a CHOICE is a tag named by its parent (`<Seq-id_local>`), the elements of a list of
//! anonymous types are `<Parent_E>`, and an ENUMERATED or BOOLEAN is the attribute `value="..."`
//! of its tag.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2069-2093
//! (`CheckStdXml`: every class of the table is `EDataSpec::eASN` of `CodeVersion` 22400, so
//! `m_StdXml` stays false and `m_SkipNextTag` (set only by an alias with `SetFullAlias` or in
//! standard XML) stays false; no member of an ASN.1 module has `notag`, an attribute list,
//! `nillable` or an `any` content, so the `Notag` and attribute branches of NCBI's functions are
//! not ported either.)
//!
//! The lexical layer is `xml_reader::XmlReader` (the one XML lexer of LOSAT, in its object mode);
//! this module is the grammar of NCBI's stream over the walk's steps: the names of the tags
//! from the frames of the object stack (`SkipStackTagName`), the members and variants by those
//! names, lists, named types, enumerated values and booleans as attributes. The style of the
//! stream (`XmlStyle`) is an input: `Asn` is the archive's, and `Std` (`SetEnforcedStdXml` set,
//! the search strategy's) is left to its own reader, where `Xml::file_header` refuses it.
//!
//! Where NCBI's stream throws a `CSerialException` the result is `ReadError::Serial`; the end of
//! the input is `ReadError::Eof` (a `CEofException`: the walk's `HandleEOF` turns it into a
//! `CSerialException` inside two frames or more); the number reads throw NCBI's `CUtilException`
//! (`ReadError::Exception`).

use super::read::{read_delayed_in, Cx, Frame, Front, ReadError, R};
use super::{type_info, EnumValue, Member, Storage, StringKind, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::is_space;
use crate::report::serial::asn_text_reader::{AsnSource, SliceSource, BAD_NUMBER, NUMBER_OVERFLOW};
use crate::report::serial::encoder::StreamPosition;
use crate::report::serial::object_frame::FrameKind;
use crate::report::serial::xml_reader::{TagState, XmlError, XmlReader, EOF_MESSAGE};

/// The style of an XML stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XmlStyle {
    /// ASN.1-style XML: `CObjectIStream::Open` with `eSerial_Xml` of blast_formatter's archive
    /// (`m_StdXml` false for every type of an ASN.1 module).
    Asn,
    /// Standard XML (`SetEnforcedStdXml` set, remote_blast.cpp:2278-2283: the search strategy
    /// of `-import_search_strategy`): not read here. The place for it is `Xml::file_header` and
    /// the members of `Front` below, which an `Std` stream would answer by the rules of
    /// `x_IsStdXml()`; until then it refuses the input (`ReadError::Unsupported`).
    Std,
}

impl From<XmlError> for ReadError {
    fn from(error: XmlError) -> Self {
        match error {
            XmlError::Format(message) if message == EOF_MESSAGE => ReadError::Eof,
            XmlError::Format(message) => ReadError::Serial(message),
            XmlError::Unsupported(message) => ReadError::Unsupported(message),
        }
    }
}

/// The XML front end: the lexical reader of one XML stream of documents.
pub(crate) struct Xml<S: AsnSource> {
    pub(crate) x: XmlReader<S>,
    style: XmlStyle,
}

impl<S: AsnSource> Xml<S> {
    pub(crate) fn new(src: S, style: XmlStyle) -> Self {
        Xml {
            x: XmlReader::object_stream(src),
            style,
        }
    }

    /// `EndOfData()`: whether no object follows: the end of the input, or white space, comments
    /// and processing instructions up to it, or a failure in them.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:118-129
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:588-597
    // (`CObjectIStream::EndOfData`: the end of the input; a failed stream is
    // `ObjectStream::failed`.)
    pub(crate) fn end_of_data(&mut self) -> bool {
        self.x.src.byte_at(self.x.pos).is_none() || self.x.skip_ws_and_comments().is_err()
    }

    /// The bytes replaced by `#` so far, posted at the line they were read on.
    fn drain(&mut self, cx: &mut Cx<'_>) {
        for (byte, line) in self.x.replaced.drain(..) {
            cx.post_bad_char(byte, &[], StreamPosition::Line(line as u64));
        }
    }

    /// A step's result with the messages it posted before it ended.
    fn posted<T>(&mut self, cx: &mut Cx<'_>, result: Result<T, XmlError>) -> R<T> {
        self.drain(cx);
        Ok(result?)
    }

    /// `OpenTag`: the opening tag `e` (the name given back by `UndoClassMember`, or read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1563-1573
    fn open_tag(&mut self, e: &str) -> R<()> {
        self.read_tag_name()?;
        if self.x.last_tag != e.as_bytes() {
            return Err(ReadError::Serial(format!(
                "tag '{e}' expected: {}",
                String::from_utf8_lossy(&self.x.last_tag)
            )));
        }
        Ok(())
    }

    /// The name of the next opening tag, in `last_tag`: the one given back, or read
    /// (`ReadName` of what `BeginOpeningTag` finds, or `RejectedName`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:525-535
    fn read_tag_name(&mut self) -> R<()> {
        match self.x.rejected.take() {
            Some(name) => {
                self.x.last_tag = name;
                self.x.state = TagState::InsideOpening;
            }
            None => {
                let c = self.x.begin_opening_tag()?;
                self.x.read_name(c)?;
            }
        }
        Ok(())
    }

    /// `CloseTag`: nothing after a self-closed opening tag, else the closing tag `e`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1575-1586
    fn close_tag(&mut self, e: &str) -> R<()> {
        if self.x.state == TagState::SelfClosed {
            self.x.state = TagState::Outside;
            return Ok(());
        }
        let c = self.x.begin_closing_tag()?;
        self.x.read_name(c)?;
        if self.x.last_tag != e.as_bytes() {
            return Err(ReadError::Serial(format!(
                "tag '{e}' expected: {}",
                String::from_utf8_lossy(&self.x.last_tag)
            )));
        }
        Ok(self.x.end_tag()?)
    }

    /// `OpenStackTag` at level 0 of an element of a list: the opening tag named by the stack.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1588-1602
    fn open_stack_tag(&mut self, frames: &[Frame]) -> R<()> {
        if self.x.rejected.is_none() {
            let c = self.x.begin_opening_tag()?;
            self.x.read_name(c)?;
            self.check_stack_tag(frames)
        } else {
            self.read_tag_name()
        }
    }

    /// `CloseStackTag` at level 0: nothing after a self-closed opening tag, else the closing tag
    /// named by the stack.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1604-1623
    fn close_stack_tag(&mut self, frames: &[Frame]) -> R<()> {
        if self.x.state == TagState::SelfClosed {
            self.x.state = TagState::Outside;
            return Ok(());
        }
        let c = self.x.begin_closing_tag()?;
        self.x.read_name(c)?;
        self.check_stack_tag(frames)?;
        Ok(self.x.end_tag()?)
    }

    /// The tag in `last_tag` must be the one the stack names: `SkipStackTagName` of it from
    /// level 0 leaves nothing.
    fn check_stack_tag(&self, frames: &[Frame]) -> R<()> {
        let tag = &self.x.last_tag;
        let rest = skip_stack_tag_name(frames, tag, 0)?;
        if !rest.is_empty() {
            return Err(ReadError::Serial(format!(
                "unexpected tag: {}{}",
                String::from_utf8_lossy(tag),
                String::from_utf8_lossy(rest)
            )));
        }
        Ok(())
    }

    fn util_error(&self, site: ThrowSite, what: &str) -> ReadError {
        ReadError::Exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.x.line),
        ))
    }

    /// `GetInt4`, `GetUint4` (`signed` false: no `-`) and `GetInt8` of `CIStreamBuffer`, the
    /// digits checked against `max` as NCBI checks them.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-669
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:672-700
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:703-746
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        let mut c = loop {
            let c = self.x.get_char()?;
            if !is_space(c) {
                break c;
            }
        };
        let mut negative = false;
        if signed && c == b'-' {
            negative = true;
            c = self.x.get_char()?;
        } else if c == b'+' {
            c = self.x.get_char()?;
        }
        if !c.is_ascii_digit() {
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        loop {
            let d = self.x.peek_no_eof(0);
            if !d.is_ascii_digit() {
                break;
            }
            self.x.skip(1);
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
            }
            n = n * 10 + d;
        }
        // `-Int4(n)` of n = 2^31 is the least Int4 (and likewise in 64 bits).
        Ok(if negative {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    /// `ReadInt4`: `BeginData`, `GetInt4`.
    fn int4(&mut self) -> R<i32> {
        self.x.begin_data()?;
        Ok(self.get_number(i32::MAX as u64, true)? as i32)
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2648-2675
    // (`SkipUNumber`, 2677-2703, takes only `+`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        self.x.begin_data()?;
        let mut c = self.x.skip_ws_and_comments()?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = self.x.peek(1)?;
            i = 2;
        }
        if !c.is_ascii_digit() {
            return Err(ReadError::Serial("invalid symbol in number".into()));
        }
        while self.x.peek_no_eof(i).is_ascii_digit() {
            i += 1;
        }
        self.x.skip(i);
        Ok(())
    }

    /// `ReadWord`: the word in the tag, with the messages posted.
    fn read_word(&mut self, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        let mut word = Vec::new();
        let result = self.x.read_word(&mut word);
        self.posted(cx, result)?;
        Ok(word)
    }
}

/// `SkipTagName`: `tag` after the prefix `name`.
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1508-1515
fn skip_tag_name<'t>(tag: &'t [u8], name: &[u8]) -> R<&'t [u8]> {
    if !tag.starts_with(name) {
        return Err(ReadError::Serial(format!(
            "invalid tag name: {}",
            String::from_utf8_lossy(tag)
        )));
    }
    Ok(&tag[name.len()..])
}

/// `SkipStackTagName`: what is left of `tag` after the name that the frame `level`
/// frames below the top of the stack gives it: a type's name, a member's tag (its parent's name,
/// `_` and its own), an element's (its parent's, `_E`).
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1517-1561
fn skip_stack_tag_name<'t>(frames: &[Frame], tag: &'t [u8], level: usize) -> R<&'t [u8]> {
    let Some(frame) = frames.len().checked_sub(level + 1).map(|at| &frames[at]) else {
        return Err(ReadError::Serial("illegal frame type".into()));
    };
    match frame.kind() {
        FrameKind::Named
        | FrameKind::NamedType
        | FrameKind::Array
        | FrameKind::Class
        | FrameKind::Choice => {
            let name = frame.type_name();
            if name.is_empty() {
                skip_stack_tag_name(frames, tag, level + 1)
            } else {
                skip_tag_name(tag, name.as_bytes())
            }
        }
        FrameKind::ClassMember | FrameKind::ChoiceVariant => {
            let tag = skip_stack_tag_char(frames, tag, level + 1, b'_')?;
            skip_tag_name(tag, frame.member().as_bytes())
        }
        FrameKind::ArrayElement => {
            if frames.len() > level + 1 {
                let tag = skip_stack_tag_name(frames, tag, level + 1)?;
                skip_tag_name(tag, b"_E")
            } else {
                Ok(&tag[..0])
            }
        }
    }
}

fn skip_stack_tag_char<'t>(frames: &[Frame], tag: &'t [u8], level: usize, c: u8) -> R<&'t [u8]> {
    let tag = skip_stack_tag_name(frames, tag, level)?;
    match tag.split_first() {
        Some((&first, rest)) if first == c => Ok(rest),
        _ => Err(ReadError::Serial(format!(
            "invalid tag name: {}",
            String::from_utf8_lossy(tag)
        ))),
    }
}

/// `UnexpectedMember`.
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2105-2114
fn unexpected_member(id: &[u8], items: &[Member]) -> ReadError {
    let mut message = format!(
        "\"{}\": unexpected member, should be one of: ",
        String::from_utf8_lossy(id)
    );
    for item in items {
        message.push('"');
        message.push_str(item.name());
        message.push_str("\" ");
    }
    ReadError::Serial(message)
}

/// `NStr::TruncateSpacesInPlace`: no white space at either end.
fn truncate_spaces(text: &[u8]) -> &[u8] {
    let start = text
        .iter()
        .position(|&c| !is_space(c))
        .unwrap_or(text.len());
    let end = text
        .iter()
        .rposition(|&c| !is_space(c))
        .map_or(start, |at| at + 1);
    &text[start..end]
}

/// `CEnumeratedTypeValues::FindValue`: the value of the name, or of the name with its first
/// letter in upper case.
// NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:124-138
fn find_value(values: &'static [EnumValue], name: &[u8]) -> R<i32> {
    let mut upper = name.to_vec();
    if let Some(first) = upper.first_mut() {
        *first = first.to_ascii_uppercase();
    }
    values
        .iter()
        .find(|value| value.name().as_bytes() == name)
        .or_else(|| {
            values
                .iter()
                .find(|value| value.name().as_bytes() == upper.as_slice())
        })
        .map(|value| value.value)
        .ok_or_else(|| {
            ReadError::Serial(format!(
                "invalid value of enumerated type: {}",
                String::from_utf8_lossy(name)
            ))
        })
}

impl<S: AsnSource> Front for Xml<S> {
    const NAMED_TYPES: bool = true;
    const ARRAY_FRAMES: bool = false;

    /// `WillHaveName`: the element of a list has the frame of an element when its type has no
    /// name.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1638-1648
    fn element_frame(element: u16) -> bool {
        type_info(element).name().is_empty()
    }

    // The lines are those of `END_OBJECT_FRAME()` in the functions that push the frames:
    // `CObjectIStream::Read` (objistr.cpp:959-969), `CObjectIStream::ReadClassSequential`
    // (objistr.cpp:1398-1428: the member frame ends in `ReadClassSequentialContentsEnd` at 1424),
    // `CObjectIStream::ReadChoiceSimple` (objistr.cpp:1482-1501),
    // `CObjectIStreamXml::ReadContainerContents` (objistrxml.cpp:1921-1979), `ReadNamedType`
    // (objistrxml.cpp:2041-2052) and the skip functions `SkipClassSequential` (objistr.cpp:
    // 1451-1470, the member frame at 1466), `SkipChoiceSimple` (objistr.cpp:1503-1522),
    // `SkipContainerContents` (objistrxml.cpp:1981-2000) and `SkipNamedType` (objistr.cpp:
    // 1302-1313). The container of an ASN.1 module's list has no name: no array frame
    // (`ReadContainer`, objistrxml.cpp:1860-1875, for a container with one).
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
        const XML: &str = "src/serial/objistrxml.cpp";
        const ISTR: &str = "src/serial/objistr.cpp";
        match (kind, skip) {
            (FrameKind::Named, _) => (ISTR, 968, "ncbi::CObjectIStream::Read"),
            (FrameKind::NamedType, false) => (XML, 2051, "ncbi::CObjectIStreamXml::ReadNamedType"),
            (FrameKind::NamedType, true) => (ISTR, 1311, "ncbi::CObjectIStream::SkipNamedType"),
            (FrameKind::Class, false) => (ISTR, 1427, "ncbi::CObjectIStream::ReadClassSequential"),
            (FrameKind::ClassMember, false) => {
                (ISTR, 1424, "ncbi::CObjectIStream::ReadClassSequential")
            }
            (FrameKind::Choice, false) => (ISTR, 1500, "ncbi::CObjectIStream::ReadChoiceSimple"),
            (FrameKind::ChoiceVariant, false) => {
                (ISTR, 1498, "ncbi::CObjectIStream::ReadChoiceSimple")
            }
            (FrameKind::Array, false) => (XML, 1873, "ncbi::CObjectIStreamXml::ReadContainer"),
            (FrameKind::ArrayElement, false) => {
                (XML, 1950, "ncbi::CObjectIStreamXml::ReadContainerContents")
            }
            (FrameKind::Class, true) => (ISTR, 1469, "ncbi::CObjectIStream::SkipClassSequential"),
            (FrameKind::ClassMember, true) => {
                (ISTR, 1466, "ncbi::CObjectIStream::SkipClassSequential")
            }
            (FrameKind::Choice, true) => (ISTR, 1521, "ncbi::CObjectIStream::SkipChoiceSimple"),
            (FrameKind::ChoiceVariant, true) => {
                (ISTR, 1519, "ncbi::CObjectIStream::SkipChoiceSimple")
            }
            (FrameKind::Array, true) => (XML, 1889, "ncbi::CObjectIStreamXml::SkipContainer"),
            (FrameKind::ArrayElement, true) => {
                (XML, 1993, "ncbi::CObjectIStreamXml::SkipContainerContents")
            }
        }
    }

    fn start_object(&mut self) {
        self.x.held = None;
        self.x.src.release_before(self.x.pos);
    }

    /// `ResetState` and `ReadFileHeader`: the state of the stream at its root, the byte order
    /// mark, the declaration, the `DOCTYPE`, and the name of the root's tag (given back).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:84-101
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:582-670
    // (`ResetState` returns at once when more than one frame is on the stack; at the root it
    // clears the tag state and the tag names kept. The `m_Encoding` of `ReadFileHeader` is the
    // one of the stream so far; a declaration with an `encoding` sets it. `Read` of a root
    // object starts with `SkipFileHeader`, whose first call is `ResetState`, objistr.cpp:883-900.)
    fn file_header(&mut self, _cx: &mut Cx<'_>) -> R<Vec<u8>> {
        if self.style == XmlStyle::Std {
            return Err(ReadError::Unsupported(
                "standard XML (SetEnforcedStdXml) is not read by the XML front end of the archive"
                    .into(),
            ));
        }
        self.x.reset_state();
        Ok(self.x.read_file_header()?)
    }

    /// `BeginClass`: the tag of a class with a name.
    fn begin_class(&mut self, ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        let name = type_info(ty).name();
        if !name.is_empty() {
            self.open_tag(name)?;
        }
        Ok(())
    }

    /// `BeginClassMember` of a SEQUENCE of an ASN.1 module: the member whose tag,
    /// `Class_member`, follows; the end of the members at the closing tag. A self-closed tag, a
    /// closing tag, a position past the last member or no tag ahead also ends the members.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2191-2367
    // (`Find` by the tag name and `FindDeep` by the tag name from `pos`, before the last lines,
    // look for a member whose name is the whole tag, or one of a member's class: a tag of an
    // ASN.1-style document is `Class_member`, which neither names; `FindDeep` is not ported.
    // Unknown members are not skipped, `SERIAL_SKIP_UNKNOWN_MEMBERS` being its registry
    // parameter's default no.)
    fn next_member(
        &mut self,
        members: &'static [Member],
        pos: usize,
        _block: &mut bool,
        cx: &mut Cx<'_>,
    ) -> R<Option<usize>> {
        if self.x.rejected.is_none() {
            if pos == 0 && self.x.has_attlist()? {
                self.x.read_undefined_attributes()?;
            }
            if self.x.state == TagState::SelfClosed || self.x.this_tag_is_self_closed()? {
                return Ok(None);
            }
            if self.x.next_tag_is_closing()? || pos >= members.len() || !self.x.next_is_tag()? {
                return Ok(None);
            }
        }
        self.read_tag_name()?;
        let tag = self.x.last_tag.clone();
        let id = skip_stack_tag_char(cx.frames, &tag, 1, b'_')?;
        match members[pos..]
            .iter()
            .position(|member| member.name().as_bytes() == id)
        {
            Some(offset) => Ok(Some(pos + offset)),
            None => Err(unexpected_member(id, members)),
        }
    }

    /// `EndClassMember`: the member's closing tag.
    fn end_member(&mut self, cx: &mut Cx<'_>) -> R<()> {
        self.close_stack_tag(cx.frames)
    }

    /// `EndClass`: the closing tag of a class with a name.
    fn end_class(&mut self, cx: &mut Cx<'_>) -> R<()> {
        let name = cx.frames.last().map_or("", |frame| frame.type_name());
        if !name.is_empty() {
            self.close_tag(name)?;
        }
        Ok(())
    }

    /// `BeginChoice`: the tag of a choice with a name.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2391-2409
    fn begin_choice(&mut self, ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        let name = type_info(ty).name();
        if !name.is_empty() {
            self.open_tag(name)?;
        }
        Ok(())
    }

    fn end_choice(&mut self, cx: &mut Cx<'_>) -> R<()> {
        let name = cx.frames.last().map_or("", |frame| frame.type_name());
        if !name.is_empty() {
            self.close_tag(name)?;
        }
        Ok(())
    }

    /// `BeginChoiceVariant` of a CHOICE of an ASN.1 module: the variant whose tag,
    /// `Choice_variant`, follows. A choice without a variant tag is the "choice variant id
    /// expected" of `ReadChoiceSimple`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2411-2507
    // (`FindDeep` and `Find` by the tag name come before: see `next_member`. No CHOICE of the
    // table may be empty (`MayBeEmpty`) or has a variant without a tag (`FindEmpty`).)
    fn begin_variant(
        &mut self,
        _ty: u16,
        variants: &'static [Member],
        cx: &mut Cx<'_>,
    ) -> R<usize> {
        if self.x.rejected.is_none() {
            if self.x.has_attlist()? {
                self.x.read_undefined_attributes()?;
            }
            if self.x.state == TagState::SelfClosed || self.x.next_tag_is_closing()? {
                return Err(ReadError::Serial("choice variant id expected".into()));
            }
            // `NextIsTag` looks for a variant without a tag when no tag follows (none in an
            // ASN.1 module); then `ReadName` of `BeginOpeningTag` throws for what is there.
            self.x.next_is_tag()?;
        }
        self.read_tag_name()?;
        let tag = self.x.last_tag.clone();
        let id = skip_stack_tag_char(cx.frames, &tag, 1, b'_')?;
        match variants
            .iter()
            .position(|variant| variant.name().as_bytes() == id)
        {
            Some(index) => Ok(index),
            None => Err(unexpected_member(&tag, variants)),
        }
    }

    /// `EndChoiceVariant`: the variant's closing tag.
    fn end_variant(&mut self, cx: &mut Cx<'_>) -> R<()> {
        self.close_stack_tag(cx.frames)
    }

    /// `BeginContainer` is not called (`ReadContainer` reads the container's contents).
    fn begin_list(&mut self, _set: bool, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }

    /// `HasMoreElements`, then `BeginArrayElement` for an element without a name: the next
    /// element's opening tag. The end of the input in `HasMoreElements` is no more elements (the
    /// closing tag after the list meets it again).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1765-1825
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1895-1910
    fn next_element(&mut self, _element: u16, _block: &mut bool, cx: &mut Cx<'_>) -> R<bool> {
        let no_more = match self
            .x
            .this_tag_is_self_closed()
            .and_then(|closed| Ok(closed || self.x.next_tag_is_closing()?))
        {
            Ok(no_more) => no_more,
            Err(XmlError::Format(message)) if message == EOF_MESSAGE => true,
            Err(error) => return Err(error.into()),
        };
        if no_more {
            return Ok(false);
        }
        if cx
            .frames
            .last()
            .is_some_and(|frame| frame.kind() == FrameKind::ArrayElement)
        {
            self.open_stack_tag(cx.frames)?;
        }
        Ok(true)
    }

    /// `EndArrayElement`: the closing tag of an element that had an opening tag of its own.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1912-1919
    fn end_element(&mut self, cx: &mut Cx<'_>) -> R<()> {
        if cx
            .frames
            .last()
            .is_some_and(|frame| frame.kind() == FrameKind::ArrayElement)
        {
            self.close_stack_tag(cx.frames)?;
        }
        Ok(())
    }

    /// `EndContainer` is not called (`ReadContainerContents`); an empty implicit class is not an
    /// error (no class of the table is `SetImplicitNonEmpty`).
    fn end_list(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }

    /// `BeginNamedType`: the tag of the named type (`<IUPACna>`, `<Seq-descr>`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2002-2037
    fn begin_named(&mut self, ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        self.open_tag(type_info(ty).name())
    }

    fn end_named(&mut self, cx: &mut Cx<'_>) -> R<()> {
        let name = cx.frames.last().map_or("", |frame| frame.type_name());
        self.close_tag(name)
    }

    /// `ReadInt4`, `ReadUint4`, `ReadInt8`, `ReadInt2` by the member's storage, or their skip.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:991-1025
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
    // (`ReadInt2` reads an Int4 that must fit: "integer overflow".)
    fn integer(
        &mut self,
        _ty: u16,
        storage: Storage,
        skip: bool,
        _cx: &mut Cx<'_>,
    ) -> R<Option<i64>> {
        if skip {
            self.skip_number(storage != Storage::Uint4)?;
            return Ok(None);
        }
        let value = match storage {
            Storage::Int4 => i64::from(self.int4()?),
            Storage::Uint4 => {
                self.x.begin_data()?;
                self.get_number(u64::from(u32::MAX), false)?
            }
            Storage::Int8 => {
                self.x.begin_data()?;
                self.get_number(i64::MAX as u64, true)?
            }
            Storage::Int2 => {
                let value = self.int4()?;
                if i16::try_from(value).is_err() {
                    return Err(ReadError::Serial("integer overflow".into()));
                }
                i64::from(value)
            }
        };
        Ok(Some(value))
    }

    /// `ReadEnum` of the stream: the attribute `value` of the tag (of the type's own tag when it
    /// has a name), or, for an INTEGER with named values, the number after the tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1376-1453
    fn enumerated(
        &mut self,
        ty: u16,
        values: &'static [EnumValue],
        integer: bool,
        _skip: bool,
        cx: &mut Cx<'_>,
    ) -> R<i32> {
        let name = type_info(ty).name();
        if !name.is_empty() {
            self.open_tag(name)?;
        }
        let result = self.enum_value(values, integer);
        self.drain(cx);
        let value = result?;
        if !name.is_empty() {
            self.close_tag(name)?;
        }
        Ok(value)
    }

    /// `ReadBool` / `SkipBool`: the attribute `value` of the tag, or the word in it; the tag
    /// must be empty after it.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:933-977
    fn boolean(&mut self, _ty: u16, _skip: bool, cx: &mut Cx<'_>) -> R<bool> {
        let mut text = Vec::new();
        let mut have = false;
        let result = (|| -> Result<(), XmlError> {
            while self.x.has_attlist()? {
                let c = self.x.skip_ws()?;
                self.x.read_name(c)?;
                let value = self.x.read_attribute_value()?;
                if self.x.last_tag == b"value" {
                    text = value;
                    have = true;
                }
            }
            Ok(())
        })();
        self.drain(cx);
        result?;
        if !have {
            text = self.read_word(cx)?;
        }
        let text = truncate_spaces(&text);
        let value = match text {
            b"true" | b"1" => true,
            b"false" | b"0" => false,
            _ => {
                return Err(ReadError::Serial(format!(
                    "'true' or 'false' value expected: {}",
                    String::from_utf8_lossy(text)
                )))
            }
        };
        if !self.x.end_opening_tag_self_closed()? && !self.x.next_tag_is_closing()? {
            return Err(ReadError::Serial(
                "boolean tag must have empty contents".into(),
            ));
        }
        Ok(value)
    }

    /// `ReadNull`, `SkipNull`: an empty tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1044-1048
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2725-2729
    fn null(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<()> {
        let empty =
            self.x.end_opening_tag_self_closed()? || (!skip && self.x.next_tag_is_closing()?);
        if !empty {
            return Err(ReadError::Serial("empty tag expected".into()));
        }
        Ok(())
    }

    /// `ReadDouble` (the skip, `SkipFNumber`, reads the number too): the word in the tag.
    fn real(&mut self, _ty: u16, skip: bool, cx: &mut Cx<'_>) -> R<Option<f64>> {
        let result = self.x.real();
        let value = self.posted(cx, result)?;
        Ok((!skip).then_some(value))
    }

    /// `ReadString` (a self-closed tag is the empty string) with the bytes replaced by `#` and
    /// posted (`x_VerifyChar`), or `SkipString`: the text up to the next tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2710-2723
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1233-1258
    // (`ReadString`: a self-closed tag is the empty string.)
    fn string(&mut self, _ty: u16, _kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        let mut out = Vec::new();
        if skip {
            self.x.begin_data()?;
            while self.x.read_escaped_char(b'<')?.0 >= 0 {}
            return Ok(out);
        }
        let result = (|| -> Result<(), XmlError> {
            if self.x.end_opening_tag_self_closed()? {
                return Ok(());
            }
            self.x.read_tag_data(&mut out)
        })();
        self.posted(cx, result)?;
        Ok(out)
    }

    /// `ReadBytes` of an OCTET STRING: pairs of hexadecimal digits up to the next tag (an odd
    /// last digit is the high half of a byte); `SkipByteBlock`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2524-2542
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2596-2613
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2731-2762
    fn octets(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<Option<Vec<u8>>> {
        self.x.begin_data()?;
        if skip {
            loop {
                let c = self.x.peek(0)?;
                match c {
                    b'0'..=b'9' | b'A'..=b'Z' | b'a'..=b'z' | b'+' | b'/' | b'=' => {
                        self.x.skip(1);
                    }
                    b'\r' | b'\n' => {
                        self.x.skip(1);
                        self.x.skip_end_of_line(c);
                    }
                    b'<' => return Ok(None),
                    _ => return Err(ReadError::Serial("invalid char in octet string".into())),
                }
            }
        }
        let mut bytes = Vec::new();
        loop {
            let Some(high) = self.hex_char()? else {
                return Ok(Some(bytes));
            };
            let Some(low) = self.hex_char()? else {
                bytes.push((high << 4) as u8);
                return Ok(Some(bytes));
            };
            bytes.push(((high << 4) | low) as u8);
        }
    }

    /// `StartDelayBuffer`: `BeginData`, the buffer starts here (the source keeps the bytes from
    /// here, `held`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1492-1500
    // (A member's own tag is read before its value: the rejected tag is the first of a member
    // read by its inner class's name, which no delayed member of the table is.)
    fn delay_start(&mut self, _cx: &mut Cx<'_>) -> R<usize> {
        self.x.begin_data()?;
        if self.x.rejected.is_some() {
            return Err(ReadError::Unsupported(
                "a delayed member that starts with a tag given back is not read from XML".into(),
            ));
        }
        self.x.held = Some(self.x.pos);
        Ok(self.x.pos)
    }

    /// `EndDelayBuffer`, and the buffer parsed by a new XML stream over its bytes
    /// (`CDelayBuffer::DoUpdate`: `UpdateDelayedBuffer` reads the member's type under the frames
    /// of its class and of the member, as `read_delayed_in` does): LOSAT's value of a delayed
    /// member is its parsed value, as the binary stream's.
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:490-501
    fn delay_end(&mut self, start: usize, ty: u16, cx: &mut Cx<'_>) -> R<Value> {
        self.x.held = None;
        let bytes = self.x.src.bytes(start..self.x.pos).to_vec();
        let mut buffer = Xml::new(SliceSource::new(&bytes), self.style);
        let frames = cx.frames;
        let outer = &frames[frames.len().saturating_sub(2)..];
        read_delayed_in(&mut buffer, outer, ty, cx)
    }
}

impl<S: AsnSource> Xml<S> {
    /// `ReadEnum`'s value, after the type's own tag was opened (or not).
    fn enum_value(&mut self, values: &'static [EnumValue], integer: bool) -> R<i32> {
        if self.x.state != TagState::InsideOpening {
            // outside of tag
            if !integer {
                return Err(ReadError::Serial("attribute 'value' expected".into()));
            }
            return self.int4();
        }
        let c = self.x.skip_ws()?;
        if matches!(c, b'>' | b'/') {
            // no attribute
            if !integer {
                return Err(ReadError::Serial("attribute 'value' expected".into()));
            }
            self.x.skip(1);
            self.x.state = TagState::Outside;
            return self.int4();
        }
        let mut found = false;
        while self.x.has_attlist()? {
            let c = self.x.skip_ws()?;
            self.x.read_name(c)?;
            if self.x.last_tag == b"value" {
                found = true;
                break;
            }
            self.x.read_attribute_value()?;
        }
        if !found {
            self.x.end_opening_tag_self_closed()?;
            return Err(ReadError::Serial("attribute 'value' is missing".into()));
        }
        let name = self.x.read_attribute_value()?;
        let value = find_value(values, truncate_spaces(&name))?;
        if !self.x.end_opening_tag_self_closed()? && integer {
            // read integer value
            self.x.skip_ws_and_comments()?;
            if value != self.get_number(i32::MAX as u64, true)? as i32 {
                return Err(ReadError::Serial(
                    "incompatible name and value of named integer".into(),
                ));
            }
        }
        Ok(value)
    }

    /// `GetHexChar`: the value of a digit or letter, `None` at the `<` that ends the data.
    fn hex_char(&mut self) -> R<Option<u32>> {
        let c = self.x.get_char()?;
        match c {
            b'0'..=b'9' => Ok(Some(u32::from(c - b'0'))),
            b'A'..=b'Z' => Ok(Some(u32::from(c - b'A') + 10)),
            b'a'..=b'z' => Ok(Some(u32::from(c - b'a') + 10)),
            _ => {
                self.x.skip_back(1);
                if c != b'<' {
                    return Err(ReadError::Serial("invalid char in octet string".into()));
                }
                Ok(None)
            }
        }
    }
}

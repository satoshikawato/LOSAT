//! The primitives of NCBI's XML object input stream (`CObjectIStreamXml`) as the object stack of
//! `read` calls them (`Input`): the ASN.1-style XML of `CObjectIStream::Open` with `eSerial_Xml`
//! (psiblast's and tblastn's `-in_pssm` of an XML PSSM, read with the XML manipulator;
//! blast_formatter's archive), where each type of an ASN.1 module is a tag (`<Seq-id>`), each member of a SEQUENCE
//! or variant of a CHOICE is a tag named by its parent (`<Seq-id_local>`), the elements of a list
//! of anonymous types are `<Parent_E>`, and an ENUMERATED or BOOLEAN is the attribute
//! `value="..."` of its tag.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2069-2093
//! (`CheckStdXml`: every class of the table is `EDataSpec::eASN` of `CodeVersion` 22400, so
//! `m_StdXml` stays false and `m_SkipNextTag` (set only by an alias with `SetFullAlias` or in
//! standard XML) stays false; no member of an ASN.1 module has `notag`, an attribute list,
//! `nillable` or an `any` content, so the `Notag` and attribute branches of NCBI's functions are
//! not ported either, except the `SetNotag` that `FindEmpty` gives the variant of an empty choice
//! element, `XmlInput::begin_variant`.)
//!
//! The lexical layer is `xml_reader::XmlReader` (the one XML lexer of LOSAT, in its object mode);
//! this module is the grammar of NCBI's stream over the object stack's calls: the names of the
//! tags from the frames of the stack (`SkipStackTagName`, over the frames the stack reports to
//! its input, `Input::frame_pushed`), the members and variants by those names, lists, named
//! types, enumerated values and booleans as attributes. The XML stream's own frames (a named
//! frame for each named type that is not a SEQUENCE, CHOICE or ENUMERATED, no array frame, an
//! element frame only for an element without a name) are the stack's, by this `Input`'s consts.
//! The style of the stream (`XmlStyle`) is an input: `Asn` is the one read here, and `Std`
//! (`SetEnforcedStdXml` set, the search strategy's) is refused at the file header.
//!
//! What the stream posts while it reads is posted at once with the stack trace of that moment
//! (`x_VerifyChar` of each control byte it replaces by `#`, also in the attributes it skips, with
//! the line and an empty string): `drain_posts` gives the posts to the stack in order. Where
//! NCBI's stream throws a `CSerialException` the result is `Fail::Serial` (NCBI's programs that
//! read XML report none of its text: `-in_pssm`'s "Unrecognized format", blast_formatter's
//! "Invalid input format"); the end of the input is `Fail::Eof` (a `CEofException`: the stack's
//! `HandleEOF` turns it into a `CSerialException` inside two frames or more); the numbers throw
//! the `CUtilException` of NCBI's input buffer (`Fail::Exception`), which the stack's frames
//! report line by line.

use super::read::{Post, NO_TYPE};
use super::{type_info, EnumValue, Kind, Member};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::is_space;
use crate::report::serial::asn_input::{
    invalid_enum, EnumItems, Fail, Input, Items, Names, StringKind, StringPost, R,
};
use crate::report::serial::asn_text_reader::{
    AsnSource, NumberStop, PointerType, BAD_NUMBER, FIND_VALUE, NUMBER_OVERFLOW,
};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::report::serial::object_frame::FrameKind;
use crate::report::serial::xml_reader::{TagState, XmlError, XmlReader, EOF_MESSAGE};

/// The style of an XML stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XmlStyle {
    /// ASN.1-style XML: `CObjectIStream::Open` with `eSerial_Xml` (`m_StdXml` false for every type
    /// of an ASN.1 module).
    Asn,
    /// Standard XML (`SetEnforcedStdXml` set, remote_blast.cpp:2278-2283: the search strategy of
    /// `-import_search_strategy`): not read here. The place for it is `read_file_header` and the
    /// functions below, which a `Std` stream would answer by the rules of `x_IsStdXml()`; until
    /// then it refuses the input.
    #[allow(dead_code)]
    Std,
}

/// A frame of the object stack as the XML stream names its tags by it (`Input::frame_pushed`).
#[derive(Clone, Copy, Debug)]
struct XmlFrame {
    kind: FrameKind,
    /// The type of a named, class or choice frame, the element type of an element frame.
    ty: u16,
    /// The member or variant of a member or variant frame (`SetTopMemberId`).
    member: Option<&'static Member>,
    /// The variant has no tag of its own (`SetNotag` of its frame: `FindEmpty`'s variant of an
    /// empty choice element), so `EndChoiceVariant` reads no closing tag.
    notag: bool,
}

impl XmlFrame {
    /// The name of the frame's type (empty for an anonymous one and for the frames that name
    /// none).
    fn type_name(&self) -> &'static str {
        match self.kind {
            FrameKind::Named | FrameKind::NamedType | FrameKind::Class | FrameKind::Choice
                if self.ty != NO_TYPE =>
            {
                type_info(self.ty).name()
            }
            _ => "",
        }
    }

    fn member_name(&self) -> &'static str {
        self.member.map_or("", Member::name)
    }
}

/// The XML stream over one source: the lexer, the frames of the stack, what it posted.
pub(crate) struct XmlInput<S: AsnSource> {
    x: XmlReader<S>,
    style: XmlStyle,
    /// What LOSAT's refusals call the document and the reader.
    names: Names,
    frames: Vec<XmlFrame>,
    posts: Vec<Post>,
}

/// A `ThrowError` of `CObjectIStreamXml` at `line` of `function` (`eFormatError` unless `kind`).
const fn xml_site(line: u32, kind: &'static str, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/serial/objistrxml.cpp",
        line,
        kind,
        function: Some(function),
    }
}

const FORMAT_ERROR: &str = "CSerialException::eFormatError";
/// `CSerialException::eMissingValue` (`fMissingValue`): `ReadWithSetFlagMember` takes it for an
/// OPTIONAL or DEFAULT member (`read::TypedReader`).
pub(crate) const MISSING_VALUE: &str = "CSerialException::eMissingValue";

// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1376-1453
const READ_ENUM: &str = "ncbi::CObjectIStreamXml::ReadEnum";
const ATTRIBUTE_VALUE_EXPECTED: ThrowSite = xml_site(1401, FORMAT_ERROR, READ_ENUM);
const ATTRIBUTE_VALUE_EXPECTED_OUTSIDE: ThrowSite = xml_site(1444, FORMAT_ERROR, READ_ENUM);
const ATTRIBUTE_VALUE_MISSING: ThrowSite = xml_site(1425, MISSING_VALUE, READ_ENUM);
const INCOMPATIBLE_NAMED_INTEGER: ThrowSite =
    xml_site(1435, "CSerialException::eInvalidData", READ_ENUM);
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1508-1623
const SKIP_TAG_NAME: ThrowSite =
    xml_site(1513, FORMAT_ERROR, "ncbi::CObjectIStreamXml::SkipTagName");
const SKIP_STACK_TAG_NAME: ThrowSite = xml_site(
    1550,
    "CSerialException::eIllegalCall",
    "ncbi::CObjectIStreamXml::SkipStackTagName",
);
const SKIP_STACK_TAG_CHAR: ThrowSite = xml_site(
    1559,
    FORMAT_ERROR,
    "ncbi::CObjectIStreamXml::SkipStackTagName",
);
const OPEN_TAG: ThrowSite = xml_site(1572, FORMAT_ERROR, "ncbi::CObjectIStreamXml::OpenTag");
const CLOSE_TAG: ThrowSite = xml_site(1583, FORMAT_ERROR, "ncbi::CObjectIStreamXml::CloseTag");
const OPEN_STACK_TAG: ThrowSite =
    xml_site(1596, FORMAT_ERROR, "ncbi::CObjectIStreamXml::OpenStackTag");
const CLOSE_STACK_TAG: ThrowSite =
    xml_site(1617, FORMAT_ERROR, "ncbi::CObjectIStreamXml::CloseStackTag");
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2105-2114
const UNEXPECTED_MEMBER: ThrowSite = xml_site(
    2113,
    FORMAT_ERROR,
    "ncbi::CObjectIStreamXml::UnexpectedMember",
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1501 (`ReadChoiceSimple`: "choice
// variant id expected" when `BeginChoiceVariant` finds none)
const CHOICE_VARIANT_EXPECTED: ThrowSite = ThrowSite {
    file: "src/serial/objistr.cpp",
    line: 1488,
    kind: FORMAT_ERROR,
    function: Some("ncbi::CObjectIStream::ReadChoiceSimple"),
};
const READ_BOOL: &str = "ncbi::CObjectIStreamXml::ReadBool";
const BOOL_VALUE_EXPECTED: ThrowSite = xml_site(969, FORMAT_ERROR, READ_BOOL);
const BOOL_NOT_EMPTY: ThrowSite = xml_site(975, FORMAT_ERROR, READ_BOOL);
const READ_NULL: ThrowSite = xml_site(1047, FORMAT_ERROR, "ncbi::CObjectIStreamXml::ReadNull");
const SKIP_NULL: ThrowSite = xml_site(2728, FORMAT_ERROR, "ncbi::CObjectIStreamXml::SkipNull");
const READ_DOUBLE: ThrowSite = xml_site(1040, FORMAT_ERROR, "ncbi::CObjectIStreamXml::ReadDouble");
const SKIP_S_NUMBER: ThrowSite =
    xml_site(2669, FORMAT_ERROR, "ncbi::CObjectIStreamXml::SkipSNumber");
const SKIP_U_NUMBER: ThrowSite =
    xml_site(2697, FORMAT_ERROR, "ncbi::CObjectIStreamXml::SkipUNumber");
const READ_BYTES: ThrowSite = xml_site(2539, FORMAT_ERROR, "ncbi::CObjectIStreamXml::ReadBytes");
const SKIP_BYTE_BLOCK: ThrowSite =
    xml_site(2759, FORMAT_ERROR, "ncbi::CObjectIStreamXml::SkipByteBlock");
const READ_OBJECT_POINTER: ThrowSite = xml_site(
    1473,
    "CSerialException::eNotImplemented",
    "ncbi::CObjectIStreamXml::ReadObjectPointer",
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704 (`ReadInt2`, cited at
// `object_frame::READ_INT2`, whose text is the same)

/// The `ThrowError` of the lexer's failure `message` (its text without what NCBI's message adds
/// after it): the function of NCBI's stream that throws it. The text is not reported by the
/// programs that read XML (see the module comment); the site is kept for the tests.
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:243-440
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:545-800
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:884-912
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1027-1042
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1271-1295
fn lexer_site(message: &str) -> ThrowSite {
    macro_rules! at {
        ($line:literal, $function:literal) => {
            xml_site(
                $line,
                FORMAT_ERROR,
                concat!("ncbi::CObjectIStreamXml::", $function),
            )
        };
    }
    match message {
        "'>' expected" => at!(335, "EndTag"),
        "end of tag expected" => at!(370, "EndOpeningTagSelfClosed"),
        "'<' expected" => at!(411, "BeginOpeningTag"),
        "unexpected '</'" => at!(414, "BeginOpeningTag"),
        "'</' expected" => at!(426, "BeginClosingTag"),
        "Name begins with an invalid character" => at!(436, "ReadName"),
        "unsupported encoding" => at!(564, "SkipQDecl"),
        "unknown tag in file header" => at!(634, "ReadFileHeader"),
        "unknown DOCTYPE" => at!(665, "ReadFileHeader"),
        "entity reference is too long" => at!(736, "ReadEscapedChar"),
        "invalid entity reference" => at!(740, "ReadEscapedChar"),
        "invalid char reference" => at!(746, "ReadEscapedChar"),
        "invalid symbol in char reference" => at!(761, "ReadEscapedChar"),
        "unknown entity name" => at!(792, "ReadEscapedChar"),
        "'=' expected" => at!(894, "ReadAttributeValue"),
        "attribute value must start with ' or \"" => at!(898, "ReadAttributeValue"),
        "invalid float number" => at!(1040, "ReadDouble"),
        "CDATA section expected" => at!(1283, "ReadCDSection"),
        _ => at!(275, "SkipWSAndComments"),
    }
}

/// `SkipTagName`: `tag` after the prefix `name`.
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1508-1515
fn skip_tag_name<'t>(tag: &'t [u8], name: &[u8]) -> Result<&'t [u8], (ThrowSite, String)> {
    if !tag.starts_with(name) {
        return Err((
            SKIP_TAG_NAME,
            format!("invalid tag name: {}", String::from_utf8_lossy(tag)),
        ));
    }
    Ok(&tag[name.len()..])
}

/// `SkipStackTagName`: what is left of `tag` after the name that the frame `level` frames below
/// the top of the stack gives it: a type's name, a member's tag (its parent's name, `_` and its
/// own), an element's (its parent's, `_E`).
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1517-1561
fn skip_stack_tag_name<'t>(
    frames: &[XmlFrame],
    tag: &'t [u8],
    level: usize,
) -> Result<&'t [u8], (ThrowSite, String)> {
    let Some(frame) = frames.len().checked_sub(level + 1).map(|at| &frames[at]) else {
        return Err((SKIP_STACK_TAG_NAME, "illegal frame type".into()));
    };
    match frame.kind {
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
            skip_tag_name(tag, frame.member_name().as_bytes())
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

fn skip_stack_tag_char<'t>(
    frames: &[XmlFrame],
    tag: &'t [u8],
    level: usize,
    c: u8,
) -> Result<&'t [u8], (ThrowSite, String)> {
    let tag = skip_stack_tag_name(frames, tag, level)?;
    match tag.split_first() {
        Some((&first, rest)) if first == c => Ok(rest),
        _ => Err((
            SKIP_STACK_TAG_CHAR,
            format!("invalid tag name: {}", String::from_utf8_lossy(tag)),
        )),
    }
}

/// `FindEmpty`: the first variant whose type is a container (`eTypeFamilyContainer`, through the
/// pointers of aliases), which an empty choice element holds empty: an anonymous SEQUENCE OF or
/// SET OF (a named one is a class of its own, `CSeq_descr` and the like).
// NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:406-430
pub(super) fn find_empty(variants: &[Member]) -> Option<usize> {
    variants.iter().position(|variant| {
        let mut ty = variant.ty;
        loop {
            let info = type_info(ty);
            match info.kind {
                Kind::Alias(target) => ty = target,
                Kind::List { .. } => return info.name().is_empty(),
                _ => return false,
            }
        }
    })
}

/// The name of the type of the ENUMERATED or INTEGER with named values whose values are `values`
/// (`CEnumeratedTypeValues::GetName`; each type of the table has a span of its own).
pub(super) fn enum_type_name(values: &'static [EnumValue]) -> &'static str {
    super::table::TYPES
        .iter()
        .find(|info| match info.kind {
            Kind::Enumerated { values: span, .. } => std::ptr::eq(span.values(), values),
            _ => false,
        })
        .map_or("", |info| info.name())
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

/// The members of the SEQUENCE or the variants of the CHOICE `ty`.
fn items_of(ty: u16) -> &'static [Member] {
    match type_info(ty).kind {
        Kind::Sequence(span) | Kind::Choice(span) => span.members(),
        _ => &[],
    }
}

impl<S: AsnSource> XmlInput<S> {
    pub(crate) fn new(src: S, style: XmlStyle, names: Names) -> Self {
        XmlInput {
            x: XmlReader::object_stream(src),
            style,
            names,
            frames: Vec::new(),
            posts: Vec::new(),
        }
    }

    /// `GetStackTrace()` of the frames (`encoder::object_stack_trace`).
    fn stack_trace(&self) -> String {
        object_stack_trace(
            self.frames.first().map_or("?", XmlFrame::type_name),
            self.frames
                .iter()
                .skip(1)
                .filter_map(|frame| match frame.kind {
                    FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                        frame.member.map(|member| StackStep::Member(member.name()))
                    }
                    FrameKind::ArrayElement => Some(StackStep::Element),
                    _ => None,
                }),
        )
    }

    /// The bytes the lexer replaced by `#` in a call of the stack, posted with the stack trace of
    /// the call (the frames do not change inside one), the line each was read at and no string
    /// (`kEmptyStr`).
    fn post_replaced(&mut self) {
        if self.x.replaced.is_empty() {
            return;
        }
        let trace = self.stack_trace();
        for (byte, line) in std::mem::take(&mut self.x.replaced) {
            self.posts.push(Post::BadChar(ReplacedChar {
                trace: trace.clone(),
                byte,
                position: StreamPosition::Line(line as u64),
                string: std::sync::Arc::from(&[][..]),
            }));
        }
    }

    /// A call's result, after what it posted (also when it fails).
    #[inline]
    fn posted<T>(&mut self, result: R<T>) -> R<T> {
        self.post_replaced();
        result
    }

    /// `ThrowError` at `site` with the stream's position.
    #[cold]
    #[inline(never)]
    fn throw(&self, site: ThrowSite, message: impl AsRef<[u8]>) -> Fail {
        crate::report::serial::asn_input::throw_error(site, self.position(), message)
    }

    /// The failure of the lexer: the end of the input is NCBI's `CEofException`, a refusal names
    /// the document and the reader, the rest is the `CSerialException` of the function of NCBI's
    /// stream that throws it.
    #[cold]
    #[inline(never)]
    fn fail(&self, error: XmlError) -> Fail {
        match error {
            XmlError::Format(message) if message == EOF_MESSAGE => Fail::Eof,
            XmlError::Format(message) => self.throw(lexer_site(&message), message),
            XmlError::NotBool => Fail::exception(NcbiException::with_pos(
                crate::ncbi_str::STRING_TO_BOOL,
                "String cannot be converted to bool",
                0,
            )),
            XmlError::Unsupported(what) => Fail::unsupported(format!(
                "{what} in an XML {} is not read by LOSAT's {} ({})",
                self.names.document,
                self.names.reader,
                self.position()
            )),
        }
    }

    /// The lexer's `result` as the stream's.
    #[inline]
    fn lex<T>(&self, result: Result<T, XmlError>) -> R<T> {
        result.map_err(|error| self.fail(error))
    }

    /// A `SkipStackTagName` failure as the stream's.
    #[cold]
    fn tag_error(&self, (site, message): (ThrowSite, String)) -> Fail {
        self.throw(site, message)
    }

    /// The top frame.
    fn top(&self) -> Option<&XmlFrame> {
        self.frames.last()
    }

    /// The name of the top frame's type.
    fn top_name(&self) -> &'static str {
        self.top().map_or("", XmlFrame::type_name)
    }

    /// The name of the next opening tag, in `last_tag`: the one given back, or read
    /// (`ReadName` of what `BeginOpeningTag` finds, or `RejectedName`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:525-535
    fn read_tag_name(&mut self) -> R<()> {
        match self.x.rejected.take() {
            Some(name) => {
                self.x.last_tag = name;
                self.x.state = TagState::InsideOpening;
                Ok(())
            }
            None => {
                let c = self.x.begin_opening_tag();
                let c = self.lex(c)?;
                let named = self.x.read_name(c);
                self.lex(named)
            }
        }
    }

    /// `OpenTag`: the opening tag `e` (the name given back by `UndoClassMember`, or read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1563-1573
    fn open_tag(&mut self, e: &str) -> R<()> {
        self.read_tag_name()?;
        if self.x.last_tag != e.as_bytes() {
            return Err(self.throw(
                OPEN_TAG,
                format!(
                    "tag '{e}' expected: {}",
                    String::from_utf8_lossy(&self.x.last_tag)
                ),
            ));
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
        let c = self.x.begin_closing_tag();
        let c = self.lex(c)?;
        let named = self.x.read_name(c);
        self.lex(named)?;
        if self.x.last_tag != e.as_bytes() {
            return Err(self.throw(
                CLOSE_TAG,
                format!(
                    "tag '{e}' expected: {}",
                    String::from_utf8_lossy(&self.x.last_tag)
                ),
            ));
        }
        let end = self.x.end_tag();
        self.lex(end)
    }

    /// `OpenStackTag` at level 0 of an element of a list: the opening tag named by the stack.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1588-1602
    fn open_stack_tag(&mut self) -> R<()> {
        if self.x.rejected.is_some() {
            return self.read_tag_name();
        }
        let c = self.x.begin_opening_tag();
        let c = self.lex(c)?;
        let named = self.x.read_name(c);
        self.lex(named)?;
        self.check_stack_tag(OPEN_STACK_TAG)
    }

    /// `CloseStackTag` at level 0: nothing after a self-closed opening tag, else the closing tag
    /// named by the stack.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1604-1623
    fn close_stack_tag(&mut self) -> R<()> {
        if self.x.state == TagState::SelfClosed {
            self.x.state = TagState::Outside;
            return Ok(());
        }
        let c = self.x.begin_closing_tag();
        let c = self.lex(c)?;
        let named = self.x.read_name(c);
        self.lex(named)?;
        self.check_stack_tag(CLOSE_STACK_TAG)?;
        let end = self.x.end_tag();
        self.lex(end)
    }

    /// The tag in `last_tag` must be the one the stack names: `SkipStackTagName` of it from
    /// level 0 leaves nothing (else "unexpected tag" at `site`).
    fn check_stack_tag(&self, site: ThrowSite) -> R<()> {
        let tag = &self.x.last_tag;
        let rest = skip_stack_tag_name(&self.frames, tag, 0).map_err(|e| self.tag_error(e))?;
        if !rest.is_empty() {
            return Err(self.throw(
                site,
                format!(
                    "unexpected tag: {}{}",
                    String::from_utf8_lossy(tag),
                    String::from_utf8_lossy(rest)
                ),
            ));
        }
        Ok(())
    }

    /// The member or variant of `items` that the tag in `last_tag` names (`Parent_member`), from
    /// `pos` on.
    fn find_item(&self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>> {
        let id = skip_stack_tag_char(&self.frames, &self.x.last_tag, 1, b'_')
            .map_err(|e| self.tag_error(e))?;
        Ok((pos..items.count()).find(|&index| items.name(index).as_bytes() == id))
    }

    /// `UnexpectedMember` of `id` among `items`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2105-2114
    #[cold]
    fn unexpected_member(&self, id: &[u8], items: &(impl Items + ?Sized)) -> Fail {
        let mut message = format!(
            "\"{}\": unexpected member, should be one of: ",
            String::from_utf8_lossy(id)
        );
        for index in 0..items.count() {
            message.push('"');
            message.push_str(items.name(index));
            message.push_str("\" ");
        }
        self.throw(UNEXPECTED_MEMBER, message)
    }

    /// `GetInt4`, `GetUint4` (`signed` false) and `GetInt8` of NCBI's input buffer after
    /// `BeginData` (`asn_text_reader::buffer_number`), whose errors are its `CUtilException` with
    /// the line.
    fn number(&mut self, max: u64, signed: bool) -> R<i64> {
        let data = self.x.begin_data();
        self.lex(data)?;
        let number = self.x.get_number(max, signed);
        match self.lex(number)? {
            Ok(n) => Ok(n),
            Err(stop) => Err(self.number_error(stop)),
        }
    }

    /// The `CUtilException` of NCBI's input buffer at the lexer's line.
    #[cold]
    #[inline(never)]
    fn number_error(&self, stop: NumberStop) -> Fail {
        let (site, what) = match stop {
            NumberStop::Bad => (BAD_NUMBER, "bad number"),
            _ => (NUMBER_OVERFLOW, "number overflow"),
        };
        Fail::exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.x.line),
        ))
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2648-2675
    // (`SkipUNumber`, 2677-2703, takes only `+`.)
    fn skip_digits(&mut self, signed: bool) -> R<()> {
        let data = self.x.begin_data();
        self.lex(data)?;
        let c = self.x.skip_ws_and_comments();
        let mut c = self.lex(c)?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            let next = self.x.peek(1);
            c = self.lex(next)?;
            i = 2;
        }
        if !c.is_ascii_digit() {
            let site = if signed { SKIP_S_NUMBER } else { SKIP_U_NUMBER };
            return Err(self.throw(site, "invalid symbol in number"));
        }
        while self.x.peek_no_eof(i).is_ascii_digit() {
            i += 1;
        }
        self.x.skip(i);
        Ok(())
    }

    /// `ReadEnum`'s value, after the type's own tag was opened (or not).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1376-1453
    fn enum_value(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32> {
        if self.x.state != TagState::InsideOpening {
            // outside of tag
            if !values.integer() {
                return Err(self.throw(
                    ATTRIBUTE_VALUE_EXPECTED_OUTSIDE,
                    "attribute 'value' expected",
                ));
            }
            return Ok(self.number(i32::MAX as u64, true)? as i32);
        }
        let c = self.x.skip_ws();
        let c = self.lex(c)?;
        if matches!(c, b'>' | b'/') {
            // no attribute
            if !values.integer() {
                return Err(self.throw(ATTRIBUTE_VALUE_EXPECTED, "attribute 'value' expected"));
            }
            self.x.skip(1);
            self.x.state = TagState::Outside;
            return Ok(self.number(i32::MAX as u64, true)? as i32);
        }
        let mut found = false;
        loop {
            let more = self.x.has_attlist();
            if !self.lex(more)? {
                break;
            }
            let c = self.x.skip_ws();
            let c = self.lex(c)?;
            let named = self.x.read_name(c);
            self.lex(named)?;
            if self.x.last_tag == b"value" {
                found = true;
                break;
            }
            let value = self.x.read_attribute_value();
            self.lex(value)?;
        }
        if !found {
            // `fMissingValue`: the member that holds the value may take it
            // (`ReadWithSetFlagMember`, `read::TypedReader`).
            let end = self.x.end_opening_tag_self_closed();
            self.lex(end)?;
            return Err(self.throw(ATTRIBUTE_VALUE_MISSING, "attribute 'value' is missing"));
        }
        let name = self.x.read_attribute_value();
        let name = self.lex(name)?;
        let name = truncate_spaces(&name);
        let Some(value) = values.value_of(name) else {
            return Err(invalid_enum(
                FIND_VALUE,
                String::from_utf8_lossy(name).into_owned(),
            ));
        };
        let closed = self.x.end_opening_tag_self_closed();
        if !self.lex(closed)? && values.integer() {
            // read integer value
            let ws = self.x.skip_ws_and_comments();
            self.lex(ws)?;
            if value != self.number(i32::MAX as u64, true)? as i32 {
                return Err(self.throw(
                    INCOMPATIBLE_NAMED_INTEGER,
                    "incompatible name and value of named integer",
                ));
            }
        }
        Ok(value)
    }

    /// `GetHexChar`: the value of a digit or letter, `None` at the `<` that ends the data.
    fn hex_char(&mut self) -> R<Option<u32>> {
        let c = self.x.get_char();
        let c = self.lex(c)?;
        match c {
            b'0'..=b'9' => Ok(Some(u32::from(c - b'0'))),
            b'A'..=b'Z' => Ok(Some(u32::from(c - b'A') + 10)),
            b'a'..=b'z' => Ok(Some(u32::from(c - b'a') + 10)),
            _ => {
                self.x.skip_back(1);
                if c != b'<' {
                    return Err(self.throw(READ_BYTES, "invalid char in octet string"));
                }
                Ok(None)
            }
        }
    }

    fn file_header(&mut self) -> R<Option<Vec<u8>>> {
        if self.style == XmlStyle::Std {
            return Err(Fail::unsupported(
                "standard XML (SetEnforcedStdXml) is not read by LOSAT's XML reader".into(),
            ));
        }
        self.x.reset_state();
        let name = self.x.read_file_header();
        Ok(Some(self.lex(name)?))
    }

    fn class_member(&mut self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>> {
        if self.x.rejected.is_none() {
            if pos == 0 {
                let more = self.x.has_attlist();
                if self.lex(more)? {
                    let attributes = self.x.read_undefined_attributes();
                    self.lex(attributes)?;
                }
            }
            let closed = self.x.state == TagState::SelfClosed || {
                let closed = self.x.this_tag_is_self_closed();
                self.lex(closed)?
            };
            if closed {
                return Ok(None);
            }
            let closing = self.x.next_tag_is_closing();
            if self.lex(closing)? || pos >= items.count() {
                return Ok(None);
            }
            let tag = self.x.next_is_tag();
            if !self.lex(tag)? {
                return Ok(None);
            }
        }
        self.read_tag_name()?;
        match self.find_item(items, pos)? {
            Some(index) => Ok(Some(index)),
            None => {
                let id = skip_stack_tag_char(&self.frames, &self.x.last_tag, 1, b'_')
                    .map_err(|e| self.tag_error(e))?
                    .to_vec();
                Err(self.unexpected_member(&id, items))
            }
        }
    }

    fn choice_variant(&mut self, items: &(impl Items + ?Sized)) -> R<usize> {
        if self.x.rejected.is_none() {
            let more = self.x.has_attlist();
            if self.lex(more)? {
                let attributes = self.x.read_undefined_attributes();
                self.lex(attributes)?;
            }
            if self.x.state == TagState::SelfClosed {
                return Err(self.throw(CHOICE_VARIANT_EXPECTED, "choice variant id expected"));
            }
            let closing = self.x.next_tag_is_closing();
            if self.lex(closing)? {
                // The variants of the choice frame below the variant's.
                let choice = self
                    .frames
                    .len()
                    .checked_sub(2)
                    .map_or(NO_TYPE, |at| self.frames[at].ty);
                let Some(index) = (choice != NO_TYPE)
                    .then(|| find_empty(items_of(choice)))
                    .flatten()
                else {
                    return Err(self.throw(CHOICE_VARIANT_EXPECTED, "choice variant id expected"));
                };
                if let Some(top) = self.frames.last_mut() {
                    top.notag = true;
                }
                return Ok(index);
            }
            // `NextIsTag` looks for a variant without a tag when no tag follows (none in an
            // ASN.1 module); then `ReadName` of `BeginOpeningTag` throws for what is there.
            let tag = self.x.next_is_tag();
            self.lex(tag)?;
        }
        self.read_tag_name()?;
        match self.find_item(items, 0)? {
            Some(index) => Ok(index),
            None => {
                let tag = self.x.last_tag.clone();
                Err(self.unexpected_member(&tag, items))
            }
        }
    }

    fn container_element(&mut self) -> R<bool> {
        let no_more = match self
            .x
            .this_tag_is_self_closed()
            .and_then(|closed| Ok(closed || self.x.next_tag_is_closing()?))
        {
            Ok(no_more) => no_more,
            Err(XmlError::Format(message)) if message == EOF_MESSAGE => true,
            Err(error) => return Err(self.fail(error)),
        };
        if no_more {
            return Ok(false);
        }
        if self
            .top()
            .is_some_and(|frame| frame.kind == FrameKind::ArrayElement)
        {
            self.open_stack_tag()?;
        }
        Ok(true)
    }

    fn boolean(&mut self) -> R<bool> {
        let mut text = Vec::new();
        let mut have = false;
        loop {
            let more = self.x.has_attlist();
            if !self.lex(more)? {
                break;
            }
            let c = self.x.skip_ws();
            let c = self.lex(c)?;
            let named = self.x.read_name(c);
            self.lex(named)?;
            let is_value = self.x.last_tag == b"value";
            let value = self.x.read_attribute_value();
            let value = self.lex(value)?;
            if is_value {
                text = value;
                have = true;
            }
        }
        if !have {
            let word = self.x.read_word(&mut text);
            self.lex(word)?;
        }
        let text = truncate_spaces(&text);
        let value = match text {
            b"true" | b"1" => true,
            b"false" | b"0" => false,
            _ => {
                return Err(self.throw(
                    BOOL_VALUE_EXPECTED,
                    format!(
                        "'true' or 'false' value expected: {}",
                        String::from_utf8_lossy(text)
                    ),
                ))
            }
        };
        let closed = self.x.end_opening_tag_self_closed();
        if !self.lex(closed)? {
            let closing = self.x.next_tag_is_closing();
            if !self.lex(closing)? {
                return Err(self.throw(BOOL_NOT_EMPTY, "boolean tag must have empty contents"));
            }
        }
        Ok(value)
    }

    fn skip_octet_text(&mut self) -> R<()> {
        let data = self.x.begin_data();
        self.lex(data)?;
        loop {
            let c = self.x.peek(0);
            let c = self.lex(c)?;
            match c {
                b'0'..=b'9' | b'A'..=b'Z' | b'a'..=b'z' | b'+' | b'/' | b'=' => self.x.skip(1),
                b'\r' | b'\n' => {
                    self.x.skip(1);
                    self.x.skip_end_of_line(c);
                }
                b'<' => return Ok(()),
                _ => return Err(self.throw(SKIP_BYTE_BLOCK, "invalid char in octet string")),
            }
        }
    }

    fn octet_text(&mut self) -> R<Vec<u8>> {
        let data = self.x.begin_data();
        self.lex(data)?;
        let mut bytes = Vec::new();
        loop {
            let Some(high) = self.hex_char()? else {
                return Ok(bytes);
            };
            let Some(low) = self.hex_char()? else {
                bytes.push((high << 4) as u8);
                return Ok(bytes);
            };
            bytes.push(((high << 4) | low) as u8);
        }
    }
}

/// The body of a call of the stack, then what it posted (also when it fails).
macro_rules! step {
    ($self:ident, $body:expr) => {{
        let result = $body;
        $self.posted(result)
    }};
}

impl<S: AsnSource> Input for XmlInput<S> {
    const ENCODING: crate::report::serial::asn_input::Encoding =
        crate::report::serial::asn_input::Encoding::Xml;
    const NAMED_TYPES: bool = true;
    const ARRAY_FRAMES: bool = false;
    const MISSING_VALUES: bool = true;

    /// `WillHaveName`: the element of a list has the frame of an element when its type has no
    /// name (the table's element of a list of `CRef<>`s is the pointed type).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1638-1648
    #[inline(always)]
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

    fn frame_pushed(&mut self, kind: FrameKind, ty: u16, depth: usize) {
        self.frames.push(XmlFrame {
            kind,
            ty,
            member: None,
            notag: false,
        });
        debug_assert_eq!(self.frames.len(), depth);
    }

    fn frame_member(&mut self, member: &'static Member) {
        if let Some(top) = self.frames.last_mut() {
            top.member = Some(member);
        }
    }

    fn frame_popped(&mut self, depth: usize) {
        self.frames.truncate(depth);
    }

    fn drain_posts(&mut self, posts: &mut Vec<Post>) {
        posts.append(&mut self.posts);
    }

    /// `GetPosition()`: the line of NCBI's input buffer (which only `SkipEndOfLine` counts).
    fn position(&self) -> StreamPosition {
        StreamPosition::Line(self.x.line as u64)
    }

    /// Each object of a stream is read by a stream of its own in NCBI (its line count starts
    /// again); no frame, nothing posted.
    fn new_object(&mut self) {
        self.x.line = 1;
        self.frames.clear();
        self.posts.clear();
        self.x.replaced.clear();
        self.x.src.release_before(self.x.pos);
    }

    /// `EndOfData()`: whether no object follows: the end of the input, or white space, comments
    /// and processing instructions up to it, or a failure in them.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:118-129
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:588-597
    #[cfg(test)]
    fn end_of_data(&mut self) -> bool {
        let pos = self.x.pos;
        self.x.src.byte_at(pos).is_none() || self.x.skip_ws_and_comments().is_err()
    }

    /// `ResetState` and `ReadFileHeader`: the state of the stream at its root, the byte order
    /// mark, the declaration, the `DOCTYPE`, and the name of the root's tag (given back).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:84-101
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:582-670
    // (`ResetState` returns at once when more than one frame is on the stack; at the root it
    // clears the tag state and the tag names kept. The `m_Encoding` of `ReadFileHeader` is the
    // one of the stream so far; a declaration with an `encoding` sets it. `Read` of a root
    // object starts with `SkipFileHeader`, whose first call is `ResetState`, objistr.cpp:883-900.)
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        step!(self, self.file_header())
    }

    /// `BeginClass`: the tag of a class with a name.
    fn begin_class(&mut self) -> R<()> {
        let name = self.top_name();
        step!(
            self,
            if name.is_empty() {
                Ok(())
            } else {
                self.open_tag(name)
            }
        )
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
    fn begin_member(&mut self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>> {
        step!(self, self.class_member(items, pos))
    }

    /// `EndClassMember`: the member's closing tag.
    fn end_member(&mut self) -> R<()> {
        step!(self, self.close_stack_tag())
    }

    /// `EndClass`: the closing tag of a class with a name.
    fn end_class(&mut self) -> R<()> {
        let name = self.top_name();
        step!(
            self,
            if name.is_empty() {
                Ok(())
            } else {
                self.close_tag(name)
            }
        )
    }

    /// `BeginChoice`: the tag of a choice with a name.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2391-2409
    fn begin_choice(&mut self) -> R<()> {
        let name = self.top_name();
        step!(
            self,
            if name.is_empty() {
                Ok(())
            } else {
                self.open_tag(name)
            }
        )
    }

    /// `BeginChoiceVariant` of a CHOICE of an ASN.1 module: the variant whose tag,
    /// `Choice_variant`, follows. A choice element that the stream already ended self-closed is
    /// the "choice variant id expected" of `ReadChoiceSimple` (one whose opening tag is still open
    /// is the "'>' expected" of `EndTag`, which `NextTagIsClosing` calls, since `/` is no
    /// attribute); an element whose closing tag follows holds the variant
    /// that `FindEmpty` finds (`find_empty`: `Seqdesc.modif` among them), read without a tag of
    /// its own (`SetNotag`) as an empty list, or, without one, is that error too.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2411-2507
    // NCBI reference (598d8ae6): c++/src/serial/choice.cpp:370-416
    // (`FindDeep` and `Find` by the tag name come before: see `begin_member`. No CHOICE of the
    // table may be empty (`MayBeEmpty`). NCBI 2.17.0 probe x1 of review D of line RG-F3:
    // `<Seqdesc></Seqdesc>` is read as `modif {}`, rc 0; x1b, `<Seqdesc/>`, is the error, rc 1.)
    fn begin_variant(&mut self, items: &(impl Items + ?Sized)) -> R<usize> {
        step!(self, self.choice_variant(items))
    }

    /// `EndChoiceVariant`: the variant's closing tag, unless the variant had none.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2509-2517
    fn end_variant(&mut self) -> R<()> {
        if let Some(top) = self.frames.last_mut() {
            if std::mem::take(&mut top.notag) {
                return Ok(());
            }
        }
        step!(self, self.close_stack_tag())
    }

    fn end_choice(&mut self) -> R<()> {
        let name = self.top_name();
        step!(
            self,
            if name.is_empty() {
                Ok(())
            } else {
                self.close_tag(name)
            }
        )
    }

    /// `ReadContainer` of a container without a name calls no `BeginContainer`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1860-1875
    fn begin_container(&mut self, _set: bool) -> R<()> {
        Ok(())
    }

    /// `HasMoreElements`, then `BeginArrayElement` for an element in an element frame: the next
    /// element's opening tag. The end of the input in `HasMoreElements` is no more elements (the
    /// closing tag after the list meets it again).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1765-1825
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1895-1910
    fn begin_element(&mut self) -> R<bool> {
        step!(self, self.container_element())
    }

    /// `EndArrayElement`: the closing tag of an element that had an opening tag of its own.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1912-1919
    fn end_element(&mut self) -> R<()> {
        if !self
            .top()
            .is_some_and(|frame| frame.kind == FrameKind::ArrayElement)
        {
            return Ok(());
        }
        step!(self, self.close_stack_tag())
    }

    /// `ReadContainerContents` calls no `EndContainer`; an empty implicit class is not an error
    /// (no class of the table is `SetImplicitNonEmpty`).
    fn end_container(&mut self) -> R<()> {
        Ok(())
    }

    /// `BeginNamedType`: the tag of the named type (`<IUPACna>`, `<Seq-descr>`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2002-2037
    fn begin_named(&mut self) -> R<()> {
        let name = self.top_name();
        step!(self, self.open_tag(name))
    }

    fn end_named(&mut self) -> R<()> {
        let name = self.top_name();
        step!(self, self.close_tag(name))
    }

    /// `ReadInt4`: `BeginData`, `GetInt4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:987-1025
    fn read_i32(&mut self) -> R<i32> {
        step!(self, self.number(i32::MAX as u64, true).map(|n| n as i32))
    }

    fn read_i64(&mut self) -> R<i64> {
        step!(self, self.number(i64::MAX as u64, true))
    }

    fn read_u32(&mut self) -> R<u32> {
        step!(
            self,
            self.number(u64::from(u32::MAX), false).map(|n| n as u32)
        )
    }

    fn skip_number(&mut self, signed: bool) -> R<()> {
        step!(self, self.skip_digits(signed))
    }

    /// `ReadBool` / `SkipBool`: the attribute `value` of the tag, or the word in it; the tag
    /// must be empty after it.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:933-977
    fn read_bool(&mut self) -> R<bool> {
        step!(self, self.boolean())
    }

    /// `ReadNull`: an empty tag (self-closed, or its closing tag next).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1044-1048
    fn read_null(&mut self) -> R<()> {
        step!(self, {
            let closed = self.x.end_opening_tag_self_closed();
            let empty = self.lex(closed)? || {
                let closing = self.x.next_tag_is_closing();
                self.lex(closing)?
            };
            if empty {
                Ok(())
            } else {
                Err(self.throw(READ_NULL, "empty tag expected"))
            }
        })
    }

    /// `SkipNull`: a self-closed tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2725-2729
    fn skip_null(&mut self) -> R<()> {
        step!(self, {
            let closed = self.x.end_opening_tag_self_closed();
            if self.lex(closed)? {
                Ok(())
            } else {
                Err(self.throw(SKIP_NULL, "empty tag expected"))
            }
        })
    }

    /// `ReadDouble` (the skip, `SkipFNumber`, reads the number too): the word in the tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1027-1042
    fn read_real(&mut self) -> R<f64> {
        step!(self, {
            let value = self.x.real();
            value.map_err(|error| match error {
                XmlError::Format(message) if message == "invalid float number" => {
                    self.throw(READ_DOUBLE, message)
                }
                error => self.fail(error),
            })
        })
    }

    fn skip_real(&mut self) -> R<()> {
        self.read_real().map(|_| ())
    }

    /// `ReadString` (a self-closed tag is the empty string) with the bytes replaced by `#` and
    /// posted (`x_VerifyChar`, `drain_posts`: not in `posts`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1233-1258
    fn read_string(
        &mut self,
        _kind: StringKind,
        out: &mut Vec<u8>,
        _posts: &mut Vec<StringPost>,
    ) -> R<()> {
        step!(self, {
            let closed = self.x.end_opening_tag_self_closed();
            if self.lex(closed)? {
                Ok(())
            } else {
                let data = self.x.read_tag_data(out);
                self.lex(data)
            }
        })
    }

    /// `SkipString`: the text up to the next tag (`ReadEscapedChar`, which replaces nothing).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2710-2723
    fn skip_string(&mut self, _kind: StringKind, _posts: &mut Vec<StringPost>) -> R<()> {
        step!(self, {
            let data = self.x.begin_data();
            self.lex(data)?;
            loop {
                let c = self.x.read_escaped_char(b'<');
                if self.lex(c)?.0 < 0 {
                    break Ok(());
                }
            }
        })
    }

    /// `ReadBytes` of an OCTET STRING: pairs of hexadecimal digits up to the next tag (an odd
    /// last digit is the high half of a byte).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2524-2542
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2596-2613
    fn read_octets(&mut self) -> R<Vec<u8>> {
        step!(self, self.octet_text())
    }

    /// `SkipByteBlock`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2731-2762
    fn skip_octets(&mut self) -> R<()> {
        step!(self, self.skip_octet_text())
    }

    /// `ReadEnum` of the stream: the attribute `value` of the tag (of the type's own tag when it
    /// has a name), or, for an INTEGER with named values, the number after the tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1376-1453
    fn read_enum(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32> {
        let name = values.type_name();
        step!(self, {
            if !name.is_empty() {
                self.open_tag(name)?;
            }
            let value = self.enum_value(values);
            // The value's posts come before an error of the closing tag.
            self.post_replaced();
            let value = value?;
            if !name.is_empty() {
                self.close_tag(name)?;
            }
            Ok(value)
        })
    }

    /// `ReadPointerType`: a self-closed tag (with no attribute) is `NULL`, anything else the
    /// value itself (the XML stream reads no object index or other type).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1455-1469
    fn read_pointer_type(&mut self) -> R<PointerType> {
        step!(self, {
            let more = self.x.has_attlist();
            if self.lex(more)? {
                return Ok(PointerType::This);
            }
            let null = match self.x.state {
                TagState::InsideOpening => {
                    let closed = self.x.end_opening_tag_self_closed();
                    self.lex(closed)?
                }
                state => state == TagState::SelfClosed,
            };
            Ok(if null {
                PointerType::Null
            } else {
                PointerType::This
            })
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1471-1475
    fn read_object_pointer(&mut self) -> R<()> {
        Err(self.throw(READ_OBJECT_POINTER, "Not Implemented"))
    }

    fn read_other_pointer(&mut self) -> R<Vec<u8>> {
        unreachable!("the XML stream reads no pointer to another type")
    }

    fn read_other_pointer_end(&mut self) -> R<()> {
        Ok(())
    }

    /// `StartDelayBuffer`: `BeginData`, the buffer starts here (the source keeps the bytes from
    /// here).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1492-1500
    // (A member's own tag is read before its value: the rejected tag is the first of a member
    // read by its inner class's name, which no delayed member of the table is.)
    fn delay_start(&mut self) -> R<usize> {
        step!(self, {
            let data = self.x.begin_data();
            self.lex(data)?;
            if self.x.rejected.is_some() {
                return Err(Fail::unsupported(format!(
                    "a delayed member that starts with a tag given back ({})",
                    self.position()
                )));
            }
            self.x.src.hold(Some(self.x.pos));
            Ok(self.x.pos)
        })
    }

    /// `EndDelayBuffer`: the bytes skipped, kept unparsed.
    fn delayed(&mut self, start: usize) -> Vec<u8> {
        let bytes = self.x.src.bytes(start..self.x.pos).to_vec();
        self.x.src.hold(None);
        bytes
    }
}

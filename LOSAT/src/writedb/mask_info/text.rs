//! ASN.1 text (`MSerial_AsnText`): the mask reader's [`Input`] over the shared lexical layer of
//! NCBI's ASN.1 text object stream (`report::serial::asn_text_reader`, `CObjectIStreamAsn` over
//! `CIStreamBuffer`). What is particular to the mask reader stays here: the members and variants
//! found by name or by tag in its schema, `m_BlockStart`, and the decision of the frames on the
//! end of the text (`Fail::Eof`).
//!
//! The text is read from a stream in blocks (`StreamSource`): one reader for all objects of a
//! file, as NCBI's object stream gives back the bytes it did not use when an object ends
//! (`CIStreamBuffer::Close`, strbuffer.cpp:125-139), so that the next object starts at the byte
//! after the last one read. Each object gets a new `CObjectIStream` in NCBI (`s_ReadObject`), so
//! its line count starts at 1 again (checked with makeblastdb 2.17.0: a string on the first line
//! of a second object is "at line 1").

use std::io::Read;

use super::reader::{throw_error, Fail, Input, StringPost, R};
use super::schema::{ChoiceDef, ClassDef, EnumDef};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::string_to_int;
use crate::report::serial::asn_text_reader::{
    AsnSource, AsnTextError, AsnTextErrorKind, MemberId, NcbiThrow, PointerType, Reader,
    StreamSource, FILL_BUFFER_EOF, SKIP_WHITE_SPACE_EOF, UNEXPECTED_MEMBER,
};
use crate::report::serial::encoder::StreamPosition;

/// The text of a mask file.
pub(super) type TextSource = StreamSource<Box<dyn Read>>;

pub(super) struct TextIn<S: AsnSource> {
    r: Reader<S>,
    /// `m_BlockStart`
    block_start: bool,
    /// `GetStackDepth()` (`Input::set_stack_depth`).
    depth: usize,
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1455-1463 (`BeginChoiceVariant`,
// cited at `begin_variant`)
const BEGIN_CHOICE_VARIANT: ThrowSite = ThrowSite {
    file: "src/serial/objistrasn.cpp",
    line: 1463,
    kind: "CSerialException::eFormatError",
    function: Some("ncbi::CObjectIStreamAsn::BeginChoiceVariant"),
};
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1488-1491 (`ReadChoiceSimple`, cited at
// `Rd::read_choice`)
const READ_CHOICE_SIMPLE: ThrowSite = ThrowSite {
    file: "src/serial/objistr.cpp",
    line: 1490,
    kind: "CSerialException::eFormatError",
    function: Some("ncbi::CObjectIStream::ReadChoiceSimple"),
};

impl<S: AsnSource> TextIn<S> {
    pub(super) fn new(src: S) -> Self {
        TextIn {
            r: Reader::new(src, "mask data", "mask data reader"),
            block_start: false,
            depth: 0,
        }
    }

    /// What NCBI's text stream throws for `error` of the lexical reader (`AsnTextError::ncbi`):
    /// the end of the text is a `CEofException`, except inside `SkipWhiteSpace` at a stack depth
    /// over 2, where it becomes the `CSerialException` that quotes it
    /// (`asn_text_reader::SKIP_WHITE_SPACE_EOF`).
    #[cold]
    fn fail(&self, error: AsnTextError) -> Fail {
        let line = error.line;
        match (error.kind, error.ncbi) {
            (AsnTextErrorKind::Eof, Some(NcbiThrow::WhiteSpaceEof)) if self.depth > 2 => {
                let eof = NcbiException::new(FILL_BUFFER_EOF, "end of file").what();
                Fail::serial(NcbiException::new(
                    SKIP_WHITE_SPACE_EOF,
                    [format!("line {line}: ").as_bytes(), &eof].concat(),
                ))
            }
            (AsnTextErrorKind::Eof, _) => Fail::Eof,
            (AsnTextErrorKind::Format, Some(NcbiThrow::Stream(site))) => Fail::serial(
                NcbiException::new(site, format!("line {line}: {}", error.message)),
            ),
            (AsnTextErrorKind::Format, Some(NcbiThrow::Plain(site))) => {
                Fail::serial(NcbiException::new(site, error.message))
            }
            (AsnTextErrorKind::Format, Some(NcbiThrow::Buffer(site))) => Fail::exception(
                NcbiException::new(site, format!("{} in line {line}", error.message)),
            ),
            // Every error of the lexical reader that this reader meets names NCBI's throw; one
            // that did not would be refused, not reported in LOSAT's words.
            _ => Fail::unsupported(error.message),
        }
    }

    /// The state of a new `CObjectIStreamAsn` at the next byte: line 1, no block started.
    pub(super) fn new_object(&mut self) {
        self.r.line = 1;
        self.block_start = false;
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:389-424
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:432-439,462-473
    // (The members and variants of these types have the automatic tags 0, 1, ... in module
    // order. Their names start with a lower case letter, so `GetAltItemIndex`, which looks up the
    // name with an upper case first letter, finds nothing.)
    /// `ReadMemberId` of what `SkipWhiteSpace` returns, and the member or variant it names among
    /// `names`, searched (`count` of them, by position) from `pos` on; one that is not there is an
    /// error (`UnexpectedMember`, which lists them all). `choice`: an empty name is
    /// `BeginChoiceVariant`'s error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1305-1314
    fn find_member(
        &mut self,
        count: usize,
        name: impl Fn(usize) -> &'static str,
        pos: usize,
        choice: bool,
    ) -> R<usize> {
        let member = self.r.member_or_tag().map_err(|e| self.fail(e))?;
        let (id, found) = match member {
            MemberId::Tag(digits) => {
                // `NStr::StringToInt` of `id`, the tag's digits: a number that does not fit is its
                // `CStringException`.
                let id = self.r.bytes(digits);
                let tag = string_to_int(id, 0).map_err(Fail::exception)?;
                let found = usize::try_from(tag)
                    .ok()
                    .filter(|tag| (pos..count).contains(tag));
                (id, found)
            }
            MemberId::Name(id) => {
                let id = self.r.bytes(id);
                if choice && id.is_empty() {
                    return Err(throw_error(
                        BEGIN_CHOICE_VARIANT,
                        &self.position(),
                        "choice variant id expected",
                    ));
                }
                let found = (pos..count).find(|&i| name(i).as_bytes() == id);
                (id, found)
            }
        };
        match found {
            Some(index) => Ok(index),
            None => Err(self.unexpected_member(id, (0..count).map(name))),
        }
    }

    /// `UnexpectedMember` of `id` and the items, of the members or variants `names`.
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, id: &[u8], names: impl Iterator<Item = &'static str>) -> Fail {
        let mut message = [b"\"", id, b"\": unexpected member, should be one of: "].concat();
        for name in names {
            message.extend_from_slice(format!("\"{name}\" ").as_bytes());
        }
        throw_error(UNEXPECTED_MEMBER, &self.position(), message)
    }
}

#[cfg(test)]
impl TextIn<TextSource> {
    /// The bytes the stream keeps.
    pub(super) fn buffered(&self) -> usize {
        self.r.src.buffered()
    }
}

impl<S: AsnSource> Input for TextIn<S> {
    /// `GetPosition()` of the text stream (objistr.cpp:807-810): `line N`.
    fn position(&self) -> String {
        format!("line {}", self.r.line)
    }

    #[inline]
    fn set_stack_depth(&mut self, depth: usize) {
        self.depth = depth;
    }

    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        let name = self.r.file_header().map_err(|e| self.fail(e))?;
        Ok(Some(self.r.bytes(name).to_vec()))
    }

    fn begin_class(&mut self) -> R<()> {
        self.block_start = self.r.start_block().map_err(|e| self.fail(e))?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1329-1362
    // (`CanSkipUnknownMembers()` is false by default, `m_TypeAlias` is null.)
    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        let more = self.r.next_element(&mut self.block_start);
        if !more.map_err(|e| self.fail(e))? {
            return Ok(None);
        }
        Ok(Some(self.find_member(
            def.members.len(),
            |i| def.members[i].name,
            pos,
            false,
        )?))
    }

    fn end_member(&mut self) -> R<()> {
        Ok(())
    }

    fn end_class(&mut self) -> R<()> {
        self.r.end_block().map_err(|e| self.fail(e))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1478
    // (BeginChoice has just set `m_BlockStart`, so `skipname` is false. A `}` there is no
    // variant: `ReadChoiceSimple`'s "choice variant id expected".)
    fn begin_choice(&mut self) -> R<()> {
        self.block_start = true;
        Ok(())
    }

    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize> {
        let more = self.r.next_element(&mut self.block_start);
        if !more.map_err(|e| self.fail(e))? {
            return Err(throw_error(
                READ_CHOICE_SIMPLE,
                &self.position(),
                "choice variant id expected",
            ));
        }
        self.find_member(def.variants.len(), |i| def.variants[i].name, 0, true)
    }

    fn end_variant(&mut self) -> R<()> {
        Ok(())
    }

    fn end_choice(&mut self) -> R<()> {
        self.block_start = false;
        Ok(())
    }

    fn begin_container(&mut self, _set: bool) -> R<()> {
        self.begin_class()
    }

    fn begin_element(&mut self) -> R<bool> {
        let more = self.r.next_element(&mut self.block_start);
        more.map_err(|e| self.fail(e))
    }

    fn end_container(&mut self) -> R<()> {
        self.end_class()
    }

    fn read_i32(&mut self) -> R<i32> {
        let value = self.r.int4();
        value.map_err(|e| self.fail(e))
    }

    fn read_i64(&mut self) -> R<i64> {
        let value = self.r.int8();
        value.map_err(|e| self.fail(e))
    }

    fn read_u32(&mut self) -> R<u32> {
        let value = self.r.uint4();
        value.map_err(|e| self.fail(e))
    }

    fn read_bool(&mut self) -> R<bool> {
        let value = self.r.boolean();
        value.map_err(|e| self.fail(e))
    }

    fn read_null(&mut self) -> R<()> {
        let value = self.r.null();
        value.map_err(|e| self.fail(e))
    }

    /// `ReadString` of a VisibleString (`Reader::visible_string`): each bad byte posted at its
    /// line with the piece being appended.
    fn read_string(&mut self, posts: &mut Vec<StringPost>) -> R<Vec<u8>> {
        let mut out = Vec::new();
        let result = self
            .r
            .visible_string(&mut out, true, false, &mut |line, byte, piece| {
                posts.push(StringPost::BadChar {
                    position: StreamPosition::Line(line as u64),
                    byte,
                    string: piece.to_vec(),
                })
            });
        result.map_err(|e| self.fail(e))?;
        Ok(out)
    }

    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        let value = if def.integer {
            self.r.named_integer(def.values)
        } else {
            self.r.enumerated(def.values).map(|(_, value)| value)
        };
        value.map_err(|e| self.fail(e))
    }

    /// `ReadPointerType` (`Reader::pointer_type`).
    fn read_pointer_type(&mut self) -> R<PointerType> {
        let value = self.r.pointer_type();
        value.map_err(|e| self.fail(e))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1631-1640
    // (`TObjectIndex` is `size_t`, objistr.hpp:1029.)
    fn read_object_pointer(&mut self) -> R<()> {
        let value = self.r.int8();
        value.map(|_| ()).map_err(|e| self.fail(e))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1643-1646
    /// `ReadOtherPointer` (`Reader::type_id`). The name of `[]` is the `std::string` of the largest
    /// `size_t` count of bytes (-1 cast to `size_t`) that `CTempString::operator string` asks for:
    /// `std::length_error`.
    fn read_other_pointer(&mut self) -> R<Vec<u8>> {
        let id = self.r.type_id();
        match id.map_err(|e| self.fail(e))? {
            Some(id) => Ok(self.r.bytes(id).to_vec()),
            None => Err(Fail::LengthError),
        }
    }

    /// `ReadOtherPointerEnd` of the text stream: `CObjectIStream`'s, which reads nothing
    /// (objistr.cpp:1223-1225).
    fn read_other_pointer_end(&mut self) -> R<()> {
        Ok(())
    }
}

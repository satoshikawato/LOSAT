//! NCBI's ASN.1 text object stream (`CObjectOStreamAsn`, outfmt 8; the Seq-align part of
//! outfmt 11) over its output buffer (`COStreamBuffer`), for the calls of `Encoder`.
//!
//! Each document is written by a new stream, so the indentation and the line length start
//! at 0 for every document (and NCBI's line count, `GetPosition`, at 1).
//!
//! The encoder has two forms: the plain one (`AsnTextEncoder::new`), which writes the bytes of a
//! document and reports a VisibleString byte outside 0x20..0x7E as
//! `SerialError::NonVisibleChar`, and NCBI's (`AsnTextEncoder::ncbi`), which a writer uses through
//! `encoder::Traced` to write such a document again: it writes `#` for the byte and records it
//! (`encoder::Replaced`) at NCBI's line, with the bytes that NCBI's output buffer
//! (`stream_buffer::StreamBuffer`, one reservation per call of NCBI's stream) has handed out then.
//!
//! NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:816-827
//! ```c++
//! CNcbiOstream& WriteObject(CNcbiOstream& os, TConstObjectPtr ptr, TTypeInfo info)
//! {
//!     unique_ptr<CObjectOStream> ostr( CObjectOStream::Open( s_FlagsToFormat(os), os) );
//!     ...
//!     ostr->Write(ptr,info);
//!     return os;
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:64-72
//! ```c++
//! CObjectOStreamAsn::CObjectOStreamAsn(CNcbiOstream& out,
//!                                      EFixNonPrint how)
//!     : CObjectOStream(eSerial_AsnText, out)
//! {
//!     FixNonPrint(how);
//!     m_Output.SetBackLimit(80);
//!     SetSeparator("\n");
//!     SetAutoSeparator(true);
//! }
//! ```
//! The back limit keeps at least the last 80 bytes in the buffer when it is flushed. The string
//! wrap below looks back to the start of the line or NCBI's buffer, whichever comes first: in the
//! NCBI form the buffer's start is NCBI's, in the plain form the start of the bytes the encoder
//! holds (`flush_to` keeps the line being written and the end of line before it). The two agree
//! but where a line is longer than 80 bytes when a string's first byte comes (its look back can
//! then end at NCBI's buffer start): no document that LOSAT writes has such a line (the member
//! names and the indentation before a string stay far below it).

use super::encoder::{Encoder, Replaced, SerialError};
use super::real::write_asn_text_double;
use super::stream_buffer::{StreamBuffer, ASN_TEXT_BACK_LIMIT};
use super::visible::{good_visible_char, Position};

/// The column after which a VisibleString is wrapped (`WrapAt(78, true)`).
const STRING_WRAP: usize = 78;

/// The line length from which an OCTET STRING starts a new line (`WrapAt(78, false)`).
const OCTET_WRAP: usize = 78;

/// The ASN.1 text encoder, appending one document to `out`: the plain form (`MODEL` false) or
/// NCBI's (`MODEL` true: replaced bytes, NCBI's line count and buffer).
pub struct AsnTextEncoder<'b, const MODEL: bool = false> {
    out: &'b mut Vec<u8>,
    /// Where the document starts in `out` (the start of NCBI's buffer, `m_Buffer`).
    start: usize,
    /// `COStreamBuffer::m_IndentLevel`.
    indent: usize,
    /// `COStreamBuffer::m_LineLength`: the bytes written since the last end of line.
    line_length: usize,
    /// `CObjectOStreamAsn::m_BlockStart`.
    block_start: bool,
    error: Option<SerialError>,
    /// NCBI's form: the bytes written as `#` (`replace`).
    replaced: Vec<Replaced>,
    /// NCBI's form: `COStreamBuffer::m_Line`, 1 and one more for each end of line (`GetPosition`).
    line: usize,
    /// NCBI's form: NCBI's output buffer, which decides the bytes NCBI has handed out when it posts
    /// a message (`Replaced::after`), and the bytes of the document handed to a sink (`flush_to`,
    /// `flush_through`).
    ncbi: StreamBuffer,
    drained: usize,
}

impl<'b> AsnTextEncoder<'b, false> {
    /// The plain encoder, appending a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        Self::with_out(out)
    }
}

impl<'b> AsnTextEncoder<'b, true> {
    /// NCBI's form of the encoder, appending a document to `out` (`encoder::Traced`).
    pub fn ncbi(out: &'b mut Vec<u8>) -> Self {
        Self::with_out(out)
    }
}

impl<'b, const MODEL: bool> AsnTextEncoder<'b, MODEL> {
    fn with_out(out: &'b mut Vec<u8>) -> Self {
        let start = out.len();
        Self {
            out,
            start,
            indent: 0,
            line_length: 0,
            block_start: false,
            error: None,
            replaced: Vec::new(),
            line: 1,
            ncbi: StreamBuffer::new(ASN_TEXT_BACK_LIMIT),
            drained: 0,
        }
    }

    /// `ReplaceVisibleChar` of the byte at `index` of the string being written, at the line
    /// before the byte's wrap, after the bytes NCBI's buffer has handed out (`Replaced`,
    /// `encoder::Traced`); `write_string` then writes `#` as it writes a byte.
    #[cold]
    #[inline(never)]
    fn replace(&mut self, index: usize) {
        self.replaced.push(Replaced {
            index,
            position: Position::Line(self.line),
            after: self.ncbi.flushed(),
            with_string: true,
        });
    }

    /// One reservation of `count` bytes in NCBI's buffer (`COStreamBuffer::Skip`, or `Reserve`
    /// and the bytes after it), in NCBI's form only.
    #[inline(always)]
    fn reserve(&mut self, count: usize) {
        if MODEL {
            self.ncbi.skip(count);
        }
    }

    /// One more end of line in NCBI's line count (NCBI's form).
    #[inline(always)]
    fn count_line(&mut self) {
        if MODEL {
            self.line += 1;
        }
    }

    /// `COStreamBuffer::PutString` / `PutChar`: bytes without an end of line, one call of NCBI's
    /// stream.
    #[inline(always)]
    fn put(&mut self, bytes: &[u8]) {
        self.reserve(bytes.len());
        self.out.extend_from_slice(bytes);
        self.line_length += bytes.len();
    }

    /// Where the string wrap stops looking back in `out`: NCBI's buffer start (`m_Buffer`) in
    /// NCBI's form, the start of the bytes held in the plain form.
    #[inline(always)]
    fn buffer_start(&self) -> usize {
        if MODEL {
            self.start + (self.ncbi.flushed() - self.drained.min(self.ncbi.flushed()))
        } else {
            self.start
        }
    }

    /// `COStreamBuffer::PutEol(true)`: an end of line and the indentation.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:356-379
    /// ```c++
    /// inline
    /// void COStreamBuffer::PutIndent(void)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     if (GetUseIndentation()) {
    ///         size_t count = m_IndentLevel;
    ///         memset(Skip(count), ' ', count);
    ///     }
    /// }
    ///
    /// inline
    /// void COStreamBuffer::PutEol(bool indent)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     if (!GetUseEol()) {
    ///         return;
    ///     }
    ///     char* pos = Reserve(1);
    ///     *pos = '\n';
    ///     m_CurrentPos = pos + 1;
    ///     ++m_Line;
    ///     m_LineLength = 0;
    ///     if ( indent )
    ///         PutIndent();
    /// }
    /// ```
    fn eol(&mut self) {
        self.reserve(1);
        self.out.push(b'\n');
        self.count_line();
        self.line_length = 0;
        let indent = self.indent;
        self.reserve(indent);
        self.out.resize(self.out.len() + indent, b' ');
        self.line_length += indent;
    }

    /// `CObjectOStreamAsn::StartBlock`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:539-544
    /// ```c++
    /// void CObjectOStreamAsn::StartBlock(void)
    /// {
    ///     m_Output.PutChar('{');
    ///     m_Output.IncIndentLevel();
    ///     m_BlockStart = true;
    /// }
    /// ```
    /// (`IncIndentLevel(size_t step = 2)`, include/util/strbuffer.hpp.)
    fn start_block(&mut self) {
        self.put(b"{");
        self.indent += 2;
        self.block_start = true;
    }

    /// `CObjectOStreamAsn::EndBlock`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:546-552
    /// ```c++
    /// void CObjectOStreamAsn::EndBlock(void)
    /// {
    ///     m_Output.DecIndentLevel();
    ///     m_Output.PutEol();
    ///     m_Output.PutChar('}');
    ///     m_BlockStart = false;
    /// }
    /// ```
    fn end_block(&mut self) {
        self.indent -= 2;
        self.eol();
        self.put(b"}");
        self.block_start = false;
    }

    /// `CObjectOStreamAsn::NextElement`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:554-561
    /// ```c++
    /// void CObjectOStreamAsn::NextElement(void)
    /// {
    ///     if ( m_BlockStart )
    ///         m_BlockStart = false;
    ///     else
    ///         m_Output.PutChar(',');
    ///     m_Output.PutEol();
    /// }
    /// ```
    fn next_element(&mut self) {
        if self.block_start {
            self.block_start = false;
        } else {
            self.put(b",");
        }
        self.eol();
    }

    /// `COStreamBuffer::WrapAt(lineLength, false)`: an end of line without indentation when the
    /// line holds `line_length` bytes or more.
    fn wrap_at_line_length(&mut self, line_length: usize) {
        if self.line_length >= line_length {
            self.reserve(1);
            self.out.push(b'\n');
            self.count_line();
            self.line_length = 0;
        }
    }

    /// `COStreamBuffer::PutEolAtWordEnd(lineLength)`: inserts an end of line (without
    /// indentation) before the last white space or `'` of the current line at a column up to
    /// `line_length`, looking back no further than a `"` or an end of line; without one, at
    /// column `line_length`, moved back before any `"` that precedes it (a doubled `""` is
    /// not split). The bytes after the inserted end of line start the next line.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1084-1137
    /// ```c++
    /// void COStreamBuffer::PutEolAtWordEnd(size_t lineLength)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     if (!GetUseEol()) {
    ///         return;
    ///     }
    ///     Reserve(1);
    ///     size_t linePos = m_LineLength;
    ///     char* pos = m_CurrentPos;
    ///     bool goodPlace = false;
    ///     while ( pos > m_Buffer && linePos > 0 ) {
    ///         --pos;
    ///         --linePos;
    ///         if ( linePos <= lineLength && (isspace((unsigned char) (*pos)) ||
    ///                                        *pos == '\'') ) {
    ///             goodPlace = true;
    ///             break;
    ///         }
    ///         else if ( *pos == '\n' || *pos == '"' ) {
    ///             // no suitable space found
    ///             break;
    ///         }
    ///     }
    ///
    ///     // Prevent insertion of more than one '\n'
    ///     if (pos > m_Buffer  &&  *(pos-1) == '\n') {
    ///         goodPlace = false;
    ///     }
    ///
    ///     if ( !goodPlace ) {
    ///         // no suitable space found
    ///         if ( linePos < lineLength ) {
    ///             pos += lineLength - linePos;
    ///             linePos = lineLength;
    ///         }
    ///         // assure we will not break double ""
    ///         while ( pos > m_Buffer && *(pos - 1) == '"' ) {
    ///             --pos;
    ///             --linePos;
    ///         }
    ///         if ( pos == m_Buffer ) {
    ///             // it's possible that we already put some " before...
    ///             while ( pos < m_CurrentPos && *pos == '"' ) {
    ///                 ++pos;
    ///                 ++linePos;
    ///             }
    ///         }
    ///     }
    ///     // split there
    ///     // insert '\n'
    ///     size_t count = m_CurrentPos - pos;
    ///     memmove(pos + 1, pos, count);
    ///     m_LineLength = count;
    ///     ++m_CurrentPos;
    ///     *pos = '\n';
    ///     ++m_Line;
    /// }
    /// ```
    /// `isspace` is that of the C locale; the bytes of the stream are visible characters,
    /// spaces and ends of line.
    // Inlined into both forms of `write_string`, as into its one caller before the replacing form.
    #[inline(always)]
    fn put_eol_at_word_end(&mut self, line_length: usize) {
        // `Reserve(1)`, then the look back stops at NCBI's buffer start (`m_Buffer`).
        self.reserve(1);
        let buffer = self.buffer_start();
        let mut line_pos = self.line_length;
        let mut pos = self.out.len();
        let mut good_place = false;
        while pos > buffer && line_pos > 0 {
            pos -= 1;
            line_pos -= 1;
            let byte = self.out[pos];
            if line_pos <= line_length
                && (matches!(byte, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
                    || byte == b'\'')
            {
                good_place = true;
                break;
            } else if byte == b'\n' || byte == b'"' {
                break;
            }
        }
        if pos > buffer && self.out[pos - 1] == b'\n' {
            good_place = false;
        }
        if !good_place {
            if line_pos < line_length {
                pos += line_length - line_pos;
            }
            while pos > buffer && self.out[pos - 1] == b'"' {
                pos -= 1;
            }
            if pos == buffer {
                while pos < self.out.len() && self.out[pos] == b'"' {
                    pos += 1;
                }
            }
        }
        let count = self.out.len() - pos;
        self.out.insert(pos, b'\n');
        self.count_line();
        self.line_length = count;
    }
}

impl<const MODEL: bool> Encoder for AsnTextEncoder<'_, MODEL> {
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:105-111
    /// ```c++
    /// void CObjectOStreamAsn::WriteFileHeader(TTypeInfo type)
    /// {
    ///     if ( true || m_Output.ZeroIndentLevel() ) {
    ///         WriteId(type->GetName());
    ///         m_Output.PutString(" ::= ");
    ///     }
    /// }
    /// ```
    /// (`WriteId` brackets a name with a space, `<` or `:`, objostrasn.cpp:489-504; a type
    /// name has none.)
    fn begin_document(&mut self, type_name: &'static str) {
        self.put(type_name.as_bytes());
        self.put(b" ::= ");
    }

    /// The separator `"\n"` that the auto separator writes after the object
    /// (`CObjectOStream::Write`, objostr.cpp:623-624).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:941-945
    /// ```c++
    /// void CObjectOStreamAsn::WriteSeparator(void)
    /// {
    ///     m_Output.PutString(GetSeparator());
    ///     FlushBuffer();
    /// }
    /// ```
    fn end_document(&mut self) {
        self.put(b"\n");
        // `WriteSeparator`: `FlushBuffer()`, every byte to the output.
        if MODEL {
            self.ncbi.flush_all();
        }
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:669-672
    /// ```c++
    /// void CObjectOStreamAsn::BeginClass(const CClassTypeInfo* /*classInfo*/)
    /// {
    ///     StartBlock();
    /// }
    /// ```
    fn begin_sequence(&mut self, _type_name: &str) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:674-677
    /// ```c++
    /// void CObjectOStreamAsn::EndClass(void)
    /// {
    ///     EndBlock();
    /// }
    /// ```
    fn end_sequence(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:679-684
    /// ```c++
    /// void CObjectOStreamAsn::BeginClassMember(const CMemberId& id)
    /// {
    ///     NextElement();
    ///
    ///     WriteMemberId(id);
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:646-667
    /// ```c++
    /// void CObjectOStreamAsn::WriteMemberId(const CMemberId& id)
    /// {
    ///     ...
    ///     if ( !name.empty() ) {
    ///         if (id.HaveNoPrefix() && isupper((unsigned char)name[0])) {
    ///             ...
    ///         } else {
    ///             m_Output.PutString(name);
    ///         }
    ///         m_Output.PutChar(' ');
    ///     }
    /// ```
    /// (The member names of the ASN.1 modules start with a lower-case letter.)
    fn member(&mut self, name: &'static str, _tag: u8) {
        self.next_element();
        self.put(name.as_bytes());
        self.put(b" ");
    }

    fn end_member(&mut self) {}

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:563-566
    /// ```c++
    /// void CObjectOStreamAsn::BeginContainer(const CContainerTypeInfo* /*cType*/)
    /// {
    ///     StartBlock();
    /// }
    /// ```
    fn begin_list(&mut self, _set: bool) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:573-576
    /// ```c++
    /// void CObjectOStreamAsn::BeginContainerElement(TTypeInfo /*elementType*/)
    /// {
    ///     NextElement();
    /// }
    /// ```
    fn element(&mut self) {
        self.next_element();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:568-571
    /// ```c++
    /// void CObjectOStreamAsn::EndContainer(void)
    /// {
    ///     EndBlock();
    /// }
    /// ```
    fn end_list(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:834-868
    /// ```c++
    /// void CObjectOStreamAsn::BeginChoice(const CChoiceTypeInfo* choiceType)
    /// {
    ///     if (choiceType->GetVariantInfo(kFirstMemberIndex)->GetId().IsAttlist()) {
    ///         ...
    ///     }
    ///     m_BlockStart = true;
    /// }
    /// ...
    /// void CObjectOStreamAsn::BeginChoiceVariant(const CChoiceTypeInfo* choiceType,
    ///                                            const CMemberId& id)
    /// {
    ///     if ( m_BlockStart ) {
    ///         m_BlockStart = false;
    ///     } else {
    ///         ...
    ///     }
    ///     WriteMemberId(id);
    /// }
    /// ```
    /// The variant's name and a space, on the line of what precedes it.
    fn variant(&mut self, _type_name: &'static str, name: &'static str, _tag: u8) {
        self.block_start = false;
        self.put(name.as_bytes());
        self.put(b" ");
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:842-848
    /// ```c++
    /// void CObjectOStreamAsn::EndChoice(void)
    /// {
    ///     if (TopFrame().GetNotag()) {
    ///         EndBlock();
    ///     }
    ///     m_BlockStart = false;
    /// }
    /// ```
    fn end_variant(&mut self) {
        self.block_start = false;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:155-158
    /// ```c++
    /// void CObjectOStreamAsn::WriteInt4(Int4 data)
    /// {
    ///     m_Output.PutInt4(data);
    /// }
    /// ```
    fn int(&mut self, value: i32) {
        let mut digits = itoa_i32(value);
        let text = digits.text();
        self.put(text);
    }

    /// `WriteInt8`: the decimal digits, as `WriteInt4`.
    fn big_int(&mut self, value: i64) {
        self.put(value.to_string().as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:113-127
    /// ```c++
    /// void CObjectOStreamAsn::WriteEnum(
    ///     const CEnumeratedTypeValues& values,
    ///     TEnumValueType value, const string& valueName)
    /// {
    ///     if (valueName.empty() || (m_WriteNamedIntegersByValue && values.IsInteger())) {
    ///         m_Output.PutInt4(value);
    ///     } else {
    /// // names coming from XML can begin with uppercase char;
    /// // in ASN, they must begin with lowercase one.
    ///         m_Output.PutChar((char)tolower((unsigned char)valueName[0]));
    ///         m_Output.PutString(valueName.data()+1, valueName.size()-1);
    ///     }
    /// }
    /// ```
    fn enumerated(&mut self, name: &str, _value: i32) {
        if MODEL {
            // `PutChar` of the first letter, then `PutString` of the rest.
            let (first, rest) = name.as_bytes().split_at(1.min(name.len()));
            self.put(first);
            self.put(rest);
        } else {
            self.put(name.as_bytes());
        }
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:127-131
    /// ```c++
    /// void CObjectOStreamAsn::WriteEnum(const CEnumeratedTypeValues& values,
    ///                                   TEnumValueType value)
    /// {
    ///     WriteEnum(values,value,values.FindNameEx(value, values.IsInteger()));
    /// }
    /// ```
    fn named_int(&mut self, name: &str, value: i32) {
        self.enumerated(name, value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:140-146
    /// ```c++
    /// void CObjectOStreamAsn::WriteBool(bool data)
    /// {
    ///     if ( data )
    ///         m_Output.PutString("TRUE");
    ///     else
    ///         m_Output.PutString("FALSE");
    /// }
    /// ```
    fn boolean(&mut self, value: bool) {
        self.put(if value { b"TRUE" } else { b"FALSE" });
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
    /// ```c++
    /// void CObjectOStreamAsn::WriteString(const char* ptr, size_t length)
    /// {
    ///     m_Output.PutChar('"');
    ///     ...
    ///     while ( length > 0 ) {
    ///         char c = *ptr;
    ///         if ( x_FixCharsMethod() != eFNP_Allow ) {
    ///             if ( !GoodVisibleChar(c) ) {
    ///                 ...
    ///                         c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, original, x_FixCharsSubst());
    ///         ...
    ///         m_Output.WrapAt(78, true);
    ///         m_Output.PutChar(c);
    ///         if ( c == '"' ) {
    ///             m_Output.PutChar('"');
    ///         }
    ///     }
    ///     m_Output.PutChar('"');
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:382-393
    /// ```c++
    /// inline
    /// void COStreamBuffer::WrapAt(size_t lineLength, bool keepWord)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     if ( keepWord ) {
    ///         if ( GetCurrentLineLength() > lineLength )
    ///             PutEolAtWordEnd(lineLength);
    ///     }
    /// ```
    /// A byte outside 0x20..0x7E (`GoodVisibleChar`, include/serial/impl/objistr.inl:456-459)
    /// is NCBI's replacement with an error message: here an error (`write_string`).
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string::<false>(value);
    }

    /// `visible_str` that writes NCBI's `#` for a byte outside 0x20..0x7E and records it
    /// (`replace`).
    fn visible_str_replacing(&mut self, value: &[u8]) {
        // Only NCBI's form replaces (`encoder::Traced` writes through it); the plain form
        // reports the byte.
        if MODEL {
            self.write_string::<true>(value);
        } else {
            self.write_string::<false>(value);
        }
    }

    /// `string_store` (`visible_str`) that writes NCBI's `#` (`visible_str_replacing`).
    fn string_store_replacing(&mut self, value: &[u8]) {
        // Only NCBI's form replaces (`encoder::Traced` writes through it); the plain form
        // reports the byte.
        if MODEL {
            self.write_string::<true>(value);
        } else {
            self.write_string::<false>(value);
        }
    }

    /// `WriteDouble2(data, DBL_DIG)` (`real::write_asn_text_double`).
    fn real(&mut self, value: f64) {
        let before = self.out.len();
        write_asn_text_double(value, self.out);
        self.line_length += self.out.len() - before;
        if MODEL {
            self.reserve_real(before);
        }
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:868-895
    /// ```c++
    /// void CObjectOStreamAsn::BeginBytes(const ByteBlock& )
    /// {
    ///     m_Output.PutChar('\'');
    /// }
    ///
    /// static const char HEX[] = "0123456789ABCDEF";
    ///
    /// void CObjectOStreamAsn::WriteBytes(const char* bytes, size_t length)
    /// {
    ///     while ( length-- > 0 ) {
    ///         char c = *bytes++;
    ///         m_Output.WrapAt(78, false);
    ///         m_Output.PutChar(HEX[(c >> 4) & 0xf]);
    ///         m_Output.PutChar(HEX[c & 0xf]);
    ///     }
    /// }
    ///
    /// void CObjectOStreamAsn::EndBytes(const ByteBlock& )
    /// {
    ///     m_Output.WrapAt(78, false);
    ///     m_Output.PutString("\'H");
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:382-393
    /// ```c++
    /// inline
    /// void COStreamBuffer::WrapAt(size_t lineLength, bool keepWord)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     if ( keepWord ) {
    ///         if ( GetCurrentLineLength() > lineLength )
    ///             PutEolAtWordEnd(lineLength);
    ///     }
    ///     else {
    ///         if ( GetCurrentLineLength() >= lineLength )
    ///             PutEol(false);
    ///     }
    /// }
    /// ```
    fn octets(&mut self, value: &[u8]) {
        const HEX: &[u8; 16] = b"0123456789ABCDEF";
        self.put(b"'");
        for &byte in value {
            self.wrap_at_line_length(OCTET_WRAP);
            if MODEL {
                // Two `PutChar`s.
                self.put(&[HEX[usize::from(byte >> 4)]]);
                self.put(&[HEX[usize::from(byte & 0xf)]]);
            } else {
                self.put(&[HEX[usize::from(byte >> 4)], HEX[usize::from(byte & 0xf)]]);
            }
        }
        self.wrap_at_line_length(OCTET_WRAP);
        self.put(b"'H");
    }

    /// `WriteStringStore` writes the bytes as `WriteString` does.
    fn string_store(&mut self, value: &[u8]) {
        self.visible_str(value);
    }

    fn buffered(&self) -> usize {
        self.out.len() - self.start
    }

    /// Keeps the last end of line and the line after it: `put_eol_at_word_end` looks back to
    /// the start of the current line and to the byte before it (an end of line stops it, and
    /// an end of line right before a place stops it from adding a second), and the buffer
    /// start of NCBI's stream only matters when the line is the first one.
    /// The plain form keeps the line being written and the end of line before it (the string
    /// wrap looks back no further); NCBI's form hands out what NCBI's buffer has handed out.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        if all {
            let through = self.drained + (self.out.len() - self.start);
            return self.drain_through(sink, through);
        }
        if MODEL {
            return self.drain_through(sink, self.ncbi.flushed());
        }
        self.flush_keeping(sink, 0)
    }

    fn flush_keeping(&mut self, sink: &mut dyn std::io::Write, keep: usize) -> std::io::Result<()> {
        if MODEL {
            return self.flush_to(sink, false);
        }
        let held = &self.out[self.start..];
        let Some(last_eol) = held.iter().rposition(|&byte| byte == b'\n') else {
            return Ok(());
        };
        let keep_from = last_eol.min(held.len().saturating_sub(keep));
        if keep_from == 0 {
            return Ok(());
        }
        self.drain_through(sink, self.drained + keep_from)
    }

    fn flush_through(
        &mut self,
        sink: &mut dyn std::io::Write,
        through: usize,
    ) -> std::io::Result<()> {
        self.drain_through(sink, through)
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }

    fn take_replaced(&mut self) -> Vec<Replaced> {
        std::mem::take(&mut self.replaced)
    }
}

impl<const MODEL: bool> AsnTextEncoder<'_, MODEL> {
    /// Hands the document's bytes before its byte `through` to `sink` and drops them.
    fn drain_through(
        &mut self,
        sink: &mut dyn std::io::Write,
        through: usize,
    ) -> std::io::Result<()> {
        if through <= self.drained {
            return Ok(());
        }
        let end = self.start + (through - self.drained);
        sink.write_all(&self.out[self.start..end])?;
        self.out.drain(self.start..end);
        self.drained = through;
        Ok(())
    }

    /// NCBI's reservations of a REAL written from `before` in `out` (`WriteDouble2`,
    /// objostrasn.cpp:176-226, `m_FastWriteDouble` true, objostr.cpp:79): one `PutString` of a
    /// constant, or `{ `, a `-`, the digits, `, 10, `, `PutInt4` of the exponent and ` }`.
    #[cold]
    fn reserve_real(&mut self, before: usize) {
        let text = &self.out[before..];
        let written = text.len();
        let digits_end = text.windows(6).position(|w| w == b", 10, ");
        match digits_end {
            Some(comma)
                if text.starts_with(b"{ ")
                    && !text.starts_with(b"{ 0,")
                    && !text.starts_with(b"{ -0,") =>
            {
                let sign = usize::from(text[2] == b'-');
                let exponent = written - (comma + 6) - 2;
                let digits = comma - 2 - sign;
                self.reserve(2);
                if sign == 1 {
                    self.reserve(1);
                }
                self.reserve(digits);
                self.reserve(6);
                self.reserve(exponent);
                self.reserve(2);
            }
            _ => self.reserve(written),
        }
    }

    /// The loop of `WriteString` (quoted at `visible_str`): without `REPLACE` a byte outside
    /// 0x20..0x7E is the error `NonVisibleChar` and is not written, with it the byte is
    /// recorded at its line (before its wrap, `replace`) and `#` is written in its place.
    #[inline(always)]
    fn write_string<const REPLACE: bool>(&mut self, value: &[u8]) {
        self.put(b"\"");
        for (index, &byte) in value.iter().enumerate() {
            let byte = if good_visible_char(byte) {
                byte
            } else if REPLACE {
                self.replace(index);
                b'#'
            } else {
                if self.error.is_none() {
                    self.error = Some(SerialError::NonVisibleChar(byte));
                }
                continue;
            };
            if self.line_length > STRING_WRAP {
                self.put_eol_at_word_end(STRING_WRAP);
            }
            self.put(&[byte]);
            if byte == b'"' {
                self.put(b"\"");
            }
        }
        self.put(b"\"");
    }
}

/// The decimal digits of an `Int4` (`COStreamBuffer::PutInt4`, NStr's `IntToString`).
pub(crate) struct DecimalI32 {
    buffer: [u8; 11],
    start: usize,
}

impl DecimalI32 {
    pub(crate) fn text(&mut self) -> &[u8] {
        &self.buffer[self.start..]
    }
}

/// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:951-975
/// ```c++
/// void COStreamBuffer::PutInt4(Int4 v)
///     THROWS1((CIOException, bad_alloc))
/// {
///     const size_t BSIZE = (sizeof(v)*CHAR_BIT) / 3 + 2;
///     char b[BSIZE];
///     Int4 n = v;
///     if ( n < 0 ) {
///         n = -n;
///     }
///     char* pos = b + BSIZE;
///     do {
///         Int4 a = '0'+n;
///         *--pos = char(a-10*(n/=Uint4(10)));
///     } while ( n );
///     if ( v < 0 ) {
///         *--pos = '-';
///     }
/// ```
/// `i32::MIN` is written as its digits. NCBI's `-n` overflows there and its build writes `-8`
/// (the archive's `WordThreshold` of a `-threshold` past `Int4`, NCBI BLAST+ 2.17.0); LOSAT does
/// not reproduce that overflow. Positions, lengths, counts and scores of BLAST do not reach it.
pub(crate) fn itoa_i32(value: i32) -> DecimalI32 {
    let mut buffer = [0u8; 11];
    let mut pos = buffer.len();
    let mut n = value.unsigned_abs();
    loop {
        pos -= 1;
        buffer[pos] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if value < 0 {
        pos -= 1;
        buffer[pos] = b'-';
    }
    DecimalI32 { buffer, start: pos }
}

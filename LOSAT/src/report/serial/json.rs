//! NCBI's JSON object stream (`CObjectOStreamJson`, outfmt 12) for the calls of `Encoder`.
//!
//! Each document is written by a new stream (`serialobject.cpp` `WriteObject`, see
//! `asn_text`), so the indentation starts at 0 for every document, and a document follows
//! the previous one without a separator.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:68-80
//! ```c++
//! CObjectOStreamJson::CObjectOStreamJson(CNcbiOstream& out, bool deleteOut)
//!     : CObjectOStream(eSerial_Json, out, deleteOut ? eTakeOwnership : eNoOwnership),
//!     m_FileHeader(false),
//!     m_BlockStart(false),
//!     m_ExpectValue(false),
//!     m_PreserveKeys(false),
//!     m_StringEncoding(eEncoding_UTF8),
//!     m_BinaryFormat(eDefault),
//!     m_WrapAt(0)
//! {
//!     m_WriteNamedIntegersByValue = true;
//!     m_EnforceWritingDefaults = true;
//! }
//! ```
//! The members that BLAST sets are always written, DEFAULT ones too
//! (`m_EnforceWritingDefaults`), as the type writers do.

use super::asn_text::itoa_i32;
use super::encoder::{Encoder, SerialError};
use super::real::write_json_double;

/// The JSON encoder, appending one document to `out`.
pub struct JsonEncoder<'b> {
    out: &'b mut Vec<u8>,
    /// `COStreamBuffer::m_IndentLevel`.
    indent: usize,
    /// `CObjectOStreamJson::m_BlockStart`.
    block_start: bool,
    /// A report written into the array of a 15 document: no wrapper object, the indentation
    /// one level in (`array_element`).
    array_element: bool,
    error: Option<SerialError>,
}

impl<'b> JsonEncoder<'b> {
    /// An encoder that appends a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        Self {
            out,
            indent: 0,
            block_start: false,
            array_element: false,
            error: None,
        }
    }

    /// An encoder for one report of a 15 document (`CBlastOStreamJson`,
    /// blastxml2_format.cpp:816-831): `WriteFileHeader` is `IncIndentLevel` alone and `EndOfWrite`
    /// is `DecIndentLevel` and `PutEol`, so the report is the object `{` (members one level in)
    /// `\n  }` and an end of line, with no wrapper object around it.
    pub fn array_element(out: &'b mut Vec<u8>) -> Self {
        let mut encoder = Self::new(out);
        encoder.indent = 2;
        encoder.array_element = true;
        encoder
    }

    /// `COStreamBuffer::PutEol(true)` (`asn_text::AsnTextEncoder::eol`).
    fn eol(&mut self) {
        self.out.push(b'\n');
        let indent = self.indent;
        self.out.resize(self.out.len() + indent, b' ');
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:811-818
    /// ```c++
    /// void CObjectOStreamJson::StartBlock(void)
    /// {
    ///     BeginValue();
    ///     m_Output.PutChar('{');
    ///     m_Output.IncIndentLevel();
    ///     m_BlockStart = true;
    ///     m_ExpectValue = false;
    /// }
    /// ```
    /// (`BeginValue` writes a skipped member's key, which only XML-like types have.)
    fn start_block(&mut self) {
        self.out.push(b'{');
        self.indent += 2;
        self.block_start = true;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:820-827
    /// ```c++
    /// void CObjectOStreamJson::EndBlock(void)
    /// {
    ///     m_Output.DecIndentLevel();
    ///     m_Output.PutEol();
    ///     m_Output.PutChar('}');
    ///     m_BlockStart = false;
    ///     m_ExpectValue = false;
    /// }
    /// ```
    fn end_block(&mut self) {
        self.indent -= 2;
        self.eol();
        self.out.push(b'}');
        self.block_start = false;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:829-838
    /// ```c++
    /// void CObjectOStreamJson::NextElement(void)
    /// {
    ///     if ( m_BlockStart ) {
    ///         m_BlockStart = false;
    ///     } else {
    ///         m_Output.PutChar(',');
    ///     }
    ///     m_Output.PutEol();
    ///     m_ExpectValue = true;
    /// }
    /// ```
    fn next_element(&mut self) {
        if self.block_start {
            self.block_start = false;
        } else {
            self.out.push(b',');
        }
        self.eol();
    }

    /// `x_WriteString`: the string in quotes, each byte through `WriteEncodedChar`, up to
    /// the first NUL byte (`*src`); a VisibleString is read as UTF-8 (`m_StringEncoding`),
    /// so a byte of 0x80 or more is copied.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:771-778
    /// ```c++
    /// void CObjectOStreamJson::x_WriteString(const string& value, EStringType type)
    /// {
    ///     m_Output.PutChar('\"');
    ///     for (const char* src = value.c_str(); *src; ++src) {
    ///         WriteEncodedChar(src,type);
    ///     }
    ///     m_Output.PutChar('\"');
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:756-769
    /// ```c++
    /// void CObjectOStreamJson::WriteEncodedChar(const char*& src, EStringType type)
    /// {
    ///     EEncoding enc_in( type == eStringTypeUTF8 ? eEncoding_UTF8 : m_StringEncoding);
    ///     EEncoding enc_out(eEncoding_UTF8);
    ///
    ///     if (enc_in == enc_out || enc_in == eEncoding_Unknown || (*src & 0x80) == 0) {
    ///         WriteEscapedChar(*src, enc_in);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:731-754
    /// ```c++
    /// void CObjectOStreamJson::WriteEscapedChar(char c, EEncoding enc_in)
    /// {
    ///     switch ( c ) {
    ///     case '"':
    ///         m_Output.PutString("\\\"");
    ///         break;
    ///     case '\\':
    ///         m_Output.PutString("\\\\");
    ///         break;
    ///     default:
    ///         if ( (unsigned int)c <  0x20 ||
    ///             ((unsigned int)c >= 0x80 && enc_in != eEncoding_UTF8) ) {
    ///             m_Output.PutString("\\u00");
    ///             Uint1 ch = c;
    ///             unsigned hi = ch >> 4;
    ///             unsigned lo = ch & 0xF;
    ///             m_Output.PutChar(HEX[hi]);
    ///             m_Output.PutChar(HEX[lo]);
    ///         } else {
    ///             m_Output.PutChar(c);
    ///         }
    ///         break;
    ///     }
    /// }
    /// ```
    /// (`static const char* const HEX = "0123456789ABCDEF";`, objostrjson.cpp:570.)
    fn write_string(&mut self, value: &[u8]) {
        const HEX: &[u8; 16] = b"0123456789ABCDEF";
        self.out.push(b'"');
        for &byte in value {
            match byte {
                0 => break,
                b'"' => self.out.extend_from_slice(b"\\\""),
                b'\\' => self.out.extend_from_slice(b"\\\\"),
                byte if byte < 0x20 => {
                    self.out.extend_from_slice(b"\\u00");
                    self.out.push(HEX[usize::from(byte >> 4)]);
                    self.out.push(HEX[usize::from(byte & 0xf)]);
                }
                byte => self.out.push(byte),
            }
        }
        self.out.push(b'"');
    }

    /// `WriteKey`: the key with every `-` replaced by `_` (`m_PreserveKeys` is false), then
    /// the name separator.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:780-788
    /// ```c++
    /// void CObjectOStreamJson::WriteKey(const string& key)
    /// {
    ///     string s(key);
    ///     if (!m_PreserveKeys) {
    ///         NStr::ReplaceInPlace(s,"-","_");
    ///     }
    ///     x_WriteString(s);
    ///     NameSeparator();
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:858-865
    /// ```c++
    /// void CObjectOStreamJson::NameSeparator(void)
    /// {
    ///     m_Output.PutChar(':');
    ///     if (m_Output.GetUseIndentation()) {
    ///         m_Output.PutChar(' ');
    ///     }
    ///     m_ExpectValue = true;
    /// }
    /// ```
    fn write_key(&mut self, key: &str) {
        let key: Vec<u8> = key
            .bytes()
            .map(|byte| if byte == b'-' { b'_' } else { byte })
            .collect();
        self.write_string(&key);
        self.out.extend_from_slice(b": ");
    }
}

impl Encoder for JsonEncoder<'_> {
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:137-150
    /// ```c++
    /// void CObjectOStreamJson::WriteFileHeader(TTypeInfo type)
    /// {
    ///     if (!m_JsonpPrefix.empty() || !m_JsonpSuffix.empty()) {
    ///         m_Output.PutString(m_JsonpPrefix);
    ///     }
    ///     if (type->GetDataSpec() != EDataSpec::eJSON) {
    ///         m_FileHeader = true;
    ///         StartBlock();
    ///         if (!type->GetName().empty()) {
    ///             m_Output.PutEol();
    ///             WriteKey(type->GetName());
    ///         }
    ///     }
    /// }
    /// ```
    fn begin_document(&mut self, type_name: &str) {
        if self.array_element {
            return;
        }
        self.start_block();
        self.eol();
        self.write_key(type_name);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:152-166
    /// ```c++
    /// void CObjectOStreamJson::EndOfWrite(void)
    /// {
    ///     if (m_FileHeader) {
    ///         EndBlock();
    ///         m_FileHeader = false;
    ///     } else {
    ///         ...
    ///     }
    ///     if (!m_JsonpPrefix.empty() || !m_JsonpSuffix.empty()) {
    ///         m_Output.PutString(m_JsonpSuffix);
    ///     }
    ///     m_Output.PutEol();
    ///     CObjectOStream::EndOfWrite();
    /// }
    /// ```
    fn end_document(&mut self) {
        if self.array_element {
            self.indent -= 2;
            self.eol();
            return;
        }
        self.end_block();
        self.eol();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:473-479
    /// ```c++
    /// void CObjectOStreamJson::BeginClass(const CClassTypeInfo* /*classInfo*/)
    /// {
    ///     if (GetStackDepth() > 1 && FetchFrameFromTop(1).GetNotag()) {
    ///         return;
    ///     }
    ///     StartBlock();
    /// }
    /// ```
    fn begin_sequence(&mut self, _type_name: &str) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:482-488
    /// ```c++
    /// void CObjectOStreamJson::EndClass(void)
    /// {
    ///     if (GetStackDepth() > 1 && FetchFrameFromTop(1).GetNotag()) {
    ///         return;
    ///     }
    ///     EndBlock();
    /// }
    /// ```
    fn end_sequence(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:517-521
    /// ```c++
    ///     if (id.HasAnyContent()) {
    ///         return;
    ///     }
    ///     NextElement();
    ///     WriteMemberId(id);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:715-719
    /// ```c++
    /// void CObjectOStreamJson::WriteMemberId(const CMemberId& id)
    /// {
    ///     WriteKey(id.GetName());
    ///     m_SkippedMemberId.erase();
    /// }
    /// ```
    fn member(&mut self, name: &str, _tag: u8) {
        self.next_element();
        self.write_key(name);
    }

    fn end_member(&mut self) {}

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:442-452
    /// ```c++
    /// void CObjectOStreamJson::BeginContainer(const CContainerTypeInfo* containerType)
    /// {
    ///     CObjectTypeInfo type(GetRealTypeInfo(containerType->GetElementType()));
    ///     if (type.GetTypeFamily() == eTypeFamilyPrimitive && type.GetPrimitiveValueType() == ePrimitiveValueAny) {
    ///         ...
    ///     }
    ///     BeginArray();
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:840-847
    /// ```c++
    /// void CObjectOStreamJson::BeginArray(void)
    /// {
    ///     BeginValue();
    ///     m_Output.PutChar('[');
    ///     m_Output.IncIndentLevel();
    ///     m_BlockStart = true;
    ///     m_ExpectValue = false;
    /// }
    /// ```
    fn begin_list(&mut self, _set: bool) {
        self.out.push(b'[');
        self.indent += 2;
        self.block_start = true;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:463-466
    /// ```c++
    /// void CObjectOStreamJson::BeginContainerElement(TTypeInfo /*elementType*/)
    /// {
    ///     NextElement();
    /// }
    /// ```
    fn element(&mut self) {
        self.next_element();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:849-856
    /// ```c++
    /// void CObjectOStreamJson::EndArray(void)
    /// {
    ///     m_Output.DecIndentLevel();
    ///     m_Output.PutEol();
    ///     m_Output.PutChar(']');
    ///     m_BlockStart = false;
    ///     m_ExpectValue = false;
    /// }
    /// ```
    fn end_list(&mut self) {
        self.indent -= 2;
        self.eol();
        self.out.push(b']');
        self.block_start = false;
    }

    /// A CHOICE is an object with one key, the variant's name.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:533-559
    /// ```c++
    /// void CObjectOStreamJson::BeginChoice(const CChoiceTypeInfo* /*choiceType*/)
    /// {
    ///     if (GetStackDepth() > 1 && FetchFrameFromTop(1).GetNotag()) {
    ///         return;
    ///     }
    ///     StartBlock();
    /// }
    /// ...
    /// void CObjectOStreamJson::BeginChoiceVariant(const CChoiceTypeInfo* /*choiceType*/,
    ///                                             const CMemberId& id)
    /// {
    ///     if (id.HasNotag() || id.IsAttlist()) {
    ///         ...
    ///     }
    ///     NextElement();
    ///     WriteMemberId(id);
    /// }
    /// ```
    fn variant(&mut self, _type_name: &str, name: &str, _tag: u8) {
        self.start_block();
        self.next_element();
        self.write_key(name);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:541-547
    /// ```c++
    /// void CObjectOStreamJson::EndChoice(void)
    /// {
    ///     if (GetStackDepth() > 1 && FetchFrameFromTop(1).GetNotag()) {
    ///         return;
    ///     }
    ///     EndBlock();
    /// }
    /// ```
    fn end_variant(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:180-183
    /// ```c++
    /// void CObjectOStreamJson::WriteInt4(Int4 data)
    /// {
    ///     WriteKeywordValue(NStr::IntToString(data));
    /// }
    /// ```
    fn int(&mut self, value: i32) {
        let mut digits = itoa_i32(value);
        self.out.extend_from_slice(digits.text());
    }

    /// `WriteInt8`: the decimal digits, as `WriteInt4`.
    fn big_int(&mut self, value: i64) {
        self.out.extend_from_slice(value.to_string().as_bytes());
    }

    /// An ENUMERATED value is its name as a string (not an INTEGER with named values).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:382-397
    /// ```c++
    /// void CObjectOStreamJson::WriteEnum(const CEnumeratedTypeValues& values,
    ///                         TEnumValueType value)
    /// {
    ///     string value_str;
    ///     if (values.IsInteger()) {
    ///         ...
    ///     } else {
    ///         value_str = values.FindNameEx(value, values.IsInteger());
    ///         WriteValue(value_str);
    ///     }
    /// }
    /// ```
    fn enumerated(&mut self, name: &str, _value: i32) {
        self.write_string(name.as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:168-171
    /// ```c++
    /// void CObjectOStreamJson::WriteBool(bool data)
    /// {
    ///     WriteKeywordValue( data ? "true" : "false");
    /// }
    /// ```
    fn boolean(&mut self, value: bool) {
        self.out
            .extend_from_slice(if value { b"true" } else { b"false" });
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:232-236
    /// ```c++
    /// void CObjectOStreamJson::WriteString(const string& str,
    ///                             EStringType type)
    /// {
    ///     WriteValue(str,type);
    /// }
    /// ```
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string(value);
    }

    /// `WriteDouble2(data, DBL_DIG)` (`real::write_json_double`), which raises an error for
    /// a NaN or an infinity.
    fn real(&mut self, value: f64) {
        if let Err(error) = write_json_double(value, self.out) {
            if self.error.is_none() {
                self.error = Some(SerialError::JsonDouble(error));
            }
        }
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }
}

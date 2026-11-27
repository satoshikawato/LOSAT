//! NCBI's JSON object stream (`CObjectOStreamJson`, outfmt 12) for the calls of `Encoder`.
//!
//! Each document is written by a new stream (`serialobject.cpp` `WriteObject`, see
//! `asn_text`), so the indentation starts at 0 for every document, and a document follows
//! the previous one without a separator.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:68-80
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
    /// The error of a value that the JSON documents do not hold (`Encoder::null`, `utf8_str`,
    /// `delayed`).
    fn unwritten(&mut self, what: &'static str) {
        if self.error.is_none() {
            self.error = Some(SerialError::Unwritten(what));
        }
    }

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

    /// The length of the output written so far (`XmlEncoder::written`).
    pub fn written(&self) -> usize {
        self.out.len()
    }

    /// `COStreamBuffer::PutEol` with its flag true (`asn_text::AsnTextEncoder::eol`).
    fn eol(&mut self) {
        self.out.push(b'\n');
        let indent = self.indent;
        self.out.resize(self.out.len() + indent, b' ');
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:811-818
    /// (`BeginValue` writes a skipped member's key, which only XML-like types have.)
    fn start_block(&mut self) {
        self.out.push(b'{');
        self.indent += 2;
        self.block_start = true;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:820-827
    fn end_block(&mut self) {
        self.indent -= 2;
        self.eol();
        self.out.push(b'}');
        self.block_start = false;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:829-838
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:756-769
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:731-754
    /// (The table of upper-case hex digits, objostrjson.cpp:570.)
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:858-865
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
    const STREAM: Option<super::encoder::Stream> = Some(super::encoder::Stream::Json);

    /// `m_EnforceWritingDefaults` is set to `true` (the constructor, cited in the module).
    fn writes_defaults(&self) -> bool {
        true
    }

    /// The document so far (`Traced` follows NCBI's buffer over it).
    fn held(&self) -> &[u8] {
        self.out
    }

    fn hand_on(&mut self, sink: &mut dyn std::io::Write, count: usize) -> std::io::Result<()> {
        sink.write_all(&self.out[..count])?;
        self.out.drain(..count);
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:137-150
    fn begin_document(&mut self, type_name: &str) {
        if self.array_element {
            return;
        }
        self.start_block();
        self.eol();
        self.write_key(type_name);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:152-166
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
    fn begin_sequence(&mut self, _type_name: &str) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:482-488
    fn end_sequence(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:517-521
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:715-719
    fn member(&mut self, name: &str, _tag: u8) {
        self.next_element();
        self.write_key(name);
    }

    fn end_member(&mut self) {}

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:442-452
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:840-847
    fn begin_list(&mut self, _set: bool) {
        self.out.push(b'[');
        self.indent += 2;
        self.block_start = true;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:463-466
    fn element(&mut self) {
        self.next_element();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:849-856
    fn end_list(&mut self) {
        self.indent -= 2;
        self.eol();
        self.out.push(b']');
        self.block_start = false;
    }

    /// A CHOICE is an object with one key, the variant's name.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:533-559
    fn variant(&mut self, _type_name: &str, name: &str, _tag: u8) {
        self.start_block();
        self.next_element();
        self.write_key(name);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:541-547
    fn end_variant(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:180-183
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
    fn enumerated(&mut self, name: &str, _value: i32) {
        self.write_string(name.as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:382-393
    fn named_int(&mut self, name: &str, _value: i32) {
        self.write_string(name.as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:168-171
    fn boolean(&mut self, value: bool) {
        self.out
            .extend_from_slice(if value { b"true" } else { b"false" });
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:232-236
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string(value);
    }

    /// The JSON documents (outfmt 12, 13, 15) have no NULL, UTF8String or delayed member (the
    /// values of a PSSM file's members, `asn_types::write`).
    fn null(&mut self) {
        self.unwritten("a NULL value");
    }
    fn null_element(&mut self, _element_type: &str) {
        if self.error.is_none() {
            self.error = Some(SerialError::NullElement);
        }
    }
    fn utf8_str(&mut self, _value: &[u8]) {
        self.unwritten("a UTF8String");
    }
    /// The JSON stream copies no delay buffer (`data_format` is none): `asn_types::write`
    /// writes the member by its values.
    fn delayed(&mut self, _bytes: &[u8]) {
        unreachable!("no delay buffer is of the JSON stream's format")
    }
    fn unread_delayed(&mut self) {
        self.unwritten("a member read into a delay buffer");
    }

    /// The JSON documents (outfmt 12, 13, 15) have no OCTET STRING.
    fn octets(&mut self, _value: &[u8]) {
        if self.error.is_none() {
            self.error = Some(SerialError::OctetString);
        }
    }

    /// `WriteDouble2` of `data` with `DBL_DIG` (`real::write_json_double`), which raises an error
    /// for a NaN or an infinity.
    fn real(&mut self, value: f64) {
        if let Err(error) = write_json_double(value, self.out) {
            if self.error.is_none() {
                self.error = Some(SerialError::JsonDouble(error));
            }
        }
    }

    fn buffered(&self) -> usize {
        self.out.len()
    }

    /// The JSON stream never looks back at what it wrote: every byte is final, and all of them
    /// go (with what preceded the document in `out`).
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, _all: bool) -> std::io::Result<()> {
        sink.write_all(self.out)?;
        self.out.clear();
        Ok(())
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }
}

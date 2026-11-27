//! The metadata file of a database (`.pjs` / `.njs`) that makeblastdb writes after a
//! successful build: `CSeqDB::GetDBMetaData` (`c++/src/objtools/blast/seqdb_reader/seqdb.cpp`)
//! serialized as JSON by `CObjectOStreamJson` (`c++/src/serial/objostrjson.cpp`).

use std::io::Write;
use std::path::Path;

use super::build_db::VOLUME_EXTENSIONS;
use super::ofstream::{BasicStream, IoState, OfStream, StreamBuf, IOS_FAILURE_WHAT};

/// `CBlast_db_metadata` (`c++/src/objects/blastdb/blastdb.asn:79-97`) as makeblastdb fills it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DbMetadata {
    pub dbname: Vec<u8>,
    pub dbtype: &'static str,
    pub db_version: u32,
    pub description: Vec<u8>,
    pub number_of_letters: u64,
    /// `Int8` (blastdb.def gives `number-of-sequences` the type `Int8`), the count of every volume;
    /// NCBI fills it from an `int` (`GetDBMetaData`), which LOSAT does not narrow.
    pub number_of_sequences: i64,
    pub last_updated: String,
    pub number_of_volumes: i32,
    /// Set for version 5 when the database has a tax id other than 0 (`-taxid`,
    /// `-taxid_map`); never without them.
    pub number_of_taxids: Option<i32>,
    pub bytes_total: i64,
    pub bytes_to_cache: i64,
    /// The names of the files: bytes, as NCBI's `std::string` holds the name of an LMDB file
    /// that is not UTF-8 (line LIB-B).
    pub files: Vec<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:731-752
// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:755-777
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7080
// makeblastdb sets the string encoding to `eEncoding_Ascii`
// (makeblastdb.cpp:1320) and every member is a VisibleString: a byte from 0x80
// goes through `CUtf8::AsUTF8` with `eEncoding_Ascii`, which appends it as it is,
// and a control character is written as `\u00XX` with capital hexadecimal
// digits. The string ends at its first NUL (`*src`). makeblastdb 2.17.0 wrote the
// title `é` (UTF-8 C3 A9) as the same two bytes, a tab as `\u0009`.
// Each byte is its own `PutChar`, an escape its own `PutString`: these calls are the points
// where the serializer's buffer may be written to the file (`OStreamBuffer`).
/// `x_WriteString`: a JSON string of `value`.
fn write_json_string<O: JsonOut>(out: &mut O, value: &[u8]) -> Result<(), O::Error> {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    out.put_char(b'"')?;
    for &c in value.iter().take_while(|&&c| c != 0) {
        match c {
            b'"' => out.put_string(b"\\\"")?,
            b'\\' => out.put_string(b"\\\\")?,
            c if c < 0x20 => {
                out.put_string(b"\\u00")?;
                out.put_char(HEX[(c >> 4) as usize])?;
                out.put_char(HEX[(c & 0xF) as usize])?;
            }
            c => out.put_char(c)?,
        }
    }
    out.put_char(b'"')
}

/// What `CObjectOStreamJson` calls on its `COStreamBuffer` for the metadata object: the
/// serializer's buffer (`OStreamBuffer`), or the bytes alone (`JsonText`).
pub(crate) trait JsonOut {
    type Error;
    /// `PutChar` of `c`.
    fn put_char(&mut self, c: u8) -> Result<(), Self::Error>;
    /// `PutString` of `str` and its length.
    fn put_string(&mut self, s: &[u8]) -> Result<(), Self::Error>;
    /// `PutEol()`: a line feed and the indentation.
    fn put_eol(&mut self) -> Result<(), Self::Error>;
    /// `IncIndentLevel()` (two spaces).
    fn inc_indent_level(&mut self);
    /// `DecIndentLevel()`.
    fn dec_indent_level(&mut self);
}

/// The bytes of the serialized object (the tests' view of the document).
#[cfg(test)]
#[derive(Debug, Default)]
pub(crate) struct JsonText {
    pub(crate) bytes: Vec<u8>,
    indent_level: usize,
}

#[cfg(test)]
impl JsonOut for JsonText {
    type Error = std::convert::Infallible;
    fn put_char(&mut self, c: u8) -> Result<(), Self::Error> {
        self.bytes.push(c);
        Ok(())
    }
    fn put_string(&mut self, s: &[u8]) -> Result<(), Self::Error> {
        self.bytes.extend_from_slice(s);
        Ok(())
    }
    fn put_eol(&mut self) -> Result<(), Self::Error> {
        self.bytes.push(b'\n');
        self.bytes
            .resize(self.bytes.len() + self.indent_level, b' ');
        Ok(())
    }
    fn inc_indent_level(&mut self) {
        self.indent_level += 2;
    }
    fn dec_indent_level(&mut self) {
        self.indent_level -= 2;
    }
}

impl DbMetadata {
    // NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:811-838
    // NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:858-865
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1310-1318
    // The members go out in the order of blastdb.asn:79-97 with their ASN.1 names
    // (`PreserveKeyNames`), two spaces of indentation per level, and a line feed
    // after the closing brace (`NcbiEndl`). The members `version` and `db-version` are
    // set by the object's defaults ("1.2", and 5 or 4 by `SetDb_version`); the
    // compressed-file members are not set.
    // NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:490-522,780-809,840-856
    // `WriteObject` of the object is `CObjectOStream::WriteClass` (objostr.cpp:892-904): its
    // `BeginClass` is `StartBlock`, each member is `BeginClassMember`, its value, and
    // `EndClassMember`, then `EndBlock`. The numbers are `WriteKeywordValue` of their decimal
    // text; `files` is a container of strings (`BeginContainer` is `BeginArray`, each element
    // `NextElement` and `WriteValue`). `BeginValue` writes nothing here: the object has no
    // member without a tag.
    /// `WriteObject` of `json_out` for `obj_info`: the metadata object, through `out`'s calls.
    pub(crate) fn write_json<O: JsonOut>(&self, out: &mut O) -> Result<(), O::Error> {
        let mut block_start = true;
        let mut key = |out: &mut O, name: &str| -> Result<(), O::Error> {
            if block_start {
                block_start = false;
            } else {
                out.put_char(b',')?;
            }
            out.put_eol()?;
            write_json_string(out, name.as_bytes())?;
            out.put_char(b':')?;
            out.put_char(b' ')
        };
        out.put_char(b'{')?;
        out.inc_indent_level();
        key(out, "version")?;
        write_json_string(out, b"1.2")?;
        key(out, "dbname")?;
        write_json_string(out, &self.dbname)?;
        key(out, "dbtype")?;
        write_json_string(out, self.dbtype.as_bytes())?;
        key(out, "db-version")?;
        out.put_string(self.db_version.to_string().as_bytes())?;
        key(out, "description")?;
        write_json_string(out, &self.description)?;
        key(out, "number-of-letters")?;
        out.put_string(self.number_of_letters.to_string().as_bytes())?;
        key(out, "number-of-sequences")?;
        out.put_string(self.number_of_sequences.to_string().as_bytes())?;
        key(out, "last-updated")?;
        write_json_string(out, self.last_updated.as_bytes())?;
        key(out, "number-of-volumes")?;
        out.put_string(self.number_of_volumes.to_string().as_bytes())?;
        if let Some(taxids) = self.number_of_taxids {
            key(out, "number-of-taxids")?;
            out.put_string(taxids.to_string().as_bytes())?;
        }
        key(out, "bytes-total")?;
        out.put_string(self.bytes_total.to_string().as_bytes())?;
        key(out, "bytes-to-cache")?;
        out.put_string(self.bytes_to_cache.to_string().as_bytes())?;
        key(out, "files")?;
        out.put_char(b'[')?;
        out.inc_indent_level();
        for (i, file) in self.files.iter().enumerate() {
            if i > 0 {
                out.put_char(b',')?;
            }
            out.put_eol()?;
            write_json_string(out, file)?;
        }
        out.dec_indent_level();
        out.put_eol()?;
        out.put_char(b']')?;
        out.dec_indent_level();
        out.put_eol()?;
        out.put_char(b'}')
    }

    /// The bytes of the metadata file: the object and the line feed of `NcbiEndl`.
    #[cfg(test)]
    pub(crate) fn to_json(&self) -> Vec<u8> {
        let mut text = JsonText {
            bytes: Vec::with_capacity(512 + 32 * self.files.len()),
            indent_level: 0,
        };
        let Ok(()) = self.write_json(&mut text);
        text.bytes.push(b'\n');
        text.bytes
    }
}

/// `GetName` of a `CDirEntry` of `path`: the part after the last path separator.
fn file_name(path: &[u8]) -> &[u8] {
    let is_separator = |c: &u8| *c == b'/' || (cfg!(windows) && (*c == b'\\' || *c == b':'));
    match path.iter().rposition(is_separator) {
        Some(at) => &path[at + 1..],
        None => path,
    }
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:4039-4042
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3204-3216
// `IsFile` follows links. Between `Exists` and `GetLength` NCBI posts an error when
// the length cannot be read; a file whose `stat` succeeded a moment before gives its
// length here (`None` only when it is not a regular file).
/// `Exists` and `GetLength` of a `CFile` of `path`: the length of the regular file `path`.
fn file_length(path: impl AsRef<Path>) -> Option<i64> {
    let meta = std::fs::metadata(path).ok()?;
    meta.is_file().then_some(meta.len() as i64)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1610-1702
// `volumes` and `alias_files` are `FindVolumePaths` of `paths` and `alias` with the last argument
// true: each path once. `user_path` is the `-metadata_output_prefix` of makeblastdb and of
// blastdbcmd (`metadata_user_path`). makeblastdb's `lmdbs` is the one LMDB file of the version 5
// database it wrote, `<db>.pdb`; blastdbcmd's are the LMDB files of the database's volumes
// (`GetLMDBFileNames`, one per `CSeqDBLMDBEntry`). The `GetDir` of `file` followed by its `GetBase`
// is the LMDB file's path without its extension. A file the earlier database of the same name left
// (a `<db>.pjs` that the removal of a database of several volumes does not take) is listed and
// counted as NCBI lists it.
/// `x_GetDBFilesMetaData`: the bytes of the database's files, the bytes of its index
/// and sequence files, and the sorted names of its files, each after `user_path`. An LMDB
/// file's name is the bytes of the volume's header (`lmdb_file_names`), which need not be
/// UTF-8; NCBI finds, counts and lists the files of such a name (line LIB-B).
pub(crate) fn db_files_metadata(
    volumes: &[Vec<u8>],
    alias_files: &[Vec<u8>],
    lmdbs: &[std::path::PathBuf],
    is_protein: bool,
    user_path: &[u8],
) -> (i64, i64, Vec<Vec<u8>>) {
    use crate::utils::os_bytes::{path_bytes, path_from_bytes, path_of_bytes};
    let mut db_files = Vec::new();
    let mut disk_bytes = 0i64;
    let mut cached_bytes = 0i64;
    let listed = |name: &[u8]| [user_path, file_name(name)].concat();
    for alias in alias_files {
        if let Some(length) = file_length(path_of_bytes(alias)) {
            db_files.push(listed(alias));
            disk_bytes += length;
        }
    }
    let mol = if is_protein { 'p' } else { 'n' };
    for path in volumes {
        for ext in VOLUME_EXTENSIONS {
            let file = [path.as_slice(), format!(".{mol}{ext}").as_bytes()].concat();
            if let Some(length) = file_length(path_of_bytes(&file)) {
                db_files.push(listed(&file));
                disk_bytes += length;
                if ext == "in" || ext == "sq" {
                    cached_bytes += length;
                }
            }
        }
    }
    for lmdb in lmdbs {
        if let Some(length) = file_length(lmdb) {
            let lmdb = path_bytes(lmdb);
            db_files.push(listed(&lmdb[..]));
            disk_bytes += length;
            let base = lmdb
                .iter()
                .rposition(|&b| b == b'.')
                .map_or(&lmdb[..], |dot| &lmdb[..dot]);
            for ext in ["os", "ot", "tf", "to"] {
                let file = [base, format!(".{mol}{ext}").as_bytes()].concat();
                if let Some(length) = file_length(path_from_bytes(file.as_slice())) {
                    db_files.push(listed(&file[..]));
                    disk_bytes += length;
                }
            }
        }
    }
    db_files.sort();
    (disk_bytes, cached_bytes, db_files)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1301-1305
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1208-1212
/// The `user_path` of `-metadata_output_prefix`, which makeblastdb and blastdbcmd read alike: the
/// value with the path separator added when it is not empty and does not end with one, and
/// empty without the option. makeblastdb 2.17.0 writes the value before each name of `files`
/// (`x`, `/data/blast/`, a quote, a backslash, a tab, UTF-8, escaped as the JSON of every
/// string), and an empty value adds nothing (line RG item D3a probes `P_*`).
pub(crate) fn metadata_user_path(prefix: Option<&[u8]>) -> Vec<u8> {
    let mut user_path = prefix.unwrap_or_default().to_vec();
    let separator = std::path::MAIN_SEPARATOR as u8;
    if !user_path.is_empty() && user_path.last() != Some(&separator) {
        user_path.push(separator);
    }
    user_path
}

/// `CIOException::EErrCode` of the serializer's buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IoErrCode {
    /// `eWrite`.
    Write,
    /// `eFlush`.
    Flush,
}

/// A `CIOException` that `COStreamBuffer` throws: its code and message (`GetMsg()`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IoException {
    pub(crate) code: IoErrCode,
    pub(crate) msg: &'static str,
}

/// What a call of the serializer's buffer throws: its `CIOException`, or the `ios_base::failure`
/// of an output stream set to throw on `badbit` (`BasicStream::set_exceptions`, blast_formatter's
/// `-out`), thrown by the stream's write or flush before the buffer looks at the state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum BufferError {
    Io(IoException),
    /// `ios_base::failure`, whose `what()` is `IOS_FAILURE_WHAT`.
    Stream,
}

impl From<IoException> for BufferError {
    fn from(exception: IoException) -> Self {
        BufferError::Io(exception)
    }
}

impl IoException {
    /// `GetType()`, then `"::"`, then `GetErrCodeString()`.
    fn type_and_code(&self) -> &'static str {
        match self.code {
            IoErrCode::Write => "CIOException::eWrite",
            IoErrCode::Flush => "CIOException::eFlush",
        }
    }
}

// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:46,783-795
// The JSON writer sets no back limit and no cancel callback, and makeblastdb and blastdbcmd give
// it the stream without ownership (`eNoOwnership`). The line count and line length serve the
// line positions and the wrapping of other formats; they write nothing here.
/// `COStreamBuffer` over the stream the metadata goes to: makeblastdb's metadata file
/// (`OfStream`), blastdbcmd's `-out` file or standard output (`super::cout`); and over
/// blast_formatter's `-out` (`blast_formatter::archive`, a Seq-annot's bytes through
/// `write_bytes`).
pub(crate) struct OStreamBuffer<'a, B: StreamBuf> {
    /// `m_Output`.
    output: &'a mut BasicStream<B>,
    /// `m_Closed`.
    closed: bool,
    /// `m_Error`.
    error: Option<&'static str>,
    /// `m_IndentLevel`.
    indent_level: usize,
    /// `m_BufferPos`: the bytes given to the stream.
    buffer_pos: u64,
    /// `m_Buffer` up to `m_BufferEnd`.
    buffer: Vec<u8>,
    /// `m_CurrentPos` less `m_Buffer`.
    current: usize,
}

impl<'a, B: StreamBuf> OStreamBuffer<'a, B> {
    pub(crate) fn new(output: &'a mut BasicStream<B>) -> Self {
        Self {
            output,
            closed: false,
            error: None,
            indent_level: 0,
            buffer_pos: 0,
            buffer: vec![0; 4096],
            current: 0,
        }
    }

    /// `m_Output`.
    pub(crate) fn output(&mut self) -> &mut BasicStream<B> {
        self.output
    }

    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:185-229
    // `fail()` is `error.is_some()`.
    /// `GetStreamPos()`.
    fn stream_pos(&self) -> u64 {
        self.buffer_pos + self.current as u64
    }

    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:283-312
    // `DoSkip` (strbuffer.inl:273-281) is `DoReserve` and the same advance.
    /// `Skip` of `count`: the position of the `count` bytes taken.
    fn skip(&mut self, count: usize) -> Result<usize, BufferError> {
        let pos = self.reserve(count)?;
        self.current = pos + count;
        Ok(pos)
    }

    /// `Reserve(count)`: the position where `count` bytes fit.
    fn reserve(&mut self, count: usize) -> Result<usize, BufferError> {
        if self.current + count <= self.buffer.len() {
            Ok(self.current)
        } else {
            self.do_reserve(count)
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:921-948
    // `BiggerBufferSize` doubles (strbuffer.cpp:48-52).
    fn do_reserve(&mut self, count: usize) -> Result<usize, BufferError> {
        self.flush_buffer(false)?;
        let need_size = self.current + count;
        let mut buffer_size = self.buffer.len();
        if buffer_size < need_size {
            while buffer_size < need_size {
                buffer_size *= 2;
            }
            self.buffer.resize(buffer_size, 0);
        }
        Ok(self.current)
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-907
    // With no back limit, both forms write all the buffer holds. A failed write keeps it there.
    /// `FlushBuffer` with `fullBuffer`.
    fn flush_buffer(&mut self, _full_buffer: bool) -> Result<(), BufferError> {
        let count = self.current;
        if count != 0 {
            self.output.write(&self.buffer[..count]);
            if self.output.take_thrown() {
                return Err(BufferError::Stream);
            }
            if self.output.fail() {
                self.error = Some("write fault");
                return Err(IoException {
                    code: IoErrCode::Write,
                    msg: "write fault",
                }
                .into());
            }
            self.current = 0;
            self.buffer_pos += count as u64;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:813-837,910-918
    // `eFlush` sets no `m_Error`. The guard's destructor restores the state and catches the
    // `ios_base::failure` of a stream that throws on it.
    /// `Flush()`.
    pub(crate) fn flush(&mut self) -> Result<(), BufferError> {
        let state = self.output.rdstate();
        self.output.clear(IoState::GOOD);
        let result = self.flush_buffer(true).and_then(|()| {
            self.output.flush();
            if self.output.take_thrown() {
                Err(BufferError::Stream)
            } else if self.output.fail() {
                Err(IoException {
                    code: IoErrCode::Flush,
                    msg: "COStreamBuffer::Flush: failed",
                }
                .into())
            } else {
                Ok(())
            }
        });
        self.restore_state(state);
        result
    }

    /// The destructor of `STemporarilyClearStreamState`: `setstate` of the state it saved, whose
    /// `ios_base::failure` it catches.
    fn restore_state(&mut self, state: IoState) {
        self.output.setstate(state);
        let _ = self.output.take_thrown();
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1143-1161
    /// `Write` of the bytes `data`: the buffer filled up and flushed whenever it is full, without
    /// growing it.
    pub(crate) fn write_bytes(&mut self, mut data: &[u8]) -> Result<(), BufferError> {
        while !data.is_empty() {
            let mut available = self.buffer.len() - self.current;
            if available == 0 {
                self.flush_buffer(false)?;
                available = self.buffer.len() - self.current;
            }
            if available >= data.len() {
                break;
            }
            let (now, rest) = data.split_at(available);
            self.buffer[self.current..self.current + available].copy_from_slice(now);
            self.current += available;
            data = rest;
        }
        self.buffer[self.current..self.current + data.len()].copy_from_slice(data);
        self.current += data.len();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:559-564
    /// `EndOfWrite`'s `FlushBuffer`: what the buffer holds goes to the stream.
    pub(crate) fn end_of_write(&mut self) -> Result<(), BufferError> {
        self.flush_buffer(true)
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:843-862
    // `m_Output` as a condition is true while the stream's `fail()` is false.
    /// `Close()`.
    fn close(&mut self) -> Result<(), BufferError> {
        if !self.closed && !self.output.fail() {
            self.closed = true;
            let state = self.output.rdstate();
            self.output.clear(IoState::GOOD);
            let result = self.flush_buffer(true);
            self.restore_state(state);
            result?;
        }
        self.error = None;
        self.indent_level = 0;
        self.current = 0;
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:251-263,314-381
// `IncIndentLevel` and `DecIndentLevel` step by 2 (strbuffer.hpp:244-245). The strings of the
// metadata object are written byte by byte (`x_WriteString`); `PutString` gets the numbers
// and the escapes, all shorter than 1024 bytes.
impl<B: StreamBuf> JsonOut for OStreamBuffer<'_, B> {
    type Error = BufferError;
    fn put_char(&mut self, c: u8) -> Result<(), BufferError> {
        let pos = self.skip(1)?;
        self.buffer[pos] = c;
        Ok(())
    }
    fn put_string(&mut self, s: &[u8]) -> Result<(), BufferError> {
        let pos = self.skip(s.len())?;
        self.buffer[pos..pos + s.len()].copy_from_slice(s);
        Ok(())
    }
    fn put_eol(&mut self) -> Result<(), BufferError> {
        let pos = self.reserve(1)?;
        self.buffer[pos] = b'\n';
        self.current = pos + 1;
        let count = self.indent_level;
        let pos = self.skip(count)?;
        self.buffer[pos..pos + count].fill(b' ');
        Ok(())
    }
    fn inc_indent_level(&mut self) {
        self.indent_level += 2;
    }
    fn dec_indent_level(&mut self) {
        self.indent_level -= 2;
    }
}

// NCBI reference (598d8ae6): c++/include/serial/objostr.hpp:312-328
const F_WRITE_ERROR: u32 = 1 << 1;
const F_NOT_OPEN: u32 = 1 << 7;

/// The parts of `CObjectOStream` that act on its stream after the object is written: the fail
/// flags and the close of its destructor, with the posts of the application's diagnostic prefix
/// `prefix` (makeblastdb's and blastdbcmd's `[program] `; blast_formatter's posts have none).
pub(crate) struct ObjectOStream<'a, B: StreamBuf> {
    /// `m_Output`.
    pub(crate) output: OStreamBuffer<'a, B>,
    /// `m_Fail`.
    fail: u32,
    /// `SetDiagPostPrefix` of the application.
    prefix: Option<&'static str>,
}

/// What `CObjectOStream::Close` throws: the `CSerialException` of its handler, or an exception
/// that is not a `CException` (the stream's `ios_base::failure`), which its handler lets through.
enum Closed {
    Serial(SerialException),
    Stream,
}

/// The `CSerialException` that `ThrowError1` throws over the exception it was given.
struct SerialException {
    /// `byte N: <message>`.
    msg: String,
    predecessor: IoException,
}

impl<'a, B: StreamBuf> ObjectOStream<'a, B> {
    pub(crate) fn new(output: OStreamBuffer<'a, B>, prefix: Option<&'static str>) -> Self {
        ObjectOStream {
            output,
            fail: 0,
            prefix,
        }
    }

    /// The prefix of the posts: `[<prefix>] `, or nothing.
    fn tag(&self) -> String {
        self.prefix.map(|p| format!("[{p}] ")).unwrap_or_default()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384,495-513
    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-84
    // NCBI reference (598d8ae6): c++/include/serial/error_codes.hpp:46
    // Every call here comes from the close in the destructor: the frames of the object are
    // gone (an exception pops them on its way out), and in a base class's destructor the
    // virtual `GetPosition` is `CObjectOStream`'s, the byte position, not the line position
    // of `CObjectOStreamJson`. The diagnostics are "Severity: (code.subcode) [program] text"
    // (blast_formatter's have no prefix).
    /// `SetFailFlags` with `flags` and a message.
    fn set_fail_flags(&mut self, flags: u32, message: &str, err: &mut dyn Write) {
        let old = self.fail;
        self.fail |= flags;
        if old == 0 && flags != 0 {
            let header = format!("Error: (802.5) {}", self.tag());
            let line = format!(
                "CObjectOStream: error at byte {}: stack is empty: {message}",
                self.output.stream_pos()
            );
            self.post(err, &header, &line, false);
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:386-400
    fn in_good_state(&mut self, err: &mut dyn Write) -> bool {
        if self.fail != 0 {
            false
        } else if let Some(error) = self.output.error {
            self.set_fail_flags(F_WRITE_ERROR, error, err);
            self.output.error = None;
            false
        } else {
            true
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:522-556
    // `RethrowError` with a flag, a message and `exc` is `ThrowError1` with `DIAG_COMPILE_INFO`,
    // the flag, the message and the address of `exc`
    // (objostr.hpp:578-579). `DefaultFlush` is `Flush` (no `fFlagNoAutoFlush`).
    /// `RethrowError` with `fWriteError`, a message and `exc`.
    fn rethrow_write_error(
        &mut self,
        message: &str,
        exc: IoException,
        err: &mut dyn Write,
    ) -> SerialException {
        let _ = self.output.flush();
        self.set_fail_flags(F_WRITE_ERROR, message, err);
        SerialException {
            msg: format!("byte {}: {message}", self.output.stream_pos()),
            predecessor: exc,
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:327-355
    // The handler takes `CException`s only: the stream's `ios_base::failure` leaves `Close`.
    /// `Close()`.
    fn close(&mut self, err: &mut dyn Write) -> Result<(), Closed> {
        if self.fail != F_NOT_OPEN {
            let result = self.output.flush().and_then(|()| {
                self.fail = F_NOT_OPEN;
                self.output.close()
            });
            match result {
                Err(BufferError::Io(exc)) if self.in_good_state(err) => {
                    let exc = self.rethrow_write_error("cannot close output stream", exc, err);
                    return Err(Closed::Serial(exc));
                }
                Err(BufferError::Stream) => return Err(Closed::Stream),
                _ => {}
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:317-324
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:676-700
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7873-7893
    // `CNcbiDiag::x_Put` (ncbidiag.cpp:7894-7939) posts each exception of the chain, the oldest
    // first, as a message of its own whose code is its type and error code, before the message
    // of the report itself ("title" and the newest exception). The predecessor of the
    // `CSerialException` is the `CIOException` of the failed flush, which has no predecessor
    // of its own; neither has extra text (`ReportExtra`).
    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:626-633,653-657
    // `NCBI_CATCH_X` reports any other `std::exception` as "[title] Exception: what()".
    /// `~CObjectOStream`.
    pub(crate) fn destroy(&mut self, err: &mut dyn Write) {
        let tag = self.tag();
        match self.close(err) {
            Err(Closed::Serial(exc)) => {
                let (kind, msg) = (exc.predecessor.type_and_code(), exc.predecessor.msg);
                self.post(err, &format!("Error: ({kind}) {tag}"), msg, true);
                let header = format!("Error: ({SERIAL_IO_ERROR}) {tag}");
                self.post(err, &header, &exc.msg, true);
                let line = format!(
                    "Cannot close serializing output stream ({SERIAL_IO_ERROR}) {}",
                    exc.msg
                );
                self.post(err, &format!("Error: (802.4) {tag}"), &line, true);
            }
            Err(Closed::Stream) => {
                let line = format!(
                    "[Cannot close serializing output stream] Exception: {IOS_FAILURE_WHAT}"
                );
                self.post(err, &format!("Error: (802.4) {tag}"), &line, true);
            }
            Ok(()) => {}
        }
        self.destroy_buffer(err);
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:797-811
    // NCBI reference (598d8ae6): c++/include/util/error_codes.hpp:48
    // After `CObjectOStream::Close` the buffer is closed, or the stream failed: this `Close`
    // writes nothing in makeblastdb's and blastdbcmd's sequences; its report is that of
    // `destroy` for an exception without predecessor.
    /// `~COStreamBuffer` (a member of the object stream, destroyed after its body).
    fn destroy_buffer(&mut self, err: &mut dyn Write) {
        let tag = self.tag();
        match self.output.close() {
            Err(BufferError::Io(exc)) => {
                let (kind, msg) = (exc.type_and_code(), exc.msg);
                self.post(err, &format!("Error: ({kind}) {tag}"), msg, true);
                let line = format!("~COStreamBuffer: exception while closing ({kind}) {msg}");
                self.post(err, &format!("Error: (205.2) {tag}"), &line, true);
            }
            Err(BufferError::Stream) => {
                let line = format!(
                    "[~COStreamBuffer: exception while closing] Exception: {IOS_FAILURE_WHAT}"
                );
                self.post(err, &format!("Error: (205.2) {tag}"), &line, true);
            }
            Ok(()) => {}
        }
    }

    /// Writes the diagnostic line of `header` and `message` (the header alone after the main
    /// thread's diagnostic stream failed, `ncbi_exception::post_text`), in LOSAT's words when
    /// `words` (`crate::wording::losat_text`: without the exception's class and code,
    /// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12). Nothing after a write found a closed
    /// pipe: SIGPIPE has ended NCBI's process there (`StreamBuf::closed_pipe`).
    fn post(&mut self, err: &mut dyn Write, header: &str, message: &str, words: bool) {
        if self.output.output().rdbuf().closed_pipe() {
            return;
        }
        let line = crate::ncbi_exception::post_text(header.as_bytes(), message.as_bytes());
        let line = String::from_utf8(line).expect("the parts are UTF-8");
        let line = if words {
            crate::wording::losat_text(&line).into_owned()
        } else {
            line
        };
        let _ = err.write_all(format!("{line}\n").as_bytes());
    }
}

/// `GetType()::GetErrCodeString()` of the serial stream's error of a failed close.
const SERIAL_IO_ERROR: &str = "CSerialException::eIoError";

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1306-1318
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
// `SeqDB_GetMetadataFileExtension` gives "pjs" or "njs". The `ofstream` and the serializer
// are those of makeblastdb (`super::ofstream`, `OStreamBuffer`): a file that cannot be opened
// or written is not refused, and the serializer reports what it finds when it is destroyed.
// The serializer's buffer holds 4096 bytes and is written when a call does not fit: a file of
// at most 4096 bytes reaches the stream at `json_out->Flush()`, a larger one first in parts of
// up to 4096 bytes during `WriteObject`. An exception (`CIOException`) leaves the block, the
// serializer is destroyed (its `Close` flushes again and reports), then the `ofstream`, and
// `CATCH_ALL` takes the exception: `eFlush` (the stream's `flush()` failed, as on a full
// device) ends makeblastdb with exit status 6, `eWrite` (a write of the buffer failed: a file
// that cannot be opened, or a full device past the stream's buffer) without a message, exit
// status 0.
/// Writes the metadata file `<db>.pjs` / `<db>.njs`, with the serializer's diagnostics to
/// `err`; `Some(message)` for the `eFlush` exception that ends makeblastdb ("BLAST failed to
/// write output: message").
pub(crate) fn write_metadata_file(
    new_db: &[u8],
    is_protein: bool,
    metadata: &DbMetadata,
    err: &mut dyn Write,
) -> Option<&'static str> {
    let extension: &[u8] = if is_protein { b".pjs" } else { b".njs" };
    let name = [new_db, extension].concat();
    let mut out = OfStream::create(&crate::utils::os_bytes::path_of_bytes(&name));
    write_metadata(&mut out, metadata, "makeblastdb", err)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1213-1220
// blastdbcmd writes the object to its `-out` stream as makeblastdb writes it to its file: the
// same serializer, `json_out->Flush()`, `out.flush()` and `NcbiEndl`, then the serializer's
// destructor; its `CATCH_ALL` is the same macro.
/// `metadata` serialized to `out` as `CObjectOStreamJson` writes it, then the line feed of
/// `NcbiEndl`, with the serializer's diagnostics (naming `program`) to `err`; `Some(message)` for
/// the `eFlush` exception that ends the program ("BLAST failed to write output: message", exit
/// status 6). After a write found a closed pipe (`StreamBuf::closed_pipe`, where SIGPIPE ends
/// NCBI's process) nothing is posted and `None` is returned.
pub(crate) fn write_metadata<B: StreamBuf>(
    out: &mut BasicStream<B>,
    metadata: &DbMetadata,
    program: &'static str,
    err: &mut dyn Write,
) -> Option<&'static str> {
    let mut json_out = ObjectOStream::new(OStreamBuffer::new(out), Some(program));
    let result = metadata
        .write_json(&mut json_out.output)
        .and_then(|()| json_out.output.flush())
        .map(|()| {
            let out = json_out.output.output();
            out.flush();
            out.endl();
        });
    json_out.destroy(err);
    if json_out.output.output().rdbuf().closed_pipe() {
        return None;
    }
    result.err().and_then(|exception| match exception {
        BufferError::Io(exception) if exception.code == IoErrCode::Flush => Some(exception.msg),
        _ => None,
    })
}

#[cfg(test)]
#[path = "metadata_tests.rs"]
mod tests;

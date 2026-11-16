//! The metadata file of a database (`.pjs` / `.njs`) that makeblastdb writes after a
//! successful build: `CSeqDB::GetDBMetaData` (`c++/src/objtools/blast/seqdb_reader/seqdb.cpp`)
//! serialized as JSON by `CObjectOStreamJson` (`c++/src/serial/objostrjson.cpp`).

use std::io::Write;
use std::path::Path;

use super::build_db::VOLUME_EXTENSIONS;
use super::ofstream::OfStream;

/// `CBlast_db_metadata` (`c++/src/objects/blastdb/blastdb.asn:79-97`) as makeblastdb fills it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DbMetadata {
    pub dbname: Vec<u8>,
    pub dbtype: &'static str,
    pub db_version: u32,
    pub description: Vec<u8>,
    pub number_of_letters: u64,
    /// `Int8` (blastdb.def: `number-of-sequences._type = Int8`), the count of every volume;
    /// NCBI fills it from an `int` (`GetDBMetaData`), which LOSAT does not narrow.
    pub number_of_sequences: i64,
    pub last_updated: String,
    pub number_of_volumes: i32,
    /// Set for version 5 when the database has a tax id other than 0 (`-taxid`,
    /// `-taxid_map`); never without them.
    pub number_of_taxids: Option<i32>,
    pub bytes_total: i64,
    pub bytes_to_cache: i64,
    pub files: Vec<String>,
}

// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:731-752
// ```c++
// void CObjectOStreamJson::WriteEscapedChar(char c, EEncoding enc_in)
// {
//     switch ( c ) {
//     case '"':
//         m_Output.PutString("\\\"");
//         break;
//     case '\\':
//         m_Output.PutString("\\\\");
//         break;
//     default:
//         if ( (unsigned int)c <  0x20 ||
//             ((unsigned int)c >= 0x80 && enc_in != eEncoding_UTF8) ) {
//             m_Output.PutString("\\u00");
//             Uint1 ch = c;
//             unsigned hi = ch >> 4;
//             unsigned lo = ch & 0xF;
//             m_Output.PutChar(HEX[hi]);
//             m_Output.PutChar(HEX[lo]);
//         } else {
//             m_Output.PutChar(c);
//         }
//         break;
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:755-777
// ```c++
// void CObjectOStreamJson::WriteEncodedChar(const char*& src, EStringType type)
// {
//     EEncoding enc_in( type == eStringTypeUTF8 ? eEncoding_UTF8 : m_StringEncoding);
//     EEncoding enc_out(eEncoding_UTF8);
//
//     if (enc_in == enc_out || enc_in == eEncoding_Unknown || (*src & 0x80) == 0) {
//         WriteEscapedChar(*src, enc_in);
//     } else {
//         CStringUTF8 tmp( CUtf8::AsUTF8( CTempString(src,1),enc_in));
//         for ( string::const_iterator t = tmp.begin(); t != tmp.end(); ++t ) {
//             m_Output.PutChar(*t);
//         }
//     }
// }
//
// void CObjectOStreamJson::x_WriteString(const string& value, EStringType type)
// {
//     m_Output.PutChar('\"');
//     for (const char* src = value.c_str(); *src; ++src) {
//         WriteEncodedChar(src,type);
//     }
//     m_Output.PutChar('\"');
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7080
// ```c++
// CStringUTF8& CUtf8::x_Append( CStringUTF8& self, const CTempString& src,
//     EEncoding encoding, EValidate validate)
// {
//     ...
//     if (encoding == eEncoding_UTF8 || encoding == eEncoding_Ascii) {
//         self.append(src);
//         return self;
//     }
// ```
// makeblastdb sets the string encoding to `eEncoding_Ascii`
// (makeblastdb.cpp:1320) and every member is a VisibleString: a byte from 0x80
// goes through `CUtf8::AsUTF8(..., eEncoding_Ascii)`, which appends it as it is,
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
    /// `PutChar(c)`.
    fn put_char(&mut self, c: u8) -> Result<(), Self::Error>;
    /// `PutString(str, length)`.
    fn put_string(&mut self, s: &[u8]) -> Result<(), Self::Error>;
    /// `PutEol()`: a line feed and the indentation.
    fn put_eol(&mut self) -> Result<(), Self::Error>;
    /// `IncIndentLevel()` (two spaces).
    fn inc_indent_level(&mut self);
    /// `DecIndentLevel()`.
    fn dec_indent_level(&mut self);
}

/// The bytes of the serialized object.
#[derive(Debug, Default)]
pub(crate) struct JsonText {
    pub(crate) bytes: Vec<u8>,
    indent_level: usize,
}

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
    // ```c++
    // void CObjectOStreamJson::StartBlock(void)
    // {
    //     BeginValue();
    //     m_Output.PutChar('{');
    //     m_Output.IncIndentLevel();
    //     m_BlockStart = true;
    //     m_ExpectValue = false;
    // }
    //
    // void CObjectOStreamJson::EndBlock(void)
    // {
    //     m_Output.DecIndentLevel();
    //     m_Output.PutEol();
    //     m_Output.PutChar('}');
    //     m_BlockStart = false;
    //     m_ExpectValue = false;
    // }
    //
    // void CObjectOStreamJson::NextElement(void)
    // {
    //     if ( m_BlockStart ) {
    //         m_BlockStart = false;
    //     } else {
    //         m_Output.PutChar(',');
    //     }
    //     m_Output.PutEol();
    //     m_ExpectValue = true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:858-865
    // ```c++
    // void CObjectOStreamJson::NameSeparator(void)
    // {
    //     m_Output.PutChar(':');
    //     if (m_Output.GetUseIndentation()) {
    //         m_Output.PutChar(' ');
    //     }
    //     m_ExpectValue = true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1310-1318
    // ```c++
    //         ofstream out(metadata_filename.c_str());
    //         unique_ptr<CObjectOStreamJson> json_out(new CObjectOStreamJson(out, eNoOwnership));
    //         json_out->SetDefaultStringEncoding(eEncoding_Ascii);
    //         json_out->PreserveKeyNames();
    //         CConstObjectInfo obj_info(m, m->GetTypeInfo());
    //         json_out->WriteObject(obj_info);
    //         json_out->Flush();
    //         out.flush();
    //         out << NcbiEndl;
    // ```
    // The members go out in the order of blastdb.asn:79-97 with their ASN.1 names
    // (`PreserveKeyNames`), two spaces of indentation per level, and a line feed
    // after the closing brace (`NcbiEndl`). The members `version` and `db-version` are
    // set by the object's defaults ("1.2", and 5 or 4 by `SetDb_version`); the
    // compressed-file members are not set.
    // NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:490-522,780-809,840-856
    // ```c++
    // void CObjectOStreamJson::BeginClassMember(const CMemberId& id)
    // {
    //     if (m_ExpectValue) {
    //         return;
    //     }
    //     ...
    //     NextElement();
    //     WriteMemberId(id);
    // }
    // ...
    // void CObjectOStreamJson::WriteKey(const string& key)
    // {
    //     string s(key);
    //     if (!m_PreserveKeys) {
    //         NStr::ReplaceInPlace(s,"-","_");
    //     }
    //     x_WriteString(s);
    //     NameSeparator();
    // }
    // ...
    // void CObjectOStreamJson::WriteKeywordValue(const string& value)
    // {
    //     BeginValue();
    //     m_Output.PutString(value);
    //     m_ExpectValue = false;
    // }
    // ...
    // void CObjectOStreamJson::BeginArray(void)
    // {
    //     BeginValue();
    //     m_Output.PutChar('[');
    //     m_Output.IncIndentLevel();
    //     m_BlockStart = true;
    //     m_ExpectValue = false;
    // }
    //
    // void CObjectOStreamJson::EndArray(void)
    // {
    //     m_Output.DecIndentLevel();
    //     m_Output.PutEol();
    //     m_Output.PutChar(']');
    //     m_BlockStart = false;
    //     m_ExpectValue = false;
    // }
    // ```
    // `WriteObject` of the object is `CObjectOStream::WriteClass` (objostr.cpp:892-904): its
    // `BeginClass` is `StartBlock`, each member is `BeginClassMember`, its value, and
    // `EndClassMember`, then `EndBlock`. The numbers are `WriteKeywordValue` of their decimal
    // text; `files` is a container of strings (`BeginContainer` is `BeginArray`, each element
    // `NextElement` and `WriteValue`). `BeginValue` writes nothing here: the object has no
    // member without a tag.
    /// `json_out->WriteObject(obj_info)`: the metadata object, through `out`'s calls.
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
            write_json_string(out, file.as_bytes())?;
        }
        out.dec_indent_level();
        out.put_eol()?;
        out.put_char(b']')?;
        out.dec_indent_level();
        out.put_eol()?;
        out.put_char(b'}')
    }

    /// The bytes of the metadata file: the object and the line feed of `NcbiEndl`.
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

/// `CDirEntry(path).GetName()`: the part after the last path separator.
fn file_name(path: &str) -> &str {
    let is_separator = |c: char| c == '/' || (cfg!(windows) && (c == '\\' || c == ':'));
    match path.rfind(is_separator) {
        Some(at) => &path[at + 1..],
        None => path,
    }
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:4039-4042
// ```c++
// bool CFile::Exists(void) const
// {
//     return IsFile();
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3204-3216
// ```c++
// Int8 CFile::GetLength(void) const
// {
//     TNcbiSys_stat st;
//     if (NcbiSys_stat(_T_XCSTRING(GetPath()), &st) != 0) {
//         LOG_ERROR_ERRNO(40, "CFile:GetLength(): stat() failed for: " + GetPath());
//         return -1L;
//     }
//     if ( GetType(st) != eFile ) {
//         LOG_ERROR_NCBI(40, "CFile:GetLength(): Not a file: " + GetPath(), CNcbiError::eOperationNotPermitted);
//         return -1L;
//     }
//     return st.st_size;
// }
// ```
// `IsFile` follows links. Between `Exists` and `GetLength` NCBI posts an error when
// the length cannot be read; a file whose `stat` succeeded a moment before gives its
// length here (`None` only when it is not a regular file).
/// `CFile(path).Exists()` and `GetLength()`: the length of the regular file `path`.
fn file_length(path: &str) -> Option<i64> {
    let meta = std::fs::metadata(path).ok()?;
    meta.is_file().then_some(meta.len() as i64)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1611-1702
// ```c++
// void CSeqDB::x_GetDBFilesMetaData(Int8 & disk_bytes, Int8 & cached_bytes, vector<string> & db_files, const string & user_path) const
// {
//     vector<string> paths;
//     vector<string> alias;
//     m_Impl->FindVolumePaths(paths, alias, true);
//     ...
//     ITERATE(vector<string>, a, alias) {
//    		CFile af(*a);
//    		if (af.Exists()) {
//    			string fn = user_path + af.GetName();
//    			db_files.push_back(fn);
//    			Int8 afl = af.GetLength();
//    			if (afl != -1) {
//    				disk_bytes += afl;
//    			...
//     vector<string> extn;
//     const bool is_protein(GetSequenceType() == CSeqDB::eProtein);
//     SeqDB_GetFileExtensions(is_protein, extn, EBlastDbVersion::eBDB_Version4);
//
//     const string kExtnMol(1, is_protein ? 'p' : 'n');
//     const string index_ext = kExtnMol + "in";
//     const string seq_ext = kExtnMol + "sq";
//
//     ITERATE(vector<string>, path, paths) {
//         ITERATE(vector<string>, ext, extn) {
//             CFile file(*path + "." + *ext);
//             if (file.Exists()) {
//             	string f = user_path + file.GetName();
//             	db_files.push_back(f);
//                 Int8 length = file.GetLength();
//                 if (length != -1) {
//                     disk_bytes += length;
//                     if((*ext == index_ext) || (*ext == seq_ext)) {
//                     	cached_bytes += length;
//                     }
//     ...
//     if (GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) {
//     	vector<string> lmdb_list;
//     	m_Impl->GetLMDBFileNames(lmdb_list);
//
//      	ITERATE(vector<string>, l, lmdb_list) {
//      		CFile file(*l);
//      		if (file.Exists()) {
//      		   string f = user_path + file.GetName();
//       	 	   db_files.push_back(f);
//       	 	   Int8 length = file.GetLength();
//        	   	   if (length != -1) {
//        	   		   disk_bytes += length;
//        	   	   ...
//        	   	   static const char * v5_exts[]={"os", "ot", "tf", "to", NULL};
//        	   	   for(const char ** p=v5_exts; *p != NULL; p++) {
//        	    		CFile v(file.GetDir() + file.GetBase() + "." + kExtnMol + (*p));
//        	    		if (v.Exists()) {
//        	    			string vf = user_path + v.GetName();
//        	    			db_files.push_back(vf);
//        	    			Int8 vl = v.GetLength();
//        	    			if (vl != -1) {
//        	    				disk_bytes += vl;
//        	    			...
//     sort(db_files.begin(), db_files.end());
// }
// ```
// `volumes` and `alias_files` are `FindVolumePaths(paths, alias, true)`: each path once.
// makeblastdb gives no `-metadata_output_prefix` (LOSAT refuses the option), so its
// `user_path` is empty, and its `lmdbs` is the one LMDB file of the version 5 database it
// wrote, `<db>.pdb`; blastdbcmd's are those of its `-metadata_output_prefix` and the LMDB
// files of the database's volumes (`GetLMDBFileNames`, one per `CSeqDBLMDBEntry`).
// `file.GetDir() + file.GetBase()` is the LMDB file's path without its extension. A file
// the earlier database of the same name left (a `<db>.pjs` that the removal of a database
// of several volumes does not take) is listed and counted as NCBI lists it.
/// `x_GetDBFilesMetaData`: the bytes of the database's files, the bytes of its index
/// and sequence files, and the sorted names of its files, each after `user_path`.
pub(crate) fn db_files_metadata(
    volumes: &[String],
    alias_files: &[String],
    lmdbs: &[String],
    is_protein: bool,
    user_path: &str,
) -> (i64, i64, Vec<String>) {
    let mut db_files = Vec::new();
    let mut disk_bytes = 0i64;
    let mut cached_bytes = 0i64;
    let listed = |name: &str| format!("{user_path}{}", file_name(name));
    for alias in alias_files {
        if let Some(length) = file_length(alias) {
            db_files.push(listed(alias));
            disk_bytes += length;
        }
    }
    let mol = if is_protein { 'p' } else { 'n' };
    for path in volumes {
        for ext in VOLUME_EXTENSIONS {
            let file = format!("{path}.{mol}{ext}");
            if let Some(length) = file_length(&file) {
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
            db_files.push(listed(lmdb));
            disk_bytes += length;
            let base = lmdb.rfind('.').map_or(lmdb.as_str(), |dot| &lmdb[..dot]);
            for ext in ["os", "ot", "tf", "to"] {
                let file = format!("{base}.{mol}{ext}");
                if let Some(length) = file_length(&file) {
                    db_files.push(listed(&file));
                    disk_bytes += length;
                }
            }
        }
    }
    db_files.sort();
    (disk_bytes, cached_bytes, db_files)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1730-1734
// ```c++
// 	CTimeFormat timeFmt = CTimeFormat::GetPredefined(CTimeFormat::eISO8601_DateTimeSec);
// 	CTime date(GetDate(), kBlastDbDateFormat);
// 	m->SetLast_updated(date.AsString(timeFmt));
// ```
// `kBlastDbDateFormat` is "b d, Y  H:m P" (seqdbcommon.hpp), the form of the date the
// writer stores, and the ISO 8601 form is "Y-M-DTh:m:s". The stored date has no
// seconds, so they are 00. The volumes of one build share their date.
/// `last-updated`: the date of the index files in ISO 8601 form.
pub(crate) fn last_updated(now: &super::time::LocalTime) -> String {
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:00",
        now.year, now.month, now.day, now.hour, now.minute
    )
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

impl IoException {
    /// `GetType() + "::" + GetErrCodeString()`.
    fn type_and_code(&self) -> &'static str {
        match self.code {
            IoErrCode::Write => "CIOException::eWrite",
            IoErrCode::Flush => "CIOException::eFlush",
        }
    }
}

// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:46,783-795
// ```c++
// static const size_t KInitialBufferSize = 4096;
// ...
// COStreamBuffer::COStreamBuffer(CNcbiOstream& out, bool deleteOut)
//     THROWS1((bad_alloc))
//     : m_Output(out), m_DeleteOutput(deleteOut), m_Closed(false), m_Error(0),
//       m_IndentLevel(0), m_BufferPos(0),
//       m_Buffer(new char[KInitialBufferSize]),
//       m_CurrentPos(m_Buffer),
//       m_BufferEnd(m_Buffer + KInitialBufferSize),
//       m_Line(1), m_LineLength(0),
//       m_BackLimit(0), m_UseIndentation(true), m_UseEol(true),
//       m_CanceledCallback(0)
// ```
// The JSON writer sets no back limit and no cancel callback, and makeblastdb gives it the
// stream without ownership (`eNoOwnership`). The line count and line length serve the
// line positions and the wrapping of other formats; they write nothing here.
/// `COStreamBuffer` over the metadata file's `ofstream`.
pub(crate) struct OStreamBuffer<'a> {
    /// `m_Output`.
    output: &'a mut OfStream,
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
    /// `m_CurrentPos - m_Buffer`.
    current: usize,
}

impl<'a> OStreamBuffer<'a> {
    pub(crate) fn new(output: &'a mut OfStream) -> Self {
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
    pub(crate) fn output(&mut self) -> &mut OfStream {
        self.output
    }

    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:185-229
    // ```c++
    // bool COStreamBuffer::fail(void) const
    // {
    //     return m_Error != 0;
    // }
    // ...
    // Int8 COStreamBuffer::GetStreamPosAsInt8(void) const
    // {
    //     return m_BufferPos + (m_CurrentPos - m_Buffer);
    // }
    // ```
    // `fail()` is `error.is_some()`.
    /// `GetStreamPos()`.
    fn stream_pos(&self) -> u64 {
        self.buffer_pos + self.current as u64
    }

    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:283-312
    // ```c++
    // char* COStreamBuffer::Skip(size_t count)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     char* pos = m_CurrentPos;
    //     char* end = pos + count;
    //     if ( end <= m_BufferEnd ) {
    //         // enough space in buffer
    //         m_CurrentPos = end;
    //         m_LineLength += count;
    //         return pos;
    //     }
    //     else {
    //         return DoSkip(count);
    //     }
    // }
    //
    // char* COStreamBuffer::Reserve(size_t count)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     char* pos = m_CurrentPos;
    //     char* end = pos + count;
    //     if ( end <= m_BufferEnd ) {
    //         // enough space in buffer
    //         return pos;
    //     }
    //     else {
    //         return DoReserve(count);
    //     }
    // }
    // ```
    // `DoSkip` (strbuffer.inl:273-281) is `DoReserve` and the same advance.
    /// `Skip(count)`: the position of the `count` bytes taken.
    fn skip(&mut self, count: usize) -> Result<usize, IoException> {
        let pos = self.reserve(count)?;
        self.current = pos + count;
        Ok(pos)
    }

    /// `Reserve(count)`: the position where `count` bytes fit.
    fn reserve(&mut self, count: usize) -> Result<usize, IoException> {
        if self.current + count <= self.buffer.len() {
            Ok(self.current)
        } else {
            self.do_reserve(count)
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:921-948
    // ```c++
    // char* COStreamBuffer::DoReserve(size_t count)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     FlushBuffer(false);
    //     size_t usedSize = m_CurrentPos - m_Buffer;
    //     size_t needSize = usedSize + count;
    //     size_t bufferSize = m_BufferEnd - m_Buffer;
    //     if ( bufferSize < needSize ) {
    //         // realloc too small buffer
    //         do {
    //             bufferSize = BiggerBufferSize(bufferSize);
    //         } while ( bufferSize < needSize );
    //         ...
    //     }
    //     return m_CurrentPos;
    // }
    // ```
    // `BiggerBufferSize` doubles (strbuffer.cpp:48-52).
    fn do_reserve(&mut self, count: usize) -> Result<usize, IoException> {
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
    // ```c++
    // void COStreamBuffer::FlushBuffer(bool fullBuffer)
    //     THROWS1((CIOException))
    // {
    //     ...
    //     size_t used = GetUsedSpace();
    //     size_t count;
    //     size_t leave;
    //     if ( fullBuffer ) {
    //         count = used;
    //         leave = 0;
    //     }
    //     else {
    //         leave = m_BackLimit;
    //         if ( used < leave )
    //             return; // none to flush
    //         count = used - leave;
    //     }
    //     if ( count != 0 ) {
    //         if ( !m_Output.write(m_Buffer, count) ) {
    //             m_Error = "write fault";
    //             NCBI_THROW(CIOException,eWrite,m_Error);
    //         }
    //         if ( leave != 0 ) {
    //             memmove(m_Buffer, m_Buffer + count, leave);
    //             m_CurrentPos -= count;
    //         }
    //         else {
    //             m_CurrentPos = m_Buffer;
    //         }
    //         m_BufferPos += CT_OFF_TYPE(count);
    //     }
    // }
    // ```
    // With no back limit, both forms write all the buffer holds. A failed write keeps it there.
    /// `FlushBuffer(fullBuffer)`.
    fn flush_buffer(&mut self, _full_buffer: bool) -> Result<(), IoException> {
        let count = self.current;
        if count != 0 {
            self.output.write(&self.buffer[..count]);
            if self.output.fail() {
                self.error = Some("write fault");
                return Err(IoException {
                    code: IoErrCode::Write,
                    msg: "write fault",
                });
            }
            self.current = 0;
            self.buffer_pos += count as u64;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:813-837,910-918
    // ```c++
    // struct STemporarilyClearStreamState
    // {
    //     STemporarilyClearStreamState(ios& stream)
    //         : m_Stream(stream),
    //           m_State(stream.rdstate())
    //         {
    //             m_Stream.clear();
    //         }
    //     ~STemporarilyClearStreamState()
    //         {
    //             try {
    //                 m_Stream.setstate(m_State);
    //             }
    //             ...
    //         }
    // ...
    // void COStreamBuffer::Flush(void)
    //     THROWS1((CIOException))
    // {
    //     STemporarilyClearStreamState state(m_Output);
    //     FlushBuffer();
    //     if ( !m_Output.flush() ) {
    //         NCBI_THROW(CIOException,eFlush,"COStreamBuffer::Flush: failed");
    //     }
    // }
    // ```
    // `eFlush` sets no `m_Error`.
    /// `Flush()`.
    pub(crate) fn flush(&mut self) -> Result<(), IoException> {
        let state = self.output.rdstate();
        self.output.clear(super::ofstream::IoState::GOOD);
        let result = self.flush_buffer(true).and_then(|()| {
            self.output.flush();
            if self.output.fail() {
                Err(IoException {
                    code: IoErrCode::Flush,
                    msg: "COStreamBuffer::Flush: failed",
                })
            } else {
                Ok(())
            }
        });
        self.output.setstate(state);
        result
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:843-862
    // ```c++
    // void COStreamBuffer::Close(void)
    // {
    //     if ( !m_Closed && m_Output ) {
    //         m_Closed = true;
    //         if ( m_DeleteOutput ) {
    //             ...
    //         }
    //         else {
    //             STemporarilyClearStreamState state(m_Output);
    //             FlushBuffer();
    //         }
    //     }
    //     m_Error = 0;
    //     m_IndentLevel = 0;
    //     m_CurrentPos = m_Buffer;
    //     m_Line = 1;
    //     m_LineLength = 0;
    // }
    // ```
    // `m_Output` as a condition is `!fail()` of the stream.
    /// `Close()`.
    fn close(&mut self) -> Result<(), IoException> {
        if !self.closed && !self.output.fail() {
            self.closed = true;
            let state = self.output.rdstate();
            self.output.clear(super::ofstream::IoState::GOOD);
            let result = self.flush_buffer(true);
            self.output.setstate(state);
            result?;
        }
        self.error = None;
        self.indent_level = 0;
        self.current = 0;
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:251-263,314-381
// ```c++
// void COStreamBuffer::IncIndentLevel(size_t step)
//     THROWS1_NONE
// {
//     m_IndentLevel += step;
// }
// ...
// void COStreamBuffer::PutChar(char c)
//     THROWS1((CIOException))
// {
//     *Skip(1) = c;
// }
// ...
// void COStreamBuffer::PutString(const char* str, size_t length)
//     THROWS1((CIOException, bad_alloc))
// {
//     if ( length < 1024 ) {
//         memcpy(Skip(length), str, length);
//     }
//     else {
//         Write(str, length);
//     }
// }
// ...
// void COStreamBuffer::PutIndent(void)
//     THROWS1((CIOException, bad_alloc))
// {
//     if (GetUseIndentation()) {
//         size_t count = m_IndentLevel;
//         memset(Skip(count), ' ', count);
//     }
// }
//
// void COStreamBuffer::PutEol(bool indent)
//     THROWS1((CIOException, bad_alloc))
// {
//     if (!GetUseEol()) {
//         return;
//     }
//     char* pos = Reserve(1);
//     *pos = '\n';
//     m_CurrentPos = pos + 1;
//     ++m_Line;
//     m_LineLength = 0;
//     if ( indent )
//         PutIndent();
// }
// ```
// `IncIndentLevel` and `DecIndentLevel` step by 2 (strbuffer.hpp:244-245). The strings of the
// metadata object are written byte by byte (`x_WriteString`); `PutString` gets the numbers
// and the escapes, all shorter than 1024 bytes.
impl JsonOut for OStreamBuffer<'_> {
    type Error = IoException;
    fn put_char(&mut self, c: u8) -> Result<(), IoException> {
        let pos = self.skip(1)?;
        self.buffer[pos] = c;
        Ok(())
    }
    fn put_string(&mut self, s: &[u8]) -> Result<(), IoException> {
        let pos = self.skip(s.len())?;
        self.buffer[pos..pos + s.len()].copy_from_slice(s);
        Ok(())
    }
    fn put_eol(&mut self) -> Result<(), IoException> {
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
// ```c++
//     enum EFailFlags {
//         /// No error
//         fNoError       = 0,             eNoError     = fNoError,
//         /// An unknown error when writing into output file
//         fWriteError    = 1 << 1,        eWriteError  = fWriteError,
//         ...
//         /// No output file
//         fNotOpen       = 1 << 7,        eNotOpen     = fNotOpen,
// ```
const F_WRITE_ERROR: u32 = 1 << 1;
const F_NOT_OPEN: u32 = 1 << 7;

/// The parts of `CObjectOStreamJson` (`CObjectOStream`) that act on its stream after the
/// object is written: the fail flags and the close of its destructor.
struct ObjectOStream<'a> {
    /// `m_Output`.
    output: OStreamBuffer<'a>,
    /// `m_Fail`.
    fail: u32,
}

/// The `CSerialException` that `ThrowError1` throws over the exception it was given.
struct SerialException {
    /// `byte N: <message>`.
    msg: String,
    predecessor: IoException,
}

impl ObjectOStream<'_> {
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384,495-513
    // ```c++
    // CObjectOStream::TFailFlags CObjectOStream::SetFailFlags(TFailFlags flags,
    //                                                         const char* message)
    // {
    //     TFailFlags old = m_Fail;
    //     m_Fail |= flags;
    //     if ( !old && flags ) {
    //         // first fail
    //         ERR_POST_X(5, "CObjectOStream: error at "<<
    //                    GetPosition()<<": "<<GetStackTrace() << ": " << message);
    //     }
    //     return old;
    // }
    // ...
    // string CObjectOStream::GetStackTrace(void) const
    // {
    //     return GetStackTraceASN();
    // }
    // ...
    // string CObjectOStream::GetPosition(void) const
    // {
    //     return "byte "+NStr::Int8ToString(NcbiStreamposToInt8(GetStreamPos()));
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-84
    // ```c++
    // string CObjectStack::GetStackTraceASN(void) const
    // {
    //     if (!GetStackDepth()) {
    //         return "stack is empty";
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/error_codes.hpp:46
    // ```c++
    // NCBI_DEFINE_ERRCODE_X(Serial_OStream,    802, 13);
    // ```
    // Every call here comes from the close in the destructor: the frames of the object are
    // gone (an exception pops them on its way out), and in a base class's destructor the
    // virtual `GetPosition` is `CObjectOStream`'s, the byte position, not the line position
    // of `CObjectOStreamJson`. makeblastdb's diagnostics are "Severity: (code.subcode)
    // [makeblastdb] text".
    /// `SetFailFlags(flags, message)`.
    fn set_fail_flags(&mut self, flags: u32, message: &str, err: &mut dyn Write) {
        let old = self.fail;
        self.fail |= flags;
        if old == 0 && flags != 0 {
            let _ = writeln!(
                err,
                "Error: (802.5) [makeblastdb] CObjectOStream: error at byte {}: stack is empty: \
                 {message}",
                self.output.stream_pos()
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:386-400
    // ```c++
    // bool CObjectOStream::InGoodState(void)
    // {
    //     if ( fail() ) {
    //         // fail flag already set
    //         return false;
    //     }
    //     else if ( m_Output.fail() ) {
    //         // IO exception thrown without setting fail flag
    //         SetFailFlags(fWriteError, m_Output.GetError());
    //         m_Output.ResetFail();
    //         return false;
    //     }
    //     else {
    //         // ok
    //         return true;
    //     }
    // }
    // ```
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

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:522-555
    // ```c++
    // void CObjectOStream::ThrowError1(const CDiagCompileInfo& diag_info,
    //                                  TFailFlags flags, const string& message,
    //                                  CException* exc)
    // {
    //     CSerialException::EErrCode err;
    //     try {
    //         DefaultFlush();
    //     } catch(...) {
    //     }
    //     string msg(message);
    //     if (flags == fUnassigned) {
    //         msg = "cannot write unassigned member "+message;
    //     }
    //     SetFailFlags(flags, msg.c_str());
    //     msg.insert(0,GetPosition()+": ");
    //     switch (flags)
    //     {
    //     ...
    //     case fWriteError:     err = CSerialException::eIoError;        break;
    //     ...
    //     }
    //     throw CSerialException(diag_info,exc,err,msg);
    // }
    // ```
    // `RethrowError(flag, mess, exc)` is `ThrowError1(DIAG_COMPILE_INFO, flag, mess, &exc)`
    // (objostr.hpp:578-579). `DefaultFlush` is `Flush` (no `fFlagNoAutoFlush`).
    /// `RethrowError(fWriteError, message, exc)`.
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
    // ```c++
    // void CObjectOStream::DefaultFlush(void)
    // {
    //     if ( GetFlags() & fFlagNoAutoFlush ) {
    //         FlushBuffer();
    //     }
    //     else {
    //         Flush();
    //     }
    // }
    //
    //
    // void CObjectOStream::Close(void)
    // {
    //     if (m_Fail != fNotOpen) {
    //         try {
    //             DefaultFlush();
    //             if ( m_Objects )
    //                 m_Objects->Clear();
    //             ClearStack();
    //             m_Fail = fNotOpen;
    //             ResetState();
    //             m_Output.Close();
    //         }
    //         catch (CException& exc) {
    //             if ( InGoodState() )
    //                 RethrowError(fWriteError, "cannot close output stream",exc);
    //         }
    //     }
    // }
    // ```
    /// `Close()`.
    fn close(&mut self, err: &mut dyn Write) -> Result<(), SerialException> {
        if self.fail != F_NOT_OPEN {
            let result = self.output.flush().and_then(|()| {
                self.fail = F_NOT_OPEN;
                self.output.close()
            });
            if let Err(exc) = result {
                if self.in_good_state(err) {
                    return Err(self.rethrow_write_error("cannot close output stream", exc, err));
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:317-324
    // ```c++
    // CObjectOStream::~CObjectOStream(void)
    // {
    //     try {
    //         Close();
    //         ResetLocalHooks();
    //     }
    //     NCBI_CATCH_X(4, "Cannot close serializing output stream")
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:676-700
    // ```c++
    // void CExceptionReporter::ReportDefaultEx(int err_code, int err_subcode,
    //     const CDiagCompileInfo& info, const string& title,const exception& ex,
    //     TDiagPostFlags flags)
    // {
    //     ...
    //     } else {
    //         CNcbiDiag d(info, cex->GetSeverity(), flags);
    //         d.SetOmitStackTrace(true);
    //         d << ErrCode(err_code, err_subcode) << title << " " <<  *cex;
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7873-7893
    // ```c++
    //     const CException* pex;
    //     const CException* main_pex = NULL;
    //     stack<const CException*> pile;
    //     // invert the order
    //     for (pex = &ex; pex; pex = pex->GetPredecessor()) {
    //         pile.push(pex);
    //         if (!main_pex  &&  pex->HasMainText())
    //             main_pex = pex;
    //     }
    //     ...
    //     if ( !IsOssEmpty(*m_Buffer.m_Stream) ) {
    //         *this << "(" << main_pex->GetType() << "::"
    //                      << main_pex->GetErrCodeString() << ") "
    //               << s_GetExceptionText(main_pex);
    //     }
    // ```
    // `CNcbiDiag::x_Put` (ncbidiag.cpp:7894-7939) posts each exception of the chain, the oldest
    // first, as a message of its own whose code is its type and error code, before the message
    // of the report itself ("title" and the newest exception). The predecessor of the
    // `CSerialException` is the `CIOException` of the failed flush, which has no predecessor
    // of its own; neither has extra text (`ReportExtra`).
    /// `~CObjectOStream`.
    fn destroy(&mut self, err: &mut dyn Write) {
        if let Err(exc) = self.close(err) {
            let _ = writeln!(
                err,
                "Error: ({}) [makeblastdb] {}",
                exc.predecessor.type_and_code(),
                exc.predecessor.msg
            );
            let _ = writeln!(
                err,
                "Error: (CSerialException::eIoError) [makeblastdb] {}",
                exc.msg
            );
            let _ = writeln!(
                err,
                "Error: (802.4) [makeblastdb] Cannot close serializing output stream \
                 (CSerialException::eIoError) {}",
                exc.msg
            );
        }
        self.destroy_buffer(err);
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:797-811
    // ```c++
    // COStreamBuffer::~COStreamBuffer(void)
    // {
    //     try {
    //         Close();
    //     }
    //     NCBI_CATCH_X(2, "~COStreamBuffer: exception while closing");
    // ```
    // NCBI reference (598d8ae6): c++/include/util/error_codes.hpp:48
    // ```c++
    // NCBI_DEFINE_ERRCODE_X(Util_Stream,      205,   2);
    // ```
    // After `CObjectOStream::Close` the buffer is closed, or the stream failed: this `Close`
    // writes nothing in makeblastdb's sequence; its report is that of `destroy` for an
    // exception without predecessor.
    /// `~COStreamBuffer` (a member of the object stream, destroyed after its body).
    fn destroy_buffer(&mut self, err: &mut dyn Write) {
        if let Err(exc) = self.output.close() {
            let _ = writeln!(
                err,
                "Error: ({}) [makeblastdb] {}",
                exc.type_and_code(),
                exc.msg
            );
            let _ = writeln!(
                err,
                "Error: (205.2) [makeblastdb] ~COStreamBuffer: exception while closing ({}) {}",
                exc.type_and_code(),
                exc.msg
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1306-1318
// ```c++
//     	CRef<CBlast_db_metadata> m = sdb.GetDBMetaData(output_prefix);
//     	string extn (kEmptyStr);
//     	SeqDB_GetMetadataFileExtension(is_protein, extn);
//     	string metadata_filename = new_db + "." + extn;
//         ofstream out(metadata_filename.c_str());
//         unique_ptr<CObjectOStreamJson> json_out(new CObjectOStreamJson(out, eNoOwnership));
//         json_out->SetDefaultStringEncoding(eEncoding_Ascii);
//         json_out->PreserveKeyNames();
//         CConstObjectInfo obj_info(m, m->GetTypeInfo());
//         json_out->WriteObject(obj_info);
//         json_out->Flush();
//         out.flush();
//         out << NcbiEndl;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
// ```c++
//     catch (const CIOException& e) {                                         \
//         if (e.GetErrCode() == CIOException::eFlush) {                       \
//             LOG_POST(Error << "BLAST failed to write output: " << e.GetMsg());\
//             exit_code = BLAST_OUTPUT_ERROR;                                 \
//         }                                                                   \
//     }                                                                       \
// ```
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
    new_db: &str,
    is_protein: bool,
    metadata: &DbMetadata,
    err: &mut dyn Write,
) -> Option<&'static str> {
    let name = format!("{new_db}.{}", if is_protein { "pjs" } else { "njs" });
    let mut out = OfStream::create(Path::new(&name));
    let escaped = {
        let mut json_out = ObjectOStream {
            output: OStreamBuffer::new(&mut out),
            fail: 0,
        };
        let result = metadata
            .write_json(&mut json_out.output)
            .and_then(|()| json_out.output.flush())
            .map(|()| {
                let out = json_out.output.output();
                out.flush();
                out.endl();
            });
        json_out.destroy(err);
        result.err()
    };
    drop(out);
    escaped
        .filter(|exception| exception.code == IoErrCode::Flush)
        .map(|exception| exception.msg)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The file makeblastdb 2.17.0 wrote for the three proteins of survey A's
    // f1_tiny.fa (`-title t -out db`, version 5).
    #[test]
    fn the_json_layout_is_ncbi_s() {
        let metadata = DbMetadata {
            dbname: b"db".to_vec(),
            dbtype: "Protein",
            db_version: 5,
            description: b"t".to_vec(),
            number_of_letters: 120,
            number_of_sequences: 3,
            last_updated: "2026-10-29T20:32:00".to_string(),
            number_of_volumes: 1,
            number_of_taxids: None,
            bytes_total: 37723,
            bytes_to_cache: 228,
            files: [
                "db.pdb", "db.phr", "db.pin", "db.pot", "db.psq", "db.ptf", "db.pto",
            ]
            .map(String::from)
            .to_vec(),
        };
        let expected =
            "{\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \
             \"db-version\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 120,\n  \
             \"number-of-sequences\": 3,\n  \"last-updated\": \"2026-10-29T20:32:00\",\n  \
             \"number-of-volumes\": 1,\n  \"bytes-total\": 37723,\n  \"bytes-to-cache\": 228,\n  \
             \"files\": [\n    \"db.pdb\",\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pot\",\n    \
             \"db.psq\",\n    \"db.ptf\",\n    \"db.pto\"\n  ]\n}\n";
        assert_eq!(String::from_utf8(metadata.to_json()).unwrap(), expected);
    }

    // makeblastdb 2.17.0 with `-title` = `q"b\s<TAB>t<0x01>c é end ` (é in UTF-8).
    #[test]
    fn strings_are_escaped_as_ncbi_s_json_writer_does() {
        let escaped = |value: &[u8]| {
            let mut out = JsonText::default();
            let Ok(()) = write_json_string(&mut out, value);
            out.bytes
        };
        assert_eq!(
            escaped("q\"b\\s\tt\u{1}c é end ".as_bytes()),
            b"\"q\\\"b\\\\s\\u0009t\\u0001c \xc3\xa9 end \"".to_vec()
        );
        assert_eq!(
            escaped(b"a\x1bb\x7f\xff\x00c"),
            b"\"a\\u001Bb\x7f\xff\"".to_vec()
        );
    }

    #[test]
    fn the_date_is_iso_8601_without_seconds() {
        let at = |hour, minute| super::super::time::LocalTime {
            year: 2026,
            month: 3,
            day: 9,
            hour,
            minute,
            second: 41,
        };
        assert_eq!(last_updated(&at(0, 5)), "2026-03-09T00:05:00");
        assert_eq!(last_updated(&at(23, 59)), "2026-03-09T23:59:00");
    }

    /// A directory of its own under the system's temporary directory.
    struct TempDir(std::path::PathBuf);

    impl TempDir {
        fn new(name: &str) -> Self {
            let dir = std::env::temp_dir().join(format!(
                "losat-metadata-{name}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            Self(dir)
        }

        fn db(&self) -> String {
            self.0.join("db").to_string_lossy().into_owned()
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// The metadata of a protein database of `volumes` volumes of four files each, with a
    /// description of `pad` bytes.
    fn metadata(volumes: usize, pad: usize) -> DbMetadata {
        let mut files: Vec<String> = (0..volumes)
            .flat_map(|v| ["phr", "pin", "pog", "psq"].map(|e| format!("db.{v:03}.{e}")))
            .collect();
        files
            .extend(["db.pal", "db.pdb", "db.pos", "db.pot", "db.ptf", "db.pto"].map(String::from));
        DbMetadata {
            dbname: b"db".to_vec(),
            dbtype: "Protein",
            db_version: 5,
            description: vec![b'x'; pad],
            number_of_letters: 7800,
            number_of_sequences: 130,
            last_updated: "2026-11-12T12:56:00".to_string(),
            number_of_volumes: volumes as i32,
            number_of_taxids: None,
            bytes_total: 412_004,
            bytes_to_cache: 61_380,
            files,
        }
    }

    /// The metadata file's text (when it is a regular file: `/dev/full` reads without end), its
    /// messages and the exception that ends makeblastdb.
    fn written(
        path: &str,
        metadata: &DbMetadata,
    ) -> (Option<Vec<u8>>, String, Option<&'static str>) {
        let mut err = Vec::new();
        let escaped = write_metadata_file(path, true, metadata, &mut err);
        let pjs = format!("{path}.pjs");
        let text = std::fs::metadata(&pjs)
            .is_ok_and(|meta| meta.is_file())
            .then(|| std::fs::read(&pjs).ok())
            .flatten();
        (text, String::from_utf8(err).unwrap(), escaped)
    }

    // The object is written through the serializer's buffer: what reaches the file is the text
    // of `to_json`, whatever its length.
    #[test]
    fn the_metadata_file_holds_the_object_and_a_line_feed() {
        for volumes in [1, 130] {
            let tmp = TempDir::new(&format!("ok{volumes}"));
            let metadata = metadata(volumes, 3);
            let (text, err, escaped) = written(&tmp.db(), &metadata);
            assert_eq!(text, Some(metadata.to_json()), "{volumes}");
            assert_eq!((err.as_str(), escaped), ("", None));
        }
    }

    // A file that cannot be opened (a directory of its name): the serializer's `Close` finds the
    // write fault at the end of what it holds. The text up to 4096 bytes is written at
    // `Flush`; a longer one fails at its first write, when the 4097th byte does not fit.
    // makeblastdb 2.17.0: `byte 479` (b_p5_pjs) and `byte 4096` (b_m_pjs), exit status 0.
    #[test]
    fn a_metadata_file_that_cannot_be_opened_reports_a_write_fault() {
        for (volumes, at) in [(1, None), (130, Some(4096))] {
            let tmp = TempDir::new(&format!("blocked{volumes}"));
            std::fs::create_dir(format!("{}.pjs", tmp.db())).unwrap();
            let metadata = metadata(volumes, 3);
            let at = at.unwrap_or(metadata.to_json().len() - 1);
            let (_, err, escaped) = written(&tmp.db(), &metadata);
            assert_eq!(
                err,
                format!(
                    "Error: (802.5) [makeblastdb] CObjectOStream: error at byte {at}: stack is \
                     empty: write fault\n"
                )
            );
            assert_eq!(escaped, None);
        }
    }

    // The buffer is written when a call does not fit: a line feed at byte 4092 fits, its
    // indentation of 4 spaces does not, so the first write holds 4093 bytes.
    #[test]
    fn the_buffer_is_written_where_a_call_does_not_fit() {
        let plain = metadata(130, 0).to_json();
        let eol = plain[..4097]
            .windows(5)
            .rposition(|w| w == b"\n    ")
            .unwrap();
        let metadata = metadata(130, 4092 - eol);
        assert_eq!(&metadata.to_json()[4092..4097], b"\n    ");
        let tmp = TempDir::new("straddle");
        std::fs::create_dir(format!("{}.pjs", tmp.db())).unwrap();
        let (_, err, _) = written(&tmp.db(), &metadata);
        assert_eq!(
            err,
            "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 4093: stack is empty: \
             write fault\n"
        );
    }

    // A device that takes no byte: a text of up to 4096 bytes sits in the stream's buffer of
    // 8191 bytes until `out.flush()` in `COStreamBuffer::Flush` fails, `eFlush` leaves
    // makeblastdb's block, and the serializer's `Close` fails again and reports (f_p5_pjs,
    // exit status 6). A longer text gets its first 4096 bytes into the stream's buffer; the
    // next 4096 do not fit the 4095 left and their write fails (f_m_pjs: `byte 8192`, exit
    // status 0).
    #[cfg(target_os = "linux")]
    #[test]
    fn a_metadata_file_on_a_full_device_cannot_be_written() {
        let tmp = TempDir::new("full");
        std::os::unix::fs::symlink("/dev/full", format!("{}.pjs", tmp.db())).unwrap();
        let small = metadata(1, 3);
        let at = small.to_json().len() - 1;
        let (_, err, escaped) = written(&tmp.db(), &small);
        assert_eq!(
            err,
            format!(
                "Error: (802.5) [makeblastdb] CObjectOStream: error at byte {at}: stack is empty: \
                 cannot close output stream\n\
                 Error: (CIOException::eFlush) [makeblastdb] COStreamBuffer::Flush: failed\n\
                 Error: (CSerialException::eIoError) [makeblastdb] byte {at}: cannot close output \
                 stream\n\
                 Error: (802.4) [makeblastdb] Cannot close serializing output stream \
                 (CSerialException::eIoError) byte {at}: cannot close output stream\n"
            )
        );
        assert_eq!(escaped, Some("COStreamBuffer::Flush: failed"));
        let (_, err, escaped) = written(&tmp.db(), &metadata(130, 3));
        assert_eq!(
            err,
            "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 8192: stack is empty: \
             write fault\n"
        );
        assert_eq!(escaped, None);
    }
}

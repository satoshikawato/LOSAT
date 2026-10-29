//! The metadata file of a database (`.pjs` / `.njs`) that makeblastdb writes after a
//! successful build: `CSeqDB::GetDBMetaData` (`c++/src/objtools/blast/seqdb_reader/seqdb.cpp`)
//! serialized as JSON by `CObjectOStreamJson` (`c++/src/serial/objostrjson.cpp`).

use std::io::Write;
use std::path::PathBuf;

use super::build_db::VOLUME_EXTENSIONS;
use super::WriteDbError;

/// `CBlast_db_metadata` (`c++/src/objects/blastdb/blastdb.asn:79-97`) as makeblastdb fills it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DbMetadata {
    pub dbname: Vec<u8>,
    pub dbtype: &'static str,
    pub db_version: u32,
    pub description: Vec<u8>,
    pub number_of_letters: u64,
    pub number_of_sequences: i32,
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
/// `x_WriteString`: a JSON string of `value`.
fn write_json_string(out: &mut Vec<u8>, value: &[u8]) {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    out.push(b'"');
    for &c in value.iter().take_while(|&&c| c != 0) {
        match c {
            b'"' => out.extend_from_slice(b"\\\""),
            b'\\' => out.extend_from_slice(b"\\\\"),
            c if c < 0x20 => {
                out.extend_from_slice(b"\\u00");
                out.push(HEX[(c >> 4) as usize]);
                out.push(HEX[(c & 0xF) as usize]);
            }
            c => out.push(c),
        }
    }
    out.push(b'"');
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1318-1328
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
    /// The bytes of the metadata file.
    pub(crate) fn to_json(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(512 + 32 * self.files.len());
        let mut first = true;
        let mut key = |out: &mut Vec<u8>, name: &str| {
            if !first {
                out.push(b',');
            }
            first = false;
            out.extend_from_slice(b"\n  ");
            write_json_string(out, name.as_bytes());
            out.extend_from_slice(b": ");
        };
        out.push(b'{');
        key(&mut out, "version");
        write_json_string(&mut out, b"1.2");
        key(&mut out, "dbname");
        write_json_string(&mut out, &self.dbname);
        key(&mut out, "dbtype");
        write_json_string(&mut out, self.dbtype.as_bytes());
        key(&mut out, "db-version");
        out.extend_from_slice(self.db_version.to_string().as_bytes());
        key(&mut out, "description");
        write_json_string(&mut out, &self.description);
        key(&mut out, "number-of-letters");
        out.extend_from_slice(self.number_of_letters.to_string().as_bytes());
        key(&mut out, "number-of-sequences");
        out.extend_from_slice(self.number_of_sequences.to_string().as_bytes());
        key(&mut out, "last-updated");
        write_json_string(&mut out, self.last_updated.as_bytes());
        key(&mut out, "number-of-volumes");
        out.extend_from_slice(self.number_of_volumes.to_string().as_bytes());
        if let Some(taxids) = self.number_of_taxids {
            key(&mut out, "number-of-taxids");
            out.extend_from_slice(taxids.to_string().as_bytes());
        }
        key(&mut out, "bytes-total");
        out.extend_from_slice(self.bytes_total.to_string().as_bytes());
        key(&mut out, "bytes-to-cache");
        out.extend_from_slice(self.bytes_to_cache.to_string().as_bytes());
        key(&mut out, "files");
        out.push(b'[');
        for (i, file) in self.files.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            out.extend_from_slice(b"\n    ");
            write_json_string(&mut out, file.as_bytes());
        }
        out.extend_from_slice(b"\n  ]");
        out.extend_from_slice(b"\n}\n");
        out
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
    meta.is_file().then(|| meta.len() as i64)
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
// makeblastdb gives no `-metadata_output_prefix` here (LOSAT refuses the option), so
// `user_path` is empty. `lmdb` is the one LMDB file of a version 5 database that
// makeblastdb wrote (`GetLMDBFileNames` of a database whose volumes name one LMDB
// file), `<db>.pdb`; `file.GetDir() + file.GetBase()` is then `<db>`. A file the
// earlier database of the same name left (a `<db>.pjs` that the removal of a database
// of several volumes does not take) is listed and counted as NCBI lists it.
/// `x_GetDBFilesMetaData`: the bytes of the database's files, the bytes of its index
/// and sequence files, and the sorted names of its files.
pub(crate) fn db_files_metadata(
    volumes: &[String],
    alias_files: &[String],
    lmdb: Option<&str>,
    is_protein: bool,
) -> (i64, i64, Vec<String>) {
    let mut db_files = Vec::new();
    let mut disk_bytes = 0i64;
    let mut cached_bytes = 0i64;
    for alias in alias_files {
        if let Some(length) = file_length(alias) {
            db_files.push(file_name(alias).to_string());
            disk_bytes += length;
        }
    }
    let mol = if is_protein { 'p' } else { 'n' };
    for path in volumes {
        for ext in VOLUME_EXTENSIONS {
            let file = format!("{path}.{mol}{ext}");
            if let Some(length) = file_length(&file) {
                db_files.push(file_name(&file).to_string());
                disk_bytes += length;
                if ext == "in" || ext == "sq" {
                    cached_bytes += length;
                }
            }
        }
    }
    if let Some(lmdb) = lmdb {
        if let Some(length) = file_length(lmdb) {
            db_files.push(file_name(lmdb).to_string());
            disk_bytes += length;
            let base = lmdb.rfind('.').map_or(lmdb, |dot| &lmdb[..dot]);
            for ext in ["os", "ot", "tf", "to"] {
                let file = format!("{base}.{mol}{ext}");
                if let Some(length) = file_length(&file) {
                    db_files.push(file_name(&file).to_string());
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

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1314-1318
// ```c++
//     	CRef<CBlast_db_metadata> m = sdb.GetDBMetaData(output_prefix);
//     	string extn (kEmptyStr);
//     	SeqDB_GetMetadataFileExtension(is_protein, extn);
//     	string metadata_filename = new_db + "." + extn;
//         ofstream out(metadata_filename.c_str());
// ```
// `SeqDB_GetMetadataFileExtension` gives "pjs" or "njs". `ofstream` replaces a file
// of that name; NCBI does not test the stream.
/// Writes the metadata file `<db>.pjs` / `<db>.njs`.
pub(crate) fn write_metadata_file(
    new_db: &str,
    is_protein: bool,
    metadata: &DbMetadata,
) -> Result<(), WriteDbError> {
    let name = format!("{new_db}.{}", if is_protein { "pjs" } else { "njs" });
    let write = || -> std::io::Result<()> {
        let mut file = std::fs::File::create(&name)?;
        file.write_all(&metadata.to_json())
    };
    write().map_err(|source| WriteDbError::Io {
        path: PathBuf::from(&name),
        source,
    })
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
        let mut out = Vec::new();
        write_json_string(&mut out, "q\"b\\s\tt\u{1}c é end ".as_bytes());
        assert_eq!(
            out,
            b"\"q\\\"b\\\\s\\u0009t\\u0001c \xc3\xa9 end \"".to_vec()
        );
        let mut out = Vec::new();
        write_json_string(&mut out, b"a\x1bb\x7f\xff\x00c");
        assert_eq!(out, b"\"a\\u001Bb\x7f\xff\"".to_vec());
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
}

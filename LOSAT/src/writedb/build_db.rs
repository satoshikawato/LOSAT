//! A database build (`CBuildDatabase`,
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp`): the log of
//! makeblastdb, the removal of an existing database of the same name, the
//! loop that hands the sequences of a FASTA input to the writer, and the end
//! of the build, which removes the files of a failed one.

use std::io::{Read, Write};
use std::path::Path;

use super::fasta_source::{FastaBioseqSource, Warn};
use super::time::LocalTime;
use super::writer::WriteDb;
use super::WriteDbError;
use crate::seqdb::{SeqDbError, SeqType};

/// `CBuildDatabase` for a FASTA input without parsed Seq-ids.
pub struct BuildDatabase<'a> {
    is_protein: bool,
    log: &'a mut dyn Write,
    output_db: WriteDb,
    output_dbname: String,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:852-927
// ```c++
//     while ( !current.empty()  ||  !tail.empty() ) {
//         std::list<string> pretail;
//         if ( !current.empty() ) {
//             NStr::Split(current, kSep, pretail);
//             current.erase();
//             if (pretail.front().empty()
//     ...
//                 ) {
//                 head.clear();
//     ...
//         if ( !head.empty() ) { // empty heads should accept anything
//             string& last = head.back();
//             if (last == DIR_CURRENT) {
//                 if (!next.empty()) {
//                     head.pop_back();
//                 }
//     ...
//             } else if (next == DIR_CURRENT) {
//                 continue;
//             } else if (next.empty()) {
//                 continue; // leave out empty components in most cases
//             } else if (next == DIR_PARENT) {
//     ...
//                 if (last.empty()) {
//                     continue;
//     ...
//                 } else if (last != DIR_PARENT) {
//                     head.pop_back();
//                     continue;
//                 }
//             }
//         }
// ```
/// `CDirEntry::NormalizePath` of an absolute Unix path, without following
/// links: `.` and empty components go, `..` takes its parent with it.
#[cfg(not(windows))]
fn normalize_absolute_path(path: &str) -> String {
    let current = path.trim_end_matches('/');
    if current.is_empty() {
        return "/".to_string();
    }
    let mut head: Vec<&str> = Vec::new();
    for next in current.split('/') {
        if let Some(&last) = head.last() {
            if last == "." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if next == "." || next.is_empty() {
                continue;
            } else if next == ".." {
                if last.is_empty() {
                    continue;
                } else if last != ".." {
                    head.pop();
                    continue;
                }
            }
        }
        head.push(next);
    }
    if head.len() == 1 && head[0].is_empty() {
        return "/".to_string();
    }
    head.join("/")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:665-670
// ```c++
// string CDirEntry::CreateAbsolutePath(const string& path, ERelativeToWhat rtw)
// {
//     if ( IsAbsolutePath(path) ) {
//         return NormalizePath(path);
//     }
//     string result;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:693-696
// ```c++
//     switch (rtw) {
//         case eRelativeToCwd:
//             result = ConcatPath(CDir::GetCwd(), path);
//             break;
// ```
// The reference platform is Unix. On Windows the path is made absolute by
// the system's own rule; under WASI, where a program sees only the
// directories it was given, the name is used as it is.
/// `CDirEntry::CreateAbsolutePath(dbname)`.
fn create_absolute_path(dbname: &str) -> String {
    #[cfg(windows)]
    {
        std::path::absolute(dbname)
            .map(|path| path.to_string_lossy().into_owned())
            .unwrap_or_else(|_| dbname.to_string())
    }
    #[cfg(target_os = "wasi")]
    {
        dbname.to_string()
    }
    #[cfg(not(any(windows, target_os = "wasi")))]
    {
        if dbname.starts_with('/') {
            return normalize_absolute_path(dbname);
        }
        match std::env::current_dir() {
            Ok(cwd) => {
                let cwd = cwd.to_string_lossy();
                normalize_absolute_path(&format!("{}/{dbname}", cwd.trim_end_matches('/')))
            }
            Err(_) => dbname.to_string(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1051-1071
// ```c++
// void CBuildDatabase::CreateDirectories(const string& dbname)
// {
//     CDirEntry dir_entry(dbname);
//     string dir_name = dir_entry.GetDir(CDirEntry::eIfEmptyPath_Empty);
//     if (dir_name.empty()) {
//         return;
//     }
//
//     CDir d(dir_name);
//     if ( !d.Exists() ) {
//         if ( !d.CreatePath() ) {
//             string msg("Failed to create directory '" + d.GetName() + "'");
//             NCBI_THROW(CMultisourceException, eOutputFileError, msg);
//         }
//     }
//     if (!d.CheckAccess(CDirEntry::fWrite)) {
//         string msg("You do not have write permissions on '" +
//                    d.GetName() + "'");
//         NCBI_THROW(CMultisourceException, eOutputFileError, msg);
//     }
// }
// ```
// `d.GetName()` is the last component of the directory. The write
// permission is not probed here: a directory that cannot be written fails
// when the first file of the database is created.
fn create_directories(dbname: &str) -> Result<(), WriteDbError> {
    let is_separator = |c: char| c == '/' || (cfg!(windows) && c == '\\');
    let Some(at) = dbname.rfind(is_separator) else {
        return Ok(());
    };
    let dir_name = &dbname[..=at];
    let dir = Path::new(dir_name);
    if !dir.is_dir() && std::fs::create_dir_all(dir).is_err() {
        let trimmed = dir_name.trim_end_matches(is_separator);
        let name = trimmed.rsplit(is_separator).next().unwrap_or(trimmed);
        return Err(WriteDbError::Other(format!(
            "Failed to create directory '{name}'"
        )));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2629-2641
// ```c++
//     extn.push_back(kExtnMol + "al");   // alias file
//     extn.push_back(kExtnMol + "in");   // index file
//     extn.push_back(kExtnMol + "hr");   // header file
//     extn.push_back(kExtnMol + "sq");   // sequence file
//     extn.push_back(kExtnMol + "ni");   // ISAM numeric index file
//     extn.push_back(kExtnMol + "nd");   // ISAM numeric data file
//     if (dbver == eBDB_Version4) {
//     	extn.push_back(kExtnMol + "si");   // ISAM string index file
//     	extn.push_back(kExtnMol + "sd");   // ISAM string data file
//     }
//     extn.push_back(kExtnMol + "pi");   // ISAM PIG index file
//     extn.push_back(kExtnMol + "pd");   // ISAM PIG data file
//     extn.push_back(kExtnMol + "js");   // BLASTDB metadata
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2647-2656
// ```c++
//     // Contain masking information
//     extn.push_back(kExtnMol + "aa");   // ISAM mask index file
//     extn.push_back(kExtnMol + "ab");   // ISAM mask data file (big-endian)
//     extn.push_back(kExtnMol + "ac");   // ISAM mask data file (little-endian)
//     extn.push_back(kExtnMol + "og");   // OID to GI file
//     extn.push_back(kExtnMol + "hi");   // ISAM sequence hash index file
//     extn.push_back(kExtnMol + "hd");   // ISAM sequence hash data file
//     extn.push_back(kExtnMol + "ti");   // ISAM trace id index file
//     extn.push_back(kExtnMol + "td");   // ISAM trace id data file
//     extn.push_back(kExtnMol + "xm");   // OID mask file to exclude models
// ```
/// `SeqDB_GetFileExtensions(is_protein, extn, eBDB_Version4)` without the
/// letter of the molecule.
const VOLUME_EXTENSIONS: [&str; 20] = [
    "al", "in", "hr", "sq", "ni", "nd", "si", "sd", "pi", "pd", "js", "aa", "ab", "ac", "og", "hi",
    "hd", "ti", "td", "xm",
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2663-2668
// ```c++
// 	static const char * ext[]={"db", "os", "ot", "tf", "to", "db-lock", "tf-lock", NULL};
// 	extn.clear();
// 	const string kExtnMol(1, db_is_protein ? 'p' : 'n');
// 	for(const char ** p=ext; *p != NULL; p++) {
// 		extn.push_back(kExtnMol + (*p));
// 	}
// ```
/// `SeqDB_GetLMDBFileExtensions` without the letter of the molecule.
const LMDB_EXTENSIONS: [&str; 7] = ["db", "os", "ot", "tf", "to", "db-lock", "tf-lock"];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1542-1593
// ```c++
// bool DeleteBlastDb(const string& dbpath, CSeqDB::ESeqType seq_type)
// {
//     int num_files_removed = 0;
//     vector<string> db_files, alias_files;
//     bool is_protein = (seq_type == CSeqDB::eProtein);
//
//     vector<string> extn;
//     SeqDB_GetFileExtensions( is_protein, extn, eBDB_Version4);
//     vector<string> lmdb_extn;
//     SeqDB_GetLMDBFileExtensions(is_protein, lmdb_extn);
//     ITERATE(vector<string>, lmdb, lmdb_extn) {
//     	CNcbiOstrstream oss;
//     	oss << dbpath << "." << *lmdb;
//         const string fname = CNcbiOstrstreamToString(oss);
//         if (CFile(fname).Remove()) {
//         	LOG_POST(Info << "Deleted " << fname);
//             num_files_removed++;
//         }
//         else {
//         	unsigned int index = 0;
//         	string vfname = dbpath + "." + NStr::IntToString(index/10) +
//         			        NStr::IntToString(index%10) + "." + *lmdb;
//         	while (CFile(vfname).Remove()) {
//         		index++;
//         		vfname = dbpath + "." + NStr::IntToString(index/10) +
//         	    	     NStr::IntToString(index%10) + "." + *lmdb;
//
//         	}
//         }
//     }
//
//     try { CSeqDB::FindVolumePaths(dbpath, seq_type, db_files, &alias_files); }
//     catch (...) {}    // ignore any errors from the invocation above
//     ITERATE(vector<string>, f, db_files) {
//         ITERATE(vector<string>, e, extn) {
//             CNcbiOstrstream oss;
//             oss << *f << "." << *e;
//             const string fname = CNcbiOstrstreamToString(oss);
//             if (CFile(fname).Remove()) {
//                 LOG_POST(Info << "Deleted " << fname);
//                 num_files_removed++;
//             }
//         }
//     }
//     ITERATE(vector<string>, f, alias_files) {
//         if (CFile(*f).Remove()) {
//             LOG_POST(Info << "Deleted " << *f);
//             num_files_removed++;
//         }
//     }
//     return static_cast<bool>(num_files_removed != 0);
// }
// ```
// NCBI finds the volumes of the existing database with its reader and
// ignores a database it cannot read. LOSAT's reader refuses some alias files
// that NCBI reads (filters, replaced totals, group alias files); leaving
// such an alias file in place would hide the new database behind it, so an
// existing database of that kind stops the build before anything is removed.
/// `DeleteBlastDb`: removes the files of the database `dbpath`; `true` when
/// a file was removed.
fn delete_blast_db(dbpath: &str, is_protein: bool) -> Result<bool, WriteDbError> {
    let letter = if is_protein { 'p' } else { 'n' };
    let seqtype = if is_protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let (db_files, alias_files) = match crate::seqdb::alias::find_volume_paths(dbpath, seqtype) {
        Ok(found) => found,
        Err(SeqDbError::Unsupported { path, message }) => {
            return Err(WriteDbError::Other(format!(
                "an existing database is in the way and LOSAT cannot list its files ({}: \
                 {message}); remove it before building {dbpath}",
                path.display()
            )));
        }
        Err(_) => (Vec::new(), Vec::new()),
    };

    let mut num_files_removed = 0;
    let mut remove = |name: &str| {
        let removed = std::fs::remove_file(name).is_ok();
        num_files_removed += u32::from(removed);
        removed
    };
    for lmdb in LMDB_EXTENSIONS {
        if !remove(&format!("{dbpath}.{letter}{lmdb}")) {
            let mut index = 0u32;
            // The removals of the volume files are not counted by NCBI.
            while std::fs::remove_file(format!(
                "{dbpath}.{}{}.{letter}{lmdb}",
                index / 10,
                index % 10
            ))
            .is_ok()
            {
                index += 1;
            }
        }
    }
    for file in &db_files {
        for extension in VOLUME_EXTENSIONS {
            remove(&format!("{file}.{letter}{extension}"));
        }
    }
    for file in &alias_files {
        remove(file);
    }
    Ok(num_files_removed != 0)
}

/// `ostream << double` with the default format (six significant digits,
/// `%g`): the seconds of the "Adding sequences" line.
fn format_seconds(value: f64) -> String {
    if value == 0.0 || !value.is_finite() {
        return "0".to_string();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let strip = |text: &str| {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            text.to_string()
        }
    };
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", strip(mantissa), exponent.abs());
    }
    strip(&format!("{value:.*}", (5 - exponent) as usize))
}

impl<'a> BuildDatabase<'a> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1165-1177
    // ```c++
    //     CreateDirectories(dbname);
    //     const string output_dbname = CDirEntry::CreateAbsolutePath(dbname);
    //     m_LogFile << "\n\nBuilding a new DB, current time: "
    //               << CTime(CTime::eCurrent).AsString() << endl;
    //
    //     m_LogFile << "New DB name:   " << output_dbname << endl;
    //     m_LogFile << "New DB title:  " << title << endl;
    //     const string mol_type(is_protein ? "Protein" : "Nucleotide");
    //     m_LogFile << "Sequence type: " << mol_type << endl;
    //     if (DeleteBlastDb(output_dbname, ParseMoleculeTypeString(mol_type))) {
    //         m_LogFile << "Deleted existing " << mol_type
    //             << " BLAST database named " << output_dbname << endl;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1270-1273
    // ```c++
    // void CBuildDatabase::SetMembBits(const TLinkoutMap & membbits,
    //                                  bool                keep_mbits)
    // {
    //     m_LogFile << "Keep MBits: " << (keep_mbits ? "T" : "F") << endl;
    // ```
    /// The constructor, then `SetMembBits(no_bits, true)` as makeblastdb
    /// calls it. `now` is the time of the build.
    pub fn new(
        dbname: &str,
        title: &str,
        is_protein: bool,
        log: &'a mut dyn Write,
        now: &LocalTime,
    ) -> Result<Self, WriteDbError> {
        create_directories(dbname)?;
        let output_dbname = create_absolute_path(dbname);
        let mol_type = if is_protein { "Protein" } else { "Nucleotide" };
        let _ = writeln!(
            log,
            "\n\nBuilding a new DB, current time: {}",
            now.log_time()
        );
        let _ = writeln!(log, "New DB name:   {output_dbname}");
        let _ = writeln!(log, "New DB title:  {title}");
        let _ = writeln!(log, "Sequence type: {mol_type}");
        if delete_blast_db(&output_dbname, is_protein)? {
            let _ = writeln!(
                log,
                "Deleted existing {mol_type} BLAST database named {output_dbname}"
            );
        }
        let _ = writeln!(log, "Keep MBits: T");
        // The limit of 1000*1000*1000 that the constructor sets
        // (build_db.cpp:1200) is replaced by makeblastdb before a sequence
        // is added.
        let output_db = WriteDb::new(
            &output_dbname,
            is_protein,
            title.as_bytes(),
            now.index_date().as_bytes(),
            1000 * 1000 * 1000,
        );
        Ok(Self {
            is_protein,
            log,
            output_db,
            output_dbname,
        })
    }

    /// `GetOutputDbName`: the absolute name of the database.
    pub fn output_dbname(&self) -> &str {
        &self.output_dbname
    }

    /// The log stream (`m_LogFile`).
    pub fn log(&mut self) -> &mut dyn Write {
        self.log
    }

    /// `SetMaxFileSize`.
    pub fn set_max_file_size(&mut self, max_file_size: u64) {
        self.output_db.set_max_file_size(max_file_size);
    }

    /// Takes the diagnostics the writer posted since the last call.
    pub fn take_diagnostics(&mut self) -> Vec<u8> {
        self.output_db.take_diagnostics()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:859-876
    // ```c++
    //         if(bs->IsAa() != m_IsProtein ){
    //                     bs = src.GetNext();
    //                     continue;
    //        }
    //
    //        if (bs->GetLength() > 0x7fffffff)
    //        {
    //                    NCBI_THROW(CWriteDBException,
    //                    eArgErr,
    //                    "Sequences longer than 2,147,483,647 bases are not supported.  Offending sequence is " + bioseq_id);
    //        }
    //
    //         if ((bs->GetLength() == 0) || (!x_EditAndAddBioseq(bs, NULL, add_pig))){
    //             m_LogFile << "Ignoring sequence '" << bioseq_id
    //                       << "' as it has no sequence data" << endl;
    //             bs = src.GetNext();
    //             continue;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-487
    // ```c++
    //     CRef<CBlast_def_line_set> headers =
    //         CWriteDB::ExtractBioseqDeflines(*bs, m_ParseIDs, m_LongIDs,
    //                 m_ScanBioseq4CFastaReaderUsrObjct);
    //
    //     x_EditHeaders(headers);
    //
    //     // Add the sequence
    //     if (sv) {
    //         m_OutputDb->AddSequence(*bs, *sv);
    //     } else {
    //         bs = s_FixBioseqDeltas(bs);
    //         if(bs->GetInst().CanGetSeq_data())
    //         	m_OutputDb->AddSequence(*bs);
    //         else
    //         	return false;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:894-901
    // ```c++
    //     if (count) {
    //         double t = sw.Elapsed();
    //
    //         m_LogFile << "Adding sequences from FASTA; added "
    //                   << count << " sequences in " << t << " seconds." << endl;
    //     }
    //
    //     return found;
    // ```
    // The FASTA source applies `s_FixBioseqDeltas` before it hands a
    // sequence over; its errors arrive as errors of the source.
    /// `AddSequences(IBioseqSource &)`: `true` when a sequence was added.
    fn add_sequences<R: Read>(
        &mut self,
        source: &mut FastaBioseqSource<R>,
        warn: Warn<'_>,
    ) -> Result<bool, WriteDbError> {
        let mut found = false;
        let start = std::time::Instant::now();
        let mut count = 0u64;
        while let Some(bioseq) = source
            .next_bioseq(warn)
            .map_err(|error| WriteDbError::Creation(error.message))?
        {
            let bioseq_id = format!("lcl|{}", bioseq.local_id);
            if bioseq.is_aa != self.is_protein {
                continue;
            }
            if bioseq.length > 0x7fff_ffff {
                return Err(WriteDbError::Creation(format!(
                    "Sequences longer than 2,147,483,647 bases are not supported.  Offending \
                     sequence is {bioseq_id}"
                )));
            }
            let mut added = false;
            if bioseq.length != 0 {
                self.output_db.extract_bioseq_deflines(&bioseq);
                if bioseq.seq_data.is_some() {
                    self.output_db.add_sequence(bioseq)?;
                    added = true;
                }
            }
            // The serializer posts its diagnostics while the sequence is
            // added, before the reader reads (and warns about) the next
            // record; both go to the standard error stream in that order.
            let diagnostics = self.output_db.take_diagnostics();
            if !diagnostics.is_empty() {
                warn(diagnostics.strip_suffix(b"\n").unwrap_or(&diagnostics));
            }
            if !added {
                let _ = writeln!(
                    self.log,
                    "Ignoring sequence '{bioseq_id}' as it has no sequence data"
                );
                continue;
            }
            found = true;
            count += 1;
        }
        if count != 0 {
            let _ = writeln!(
                self.log,
                "Adding sequences from FASTA; added {count} sequences in {} seconds.",
                format_seconds(start.elapsed().as_secs_f64())
            );
        }
        Ok(found)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1403-1420
    // ```c++
    //     if (fasta_file) {
    //         CFastaBioseqSource fbs(fasta_file,
    //                                m_IsProtein,
    //                                m_ParseIDs,
    //                                m_LongIDs);
    //
    //         try {
    //             success = AddSequences(fbs);
    // 	    if (success == false)
    //             	NCBI_THROW(CWriteDBException, eFileErr, "No sequences added");
    //
    //         }
    //         catch (...) {
    //             EndBuild(true);
    //             throw;
    //         }
    //     }
    //     return success;
    // ```
    /// `AddFasta`: adds the sequences of one FASTA input. An error ends the
    /// build and removes its files.
    pub fn add_fasta<R: Read>(
        &mut self,
        fasta_file: R,
        warn: Warn<'_>,
    ) -> Result<(), WriteDbError> {
        let mut source = FastaBioseqSource::new(fasta_file, self.is_protein);
        let result = match self.add_sequences(&mut source, warn) {
            Ok(true) => Ok(()),
            Ok(false) => Err(WriteDbError::Creation("No sequences added".to_string())),
            Err(error) => Err(error),
        };
        if result.is_err() {
            // The error of the build is the one already in hand.
            let _ = self.end_build(true);
        }
        result
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1423-1475
    // ```c++
    // bool CBuildDatabase::EndBuild(bool erase)
    // {
    //     try {
    //         m_OutputDb->Close();
    //         return x_EndBuild(erase, NULL);
    //     } catch (const CException& e) {
    //         return x_EndBuild(true, erase ? NULL : &e);
    //     ...
    // bool CBuildDatabase::x_EndBuild(bool erase, const CException * close_exception)
    // {
    //     bool success = true;
    //
    //     vector<string> vols;
    //     vector<string> files;
    //
    //     m_OutputDb->ListVolumes(vols);
    //     m_OutputDb->ListFiles(files);
    //
    //     m_LogFile << endl;
    //
    //     _ASSERT(vols.empty() == files.empty());
    //
    //     if (vols.empty()) {
    //         m_LogFile << "No volumes were created."
    //                   << endl;
    //
    //         success = false;
    //     } else {
    //         ITERATE(vector<string>, iterf, files) {
    //             if (erase) {
    //                 CFile(*iterf).Remove();
    //                 _TRACE("Deleting " << *iterf);
    //             }
    //         }
    //     }
    //
    //     m_LogFile << endl;
    //
    //     if (close_exception) {
    //         NCBI_RETHROW(*close_exception, CWriteDBException, eArgErr,
    //                      close_exception->GetMsg());
    //     }
    //
    //     return success;
    // }
    // ```
    /// `EndBuild`: closes the database; with `erase`, or when closing fails,
    /// removes its files. `Ok(false)` when no volume was created.
    pub fn end_build(&mut self, erase: bool) -> Result<bool, WriteDbError> {
        let close_error = self.output_db.close().err();
        let erase = erase || close_error.is_some();
        if erase {
            self.output_db.abandon();
        }
        let vols = self.output_db.list_volumes();
        let files = self.output_db.list_files();
        let _ = writeln!(self.log);
        let mut success = true;
        if vols.is_empty() {
            let _ = writeln!(self.log, "No volumes were created.");
            success = false;
        } else if erase {
            for file in &files {
                let _ = std::fs::remove_file(file);
            }
        }
        let _ = writeln!(self.log);
        match close_error {
            Some(WriteDbError::Io { path, source }) => Err(WriteDbError::Io { path, source }),
            Some(
                WriteDbError::Creation(message)
                | WriteDbError::Options(message)
                | WriteDbError::Other(message),
            ) => Err(WriteDbError::Creation(message)),
            None => Ok(success),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(not(windows))]
    #[test]
    fn absolute_paths_are_normalized_as_ncbi_does() {
        assert_eq!(normalize_absolute_path("/a//b/./c/../d/"), "/a/b/d");
        assert_eq!(normalize_absolute_path("/a/b/.."), "/a");
        assert_eq!(normalize_absolute_path("/../a"), "/a");
        assert_eq!(normalize_absolute_path("/a/./"), "/a");
        assert_eq!(normalize_absolute_path("/"), "/");
        assert_eq!(normalize_absolute_path("/a/.."), "/");
    }

    // Values makeblastdb 2.17.0 printed in "added N sequences in T seconds."
    #[test]
    fn seconds_print_as_six_significant_digits() {
        assert_eq!(format_seconds(0.000254154), "0.000254154");
        assert_eq!(format_seconds(0.00107908), "0.00107908");
        assert_eq!(format_seconds(0.0001719), "0.0001719");
        assert_eq!(format_seconds(0.000169992), "0.000169992");
        assert_eq!(format_seconds(5.5), "5.5");
        assert_eq!(format_seconds(12.3456789), "12.3457");
        assert_eq!(format_seconds(0.00001234567), "1.23457e-05");
        assert_eq!(format_seconds(1234567.0), "1.23457e+06");
        assert_eq!(format_seconds(100.0), "100");
        assert_eq!(format_seconds(0.0), "0");
    }
}

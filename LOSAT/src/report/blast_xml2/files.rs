//! The files of outfmt 13 and 14: the names of the side files and of the entries of the master
//! file, and the side files on the file system (DESIGN_OF5 1.7).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2316-2330
//! ```c++
//! static string s_GetBaseName(const string & baseFile, bool isXML, bool withPath)
//! {
//! 	string dir = kEmptyStr;
//! 	string base = kEmptyStr;
//! 	string ext = kEmptyStr;
//! 	CDirEntry::SplitPath(baseFile, withPath ? &dir:NULL, &base, &ext );
//! 	if(!((isXML && NStr::CompareNocase(ext, ".xml") == 0 ) ||
//! 	     (!isXML && NStr::CompareNocase(ext, ".json") == 0))){
//! 		 base += ext;
//! 	}
//! 	if(withPath)
//! 		return dir + base;
//!
//! 	return base;
//! }
//! ```
//! The side file of the n-th report is `s_GetBaseName(-out, is_xml, true) + "_" + n + ".xml" |
//! ".json"`; the master file lists `s_GetBaseName(-out, true, false) + "_" + n + ext` for both
//! formats (blast_format.cpp:2432, :2449: the XML rule, no directory). So the master of an
//! outfmt 13 run lists `Y.Json_1.json` for `-out d/Y.Json` (the side file is `d/Y_1.json`) and
//! `k_1.json` for `-out d/k.xml` (the side file is `d/k.xml_1.json`): names of files that do not
//! exist. DECISIONS 8 keeps this as NCBI's defined behaviour.

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use super::framing::{SideFiles, Xml2Kind};

/// The characters `CDirEntry::SplitPath` takes as the end of the directory part: `/` on Unix;
/// `:`, `/` and `\` on Windows.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:105-116, :358-377
/// ```c++
/// #  define ALL_SEPARATORS    ":/\\"          // NCBI_OS_MSWIN
/// #  define ALL_SEPARATORS    "/"             // NCBI_OS_UNIX
/// void CDirEntry::SplitPath(const string& path, string* dir,
///                           string* base, string* ext)
/// {
///     // Get file name
///     size_t pos = path.find_last_of(ALL_SEPARATORS);
///     string filename = (pos == NPOS) ? path : path.substr(pos+1);
///     // Get dir
///     if ( dir ) {
///         *dir = (pos == NPOS) ? kEmptyStr : path.substr(0, pos+1);
///     }
///     // Split file name to base and extension
///     pos = filename.rfind('.');
///     if ( base ) {
///         *base = (pos == NPOS) ? filename : filename.substr(0, pos);
///     }
///     if ( ext ) {
///         *ext = (pos == NPOS) ? kEmptyStr : filename.substr(pos);
///     }
/// }
/// ```
fn is_separator(byte: u8, windows: bool) -> bool {
    byte == b'/' || (windows && (byte == b'\\' || byte == b':'))
}

#[cfg(unix)]
fn path_bytes(path: &Path) -> Vec<u8> {
    use std::os::unix::ffi::OsStrExt;
    path.as_os_str().as_bytes().to_vec()
}

/// Off Unix a path is Unicode text to the program: a `-out` name that is not valid Unicode is
/// changed by `to_string_lossy` (and `path_from_bytes`), which no NCBI run on Windows tells apart.
#[cfg(not(unix))]
fn path_bytes(path: &Path) -> Vec<u8> {
    path.to_string_lossy().into_owned().into_bytes()
}

#[cfg(unix)]
fn path_from_bytes(bytes: Vec<u8>) -> PathBuf {
    use std::os::unix::ffi::OsStringExt;
    PathBuf::from(std::ffi::OsString::from_vec(bytes))
}

#[cfg(not(unix))]
fn path_from_bytes(bytes: Vec<u8>) -> PathBuf {
    PathBuf::from(String::from_utf8_lossy(&bytes).into_owned())
}

/// `s_GetBaseName`: the `-out` value without its extension when that is `.xml` (`xml`) or
/// `.json` (not `xml`) in any case, with its directory part (`with_dir`) or without. The
/// extension starts at the last `.` of the file name, so `.xml` alone is a base of nothing with
/// the extension `.xml`, and `a.b.json` has the extension `.json`.
pub fn base_name(out: &Path, xml: bool, with_dir: bool) -> Vec<u8> {
    base_name_of(&path_bytes(out), xml, with_dir, cfg!(windows))
}

/// `base_name` of the bytes of a path, with the separators of Windows (`windows`) or of Unix (the
/// unit tests check both on every platform).
fn base_name_of(path: &[u8], xml: bool, with_dir: bool, windows: bool) -> Vec<u8> {
    let name_at = path
        .iter()
        .rposition(|&byte| is_separator(byte, windows))
        .map_or(0, |at| at + 1);
    let (dir, file) = path.split_at(name_at);
    let (mut base, ext) = match file.iter().rposition(|&byte| byte == b'.') {
        Some(dot) => (file[..dot].to_vec(), &file[dot..]),
        None => (file.to_vec(), &file[file.len()..]),
    };
    let dropped: &[u8] = if xml { b".xml" } else { b".json" };
    if !ext.eq_ignore_ascii_case(dropped) {
        base.extend_from_slice(ext);
    }
    if with_dir {
        let mut name = dir.to_vec();
        name.extend_from_slice(&base);
        name
    } else {
        base
    }
}

/// The side files of a run written to the file system next to `-out`: the n-th report goes to
/// `<base>_<n>.xml` (outfmt 14) or `<base>_<n>.json` (outfmt 13), opened for writing and
/// truncated if it exists (`IOS_BASE::out`; an older file of another name stays as it was).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-775
/// ```c++
/// BlastXML2_FormatReport(const IBlastXML2ReportData* data, string file_name)
/// {
/// 	blastxml2::CBlastOutput2 bxmlout;
/// 		CNcbiOfstream out_stream;
/// 		out_stream.open(file_name.c_str(), IOS_BASE::out);
/// 		if(!out_stream.is_open())
/// 			 NCBI_THROW(CArgException, eInvalidArg, "Cannot open output file");
/// ```
pub struct FsSideFiles {
    /// `s_GetBaseName(out, xml, true)`: the side files' path up to the `_<n>`.
    path_base: Vec<u8>,
    /// `s_GetBaseName(out, true, false)`: the name the master lists, up to the `_<n>`.
    listed_base: Vec<u8>,
    extension: &'static str,
}

impl FsSideFiles {
    /// The side files of the run with `-out` `out` of the kind `kind` (`Files14` or `Files13`).
    pub fn new(out: &Path, kind: Xml2Kind) -> Self {
        let xml = kind != Xml2Kind::Files13;
        Self {
            path_base: base_name(out, xml, true),
            listed_base: base_name(out, true, false),
            extension: if xml { ".xml" } else { ".json" },
        }
    }

    fn name(&self, base: &[u8], n: u32) -> Vec<u8> {
        let mut name = base.to_vec();
        name.push(b'_');
        name.extend_from_slice(n.to_string().as_bytes());
        name.extend_from_slice(self.extension.as_bytes());
        name
    }

    /// The path of the n-th side file.
    pub fn path(&self, n: u32) -> PathBuf {
        path_from_bytes(self.name(&self.path_base, n))
    }
}

impl SideFiles for FsSideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        Ok(Box::new(BufWriter::new(std::fs::File::create(
            self.path(n),
        )?)))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        self.name(&self.listed_base, n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `CDirEntry::SplitPath` ends the directory part at `:`, `/` or `\` on Windows and at `/` alone
    /// elsewhere (ncbifile.cpp:105-116): checked here for both on every platform, since the Windows
    /// job of the CI would be the only one to see a wrong `cfg!(windows)`.
    #[test]
    fn the_separators_of_windows_end_the_directory_part_there_only() {
        let base = |path: &str, xml: bool, with_dir: bool, windows: bool| {
            String::from_utf8(base_name_of(path.as_bytes(), xml, with_dir, windows)).unwrap()
        };
        // (path, xml, with_dir, base on Windows, base elsewhere)
        let cases = [
            ("d\\X.XML", true, true, "d\\X", "d\\X"),
            ("d\\X.XML", true, false, "X", "d\\X"),
            ("C:o.xml", true, true, "C:o", "C:o"),
            ("C:o.xml", true, false, "o", "C:o"),
            ("a/b\\c.json", false, true, "a/b\\c", "a/b\\c"),
            ("a/b\\c.json", false, false, "c", "b\\c"),
            // the dot before a separator is not the file's: `d.xml\o` has no extension on Windows
            ("d.xml\\o", true, true, "d.xml\\o", "d.xml\\o"),
            ("d.xml\\o", true, false, "o", "d.xml\\o"),
            ("C:\\d\\.json", false, true, "C:\\d\\", "C:\\d\\"),
        ];
        for (path, xml, with_dir, on_windows, elsewhere) in cases {
            assert_eq!(
                base(path, xml, with_dir, true),
                on_windows,
                "{path} windows"
            );
            assert_eq!(base(path, xml, with_dir, false), elsewhere, "{path} unix");
        }
    }

    /// The names of the side files and the master's entries on the platform of the build.
    #[cfg(windows)]
    #[test]
    fn windows_paths_name_the_side_files_and_the_master_entries() {
        let files = FsSideFiles::new(Path::new("d\\X.XML"), Xml2Kind::Files14);
        assert_eq!(files.path(1), PathBuf::from("d\\X_1.xml"));
        assert_eq!(files.listed(1), b"X_1.xml");
        let files = FsSideFiles::new(Path::new("C:o.json"), Xml2Kind::Files13);
        assert_eq!(files.path(2), PathBuf::from("C:o_2.json"));
        assert_eq!(files.listed(2), b"o.json_2.json");
    }
}

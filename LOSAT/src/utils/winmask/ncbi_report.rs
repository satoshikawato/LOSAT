//! NCBI's exception reports for the errors of the counts files, as blastn prints them: the
//! `NCBI C++ Exception:` text of `CException::what()` with one `T0` line per exception of the
//! chain.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:370-384
//! ```c
//! string CException::ReportAll(TDiagPostFlags flags) const
//! {
//!     // invert the order
//!     stack<const CException*> pile;
//!     const CException* pex;
//!     for (pex = this; pex; pex = pex->GetPredecessor()) {
//!         pile.push(pex);
//!     }
//!     CNcbiOstrstream os;
//!     os << "NCBI C++ Exception:" << '\n';
//!     for (; !pile.empty(); pile.pop()) {
//!         //indentation
//!         os << "    ";
//!         os << pile.top()->ReportThis(flags) << '\n';
//!     }
//! ```
//! A line is `T0 "<file>", line <n>: Error: (<class>::<error code>) <function>() - <message>`
//! (`CException::ReportStd`); the file is the path of the source file in NCBI's build tree.

use super::{FormatSite, WinMaskError};

/// The build tree that NCBI's reports name: that of the official Linux x64 release of BLAST+
/// 2.17.0 (`ncbi-blast-2.17.0+-x64-linux.tar.gz`). Other builds of 2.17.0 name their own tree
/// and print the same text around it (line NP, M/DECISIONS.md 14).
pub const NCBI_BUILD_ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

/// One `T0` line of the report of an exception thrown in `src/algo/winmask/<file>`.
fn throw_line(file: &str, line: u32, code: &str, function: Option<&str>, message: &str) -> String {
    let function = function
        .map(|name| format!("{name}() - "))
        .unwrap_or_default();
    format!(
        "    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/{file}\", line {line}: Error: (Exception::{code}) {function}{message}\n"
    )
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:165-172
/// ```c
/// CSeqMaskerIstatFactory::EStatType CSeqMaskerIstatFactory::DiscoverStatType(
///         string const & name, vector< string > & md, size_t & skip ) {
///     GetMetaDataLines( name, md, skip );
///     CNcbiIfstream check( name.c_str() );
///
///     if( !check ) {
///         NCBI_THROW( Exception, eOpen, string( "could not open " ) + name );
///     }
/// ```
fn open_line(name: &str) -> String {
    throw_line(
        "seq_masker_istat_factory.cpp",
        171,
        "open failed",
        Some("ncbi::CSeqMaskerIstatFactory::DiscoverStatType"),
        &format!("could not open {name}"),
    )
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:355-366
/// ```c
///     catch( CException & e)
///     {
///         NCBI_RETHROW( e, Exception, eCreateFail,
///                       "could not create a unit counts container" );
///     }
///     catch( std::exception & e )
///     {
///         NCBI_THROW( Exception, eCreateFail,
///                     std::string( "could not create a unit counts container: " ) +
///                     e.what() );
///     }
/// ```
fn create_line(line: u32, message: &str) -> String {
    throw_line(
        "seq_masker_istat_factory.cpp",
        line,
        "creation failure",
        Some("ncbi::CSeqMaskerIstatFactory::create"),
        message,
    )
}

const CREATE_FAIL: &str = "could not create a unit counts container";

impl WinMaskError {
    /// The report (`what()`) of the exception of `DiscoverStatType` when it cannot open the
    /// file; `None` for the other errors (which `DiscoverStatType` does not throw).
    pub fn discover_report(&self) -> Option<String> {
        match self {
            Self::Open { name } => Some(format!("NCBI C++ Exception:\n{}", open_line(name))),
            _ => None,
        }
    }

    /// The report (`what()`) of the exception of `CSeqMaskerIstatFactory::create`: the
    /// reader's exception and the factory's rethrow. `None` for the errors that are not
    /// NCBI's exceptions there (the formats that LOSAT does not read, the masker's checks, a
    /// bad index of the hash table).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:336-340
    /// ```c
    ///             case eOBinary: res = new CSeqMaskerIstatOBinary(
    ///                                    name,
    /// ...
    /// 			default: NCBI_THROW(Exception, eCreateFail, "unrecognized unit counts format");
    /// ```
    pub fn create_report(&self) -> Option<String> {
        let lines = match self {
            Self::Open { name } => format!("{}{}", open_line(name), create_line(358, CREATE_FAIL)),
            Self::Format { site, .. } => {
                let (file, line, code, function) = site.throw_site();
                format!(
                    "{}{}",
                    throw_line(file, line, code, function, site.message()),
                    create_line(358, CREATE_FAIL)
                )
            }
            Self::CreateFail { message } if message.starts_with(CREATE_FAIL) => {
                create_line(363, message)
            }
            Self::CreateFail { message } => {
                format!(
                    "{}{}",
                    create_line(340, message),
                    create_line(358, CREATE_FAIL)
                )
            }
            Self::UnsupportedFormat(_)
            | Self::WindowSize { .. }
            | Self::EmptyCacheBitArray { .. }
            | Self::BadIndex { .. } => return None,
        };
        Some(format!("NCBI C++ Exception:\n{lines}"))
    }

    /// The warnings that NCBI posts while it reads the file, before the error.
    pub fn warnings(&self) -> &[String] {
        match self {
            Self::Format { warnings, .. } => warnings,
            _ => &[],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The texts of NCBI BLAST+ 2.17.0 (the official Linux x64 release) blastn
    // `-window_masker_db` with hand-made counts files (line NP, M/STATE.md "M2 の実機確認").
    #[test]
    fn the_reports_are_ncbis() {
        let root = NCBI_BUILD_ROOT;
        assert_eq!(
            WinMaskError::Open { name: "nosuch.obinary".into() }.discover_report().unwrap(),
            format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open nosuch.obinary\n")
        );
        let factory = format!("    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n");
        let cases = [
            (FormatSite::FileTooShort, "seq_masker_istat_obinary.cpp\", line 63: Error: (Exception::format error) file too short"),
            (FormatSite::UnitSize, "seq_masker_istat_obinary.cpp\", line 114: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - unit size must be in [1,16]"),
            (FormatSite::HashKeySize, "seq_masker_istat_obinary.cpp\", line 125: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - hash key size must be in [1,2*unit_size - 1]"),
            (FormatSite::HashOffset, "seq_masker_istat_obinary.cpp\", line 129: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - offset must by in [0,32 - hash_key_size]"),
            (FormatSite::HashShift, "seq_masker_istat_obinary.cpp\", line 133: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - shift must be in [1, 32 - hash_key_size]"),
            (FormatSite::HashTable, "seq_masker_istat_obinary.cpp\", line 197: Error: (Exception::format error) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - not enough data to fill the hash table"),
            (FormatSite::ValuesTable, "seq_masker_istat_obinary.cpp\", line 208: Error: (Exception::format error) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - not enough data to fill the values table"),
            (FormatSite::BinFileSize, "seq_masker_istat_bin.cpp\", line 79: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - wrong file size"),
            (FormatSite::BinUnitSize, "seq_masker_istat_bin.cpp\", line 99: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - illegal unit size"),
            (FormatSite::BinLength, "seq_masker_istat_bin.cpp\", line 108: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - wrong length"),
        ];
        for (site, line) in cases {
            assert_eq!(
                WinMaskError::format(site).create_report().unwrap(),
                format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/{line}\n{factory}"),
                "{site:?}"
            );
        }
        let substr = "could not create a unit counts container: basic_string::substr: __pos (which is 7) > this->size() (which is 6)";
        assert_eq!(
            WinMaskError::CreateFail { message: substr.into() }.create_report().unwrap(),
            format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 363: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - {substr}\n")
        );
        assert_eq!(
            WinMaskError::BadIndex { key: 1, value: 2 }.create_report(),
            None
        );
    }
}

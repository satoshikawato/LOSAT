//! NCBI's exception reports for the errors of the counts files, as blastn prints them: the
//! exceptions of the readers and of `CSeqMaskerIstatFactory`, whose `what()`
//! ([`NcbiException::what`]) is the `NCBI C++ Exception:` text with one `T0` line per exception
//! of the chain.

use super::WinMaskError;
use crate::ncbi_exception::{NcbiException, ThrowSite};

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
fn open_exception(name: &str) -> NcbiException {
    const OPEN: ThrowSite = ThrowSite {
        file: "src/algo/winmask/seq_masker_istat_factory.cpp",
        line: 171,
        kind: "Exception::open failed",
        function: Some("ncbi::CSeqMaskerIstatFactory::DiscoverStatType"),
    };
    NcbiException::new(OPEN, format!("could not open {name}"))
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
const fn create_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/algo/winmask/seq_masker_istat_factory.cpp",
        line,
        kind: "Exception::creation failure",
        function: Some("ncbi::CSeqMaskerIstatFactory::create"),
    }
}

const CREATE_FAIL: &str = "could not create a unit counts container";

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:461-463
/// ```c
///     case eValidation:
///
///         return "validation error";
/// ```
const WINDOW_SIZE: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker.cpp",
    line: 126,
    kind: "CSeqMaskerException::validation error",
    function: Some("ncbi::CSeqMasker::CSeqMasker"),
};

/// The factory's rethrow of the exception of a reader (line 358).
fn rethrown(e: NcbiException) -> NcbiException {
    e.rethrow(create_site(358), CREATE_FAIL)
}

impl WinMaskError {
    /// The report (`what()`) of the exception of `DiscoverStatType` when it cannot open the
    /// file; `None` for the other errors (which `DiscoverStatType` does not throw).
    pub fn discover_report(&self) -> Option<String> {
        match self {
            Self::Open { name } => Some(open_exception(name).what()),
            _ => None,
        }
    }

    /// The report (`what()`) of the exception of `CSeqMaskerIstatFactory::create`: the
    /// reader's exception and the factory's rethrow, or the masker's window check. `None` for
    /// the errors that are not NCBI's exceptions there (LOSAT's checks of the masker's input, a
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
        let exception = match self {
            Self::Open { name } => rethrown(open_exception(name)),
            Self::Format { site, .. } => {
                rethrown(NcbiException::new(site.throw_site(), site.message()))
            }
            Self::CreateFail { message } if message.starts_with(CREATE_FAIL) => {
                NcbiException::new(create_site(363), message.as_str())
            }
            Self::CreateFail { message } => {
                rethrown(NcbiException::new(create_site(340), message.as_str()))
            }
            Self::Exception { exception, .. } => rethrown((**exception).clone()),
            // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:121-126
            // ```c
            //     if( window_size < ustat->UnitSize() ) {
            //         std::ostringstream os;
            //         os << "window size (" << (int)window_size << ") "
            //               "must be greater or equal to unit size (" <<
            //               (int)ustat->UnitSize() << ")";
            //         NCBI_THROW( CSeqMaskerException, eValidation, os.str() );
            //     }
            // ```
            // Thrown by the masker after the factory returned, so not rethrown by it;
            // `s_BuildSeqMasker` gives its `what()` to blastn's error (`winmask_filter.cpp:147-151`).
            Self::WindowSize { .. } => NcbiException::new(WINDOW_SIZE, self.to_string()),
            Self::UnitSize { .. }
            | Self::EmptyCacheBitArray { .. }
            | Self::SequenceTooLong { .. }
            | Self::BadIndex { .. } => return None,
        };
        Some(exception.what())
    }

    /// The warnings that NCBI posts while it reads the file, before the error.
    pub fn warnings(&self) -> &[String] {
        match self {
            Self::Format { warnings, .. }
            | Self::Exception { warnings, .. }
            | Self::WindowSize { warnings, .. } => warnings,
            _ => &[],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncbi_exception::NCBI_BUILD_ROOT;
    use crate::utils::winmask::FormatSite;

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
            // `-window_masker_taxid` with optimized ascii counts (line NP M3d).
            (FormatSite::OAsciiUnitSize, "seq_masker_istat_oascii.cpp\", line 101: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - unit size must be in [1,16]"),
            (FormatSite::OAsciiHashKeySize, "seq_masker_istat_oascii.cpp\", line 114: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - hash key size must be in [1,2*unit_size - 1]"),
            (FormatSite::OAsciiHashOffset, "seq_masker_istat_oascii.cpp\", line 119: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - offset must by in [0,32 - hash_key_size]"),
            (FormatSite::OAsciiHashShift, "seq_masker_istat_oascii.cpp\", line 124: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - shift must be in [1, 32 - hash_key_size]"),
            (FormatSite::OAsciiFileTooShort, "seq_masker_istat_oascii.cpp\", line 145: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - file too short"),
            (FormatSite::OAsciiHashTable, "seq_masker_istat_oascii.cpp\", line 176: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - not enough lines to fill the hash table"),
            (FormatSite::OAsciiValuesTable, "seq_masker_istat_oascii.cpp\", line 192: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - not enough lines to fill the values table"),
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

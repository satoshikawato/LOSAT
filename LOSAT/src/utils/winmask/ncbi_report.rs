//! NCBI's exception reports for the errors of the counts files, as blastn prints them: the
//! exceptions of the readers and of `CSeqMaskerIstatFactory`, whose `what()`
//! ([`NcbiException::what`]) is the `NCBI C++ Exception:` text with one `T0` line per exception
//! of the chain.

use super::WinMaskError;
use crate::ncbi_exception::{NcbiException, ThrowSite};

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:165-172
fn open_exception(name: &[u8]) -> NcbiException {
    const OPEN: ThrowSite = ThrowSite {
        file: "src/algo/winmask/seq_masker_istat_factory.cpp",
        line: 171,
        kind: "Exception::open failed",
        function: Some("ncbi::CSeqMaskerIstatFactory::DiscoverStatType"),
    };
    NcbiException::new(OPEN, [b"could not open ".as_slice(), name].concat())
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:355-366
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
    /// The exception of `DiscoverStatType` when it cannot open the file; `None` for the other
    /// errors (which `DiscoverStatType` does not throw).
    pub(crate) fn discover_exception(&self) -> Option<NcbiException> {
        match self {
            Self::Open { name } => Some(open_exception(name)),
            _ => None,
        }
    }

    /// The report (`what()`) of the exception of `CSeqMaskerIstatFactory::create`: the
    /// reader's exception and the factory's rethrow, or the masker's window check. `None` for
    /// the errors that are not NCBI's exceptions there (LOSAT's checks of the masker's input, a
    /// bad index of the hash table).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:332-340
    pub fn create_report(&self) -> Option<Vec<u8>> {
        self.create_exception().map(|exception| exception.what())
    }

    /// The exception whose report is [`WinMaskError::create_report`], for the windowmasker
    /// application, which reports it as its own (`CNcbiApplication`'s catch).
    pub fn create_exception(&self) -> Option<NcbiException> {
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
            // Thrown by the masker after the factory returned, so not rethrown by it;
            // `s_BuildSeqMasker` gives its `what()` to blastn's error (`winmask_filter.cpp:147-151`).
            Self::WindowSize { .. } => NcbiException::new(WINDOW_SIZE, self.to_string()),
            Self::UnitSize { .. }
            | Self::UnsetHashParameter { .. }
            | Self::EmptyCacheBitArray { .. }
            | Self::SequenceTooLong { .. }
            | Self::BadIndex { .. }
            | Self::PctLine { .. } => return None,
        };
        Some(exception)
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
#[path = "ncbi_report_tests.rs"]
mod tests;

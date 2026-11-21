//! The diagnostic build: the trace, debug, dump and profile variables that no test, script or CI
//! step sets are read only when the `diagnostics` feature is on.
//!
//! A default build ignores these variables: `var` and `var_os` return "not set" without reading the
//! environment, so the code behind them (the tracing statics and their printing paths) is dead and
//! the compiler drops it, and `main` writes one notice line to stderr for each of them that is set
//! (`ignored_notices`). Build with `cargo build --release --features diagnostics` to use them (the
//! WASI command adds `--target wasm32-wasip1 --no-default-features`). Both builds compile the same
//! code; there are no `#[cfg(feature = "diagnostics")]` blocks, only the constant `ON`.
//!
//! The rule for a new diagnostic variable: one that no test, script or CI step sets is read through
//! this module and listed in `NAMES` and in the list of `docs/agents/diagnostics.md` ("Diagnostic
//! build"); one that they set is read directly with `std::env`, because a test of a default build
//! would not see it otherwise. The tests below check the list, the docs block, and that no other
//! file reads a listed name directly.
//!
//! No NCBI counterpart: these are LOSAT debugging variables and NCBI BLAST+ has no equivalent;
//! this module only decides whether LOSAT reads them, and it does not change any value NCBI
//! computes.

use std::env::VarError;
use std::ffi::OsString;

/// Whether this build reads the diagnostic variables.
pub const ON: bool = cfg!(feature = "diagnostics");

/// The diagnostic variables, sorted. `docs/agents/diagnostics.md` lists the same names.
pub const NAMES: &[&str] = &[
    "BLEMIR_DEBUG",
    "BLEMIR_DEBUG_WINDOW",
    "LOSAT_DEBUG_BLASTP_TRACE_HSP",
    "LOSAT_DEBUG_BLASTP_WSSV",
    "LOSAT_DEBUG_CHAINING",
    "LOSAT_DEBUG_CUTOFFS_ALL",
    "LOSAT_DEBUG_EXTENSION",
    "LOSAT_DEBUG_HSP_SAVING",
    "LOSAT_DEBUG_SCAN_SOFF",
    "LOSAT_DEBUG_SCAN_WINDOW",
    "LOSAT_DUMP_TBLASTX_STAGE",
    "LOSAT_LINK_WIDE_KEYS",
    "LOSAT_STARTUP_TRACE",
    "LOSAT_TRACE_CHAIN_HSP",
    "LOSAT_TRACE_HSP",
    "LOSAT_TRACE_HSP_MASKS",
    "LOSAT_TRACE_LINK_SELECTIONS",
    "LOSAT_X_ADJMEMO_STATS",
    "LOSAT_X_DPCAPTURE",
    "LOSAT_X_DPCAPTURE_EVERY",
    "LOSAT_X_NEWTONEXACT_PROF",
];

/// `std::env::var_os(name)` in a diagnostic build, `None` otherwise.
///
/// No NCBI counterpart: reads a LOSAT debugging variable; it does not change any value NCBI
/// computes.
#[inline]
pub fn var_os(name: &str) -> Option<OsString> {
    debug_assert!(NAMES.contains(&name), "{name} is not in xdiag::NAMES");
    if ON {
        std::env::var_os(name)
    } else {
        None
    }
}

/// `std::env::var(name)` in a diagnostic build, `Err(VarError::NotPresent)` otherwise.
///
/// No NCBI counterpart: reads a LOSAT debugging variable; it does not change any value NCBI
/// computes.
#[inline]
pub fn var(name: &str) -> Result<String, VarError> {
    debug_assert!(NAMES.contains(&name), "{name} is not in xdiag::NAMES");
    if ON {
        std::env::var(name)
    } else {
        Err(VarError::NotPresent)
    }
}

/// The stderr line for a diagnostic variable that a build without the feature ignores.
///
/// No NCBI counterpart: LOSAT message about a LOSAT variable; it does not change any value NCBI
/// computes.
pub fn notice(name: &str) -> String {
    format!(
        "{name} is ignored: this losat was built without the diagnostics feature \
         (docs/agents/diagnostics.md, \"Diagnostic build\")"
    )
}

/// One notice per diagnostic variable that is set, in `NAMES` order; empty in a diagnostic build.
/// The loop asks `getenv` for each name, as a diagnostic build would, so the match rules of the
/// platform (case-insensitive names on Windows) are the same.
///
/// No NCBI counterpart: LOSAT message about LOSAT variables; it does not change any value NCBI
/// computes.
pub fn ignored_notices() -> Vec<String> {
    if ON {
        return Vec::new();
    }
    NAMES
        .iter()
        .filter(|name| std::env::var_os(name).is_some())
        .map(|name| notice(name))
        .collect()
}

#[cfg(test)]
#[path = "xdiag_tests.rs"]
mod tests;

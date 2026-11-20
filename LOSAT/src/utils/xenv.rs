//! The debug switches that blastn tests once per extension, read from the environment once per
//! process; `LOSAT_DEBUG_CUTOFFS` (tblastx and blastn, once per subject or per call) and
//! `LOSAT_DEBUG_BLASTN` (once per purge) are read the same way.
//!
//! `getenv` walks the whole environment on every call. In megablast and
//! blastn it is called a few times per ungapped and gapped extension to ask
//! whether `LOSAT_DEBUG_COORDS` / `LOSAT_DEBUG_COORDS_START` are set, which
//! shows up as 1-2 % of the run (more under WASI). Each function below returns
//! what the same expression returned the first time it was evaluated; the only
//! difference from reading the environment every time is for a process that
//! changes these variables while it is running. (Graduated 2026-10-22 from
//! `LOSAT_X_ENVCACHE`; the variable is gone.)
//!
//! No NCBI counterpart: LOSAT_DEBUG_COORDS, LOSAT_DEBUG_COORDS_START, LOSAT_DEBUG_CUTOFFS and
//! LOSAT_DEBUG_BLASTN are LOSAT debugging variables. NCBI BLAST+ has no equivalent, and this
//! module only changes how often the environment is read; it does not change any value NCBI
//! computes.

use std::ffi::OsString;
use std::sync::OnceLock;

/// `std::env::var("LOSAT_DEBUG_COORDS").is_ok()`
#[inline]
pub(crate) fn debug_coords_is_ok() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var("LOSAT_DEBUG_COORDS").is_ok())
}

/// `std::env::var_os("LOSAT_DEBUG_COORDS").is_some()`
#[inline]
pub(crate) fn debug_coords_is_some() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var_os("LOSAT_DEBUG_COORDS").is_some())
}

/// `std::env::var_os("LOSAT_DEBUG_COORDS_START")`
#[inline]
pub(crate) fn debug_coords_start() -> Option<OsString> {
    static VALUE: OnceLock<Option<OsString>> = OnceLock::new();
    VALUE
        .get_or_init(|| std::env::var_os("LOSAT_DEBUG_COORDS_START"))
        .clone()
}

/// `std::env::var_os("LOSAT_DEBUG_CUTOFFS").is_some()`
#[inline]
pub(crate) fn debug_cutoffs_is_some() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var_os("LOSAT_DEBUG_CUTOFFS").is_some())
}

/// `std::env::var("LOSAT_DEBUG_CUTOFFS").is_ok()`
#[inline]
pub(crate) fn debug_cutoffs_is_ok() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var("LOSAT_DEBUG_CUTOFFS").is_ok())
}

/// `std::env::var("LOSAT_DEBUG_BLASTN").is_ok()`
#[inline]
pub(crate) fn debug_blastn_is_ok() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var("LOSAT_DEBUG_BLASTN").is_ok())
}

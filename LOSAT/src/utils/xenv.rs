//! EXPERIMENT (LOSAT_X_ENVCACHE): the debug switches that blastn tests once
//! per extension are read from the environment once per process instead.
//!
//! `getenv` walks the whole environment on every call. In megablast and
//! blastn it is called a few times per ungapped and gapped extension to ask
//! whether `LOSAT_DEBUG_COORDS` / `LOSAT_DEBUG_COORDS_START` are set, which
//! shows up as 1-2 % of the run (more under WASI). With the switch on, each
//! function below returns what the same expression returned the first time
//! it was evaluated; the only difference from reading the environment every
//! time is for a process that changes these variables while it is running.

use std::ffi::OsString;
use std::sync::OnceLock;

fn cached() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_ENVCACHE").is_some())
}

/// `std::env::var("LOSAT_DEBUG_COORDS").is_ok()`
#[inline]
pub(crate) fn debug_coords_is_ok() -> bool {
    if cached() {
        static VALUE: OnceLock<bool> = OnceLock::new();
        *VALUE.get_or_init(|| std::env::var("LOSAT_DEBUG_COORDS").is_ok())
    } else {
        std::env::var("LOSAT_DEBUG_COORDS").is_ok()
    }
}

/// `std::env::var_os("LOSAT_DEBUG_COORDS").is_some()`
#[inline]
pub(crate) fn debug_coords_is_some() -> bool {
    if cached() {
        static VALUE: OnceLock<bool> = OnceLock::new();
        *VALUE.get_or_init(|| std::env::var_os("LOSAT_DEBUG_COORDS").is_some())
    } else {
        std::env::var_os("LOSAT_DEBUG_COORDS").is_some()
    }
}

/// `std::env::var_os("LOSAT_DEBUG_COORDS_START")`
#[inline]
pub(crate) fn debug_coords_start() -> Option<OsString> {
    if cached() {
        static VALUE: OnceLock<Option<OsString>> = OnceLock::new();
        VALUE
            .get_or_init(|| std::env::var_os("LOSAT_DEBUG_COORDS_START"))
            .clone()
    } else {
        std::env::var_os("LOSAT_DEBUG_COORDS_START")
    }
}

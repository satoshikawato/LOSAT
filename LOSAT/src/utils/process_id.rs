//! The process id in the names of the files that LOSAT makes: the inflated copy of a compressed
//! input (`blastinput::input_files`), the staging file of psiblast's ASN.1 PSSM
//! (`algorithm::psiblast::pssm`), and the records of the diagnostics `LOSAT_X_REACH` and
//! `LOSAT_X_DPCAPTURE`.
//!
//! No NCBI counterpart: these names are LOSAT's own; each maker opens its file with `create_new`
//! and tries the next name when one is taken, so the id only keeps the names of two processes
//! apart.

/// `std::process::id()`.
#[cfg(not(target_family = "wasm"))]
pub(crate) fn process_id() -> u32 {
    std::process::id()
}

/// 0 on the Wasm targets, which have no process id: Rust's `std::process::id()` panics there
/// ("unsupported", `library/std/src/sys/pal/wasip1/os.rs:131` of 1.92 on `wasm32-wasip1` and
/// `wasm32-wasip1-threads`). A Wasm command is one process in its host, and the makers' other
/// parts of the name (a time stamp, a counter) and `create_new` keep its names apart.
#[cfg(target_family = "wasm")]
pub(crate) fn process_id() -> u32 {
    0
}

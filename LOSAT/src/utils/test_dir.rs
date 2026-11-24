//! The directory of the tests' temporary files, the number in their names, and the name NCBI's
//! database reader gives such a directory (`ncbi_dir`).
//!
//! No NCBI counterpart for the first two: test plumbing; the product code that needs a temporary
//! directory handles WASI itself (`blastinput::input_files::temporary_directory`).

/// `std::env::temp_dir()` off WASI. On WASI Rust's `temp_dir()` panics ("no filesystem on wasm",
/// `library/std/src/sys/pal/wasip1/os.rs:119` of 1.92), and the panic aborts the whole test binary
/// (`panic=abort` on `wasm32`). There a host gives a module directories only by its preopens, so
/// the `TMPDIR` variable names the directory when it is set and not empty, and `/tmp` is the
/// default (`web/tools/wasi-test-runner.mjs` preopens `/`).
///
/// The tests call this instead of `std::env::temp_dir()`; `clippy.toml` refuses the latter.
pub(crate) fn temp_dir() -> std::path::PathBuf {
    #[cfg(not(target_os = "wasi"))]
    {
        #[allow(clippy::disallowed_methods)] // the one call off WASI; WASI takes the other branch
        std::env::temp_dir()
    }
    #[cfg(target_os = "wasi")]
    {
        match std::env::var_os("TMPDIR") {
            Some(directory) if !directory.is_empty() => std::path::PathBuf::from(directory),
            _ => std::path::PathBuf::from("/tmp"),
        }
    }
}

/// `std::process::id()` off WASI. On WASI Rust's `process::id()` panics ("unsupported",
/// `library/std/src/sys/pal/wasip1/os.rs:131` of 1.92) and the panic aborts the test binary. There
/// the number is made once per process from the clock, so that the tests of two runs at the same
/// time (the sweep of `web/tools/wasi-test-runner.mjs` runs) do not share their directories.
///
/// The tests call this instead of `std::process::id()`; `clippy.toml` refuses the latter.
pub(crate) fn process_id() -> u32 {
    #[cfg(not(target_os = "wasi"))]
    {
        #[allow(clippy::disallowed_methods)] // the one call off WASI; WASI takes the other branch
        std::process::id()
    }
    #[cfg(target_os = "wasi")]
    {
        static ID: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
        *ID.get_or_init(|| {
            let since = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default();
            // Seconds and nanoseconds folded into 32 bits; a clock that is not unique between two
            // starts gives the same number to both, which only a test run on a coarse clock meets.
            (since.as_secs() as u32).wrapping_mul(1_000_003) ^ since.subsec_nanos()
        })
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-972
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:325-335
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
// `NormalizePath` with `eFollowLinks` reads the symbolic link of each component on Unix only; the
// reader names the databases it opens (and with them their volumes, columns and alias files'
// lists) and the directories of its search path (`BLASTDB`) through it.
/// The directory `dir` as NCBI's database reader names it: its symbolic links followed on Unix,
/// as it is elsewhere. On macOS the temporary directory lies under `/var`, a link to
/// `/private/var`, so NCBI and LOSAT name its files under `/private/var` (Platform checks
/// 45315823856).
pub(crate) fn ncbi_dir(dir: &std::path::Path) -> std::path::PathBuf {
    if cfg!(unix) {
        dir.canonicalize().unwrap()
    } else {
        dir.to_path_buf()
    }
}

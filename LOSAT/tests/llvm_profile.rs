//! The lines of LLVM's profile runtime in the standard error of a LOSAT child process. An
//! instrumented build (cargo llvm-cov, the nightly coverage jobs) writes its profile when the
//! process exits; under a file-size limit (`prlimit --fsize`) that write fails, and the runtime
//! writes "LLVM Profile Error: Failed to write file ...: File too large". Those lines are the
//! instrumentation's, not LOSAT's: the tests that compare a limited run's standard error leave
//! them out (`without_profile_errors`).

/// `stderr` without the lines of LLVM's profile runtime ("LLVM Profile Error: ...") in an
/// instrumented build (`LLVM_PROFILE_FILE` set); `stderr` as it is otherwise.
pub(crate) fn without_profile_errors(stderr: &[u8]) -> Vec<u8> {
    if std::env::var_os("LLVM_PROFILE_FILE").is_none() {
        return stderr.to_vec();
    }
    stderr
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| !line.starts_with(b"LLVM Profile Error: "))
        .flatten()
        .copied()
        .collect()
}

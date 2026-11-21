use super::*;

/// LOSAT Web writes glibc's text where wasi-libc's differs: an EIO read of a picked file that
/// changed is NCBI's `Input/output error`, not `I/O error` (line WQA D18, 2026-11-21).
#[test]
fn the_web_writes_the_oracle_texts_of_the_wasi_codes() {
    assert_eq!(oracle_text(29), Some("Input/output error"));
    assert_eq!(oracle_text(2), Some("Permission denied"));
    assert_eq!(oracle_text(44), Some("No such file or directory"));
    // ENOTCAPABLE has no Linux code: LOSAT's own words.
    assert_eq!(oracle_text(76), None);
    let wasi: std::collections::BTreeSet<i32> = ORACLE_TEXTS.iter().map(|e| e.0).collect();
    assert_eq!(wasi.len(), ORACLE_TEXTS.len(), "one text per WASI code");
}

/// On Linux the table is the system's own text of each Linux code, so the Web's bytes are the
/// native CLI's.
#[cfg(target_os = "linux")]
#[test]
fn the_oracle_texts_are_the_linux_texts() {
    for &(_, linux, text) in &ORACLE_TEXTS {
        let error = std::io::Error::from_raw_os_error(linux);
        assert_eq!(ncbi_strerror(&error).as_deref(), Some(text), "code {linux}");
    }
}

#[test]
fn an_error_without_a_system_code_has_no_text() {
    assert_eq!(
        ncbi_strerror(&std::io::Error::from(std::io::ErrorKind::UnexpectedEof)),
        None
    );
}

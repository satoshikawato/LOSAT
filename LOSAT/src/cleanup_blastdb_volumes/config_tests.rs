use super::*;

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:47-59
#[test]
fn parses_blastdb_setting_and_ignores_empty_file() {
    assert_eq!(parse_blastdb_from_ini("", "empty"), Ok(None));
    assert_eq!(
        parse_blastdb_from_ini("[BLAST]\nBLASTDB = /blast/db/blast\n", "cfg"),
        Ok(Some(b"/blast/db/blast".to_vec()))
    );
    assert_eq!(
        parse_blastdb_from_ini("# comment\n; semi\n[BLAST]\n  blastdb : /opt/db  \n", "cfg"),
        Ok(Some(b"/opt/db".to_vec()))
    );
    // Section name is case-sensitive in Python ConfigParser.
    assert_eq!(
        parse_blastdb_from_ini("[blast]\nBLASTDB = /opt/db\n", "cfg"),
        Ok(None)
    );
    // DEFAULT applies only when [BLAST] section exists.
    assert_eq!(
        parse_blastdb_from_ini("[DEFAULT]\nBLASTDB = /opt/db\n", "cfg"),
        Ok(None)
    );
    assert_eq!(
        parse_blastdb_from_ini("[DEFAULT]\nBLASTDB = /opt/db\n[BLAST]\n", "cfg"),
        Ok(Some(b"/opt/db".to_vec()))
    );
    // Basic interpolation.
    assert_eq!(
        parse_blastdb_from_ini(
            "[DEFAULT]\nroot = /data\n[BLAST]\nBLASTDB = %(root)s/blast%%1\n",
            "cfg"
        ),
        Ok(Some(b"/data/blast%1".to_vec()))
    );
    // Trailing comment / text after `]` on a section header is allowed by Python SECTCRE.
    assert_eq!(
        parse_blastdb_from_ini("[BLAST] # comment\nBLASTDB = /opt/db\n", "cfg"),
        Ok(Some(b"/opt/db".to_vec()))
    );
    assert_eq!(
        parse_blastdb_from_ini("[BLAST]extra]\nBLASTDB = /opt/db\n", "cfg"),
        Ok(None)
    );
    // 10 interpolation steps succeed; 11 fail.
    let mut ten = String::from("[BLAST]\nk0 = /opt/db\n");
    for i in 1..=9 {
        ten.push_str(&format!("k{i} = %(k{})s\n", i - 1));
    }
    let eleven = format!("{ten}k10 = %(k9)s\nBLASTDB = %(k10)s\n");
    ten.push_str("BLASTDB = %(k9)s\n");
    assert_eq!(
        parse_blastdb_from_ini(&ten, "cfg"),
        Ok(Some(b"/opt/db".to_vec()))
    );
    assert!(parse_blastdb_from_ini(&eleven, "cfg").is_err());
}

#[test]
fn python_whitespace_and_unicode_option_names_match_configparser() {
    // U+001C–U+001F are Python whitespace: trimmed from keys/values and count as indentation.
    assert_eq!(
        parse_blastdb_from_ini(
            "[BLAST]\n\u{001c}BLASTDB = \u{001d}/opt/db\u{001e}\n",
            "cfg"
        ),
        Ok(Some(b"/opt/db".to_vec()))
    );
    // Leading U+001C after a 0-indent option makes a continuation line, not a duplicate option.
    assert_eq!(
        parse_blastdb_from_ini("[BLAST]\nBLASTDB = /a\n\u{001c}BLASTDB = /b\n", "cfg"),
        Ok(Some(b"/a\nBLASTDB = /b".to_vec()))
    );
    // Equal indent of U+001C produces DuplicateOptionError (both trim to `blastdb`).
    assert!(parse_blastdb_from_ini(
        "[BLAST]\n\u{001c}BLASTDB = /a\n\u{001f}BLASTDB = /b\n",
        "cfg"
    )
    .is_err());
    // Indent is measured in Unicode characters, so 1-char NBSP (2 UTF-8 bytes) is not > 1 space.
    assert!(parse_blastdb_from_ini("[BLAST]\n k = v1\n\u{00a0}v2\n", "cfg").is_err());
    assert_eq!(
        parse_blastdb_from_ini("[BLAST]\n\u{00a0}BLASTDB = /a\n  /b\n", "cfg"),
        Ok(Some(b"/a\n/b".to_vec()))
    );
    // Option names and interpolation keys use Unicode lowercasing.
    assert_eq!(
        parse_blastdb_from_ini(
            "[BLAST]\n\u{00c4} = /opt/db\nBLASTDB = %(\u{00e4})s\n",
            "cfg"
        ),
        Ok(Some(b"/opt/db".to_vec()))
    );
    assert!(parse_blastdb_from_ini("[BLAST]\n\u{00c4} = /a\n\u{00e4} = /b\n", "cfg").is_err());
    assert_eq!(
        parse_blastdb_from_ini(
            "[BLAST]\n\u{0130} = /opt/db\nBLASTDB = %(i\u{0307})s\n",
            "cfg"
        ),
        Ok(Some(b"/opt/db".to_vec()))
    );
}

#[test]
fn malformed_ini_returns_error() {
    assert!(parse_blastdb_from_ini("no_section = 1\n", "bad").is_err());
    assert!(parse_blastdb_from_ini("[BLAST]\nnot_an_option\n", "bad").is_err());
    assert!(parse_blastdb_from_ini("[BLAST]\nBLASTDB = a\nblastdb = b\n", "bad").is_err());
    assert!(parse_blastdb_from_ini("[BLAST]\nBLASTDB = a\n[BLAST]\n", "bad").is_err());
    assert!(parse_blastdb_from_ini("[BLAST]\nBLASTDB = %()s\n", "bad").is_err());
    assert!(parse_blastdb_from_ini("[BLAST]\na)b = x\nBLASTDB = %(a)b)s\n", "bad").is_err());
}

#[cfg(unix)]
#[test]
fn dev_null_returns_none() {
    assert_eq!(
        get_blastdb_from_ncbi_config(Path::new("/dev/null")),
        Ok(None)
    );
}

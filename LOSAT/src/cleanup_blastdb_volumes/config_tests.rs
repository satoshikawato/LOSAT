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

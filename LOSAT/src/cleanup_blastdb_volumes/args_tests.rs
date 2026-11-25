use std::ffi::OsString;

use super::*;

fn os(words: &[&str]) -> Vec<OsString> {
    words.iter().map(OsString::from).collect()
}

#[test]
fn help_and_version_flags() {
    for words in [
        &["-h"][..],
        &["--help"],
        &["-help"],
        &["-he"],
        &["-hx"],
        &["-help=1"],
        &["--h"],
        &["--he"],
        &["--hel"],
        &["-h", "-bogus"],
        &["-bogus", "-h"],
        &["-h", "-version"],
        &["-db", "x", "-dbtype", "prot", "-h"],
    ] {
        assert_eq!(
            parse(&os(words)),
            ParseOutcome::Help(help_text()),
            "{words:?}"
        );
    }

    let ver = format!(
        "cleanup-blastdb-volumes.py {}\n",
        crate::wording::LOSAT_VERSION
    );
    for words in [
        &["-version"][..],
        &["-vers"],
        &["-versi"],
        &["-versio"],
        &["-version", "-bogus"],
        &["-bogus", "-version"],
        &["-version", "-h"],
        &["-db", "x", "-dbtype", "prot", "-vers"],
    ] {
        assert_eq!(
            parse(&os(words)),
            ParseOutcome::Version(ver.clone()),
            "{words:?}"
        );
    }
}

#[test]
fn valid_options_and_abbreviations() {
    assert_eq!(
        parse(&os(&["-db", "mydb", "-dbtype", "prot"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"mydb".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
    assert_eq!(
        parse(&os(&["-dbt=nucl", "-db=-1.5", "-dr", "-verb", "-verbose",])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-1.5".to_vec(),
            dbtype: DbType::Nucl,
            dry_run: true,
            verbose: 2,
        })
    );
    assert_eq!(
        parse(&os(&["-db", "-", "-dbtype", "nucl"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-".to_vec(),
            dbtype: DbType::Nucl,
            dry_run: false,
            verbose: 0,
        })
    );
    assert_eq!(
        parse(&os(&["-db", "-1", "-dbtype", "prot"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-1".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
    assert_eq!(
        parse(&os(&["-db", "-.5", "-dbtype", "prot"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-.5".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
    assert_eq!(
        parse(&os(&["-db=-foo", "-dbtype=prot", "-db", "second"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"second".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
}

#[test]
fn error_messages_match_argparse() {
    let err = |msg: &str| ParseOutcome::Error(format!("{USAGE}{PROG}: error: {msg}\n"));

    let need_both = err("the following arguments are required: -db, -dbtype");
    for words in [&[][..], &["--version"]] {
        assert_eq!(parse(&os(words)), need_both, "{words:?}");
    }
    assert_eq!(
        parse(&os(&["-db", "x"])),
        err("the following arguments are required: -dbtype")
    );
    assert_eq!(
        parse(&os(&["-dbtype", "prot"])),
        err("the following arguments are required: -db")
    );
    assert_eq!(
        parse(&os(&["-v"])),
        err("ambiguous option: -v could match -version, -verbose")
    );
    assert_eq!(
        parse(&os(&["-ver"])),
        err("ambiguous option: -ver could match -version, -verbose")
    );
    assert_eq!(
        parse(&os(&["-d"])),
        err("ambiguous option: -d could match -db, -dbtype, -dry-run")
    );
    let need_db_arg = err("argument -db: expected one argument");
    for words in [
        &["-db"][..],
        &["-db", "-h", "-dbtype", "prot"],
        &["-db", "-1a", "-dbtype", "prot"],
        &["-db", "--", "-dbtype", "prot"],
        &["-db", "-h bar", "-dbtype", "prot"],
        &["-db", "-db=x y", "-dbtype", "prot"],
    ] {
        assert_eq!(parse(&os(words)), need_db_arg, "{words:?}");
    }
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "bad"])),
        err("argument -dbtype: invalid choice: 'bad' (choose from prot, nucl)")
    );
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "prot", "-dry-run=1"])),
        err("argument -dry-run: ignored explicit argument '1'")
    );
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "prot", "-verbose=1"])),
        err("argument -verbose: ignored explicit argument '1'")
    );
    let help_ignored_one = err("argument -h/--help: ignored explicit argument '1'");
    for words in [&["--help=1"][..], &["-h=1"], &["--=1"]] {
        assert_eq!(parse(&os(words)), help_ignored_one, "{words:?}");
    }
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "prot", "--", "extra"])),
        err("unrecognized arguments: -- extra")
    );
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "prot", "--dry-run"])),
        err("unrecognized arguments: --dry-run")
    );
    assert_eq!(
        parse(&os(&["-v=1"])),
        err("ambiguous option: -v=1 could match -version, -verbose")
    );
    assert_eq!(
        parse(&os(&["-d=1"])),
        err("ambiguous option: -d=1 could match -db, -dbtype, -dry-run")
    );
    assert_eq!(
        parse(&os(&["-h-"])),
        err("argument -h/--help: ignored explicit argument '-'")
    );
    assert_eq!(
        parse(&os(&["-h-foo"])),
        err("argument -h/--help: ignored explicit argument '-foo'")
    );
    assert_eq!(
        parse(&os(&["--="])),
        err("argument -h/--help: ignored explicit argument ''")
    );
    assert_eq!(
        parse(&os(&["-h=a'b"])),
        err("argument -h/--help: ignored explicit argument \"a'b\"")
    );
    assert_eq!(
        parse(&os(&["-db", "x", "-dbtype", "a'\"b"])),
        err("argument -dbtype: invalid choice: 'a\\'\"b' (choose from prot, nucl)")
    );
    assert_eq!(
        parse(&os(&["-db", "-foo bar", "-dbtype", "prot"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-foo bar".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
    assert_eq!(
        parse(&os(&["-db", "-- foo", "-dbtype", "prot"])),
        ParseOutcome::Run(ParsedArgs {
            db: b"-- foo".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
}

#[cfg(unix)]
#[test]
fn non_utf8_db_name_bytes_are_preserved() {
    use std::os::unix::ffi::OsStringExt;
    let words = vec![
        OsString::from("-db"),
        OsString::from_vec(b"dir_\xff/mydb".to_vec()),
        OsString::from("-dbtype"),
        OsString::from("prot"),
    ];
    assert_eq!(
        parse(&words),
        ParseOutcome::Run(ParsedArgs {
            db: b"dir_\xff/mydb".to_vec(),
            dbtype: DbType::Prot,
            dry_run: false,
            verbose: 0,
        })
    );
}

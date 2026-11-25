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
fn resolve_columns_follows_python_terminal_size_precedence() {
    assert_eq!(resolve_columns(None, None), 80);
    assert_eq!(resolve_columns(Some(""), None), 80);
    assert_eq!(resolve_columns(Some("abc"), Some(120)), 120);
    assert_eq!(resolve_columns(Some("0"), Some(120)), 120);
    assert_eq!(resolve_columns(Some("-5"), Some(120)), 120);
    assert_eq!(resolve_columns(Some("40.5"), Some(120)), 120);
    assert_eq!(resolve_columns(Some("  40  "), Some(120)), 40);
    assert_eq!(resolve_columns(Some("+60"), Some(120)), 60);
    assert_eq!(resolve_columns(Some("00050"), Some(120)), 50);
    assert_eq!(
        resolve_columns(Some("999999999999999999999"), None),
        usize::MAX
    );
    assert_eq!(resolve_columns(None, Some(0)), 80);
    assert_eq!(resolve_columns(None, Some(100)), 100);
}

#[test]
fn help_and_usage_wrap_by_terminal_columns() {
    // Default 80 columns.
    assert_eq!(
        usage_text(80),
        "usage: cleanup-blastdb-volumes.py [-h] -db DB -dbtype {prot,nucl} [-dry-run]\n\
         \x20                                 [-version] [-verbose]\n"
    );
    // 100 columns fits usage on a single line (98 chars <= 100 - 2).
    assert_eq!(
        usage_text(100),
        "usage: cleanup-blastdb-volumes.py [-h] -db DB -dbtype {prot,nucl} [-dry-run] [-version] [-verbose]\n"
    );
    assert_eq!(usage_text(200), usage_text(100));
    // 99 columns wraps the last token (98 chars > 99 - 2 = 97).
    assert_eq!(
        usage_text(99),
        "usage: cleanup-blastdb-volumes.py [-h] -db DB -dbtype {prot,nucl} [-dry-run] [-version]\n\
         \x20                                 [-verbose]\n"
    );
    // 45 columns puts usage options on the next line indented under "usage: " (7 spaces).
    let indent7 = " ".repeat(7);
    assert_eq!(
        usage_text(45),
        format!(
            "usage: cleanup-blastdb-volumes.py\n\
             {indent7}[-h] -db DB -dbtype {{prot,nucl}}\n\
             {indent7}[-dry-run] [-version] [-verbose]\n"
        )
    );
    // 46 columns keeps [-h] on the first line (indent 34).
    let indent34 = " ".repeat(34);
    assert_eq!(
        usage_text(46),
        format!(
            "usage: cleanup-blastdb-volumes.py [-h]\n\
             {indent34}-db DB\n\
             {indent34}-dbtype {{prot,nucl}}\n\
             {indent34}[-dry-run]\n\
             {indent34}[-version]\n\
             {indent34}[-verbose]\n"
        )
    );
    // 40 columns wraps usage and options column.
    assert_eq!(
        help_text_for_columns(40),
        "usage: cleanup-blastdb-volumes.py\n\
         \x20      [-h] -db DB -dbtype {prot,nucl}\n\
         \x20      [-dry-run] [-version]\n\
         \x20      [-verbose]\n\
         \n\
         Remove needless BLAST database\n\
         volumes.\n\
         \n\
         options:\n\
         \x20 -h, --help      show this help\n\
         \x20                 message and exit\n\
         \x20 -db DB          BLAST database name\n\
         \x20 -dbtype {prot,nucl}\n\
         \x20                 Molecule type\n\
         \x20 -dry-run        Do not delete any\n\
         \x20                 files, just list\n\
         \x20                 them\n\
         \x20 -version        show program's\n\
         \x20                 version number and\n\
         \x20                 exit\n\
         \x20 -verbose        Increase output\n\
         \x20                 verbosity\n"
    );
    // Extreme narrow width (COLUMNS=1): description/help minimum width is 11.
    assert_eq!(
        help_text_for_columns(1),
        "usage: cleanup-blastdb-volumes.py\n\
         \x20      [-h]\n\
         \x20      -db DB\n\
         \x20      -dbtype {prot,nucl}\n\
         \x20      [-dry-run]\n\
         \x20      [-version]\n\
         \x20      [-verbose]\n\
         \n\
         Remove\n\
         needless\n\
         BLAST\n\
         database\n\
         volumes.\n\
         \n\
         options:\n\
         \x20 -h, --help\n\
         \x20   show this\n\
         \x20   help\n\
         \x20   message and\n\
         \x20   exit\n\
         \x20 -db DB\n\
         \x20   BLAST\n\
         \x20   database\n\
         \x20   name\n\
         \x20 -dbtype {prot,nucl}\n\
         \x20   Molecule\n\
         \x20   type\n\
         \x20 -dry-run\n\
         \x20   Do not\n\
         \x20   delete any\n\
         \x20   files, just\n\
         \x20   list them\n\
         \x20 -version\n\
         \x20   show\n\
         \x20   program's\n\
         \x20   version\n\
         \x20   number and\n\
         \x20   exit\n\
         \x20 -verbose\n\
         \x20   Increase\n\
         \x20   output\n\
         \x20   verbosity\n"
    );
    // Error usage wraps with `columns` while the error message itself is unwrapped.
    assert_eq!(
        parse_with_columns(&os(&["-db", "x"]), 200),
        ParseOutcome::Error(format!(
            "{}{PROG}: error: the following arguments are required: -dbtype\n",
            usage_text(200)
        ))
    );
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
    let usage = usage_text(DEFAULT_COLUMNS);
    let err = |msg: &str| ParseOutcome::Error(format!("{usage}{PROG}: error: {msg}\n"));

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

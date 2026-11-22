use super::*;

#[test]
fn toolkit_variables_are_rejected_and_harmless_entries_accepted() {
    let reject = |name: &str, _value: &str| rejected_variable(OsStr::new(name));
    assert!(reject("DIAG_POST_LEVEL", "Error").is_some());
    assert!(reject("DIAG_OLD_POST_FORMAT", "false").is_some());
    assert!(reject("NCBI_CONFIG__DIAG__POST_FILTER", "!BLAST").is_some());
    assert!(reject("NCBI_CONFIG__BLAST__LONG_SEQID", "1").is_some());
    assert!(reject("NCBI_CONFIG__NCBI__MemorySizeLimit", "1").is_some());
    assert!(reject("NCBI_CONFIG_OVERRIDES", "/tmp/x").is_some());
    assert!(reject("NCBI_CONFIG_LONG_SEQID__BLAST", "1").is_some());
    assert!(reject("NCBI_CONFIG_PATH", "/etc").is_none());
    assert!(reject("ABORT_ON_THROW", "1").is_some());
    // NCBI's Seq-id matching parameters (c++/src/objects/seq/seq_id_tree.cpp:692-707): their own message.
    for name in ["OBJECTS_PACK_TEXTID", "OBJECTS_PACK_GENERAL"] {
        let message = reject(name, "0").unwrap();
        assert!(message.contains("Seq-ids are matched"), "{message}");
        assert!(reject(name, "1").is_some(), "{name}");
    }
    assert!(reject("NCBI_CONFIG__OBJECTS__PACK_TEXTID", "0").is_some());
    assert!(reject("NCBI_CONFIG__OBJECTS__PACK_GENERAL", "0").is_some());
    // `application_start` refuses a `BLAST_USAGE_REPORT` that is not a Boolean.
    assert!(reject("BLAST_USAGE_REPORT", "bogus").is_none());
    assert!(reject("BLAST_USAGE_REPORT", "false").is_none());
    assert!(reject("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "bogus").is_none());
    assert!(reject("NCBI_CONFIG__BLAST__BLASTDB", "/db").is_none());
    assert!(reject("NCBI_CONFIG__blast__blastdb", "/db").is_none());
    assert!(reject("BLASTDB", "/db").is_none());
    assert!(reject("NCBI", "/opt/ncbi").is_none());
    assert!(reject("HOME", "/home/x").is_none());
}

#[test]
fn harmless_entries_are_the_listed_ones() {
    let read = read_registry(
        b"# file comment\n; entry comment\n[BLAST]\nBLASTDB = /db \\\n  /more\n  \n[ncbi]\ndata=/data\n",
        Layer::SystemWide,
    )
    .unwrap();
    assert_eq!(read.error, None);
    assert_eq!(read.entries[0].value, b"/db \n/more");
    assert_eq!(
        check_registry_entries(Path::new("f"), &read.entries),
        Ok(())
    );
    assert!(!is_harmless_entry("BLAST", "LONG_SEQID"));
    assert!(!is_harmless_entry("DEBUG", "DIAG_POST_LEVEL"));
    assert!(!is_harmless_entry("NCBI", ".Inherits"));
    assert!(!is_harmless_entry("OBJECTS", "PACK_TEXTID"));
    assert!(!is_harmless_entry("OBJECTS", "PACK_GENERAL"));
    assert!(is_harmless_entry("blast", "blast_usage_report"));
    let other = read_registry(b"[BLAST]\nLONG_SEQID=1\n", Layer::SystemWide).unwrap();
    let refusal = check_registry_entries(Path::new("f"), &other.entries).unwrap_err();
    assert!(
        refusal.contains("the registry file f sets [BLAST] LONG_SEQID"),
        "{refusal}"
    );
    let packed = read_registry(b"[OBJECTS]\nPACK_TEXTID=0\n", Layer::SystemWide).unwrap();
    let refusal = check_registry_entries(Path::new("f"), &packed.entries).unwrap_err();
    assert!(
        refusal.contains("the registry file f sets [OBJECTS] PACK_TEXTID"),
        "{refusal}"
    );
}

/// No NCBI counterpart: a process with no user database entry and no program directory,
/// for the tests whose search path does not reach them.
fn no_home() -> Option<OsString> {
    None
}
fn no_dirs(_: bool) -> Result<Vec<PathBuf>, String> {
    Ok(Vec::new())
}
const NO_PROCESS: ProcessFacts<'static> = ProcessFacts {
    home_by_uid: &no_home,
    program_dirs: &no_dirs,
    files: &FileSystem,
};

// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:54-56 (`kConfigPathDelim` is ";" on
// Windows and ":;" elsewhere) and :331-396 (`CMetaRegistry::GetDefaultSearchPath`, cited at
// `registry_search_path`); the home directory comes from `CDir::GetHome` (`ncbi_home`). The
// test uses the delimiter and the variable of the platform it runs on, so it checks the same
// rule everywhere.
#[test]
fn the_search_path_follows_ncbi_config_path() {
    let env = |pairs: &'static [(&'static str, &'static str)]| {
        move |name: &str| {
            pairs
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| OsString::from(value))
        }
    };
    // The program's directories after the arguments are set, and before.
    let dirs = |arguments_set: bool| {
        Ok(match arguments_set {
            true => vec![PathBuf::from("/p/"), PathBuf::from("/q/")],
            false => vec![PathBuf::from("/q/")],
        })
    };
    let process = ProcessFacts {
        home_by_uid: &no_home,
        program_dirs: &dirs,
        files: &FileSystem,
    };
    let search = |pairs| registry_search_path("blastn", &env(pairs), &process).unwrap();
    let only = search(&[("NCBI_CONFIG_PATH", "/a;/b")]);
    assert_eq!(only, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
    // ':' splits the path outside Windows only (metareg.cpp:54-56).
    let colon = search(&[("NCBI_CONFIG_PATH", "/a:/b")]);
    if cfg!(windows) {
        assert_eq!(colon, vec![PathBuf::from("/a:/b")]);
    } else {
        assert_eq!(colon, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
    }
    // An empty NCBI_CONFIG_PATH splits into no directory (ncbistr_util.hpp:271-272): no
    // empty element, so nothing is searched.
    assert_eq!(search(&[("NCBI_CONFIG_PATH", "")]), Vec::<PathBuf>::new());
    // No NCBI counterpart: test plumbing that picks the platform's empty-entry path and home
    // variable (the slices must be constants to be 'static).
    let spliced = search(if cfg!(windows) {
        &[("NCBI_CONFIG_PATH", "/a;;/z"), ("APPDATA", "/h")]
    } else {
        &[("NCBI_CONFIG_PATH", "/a::/z"), ("HOME", "/h")]
    });
    assert_eq!(
        &spliced[..3],
        &[PathBuf::from("/a"), PathBuf::from("."), PathBuf::from("/h")]
    );
    assert_eq!(
        &spliced[spliced.len() - 2..],
        &[PathBuf::from("/q/"), PathBuf::from("/z")]
    );
    let no_local = search(&[("NCBI_DONT_USE_LOCAL_CONFIG", ""), ("NCBI", "/n")]);
    assert_eq!(&no_local[..1], &[PathBuf::from("/n")]);
    assert_eq!(no_local.last(), Some(&PathBuf::from("/q/")));
    // A false BLAST_USAGE_REPORT, or .ncbirc turned off in the environment (an empty value
    // too), makes the path after the arguments are set (blast_usage_report.cpp:195-209,
    // ncbireg.cpp:1634-1648); a true BLAST_USAGE_REPORT, before.
    for (variable, value, set) in [
        ("BLAST_USAGE_REPORT", "false", true),
        ("BLAST_USAGE_REPORT", "0", true),
        ("BLAST_USAGE_REPORT", "on", false),
        ("NCBI_DONT_USE_NCBIRC", "", true),
        ("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC", "", true),
    ] {
        let env = |name: &str| match name {
            "HOME" => Some(OsString::from("/h")),
            _ if name == variable => Some(OsString::from(value)),
            _ => None,
        };
        let path = registry_search_path("blastn", &env, &process).unwrap();
        let last = dirs(set).unwrap();
        assert_eq!(
            &path[path.len() - last.len()..],
            &last[..],
            "{variable}={value}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3657 (`CDir::GetHome` and
// `s_GetHomeByUID`, cited at `ncbi_home`).
#[test]
fn the_home_directory_is_ncbis() {
    let env = |pairs: &'static [(&'static str, &'static str)]| {
        move |name: &str| {
            pairs
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| OsString::from(value))
        }
    };
    let uid_home = || Some(OsString::from("/pw"));
    let home = |pairs, by_uid: &dyn Fn() -> Option<OsString>| ncbi_home(&env(pairs), by_uid);
    if cfg!(windows) {
        // APPDATA when it is set, an empty one too; then USERPROFILE.
        let both = &[("APPDATA", "/app"), ("USERPROFILE", "/user")];
        assert_eq!(home(both, &no_home), Ok(OsString::from("/app")));
        let empty = &[("APPDATA", ""), ("USERPROFILE", "/user")];
        assert_eq!(home(empty, &no_home), Ok(OsString::new()));
        assert_eq!(
            home(&[("USERPROFILE", "/user")], &no_home),
            Ok(OsString::from("/user"))
        );
        assert_eq!(home(&[("HOME", "/h")], &no_home), Ok(OsString::new()));
    } else if cfg!(unix) {
        // HOME when it is set, an empty one too; then the user ID's entry, an empty
        // directory too; without an entry the login name's, which LOSAT refuses.
        assert_eq!(home(&[("HOME", "/h")], &uid_home), Ok(OsString::from("/h")));
        assert_eq!(home(&[("HOME", "")], &uid_home), Ok(OsString::new()));
        assert_eq!(home(&[], &uid_home), Ok(OsString::from("/pw")));
        assert_eq!(home(&[], &|| Some(OsString::new())), Ok(OsString::new()));
        assert!(home(&[("USER", "u")], &no_home)
            .unwrap_err()
            .contains("by the login name"));
        assert_eq!(
            home(&[("APPDATA", "/app")], &uid_home),
            Ok(OsString::from("/pw"))
        );
    }
    // The search path leaves out an empty home (metareg.cpp:349-352).
    let process = ProcessFacts {
        home_by_uid: &|| Some(OsString::new()),
        program_dirs: &no_dirs,
        files: &FileSystem,
    };
    let path = registry_search_path("blastn", &env(&[]), &process).unwrap();
    assert_eq!(path[0], PathBuf::from("."));
    assert!(!path.contains(&PathBuf::new()), "{path:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1381
// (`BlastFindMatrixPath`) and c++/src/util/util_misc.cpp:139-179 (`g_FindDataFile`), quoted
// at `ApplicationSettings::matrix_file`.
#[test]
fn matrix_files_are_found_where_ncbi_looks() {
    let sep = if cfg!(windows) { "\\" } else { "/" };
    let delimiter = if cfg!(windows) { ";" } else { ":" };
    let path = |parts: &[&str]| PathBuf::from(parts.join(sep));
    let files = [
        path(&["dp2", "MYMAT"]),
        path(&["data_dir", "MYMAT"]),
        path(&["bm", "MYMAT"]),
        path(&["bm", "nt", "NTMAT"]),
        path(&["data", "LOCAL"]),
        path(&["bm", "mixed"]),
    ];
    let dirs = [PathBuf::from("bm")];
    let is_file = |p: &Path| files.iter().any(|file| file == p);
    let is_dir = |p: &Path| dirs.iter().any(|dir| dir == p);
    let search =
        |data_path: Option<&str>, data: Option<&str>, blastmat: Option<&str>| MatrixSearch {
            data_path: data_path.map(OsString::from),
            data: Ok(data.map(OsString::from)),
            blastmat: blastmat.map(OsString::from),
        };
    let find = |name: &str, search: &MatrixSearch| {
        find_matrix_file(name, search, &is_file, &is_dir).unwrap()
    };
    let none = search(None, None, None);
    // NCBI_DATA_PATH's directories (empty parts dropped), then [NCBI] Data, each with a
    // separator added, for the upper-cased name.
    let data_path = format!("{delimiter}dp1{delimiter}{delimiter}dp2{sep}");
    assert_eq!(
        find("mymat", &search(Some(&data_path), None, None)),
        Some(path(&["dp2", "MYMAT"]))
    );
    assert_eq!(
        find("MyMat", &search(None, Some("data_dir"), Some("bm"))),
        Some(path(&["data_dir", "MYMAT"]))
    );
    // BLASTMAT, a directory: its file, then its "nt" subdirectory's.
    assert_eq!(
        find("mymat", &search(None, None, Some("bm"))),
        Some(path(&["bm", "MYMAT"]))
    );
    assert_eq!(
        find("ntmat", &search(None, None, Some("bm"))),
        Some(path(&["bm", "nt", "NTMAT"]))
    );
    assert_eq!(find("mymat", &search(None, None, Some("nobm"))), None);
    // The name is upper-cased before the search (blast_setup.c:381-385).
    assert_eq!(find("mixed", &search(None, None, Some("bm"))), None);
    // The working directory's "data".
    assert_eq!(find("local", &none), Some(path(&["data", "LOCAL"])));
    assert_eq!(find("mymat", &none), None);
    // An empty [NCBI] Data is no directory; an unreadable registry value is the reason.
    assert_eq!(find("mymat", &search(None, Some(""), None)), None);
    let unreadable = MatrixSearch {
        data: Err("reason".to_string()),
        ..none.clone()
    };
    assert_eq!(
        find_matrix_file("mymat", &unreadable, &is_file, &is_dir),
        Err("reason".to_string())
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-971 (`NormalizePath`, quoted
// at `normalize_path`).
#[test]
fn paths_are_normalized_like_ncbis() {
    let lexical =
        |path: &str| String::from_utf8(normalize_path(path.as_bytes(), None).unwrap()).unwrap();
    for (path, normal) in [
        ("", ""),
        ("/", "/"),
        ("///", "/"),
        ("/usr/bin/", "/usr/bin"),
        ("/usr//./bin", "/usr/bin"),
        ("/../a", "/a"),
        ("/a/b/../../..", "/"),
        ("/a/b/../c", "/a/c"),
        ("a/../..", ".."),
        ("../../a", "../../a"),
        ("./a", "a"),
        ("./", "."),
        (".", "."),
        ("a/.", "a"),
        ("a/..", "."),
        ("./..", ".."),
        ("bin/blastn", "bin/blastn"),
    ] {
        assert_eq!(lexical(path), normal, "{path:?}");
    }
    // With links: /l -> /real/dir, /real/dir/x -> ../y (relative to its directory).
    let read_link = |name: &[u8]| match name {
        b"/l" => Some(b"/real/dir".to_vec()),
        b"/real/dir/x" => Some(b"../y".to_vec()),
        _ => None,
    };
    let follow = |path: &str| {
        String::from_utf8(normalize_path(path.as_bytes(), Some(&read_link)).unwrap()).unwrap()
    };
    assert_eq!(follow("/l/z"), "/real/dir/z");
    assert_eq!(follow("/l/x"), "/real/y");
    assert_eq!(follow("/l/../z"), "/real/z");
    // A loop of links: NCBI stops at 1024 with a warning.
    let looping = |name: &[u8]| (name == b"ncbi").then(|| b"ncbi".to_vec());
    assert_eq!(normalize_path(b"ncbi", Some(&looping)), None);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1452-1457, 1549-1595
// (`FindProgramExecutablePath`, cited at `unix_program_name`) and metareg.cpp:378-385,
// ncbienv.cpp:408-415 (`GetProgramDirname`, cited at `program_directories`).
#[test]
fn the_program_directories_are_ncbis() {
    let files = [
        "/cwd/bin/blastn",
        "/cwd/blastn",
        "/home/u/bin/blastn",
        "/opt/l/losat",
        "/cwd/ncbi",
    ];
    let is_file = |name: &[u8]| files.iter().any(|file| file.as_bytes() == name);
    // The working directory is /cwd: a relative name is tried there first.
    let in_cwd = |name: &[u8]| {
        let full = if name.starts_with(b"/") {
            name.to_vec()
        } else {
            [b"/cwd/".as_slice(), name].concat()
        };
        is_file(&normalize_path(&full, None).unwrap())
    };
    let name = |argv0: Option<&str>, path: &str| {
        String::from_utf8(unix_program_name(
            argv0.map(str::as_bytes),
            path.as_bytes(),
            b"/cwd",
            &in_cwd,
        ))
        .unwrap()
    };
    // An absolute name stays (normalized).
    assert_eq!(
        name(Some("/home/u/bin/../bin/blastn"), ""),
        "/home/u/bin/blastn"
    );
    // A relative name of a file in the working directory.
    assert_eq!(name(Some("./bin/blastn"), ""), "/cwd/bin/blastn");
    assert_eq!(name(Some("blastn"), "/home/u/bin"), "/cwd/blastn");
    // Else the first directory of PATH (empty parts dropped) with a file of the base name.
    assert_eq!(
        name(Some("bin2/blastn.x"), "::/nowhere::/home/u/bin:"),
        "/home/u/bin/blastn"
    );
    assert_eq!(name(Some("losat"), "/x:/opt/l/"), "/opt/l/losat");
    // Not found: argv[0] as it is.
    assert_eq!(name(Some("tblastn"), "/home/u/bin"), "tblastn");
    // No argv[0] or an empty one: the name "ncbi" is looked for, else none.
    assert_eq!(name(None, ""), "/cwd/ncbi");
    assert_eq!(name(Some(""), ""), "/cwd/ncbi");
    let none = |name: &[u8]| name == b"x";
    assert_eq!(unix_program_name(Some(b""), b"", b"/cwd", &none), b"");
    assert_eq!(unix_program_name(None, b"/bin", b"/cwd", &none), b"");
    // getcwd failed: "" + "/" + the name, and the name as it is when that is no file.
    let rooted = |name: &[u8]| name == b"x" || name == b"/x";
    assert_eq!(unix_program_name(Some(b"x"), b"", b"", &rooted), b"/x");
    assert_eq!(unix_program_name(Some(b"x"), b"", b"", &none), b"x");
    // A working directory of 4095 bytes or more is read as none (`getcwd_reads`), which is
    // what `program_directories` takes (line WQA's D11 round 4, 2026-11-20).
    let deep = |length: usize| {
        getcwd_reads(OsString::from(format!("/{}", "d".repeat(length - 1))))
            .map(|cwd| cwd.as_encoded_bytes().to_vec())
            .unwrap_or_default()
    };
    assert_eq!(
        unix_program_name(Some(b"x"), b"", &deep(4095), &rooted),
        b"/x"
    );
    let in_deep = |name: &[u8]| name == b"x" || name.ends_with(b"d/x");
    assert_eq!(
        unix_program_name(Some(b"x"), b"", &deep(4094), &in_deep),
        [deep(4094).as_slice(), b"/x"].concat()
    );
    // The directories: the name's, then the resolved name's when it is another one (the
    // last of '/', '\\' and ':' ends a directory on every platform).
    let dirs = |name: &str, resolved: Option<&str>| {
        program_dirs_of(name.as_bytes(), resolved.map(str::as_bytes))
            .into_iter()
            .map(|dir| String::from_utf8(dir).unwrap())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        dirs("/home/u/bin/blastn", Some("/opt/l/losat")),
        ["/home/u/bin/", "/opt/l/"]
    );
    assert_eq!(dirs("/opt/l/losat", Some("/opt/l/losat")), ["/opt/l/"]);
    assert_eq!(dirs("tblastn", Some("/opt/l/losat")), ["/opt/l/"]);
    assert_eq!(dirs("", None), Vec::<String>::new());
    assert_eq!(dirs("/a/b:c", None), ["/a/b:"]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
// (case-insensitive substrings `blastdb`, `genbank` and `none`) and
// c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
// `.ini`, then `.ncbirc`); an empty value in a registry file counts as no entry
// (ncbireg.cpp:984-991), while a variable set to an empty value is an entry
// (env_reg.cpp:157-167; audit finding B-1 of session SFc).
#[test]
fn data_loaders_follow_the_registry_layers_and_the_substring_rules() {
    let only = |value: &str| data_loaders_of(None, [None, Some(value.as_bytes())]);
    for (value, used) in [
        ("blastdb", true),
        ("genbank", true),
        ("BlastDB", true),
        ("GENBANK", true),
        ("blastdb,genbank", true),
        ("genbank blastdb", true),
        ("myblastdbs", true),
        ("none", false),
        ("NONE", false),
        ("blastdb,none", false),
        ("genbank none blastdb", false),
        ("nonexistent", false),
        ("x", false),
        ("blast db", false),
        ("gen bank", false),
    ] {
        assert_eq!(only(value), used, "{value:?}");
    }
    let layers = |env: Option<&str>, ini: Option<&str>, ncbirc: Option<&str>| {
        data_loaders_of(
            env.map(str::as_bytes),
            [ini.map(str::as_bytes), ncbirc.map(str::as_bytes)],
        )
    };
    assert!(layers(None, None, None));
    assert!(!layers(Some("none"), Some("blastdb"), Some("genbank")));
    assert!(layers(Some("genbank"), Some("none"), Some("none")));
    assert!(!layers(None, Some("none"), Some("blastdb")));
    assert!(layers(None, Some("blastdb"), Some("none")));
    assert!(!layers(None, None, Some("none")));
    // An empty value in a registry file is no entry: the next layer decides, and with
    // none left the data loaders stay on.
    assert!(!layers(None, Some(""), Some("x")));
    assert!(layers(None, Some(""), Some("genbank")));
    assert!(layers(None, None, Some("")));
    assert!(layers(None, Some(""), Some("")));
    // An empty variable is an entry: it answers first and turns both loaders off, also
    // over a registry file that turns them on.
    assert!(!layers(Some(""), None, None));
    assert!(!layers(Some(""), Some("blastdb"), Some("genbank")));
    assert!(!layers(Some(""), None, Some("genbank")));
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1636-1651 (`.ncbirc` is read
// unless NCBI_DONT_USE_NCBIRC is set or the registry has [NCBI] DONT_USE_NCBIRC, also
// from the environment) and c++/src/corelib/env_reg.cpp:366-375 (the variable of
// [BLAST] DATA_LOADERS is NCBI_CONFIG__BLAST__DATA_LOADERS).
#[test]
fn application_settings_read_the_environment_and_the_registry_files() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-data-loaders-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let ini = dir.join("blastn.ini");
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let settings = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_settings(
            "blastn",
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
        .map(|settings| settings.data_loaders)
    };
    let _ = std::fs::remove_file(&ini);
    let _ = std::fs::remove_file(&ncbirc);
    assert_eq!(settings(&[]), Ok(true));
    std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
    assert_eq!(settings(&[]), Ok(false));
    // The environment comes first; an empty variable is an entry that turns both
    // loaders off (env_reg.cpp:157-167).
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "blastdb")]),
        Ok(true)
    );
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
        Ok(false)
    );
    std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = genbank\n").unwrap();
    assert_eq!(settings(&[]), Ok(true));
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
        Ok(false)
    );
    std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS =\n").unwrap();
    assert_eq!(settings(&[]), Ok(true));
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
        Ok(false)
    );
    std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
    // Only the exact variable name is the entry.
    assert_eq!(
        settings(&[("NCBI_CONFIG__blast__data_loaders", "blastdb")]),
        Ok(false)
    );
    // The program's .ini comes before .ncbirc; its last entry counts, without quotes.
    std::fs::write(
        &ini,
        "[blast]\ndata_loaders = none\n[BLAST]\nDATA_LOADERS = \"GenBank\"\n",
    )
    .unwrap();
    assert_eq!(settings(&[]), Ok(true));
    std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = \"\"\n").unwrap();
    assert_eq!(settings(&[]), Ok(false));
    // A value NCBI reports as malformed ends the program (exit status 2): the start has no
    // settings, and its caller writes NCBI's report before it reaches them
    // (`ApplicationStart::write_config_load`; the library's calls too, `api::start`, line RG
    // item E1w-2).
    std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = no\"ne\n").unwrap();
    assert!(settings(&[]).unwrap_err().contains("exit status 2"));
    // The escapes are read: `\x6e` is `n`.
    std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = \\x6eone\n").unwrap();
    assert_eq!(settings(&[]), Ok(false));
    // [NCBI] DONT_USE_NCBIRC in the .ini or the environment turns .ncbirc off.
    std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
    assert_eq!(settings(&[]), Ok(true));
    std::fs::remove_file(&ini).unwrap();
    assert_eq!(
        settings(&[("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC", "")]),
        Ok(true)
    );
    assert_eq!(settings(&[("NCBI_DONT_USE_NCBIRC", "1")]), Ok(true));
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::remove_dir(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:739-747 (the
// registry's `[BLAST] MAX_SEQID_LENGTH`, used when it is not empty): the value of the first
// layer that has the entry (the environment, the program's .ini, .ncbirc; an empty value
// hides the layers below, as makeblastdb 2.17.0 reads them: pr2/fix/n4probe), for
// makeblastdb to refuse with -parse_seqids when NCBI acts on it; the entry stays harmless
// for the other programs.
#[test]
fn max_seqid_length_is_reported_from_every_layer() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-max-seqid-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let ini = dir.join("makeblastdb.ini");
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let _ = std::fs::remove_file(&ini);
    let _ = std::fs::remove_file(&ncbirc);
    let settings = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_settings(
            "makeblastdb",
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
        .map(|settings| settings.max_seqid_length)
    };
    let value = |v: &str| Ok(Ok(Some(v.to_string())));
    assert_eq!(settings(&[]), Ok(Ok(None)));
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "60")]),
        value("60")
    );
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "")]),
        value("")
    );
    std::fs::write(&ncbirc, "[BLAST]\nMAX_SEQID_LENGTH = 70\n").unwrap();
    assert_eq!(settings(&[]), value("70"));
    // An empty environment value hides the files' values.
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "")]),
        value("")
    );
    std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = 80\n").unwrap();
    assert_eq!(settings(&[]), value("80"));
    std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH =\n").unwrap();
    assert_eq!(settings(&[]), value(""));
    std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = \\x37\\x30\n").unwrap();
    assert_eq!(settings(&[]), value("70"));
    std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = 1\\xff\n").unwrap();
    assert!(matches!(settings(&[]), Ok(Err(_))));
    std::fs::write(&ini, "[BLAST]\nBLASTDB = /tmp\n").unwrap();
    assert_eq!(settings(&[]), value("70"));
    std::fs::write(&ncbirc, "[BLAST]\nLONG_SEQID = 1\n").unwrap();
    assert!(settings(&[]).is_err());
    std::fs::remove_file(&ini).unwrap();
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::remove_dir(&dir).unwrap();
}

/// A field of `registry.tsv`: backslash escapes, `\s` for a space at either end and
/// `\e` for an empty field.
fn unescape(field: &str) -> String {
    let mut out = String::new();
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some('s') => out.push(' '),
            Some('e') => {}
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:527-533
// (`registry.HasEntry` for section "BLAST", entry "BLASTDB", and `registry.Get`, cited at `blastdb_of`) and
// c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
// `.ini`, then `.ncbirc`).
/// `tests/seqdb/expected/search_path/registry.tsv`: the value of `[BLAST] BLASTDB` that
/// NCBI's blastp put into the database search path (it prints the path when a database
/// is not found) for a variable, a `blastp.ini` and a `.ncbirc`; `{R}` stands for the
/// directory of the run, `<unset>` for a variable or file that is absent and `<none>`
/// for a registry without the entry.
#[test]
fn the_registry_value_of_blastdb_is_ncbis() {
    let table =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/expected/search_path/registry.tsv");
    let text = std::fs::read_to_string(&table).unwrap();
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-blastdb-{}",
        crate::utils::test_dir::process_id()
    ));
    let place = |value: &str| value.replace("{R}", "/sandbox");
    let mut cases = 0;
    for line in text.lines().skip(1) {
        let row: Vec<String> = line.split('\t').map(unescape).collect();
        assert_eq!(row.len(), 8, "{line:?}");
        let (id, variable, value, ini, ncbirc, registry) =
            (&row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
        std::fs::create_dir_all(&dir).unwrap();
        for (name, content) in [("blastp.ini", ini), (SYSTEM_WIDE_REGISTRY, ncbirc)] {
            let file = dir.join(name);
            let _ = std::fs::remove_file(&file);
            if content != "<unset>" {
                std::fs::write(&file, place(content)).unwrap();
            }
        }
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        if value != "<unset>" {
            all.push((variable.clone(), place(value)));
        }
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        let settings = application_settings(
            "blastp",
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
        .unwrap_or_else(|reason| panic!("{id}: {reason}"));
        let expected = (registry != "<none>").then(|| place(registry));
        assert_eq!(settings.blastdb, Ok(expected), "{id}");
        cases += 1;
    }
    assert!(cases >= 15, "{cases} cases");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `tests/winmask/expected/taxid_path.tsv`: the directory whose `7227/wmasker.obinary`
/// NCBI's blastn `-window_masker_taxid 7227` used (told apart by the masks of three
/// different counts files) for the two variables, another variable, a `blastn.ini` and
/// a `.ncbirc` in the current directory and in the home directory; `{R}` stands for the
/// directory of the run and `<cwd>` for the current directory. The search path here is
/// the current directory, then the home directory, as NCBI's (`registry_search_path`).
#[test]
fn the_window_masker_directory_is_ncbis() {
    let table = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/winmask/expected/taxid_path.tsv");
    let text = std::fs::read_to_string(&table).unwrap();
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-window-masker-{}",
        crate::utils::test_dir::process_id()
    ));
    let (cwd, home) = (dir.join("cwd"), dir.join("home"));
    let place = |value: &str| value.replace("{R}", "/sandbox");
    let mut cases = 0;
    for line in text.lines().skip(1) {
        let row: Vec<String> = line.split('\t').map(unescape).collect();
        assert_eq!(row.len(), 9, "{line:?}");
        let id = &row[0];
        std::fs::create_dir_all(&cwd).unwrap();
        std::fs::create_dir_all(&home).unwrap();
        for (file, content) in [
            (cwd.join("blastn.ini"), &row[4]),
            (cwd.join(SYSTEM_WIDE_REGISTRY), &row[5]),
            (home.join(SYSTEM_WIDE_REGISTRY), &row[6]),
        ] {
            let _ = std::fs::remove_file(&file);
            if content != "<unset>" {
                std::fs::write(&file, place(content)).unwrap();
            }
        }
        let mut all = vec![(
            "NCBI_CONFIG_PATH".to_string(),
            format!("{};{}", cwd.display(), home.display()),
        )];
        for (name, value) in [
            ("NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH", &row[1]),
            ("WINDOW_MASKER_PATH", &row[2]),
        ] {
            if value != "<unset>" {
                all.push((name.to_string(), place(value)));
            }
        }
        if row[3] != "<unset>" {
            let (name, value) = row[3].split_once('=').unwrap();
            all.push((name.to_string(), place(value)));
        }
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        let settings = application_settings(
            "blastn",
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
        .unwrap_or_else(|reason| panic!("{id}: {reason}"));
        let used = settings
            .window_masker_registry
            .unwrap_or_else(|reason| panic!("{id}: {reason}"))
            .path
            .filter(|value| !value.is_empty());
        let expected = (row[7] != "<cwd>").then(|| place(&row[7]));
        assert_eq!(used, expected, "{id}");
        cases += 1;
    }
    assert_eq!(cases, 20);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// The entry is accepted in every program's registry files, read by NCBI's reader
/// (`read_registry_file`). A `.ncbirc` that NCBI reports as malformed (a syntax error, a bad
/// value or name in any entry) gets the application start's messages and gives the entries
/// above the line, to the registry of `-window_masker_taxid` too, which posts the error again
/// at each load with the program's prefix, also when the program's `.ini` turned `.ncbirc`
/// off for the application and when a variable gives the directory (NCBI reads the file
/// before it asks for the entry). A one-byte file whose stream fails is posted by that
/// registry's first load only where no read with its flags came before (NCBI BLAST+ 2.17.0 on
/// 2026-11-12, line NP M3c). `[NCBI] .Inherits` and a value that is not UTF-8 are refused for
/// `-window_masker_taxid`. A malformed file is posted and the start goes on, in LOSAT's library
/// too (`api::start`, line RG item E1w-2).
#[test]
fn window_masker_registry_values_are_read_or_refused_as_ncbi_stores_them() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-window-masker-kept-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let config_path = dir.display().to_string();
    let start_of = |program: &str, extra: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH", config_path.as_str())];
        all.extend_from_slice(extra);
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_start(
            program,
            None,
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
    };
    let registry_of = |program: &str, extra: &[(&str, &str)]| {
        start_of(program, extra)
            .settings
            .unwrap()
            .window_masker_registry
    };
    let registry = |extra: &[(&str, &str)]| registry_of("blastn", extra).unwrap();
    let path_of = |path: &str| WindowMaskerRegistry {
        path: Some(path.to_string()),
        ..WindowMaskerRegistry::default()
    };
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    std::fs::write(&ncbirc, "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n").unwrap();
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        assert_eq!(registry_of(program, &[]), Ok(path_of("/w")), "{program}");
    }
    // NCBI's reader takes the escapes and the quotes at the ends.
    std::fs::write(
        &ncbirc,
        "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = \"/w\\\\x\"\n",
    )
    .unwrap();
    assert_eq!(registry(&[]), path_of("/w\\x"));
    // A malformed line: the application goes on with the entries above it, with NCBI's
    // message at the start, and so does the registry of `-window_masker_taxid`, which posts
    // it again at each load, with the program's prefix.
    std::fs::write(
        &ncbirc,
        "[BLAST]\nBLASTDB = /dbs\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\nWINDOW_MASKER_PATH = a\"b\n",
    )
    .unwrap();
    for program in ["blastn", "blastp"] {
        let start = start_of(program, &[]);
        let at_start = start.usage_report().unwrap().to_string();
        assert!(
            at_start.starts_with("Critical: (110.6) CNcbiRegistry: "),
            "{program}"
        );
        let settings = start.settings.unwrap();
        assert_eq!(settings.blastdb, Ok(Some("/dbs".to_string())), "{program}");
        let registry = settings.window_masker_registry.unwrap();
        assert_eq!(registry.path.as_deref(), Some("/w"), "{program}");
        let mut posts = Posts::new();
        posts.post(
            CRITICAL,
            at_start.replacen("(110.6) ", &format!("(110.6) [{program}] "), 1),
        );
        assert_eq!(registry.load_posts, posts, "{program}");
        assert!(!registry.posts_once, "{program}");
    }
    assert!(application_settings(
        "blastn",
        [(
            OsString::from("NCBI_CONFIG_PATH"),
            OsString::from(&config_path)
        )]
        .into_iter(),
        &|name: &str| (name == "NCBI_CONFIG_PATH").then(|| OsString::from(&config_path)),
        &NO_PROCESS,
    )
    .is_ok());
    // A one-byte file whose stream fails, read before with the registry's flags (by the
    // usage report, or by the application without a `.ini`): its loads post nothing.
    std::fs::write(&ncbirc, b"\xef").unwrap();
    assert_eq!(registry(&[]), WindowMaskerRegistry::default());
    assert_eq!(
        registry(&[("BLAST_USAGE_REPORT", "0")]),
        WindowMaskerRegistry::default()
    );
    let ini = dir.join("blastn.ini");
    std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
    assert_eq!(registry(&[]), WindowMaskerRegistry::default());
    // Without such a read (a `.ini`, and no usage report): its first load posts it.
    let mut posts = Posts::new();
    posts.post(
        ERROR,
        "Error: (110.4) [blastn] Error reading the registry after line 1: \n",
    );
    assert_eq!(
        registry(&[("BLAST_USAGE_REPORT", "0")]),
        WindowMaskerRegistry {
            path: None,
            load_posts: posts,
            posts_once: true,
        }
    );
    // An entry set again while its value is not empty: a warning, kept as the stream that
    // failed is (NCBI BLAST+ 2.17.0 on 2026-11-12, line NP M3c review `d1`, `d2`).
    let twice = "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /u\nwindow_masker_path = /w\n";
    std::fs::write(&ncbirc, twice).unwrap();
    let mut posts = Posts::new();
    posts.post(
        WARNING,
        "Warning: (110.8) [blastn] Found multiple [WINDOW_MASKER] window_masker_path settings; using the one from line 3\n",
    );
    assert_eq!(
        registry(&[("BLAST_USAGE_REPORT", "0")]),
        WindowMaskerRegistry {
            path: Some("/w".to_string()),
            load_posts: posts.clone(),
            posts_once: true,
        }
    );
    assert_eq!(registry(&[]), path_of("/w"));
    // Above a malformed line: the warning, then the error, at each load (`d3`).
    std::fs::write(&ncbirc, format!("{twice}[my blast]\n")).unwrap();
    for extra in [&[][..], &[("BLAST_USAGE_REPORT", "0")][..]] {
        let registry = registry(extra);
        assert_eq!(registry.path.as_deref(), Some("/w"));
        assert!(!registry.posts_once);
        let mut expected = posts.clone();
        expected.post(
            CRITICAL,
            "Critical: (110.6) [blastn] CNcbiRegistry: Syntax error in system-wide configuration file: ",
        );
        assert!(
            registry.load_posts.text().starts_with(&expected.text()),
            "{}",
            String::from_utf8_lossy(&registry.load_posts.text())
        );
    }
    // An earlier entry with the empty value is no entry to NCBI's `HasEntry`.
    std::fs::write(
        &ncbirc,
        "[WINDOW_MASKER]\nWINDOW_MASKER_PATH =\nWINDOW_MASKER_PATH = /w\n",
    )
    .unwrap();
    assert_eq!(
        registry(&[("BLAST_USAGE_REPORT", "0")]).load_posts,
        Posts::new()
    );
    std::fs::write(
        &ncbirc,
        "[BLAST]\nLONG_SEQID = 1\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n",
    )
    .unwrap();
    assert_eq!(registry(&[]), path_of("/w"));
    assert_eq!(registry(&[("WINDOW_MASKER_PATH", "/v")]), path_of("/v"));
    // An entry before the first section is dropped by NCBI without an error, but its name
    // is checked (NCBI BLAST+ 2.17.0 on 2026-11-12: `x y = 1` there is a syntax error).
    std::fs::write(
        &ncbirc,
        "xy = 1\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n",
    )
    .unwrap();
    assert_eq!(registry(&[]), path_of("/w"));
    // Malformed files, with the `.ini` that turns `.ncbirc` off for the application: the
    // entries above the line, the posts at each load (also under a variable).
    let wm = "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n";
    for (text, path) in [
        (
            "[WINDOW_MASKER\nWINDOW_MASKER_PATH = /w\n".to_string(),
            None,
        ),
        (
            "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = a\"b\n".to_string(),
            None,
        ),
        (format!("[BLAST]\nLONG_SEQID = a\"b\n{wm}"), None),
        (format!("[my blast]\nx = 1\n{wm}"), None),
        (format!("[BLAST]\nlong seqid = 1\n{wm}"), None),
        (format!("x y = 1\n{wm}"), None),
        (format!("{wm}[my blast]\n"), Some("/w")),
    ] {
        std::fs::write(&ncbirc, &text).unwrap();
        for (extra, path) in [
            (&[][..], path),
            (&[("WINDOW_MASKER_PATH", "/v")][..], Some("/v")),
        ] {
            let registry = registry(extra);
            assert_eq!(registry.path.as_deref(), path, "{text:?} {extra:?}");
            assert!(
                registry.load_posts.text().starts_with(
                    // NCBI's `... file: NCBI C++ Exception:\n    T0 ...` in LOSAT's words
                    // (`crate::wording::losat_text`) goes on with the message at once.
                    b"Critical: (110.6) [blastn] CNcbiRegistry: Syntax error in system-wide configuration file: "
                ),
                "{text:?} {extra:?}: {}",
                String::from_utf8_lossy(&registry.load_posts.text())
            );
            assert!(!registry.posts_once);
        }
    }
    std::fs::write(&ncbirc, "[NCBI]\n.Inherits = other\n").unwrap();
    for extra in [&[][..], &[("WINDOW_MASKER_PATH", "/v")][..]] {
        let reason = registry_of("blastn", extra).unwrap_err();
        assert!(reason.contains(".Inherits"), "{reason}");
    }
    // A byte that is not UTF-8 in the value.
    std::fs::write(&ncbirc, b"[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\xff\n").unwrap();
    let reason = registry_of("blastn", &[]).unwrap_err();
    assert!(reason.contains("UTF-8"), "{reason}");
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781 (cited at
// `ncbi_registry::read_lines`): the escapes of a value are read, and a value that is not
// UTF-8 is not a path LOSAT reads.
/// A value of `[BLAST] BLASTDB` that LOSAT cannot read is no reason to reject a program
/// that opens no database: it is kept as the reason for the program that does.
#[test]
fn an_unreadable_blastdb_value_is_kept_for_the_program_that_opens_a_database() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-blastdb-unreadable-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let settings = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_settings(
            "blastp",
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
        .unwrap()
    };
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let ini = dir.join("blastp.ini");
    let _ = std::fs::remove_file(&ini);
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB = C:\\db\n").unwrap();
    assert_eq!(settings(&[]).blastdb, Ok(Some("C:db".to_string())));
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB = /d\\xffb\n").unwrap();
    let unreadable = settings(&[]);
    assert!(unreadable.data_loaders);
    let reason = unreadable.blastdb.unwrap_err();
    assert!(reason.contains("[BLAST] BLASTDB") && reason.contains("UTF-8"));
    // A layer above the file answers first.
    assert_eq!(
        settings(&[("NCBI_CONFIG__BLAST__BLASTDB", "/db")]).blastdb,
        Ok(Some("/db".to_string()))
    );
    std::fs::write(&ini, "[BLAST]\nBLASTDB = /ini\n").unwrap();
    assert_eq!(settings(&[]).blastdb, Ok(Some("/ini".to_string())));
    std::fs::write(&ini, "[BLAST]\nBLASTDB = /a\\\"b\n").unwrap();
    assert_eq!(settings(&[]).blastdb, Ok(Some("/a\"b".to_string())));
    assert_eq!(ApplicationSettings::default().blastdb, Ok(None));
    std::fs::remove_file(&ini).unwrap();
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::remove_dir(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-866 (cited at `log_file`).
#[test]
fn the_log_file_is_the_last_one_before_the_delimiter() {
    let words = |line: &[&str]| line.iter().map(OsString::from).collect::<Vec<_>>();
    assert_eq!(
        log_file(&words(&["-logfile", "a", "-logfile=b"])),
        Some(OsStr::new("b"))
    );
    assert_eq!(
        log_file(&words(&["-in", "x", "-logfile", "a", "-logfile"])),
        Some(OsStr::new("a"))
    );
    assert_eq!(log_file(&words(&["-logfile"])), None);
    assert_eq!(log_file(&words(&["--", "-logfile", "a"])), None);
    assert_eq!(log_file(&words(&["-logfilex", "a"])), None);
    assert_eq!(log_file(&words(&["-logfile="])), Some(OsStr::new("")));
    assert!(matches!(diag_target(None), DiagTarget::Stderr));
    assert!(matches!(
        diag_target(Some(OsStr::new("-"))),
        DiagTarget::Stderr
    ));
    for name in ["", "/dev/null"] {
        assert!(matches!(
            diag_target(Some(OsStr::new(name))),
            DiagTarget::Discard
        ));
    }
    assert!(matches!(
        diag_target(Some(OsStr::new("/nonexistent-dir-of-losat/log.txt"))),
        DiagTarget::Stderr
    ));
}

/// NCBI's reading of a `-logfile` name on Unix (`log_file_name`): the directory that
/// `SetLogFile` checks, the file that is opened and the directory whose free space counts
/// (line WQA's D11 round 2, 2026-11-20; NCBI BLAST+ 2.17.0 makes `logs` for `logs/`,
/// `run.log` for ` run.log`, `run.log` for `a/../run.log` where `a` is a link, and
/// `abs/run.log` for `<dir>/abs/run.log/`). Round 3: the directory of `//run.log` is `//`,
/// which `CDirEntry::Reset` empties (NCBI keeps the messages on `cerr`, also where `/` is
/// writable); without a working directory a relative name is opened with its '\' made '/'.
#[test]
fn a_log_name_is_read_as_ncbi_reads_it() {
    let read = |name: &[u8], cwd: Option<&[u8]>| {
        let file = log_file_name(PathRules::Unix, name, cwd).expect("Unix reads every name");
        (file.directory, file.path, file.space)
    };
    let paths = |directory: Option<&str>, path: &str, space: &str| {
        (directory.map(Into::into), path.into(), space.into())
    };
    let cwd = Some(&b"/w"[..]);
    for (name, read_as) in [
        ("run.log", paths(Some("."), "/w/run.log", "/w/")),
        ("logs/", paths(Some("."), "/w/logs", "/w/")),
        ("logs//", paths(Some("."), "/w/logs", "/w/")),
        ("run.log/", paths(Some("."), "/w/run.log", "/w/")),
        ("a/b/", paths(Some("a"), "/w/a/b", "/w/a/")),
        ("a/./", paths(Some("a"), "/w/a", "/w/")),
        (" run.log", paths(Some("."), "/w/run.log", "/w/")),
        ("run.log\t", paths(Some("."), "/w/run.log", "/w/")),
        (
            " logs/run.log",
            paths(Some(" logs"), "/w/logs/run.log", "/w/logs/"),
        ),
        (" /run.log", paths(Some(" "), "/w/run.log", "/w/")),
        ("a/../run.log", paths(Some("a/.."), "/w/run.log", "/w/")),
        ("../run.log", paths(Some(".."), "/run.log", "/")),
        (".", paths(Some("."), "/w", "/")),
        ("  ", paths(Some("."), "/w", "/")),
        (
            "/tmp/abs/run.log/",
            paths(Some("/tmp/abs"), "/tmp/abs/run.log", "/tmp/abs/"),
        ),
        ("/dev/fd/2", paths(Some("/dev/fd"), "/dev/fd/2", "/dev/fd/")),
        (
            "//dev/./stderr",
            paths(Some("//dev/."), "/dev/stderr", "//dev/./"),
        ),
        (
            "/dev/fd/../stderr",
            paths(Some("/dev/fd/.."), "/dev/stderr", "/dev/fd/../"),
        ),
        ("/", paths(Some("/"), "/", "/")),
        ("///", paths(None, "/", "")),
        // A directory of '/'s only is no directory (round 3, M1).
        ("//run.log", paths(Some(""), "/run.log", "//")),
        ("///run.log", paths(Some(""), "/run.log", "///")),
        ("//", paths(None, "/", "")),
        // On Unix '\' is a byte of a name, ':' too.
        ("a\\b.log", paths(Some("."), "/w/a\\b.log", "/w/")),
        ("C:run.log", paths(Some("."), "/w/C:run.log", "/w/")),
    ] {
        assert_eq!(read(name.as_bytes(), cwd), read_as, "{name:?}");
    }
    // Without a working directory (`GetCwd` gives "") the name stays relative, and is opened
    // by `ConvertToOSPath`'s relative branch: each '\' is a '/' (round 3, M2), unless the name
    // is a Windows disk's or network path (`IsAbsolutePathEx`).
    for (name, read_as) in [
        ("logs/", paths(Some("."), "logs", "./")),
        ("run.log", paths(Some("."), "run.log", "./")),
        (" run.log ", paths(Some("."), "run.log", "./")),
        ("a\\b.log", paths(Some("."), "a/b.log", "./")),
        ("..\\run.log", paths(Some("."), "../run.log", "./")),
        ("a/..\\run.log", paths(Some("a"), "run.log", "a/")),
        ("x..", paths(Some("."), "x..", "./")),
        ("C:\\run.log", paths(Some("."), "C:\\run.log", "./")),
        ("\\\\run.log", paths(Some("."), "\\\\run.log", "./")),
        ("/abs\\run.log", paths(Some("/"), "/abs\\run.log", "/")),
    ] {
        assert_eq!(read(name.as_bytes(), None), read_as, "{name:?}");
    }
    // NCBI keeps 2047 bytes of the name that it makes.
    let long = vec![b'x'; 3000];
    let (_, path, space) = read(&long, cwd);
    assert_eq!(path, [b"/w/".as_slice(), &long[..2044]].concat());
    assert_eq!(space, b"/w/");
    // On Unix the cut may split a character, whose first bytes NCBI opens (round 4).
    let japanese = "日".repeat(700);
    let (_, path, _) = read(japanese.as_bytes(), cwd);
    assert_eq!(
        path,
        [b"/w/".as_slice(), &japanese.as_bytes()[..2044]].concat()
    );
}

/// NCBI's reading of a `-logfile` name on Windows (by NCBI's source: line WQA's D11 round 3,
/// 2026-11-20; no Windows run): '\' and '/' both separate directories, a relative name is
/// joined to the working directory with '\' and normalized, a name from one separator is on
/// the working directory's disk (its first 3 bytes, as they are), and a relative name with a
/// disk's ':' is none that NCBI reads (`CreateAbsolutePath` throws).
#[test]
fn a_log_name_is_read_on_windows_as_ncbi_reads_it() {
    let read = |name: &str, cwd: Option<&str>| {
        log_file_name(PathRules::Windows, name.as_bytes(), cwd.map(str::as_bytes))
            .map(|file| (file.directory, file.path, file.space))
    };
    let paths = |directory: Option<&str>, path: &str, space: &str| {
        Some((directory.map(Into::into), path.into(), space.into()))
    };
    let cwd = Some("C:\\w");
    for (name, read_as) in [
        ("run.log", paths(Some("."), "C:\\w\\run.log", "C:\\w\\")),
        ("logs/", paths(Some("."), "C:\\w\\logs", "C:\\w\\")),
        ("logs\\", paths(Some("."), "C:\\w\\logs", "C:\\w\\")),
        (" run.log", paths(Some("."), "C:\\w\\run.log", "C:\\w\\")),
        (
            "logs\\run.log",
            paths(Some("logs"), "C:\\w\\logs\\run.log", "C:\\w\\logs\\"),
        ),
        (
            "a/../run.log",
            paths(Some("a/.."), "C:\\w\\run.log", "C:\\w\\"),
        ),
        // After the disk's name a ".." stays (`C:` and `..`).
        (
            "..\\..\\run.log",
            paths(Some("..\\.."), "C:\\..\\run.log", "C:\\..\\"),
        ),
        ("\\run.log", paths(Some("\\"), "C:\\\\run.log", "C:\\\\")),
        ("/run.log", paths(Some("/"), "C:\\\\run.log", "C:\\/")),
        (
            "C:\\logs\\run.log",
            paths(Some("C:\\logs"), "C:\\logs\\run.log", "C:\\logs\\"),
        ),
        ("C:/run.log/", paths(Some("C:/"), "C:\\run.log", "C:/")),
        ("C:\\", paths(Some("C:\\"), "C:\\", "C:\\")),
        (
            "\\\\srv\\share\\run.log",
            paths(
                Some("\\\\srv\\share"),
                "\\\\srv\\share\\run.log",
                "\\\\srv\\share\\",
            ),
        ),
        ("C:run.log", None),
        ("logs\\a:b", None),
        ("//run.log", paths(Some(""), "\\\\run.log", "//")),
    ] {
        assert_eq!(read(name, cwd), read_as, "{name:?}");
    }
    // A cut of NCBI's that splits a character leaves a string that NCBI's conversion to
    // UTF-16 does not read (round 4, L1): the 2047 bytes of the name, the 2 of a part with a
    // ':' (a 3-byte character); where the cuts end characters, the path is read.
    let japanese = "日".repeat(700);
    assert_eq!(read(&japanese, cwd), None);
    let whole = format!("ab{japanese}");
    let (_, path, _) = read(&whole, cwd).unwrap();
    assert_eq!(path, format!("C:\\w\\ab{}", "日".repeat(680)).into_bytes());
    assert_eq!(read("C:\\日:x\\..\\run.log", cwd), None);
    assert_eq!(
        read("C:\\é:x\\..\\run.log", cwd),
        paths(Some("C:\\é:x\\.."), "C:\\é\\..\\run.log", "C:\\é:x\\..\\")
    );
    // Without a working directory a name from one separator stays as it is.
    assert_eq!(
        read("\\run.log", None),
        paths(Some("\\"), "\\run.log", "\\")
    );
    assert_eq!(read("run.log", None), paths(Some("."), "run.log", ".\\"));
    // NCBI's `NormalizePath` on Windows.
    for (path, normalized) in [
        ("a/b\\c", "a\\b\\c"),
        ("a\\.\\b\\..\\c\\", "a\\c"),
        ("\\\\?\\C:\\x\\.\\y", "C:\\x\\y"),
        ("\\\\?\\UNC\\srv\\share\\x", "\\\\srv\\share\\x"),
        ("\\\\srv\\share\\..\\x", "\\\\srv\\x"),
        ("\\a\\..\\..\\b", "\\b"),
        ("C:", "C:"),
        ("C:.", "C:."),
        ("C:\\.", "C:\\"),
        ("C:\\.\\", "C:\\"),
        ("C:\\..", "C:\\.."),
        ("C:\\a:b\\..\\x", "C:\\a:\\..\\x"),
        ("C:\\x\\..", "C:"),
        ("/", "\\"),
        ("//", "\\"),
        ("..\\a", "..\\a"),
    ] {
        assert_eq!(
            PathRules::Windows.normalize(path.as_bytes(), None),
            Some(normalized.as_bytes().to_vec()),
            "{path:?}"
        );
    }
}

/// `CDir::GetCwd` reads a working directory of up to 4094 characters (`getcwd` with a buffer of
/// 4095 bytes), and NCBI's free space is its `Uint8` product, modulo 2^64: 2^64 bytes free is none
/// (line WQA's D11 round 3, 2026-11-20).
#[test]
fn the_working_directory_and_the_free_space_are_ncbi_s() {
    let directory = |length: usize| {
        let mut cwd = "/".to_string();
        cwd.push_str(&"d".repeat(length - 1));
        OsString::from(cwd)
    };
    assert_eq!(getcwd_reads(directory(4094)), Some(directory(4094)));
    assert_eq!(getcwd_reads(directory(4095)), None);
    assert_eq!(getcwd_reads(directory(4298)), None);
    assert_eq!(free_space(4096, 5), LOW_DISK_SPACE);
    assert_eq!(free_space(4096, 1 << 52), 0);
    assert_eq!(free_space(4096, (1 << 52) + 4), 4 * 4096);
    assert_eq!(free_space(u64::MAX, 2), u64::MAX - 1);
}

/// The command line opens the file that NCBI opens (`log_file_name`) and drops the messages
/// of a file whose file system reports no free space (line WQA's D11 round 2, 2026-11-20;
/// NCBI BLAST+ 2.17.0 tblastx `-logfile /dev/fd/2` and `/proc/self/fd/2` write no message).
/// The file that the test opens through procfs is one of its own, whatever its standard
/// streams are (round 3).
#[cfg(target_os = "linux")]
#[test]
fn the_command_line_opens_the_file_that_ncbi_opens() {
    use std::os::fd::AsRawFd;
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-log-name-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("deep/x")).unwrap();
    std::fs::create_dir(dir.join("logs")).unwrap();
    std::os::unix::fs::symlink(dir.join("deep/x"), dir.join("a")).unwrap();
    let target = |name: &str| {
        let name = format!("{}/{name}", dir.display());
        diag_target(Some(OsStr::new(&name)))
    };
    // A trailing '/' names the file without it; a directory stays `cerr`.
    assert!(matches!(target("run.log/"), DiagTarget::File(_)));
    assert!(dir.join("run.log").is_file());
    assert!(matches!(target("logs/"), DiagTarget::Stderr));
    assert!(matches!(target("missing/run.log"), DiagTarget::Stderr));
    // ".." is read before the link is followed: `a/..` is `dir`, not `deep`.
    assert!(matches!(target("a/../linked.log"), DiagTarget::File(_)));
    assert!(dir.join("linked.log").is_file());
    assert!(!dir.join("deep/linked.log").exists());
    // A directory of '/'s only is none (round 3, M1): no file at the root.
    let rooted = format!(
        "//losat-ncbi-environment-{}.log",
        crate::utils::test_dir::process_id()
    );
    assert!(matches!(
        diag_target(Some(OsStr::new(&rooted))),
        DiagTarget::Stderr
    ));
    assert!(!Path::new(&rooted).exists());
    // procfs reports no free space: a file opened through it gets nothing. Through a link
    // in a directory with free space, the same file gets the messages.
    let own = std::fs::File::create(dir.join("own.log")).unwrap();
    let fd = own.as_raw_fd();
    std::os::unix::fs::symlink(format!("/proc/self/fd/{fd}"), dir.join("own-link")).unwrap();
    for name in [
        format!("/dev/fd/{fd}"),
        format!("/proc/self/fd/{fd}"),
        format!("/dev/fd/../fd/{fd}"),
    ] {
        assert!(
            matches!(diag_target(Some(OsStr::new(&name))), DiagTarget::Discard),
            "{name}"
        );
    }
    assert!(matches!(target("own-link"), DiagTarget::File(_)));
    drop(own);
    // A relative name is opened in the working directory that `open_log_file` is given, the
    // one that NCBI reads (`ncbi_cwd`; round 4, L3: not this process's own here).
    assert!(matches!(
        open_log_file(OsStr::new("relative.log"), Some(dir.as_os_str())),
        DiagTarget::File(_)
    ));
    assert!(dir.join("relative.log").is_file());
    assert!(matches!(
        open_log_file(
            OsStr::new("relative.log"),
            Some(dir.join("missing").as_os_str())
        ),
        DiagTarget::Stderr
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}

/// The writes of an application start on the standard error stream and on the file of a
/// `-logfile` for the three targets, as the command line's entry and blastdbcmd make them:
/// the usage report's messages on stderr at once, the first `SetupDiag`, then the
/// configuration's.
fn writes_of(start: &ApplicationStart, target: &str) -> (String, String, Option<i32>) {
    let mut stderr = Vec::new();
    stderr.extend_from_slice(start.usage_report().unwrap().as_bytes());
    let file = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-start-{target}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_file(&file);
    let mut diag = match target {
        "stderr" => DiagTarget::Stderr,
        "file" => diag_target(Some(file.as_os_str())),
        _ => diag_target(Some(OsStr::new("/dev/null"))),
    };
    start.write_log_setup(&mut diag, &mut stderr);
    let status = start.write_config_load(&mut diag, &mut stderr).unwrap();
    drop(diag);
    let logged = std::fs::read_to_string(&file).unwrap_or_default();
    let _ = std::fs::remove_file(&file);
    (String::from_utf8(stderr).unwrap(), logged, status)
}

/// A start with the usage report's messages `u` and the configuration `config`.
fn start_with(u: &str, config: ConfigLoad) -> ApplicationStart {
    start_of("tblastx", u, config)
}

/// `start_with` of `program`.
fn start_of(program: &str, u: &str, config: ConfigLoad) -> ApplicationStart {
    ApplicationStart {
        program: program.to_string(),
        usage_report: Ok(u.to_string()),
        config: Ok(config),
        settings: Ok(ApplicationSettings::default()),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:3863-3874 (cited at
// `ApplicationStart::write_log_setup`), c++/src/corelib/ncbiapp.cpp:669-674 (cited at
// `ApplicationStart`) and c++/src/algo/blast/api/blast_aux.cpp:1181-1187 (cited at
// `replays_to_stderr`). NCBI BLAST+ 2.17.0 (g_reg probes: x_blastp, x_ini, x_iniUR,
// x_inibothbad, o1, o2, l, h_log, inibad_null, inibad_empty; probe_h of 2026-11-12 for every
// program): where the messages go. tblastx, rpsblast, makeblastdb and blastdbcmd replay the
// collected messages to the target; blastn, blastp, blastx, tblastn and psiblast to stderr.
#[test]
fn the_start_writes_where_ncbi_writes() {
    let (u, a, p) = ("U\n", "A\n", "P\n");
    let s = |text: &str| text.to_string();
    let loaded = start_with(u, ConfigLoad::Loaded(s(a)));
    assert_eq!(writes_of(&loaded, "stderr"), (s("U\nA\n"), s(""), None));
    assert_eq!(writes_of(&loaded, "file"), (s(u), s("U\nA\n"), None));
    assert_eq!(writes_of(&loaded, "nowhere"), (s(u), s(""), None));
    let failed = start_with(u, ConfigLoad::Failed(s(p)));
    assert_eq!(
        writes_of(&failed, "stderr"),
        (s("U\nP\nU\nP\n"), s(""), Some(2))
    );
    assert_eq!(
        writes_of(&failed, "file"),
        (s(u), s("U\nP\nU\nP\n"), Some(2))
    );
    assert_eq!(writes_of(&failed, "nowhere"), (s(u), s(""), Some(2)));
    // Without the usage report's messages.
    let failed = start_with("", ConfigLoad::Failed(s(p)));
    assert_eq!(writes_of(&failed, "stderr"), (s("P\nP\n"), s(""), Some(2)));
    assert_eq!(writes_of(&failed, "file"), (s(""), s("P\nP\n"), Some(2)));
    let loaded = start_with("", ConfigLoad::Loaded(s(a)));
    assert_eq!(writes_of(&loaded, "stderr"), (s(a), s(""), None));
    assert_eq!(writes_of(&loaded, "file"), (s(""), s(a), None));
    assert_eq!(writes_of(&loaded, "nowhere"), (s(""), s(""), None));
    // The programs with a `CBlastAppDiagHandler`: the report to the target, the replay to
    // stderr.
    for program in ["blastn", "blastp", "blastx", "tblastn", "psiblast"] {
        let failed = start_of(program, u, ConfigLoad::Failed(s(p)));
        assert_eq!(
            writes_of(&failed, "stderr"),
            (s("U\nP\nU\nP\n"), s(""), Some(2))
        );
        assert_eq!(
            writes_of(&failed, "file"),
            (s("U\nU\nP\n"), s("U\nP\n"), Some(2))
        );
        assert_eq!(
            writes_of(&failed, "nowhere"),
            (s("U\nU\nP\n"), s(""), Some(2))
        );
        let failed = start_of(program, "", ConfigLoad::Failed(s(p)));
        assert_eq!(writes_of(&failed, "file"), (s(p), s(p), Some(2)));
        let loaded = start_of(program, u, ConfigLoad::Loaded(s(a)));
        assert_eq!(writes_of(&loaded, "file"), (s(u), s("U\nA\n"), None));
    }
    for program in ["tblastx", "rpsblast", "makeblastdb", "blastdbcmd"] {
        let failed = start_of(program, u, ConfigLoad::Failed(s(p)));
        assert_eq!(
            writes_of(&failed, "file"),
            (s(u), s("U\nP\nU\nP\n"), Some(2))
        );
    }
    // A refusal at the configuration.
    let mut refused = start_with(u, ConfigLoad::Loaded(s(a)));
    refused.config = Err(s("no"));
    let mut stderr = Vec::new();
    assert_eq!(
        refused.write_config_load(&mut DiagTarget::Stderr, &mut stderr),
        Err(s("no"))
    );
    assert!(stderr.is_empty());
}

/// NCBI's post of a `.ncbirc` with `garbage` on line `line` (the official Linux x64 build;
/// g_reg probe `partial`), in LOSAT's words (`crate::wording::losat_text`).
fn critical(line: usize) -> String {
    crate::wording::losat_text(&format!(
        "Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = {line})\n\n",
        crate::ncbi_exception::NCBI_BUILD_ROOT
    ))
    .into_owned()
}

// NCBI's rpsbproc 0.5.1 (the 2025-05-13 package and bioconda's, env parity probe of
// 2026-11-12 on ws; cited at `has_usage_report`): no usage report, so neither its read of
// `.ncbirc` nor its `BLAST_USAGE_REPORT`; a `.ncbirc` that rpsbproc reports is posted once, at
// the configuration load, naming the 2025-05-13 package's build tree, and the entries before
// the bad line are kept (env parity probe I of 2026-11-14, `i_probe/`); a malformed
// `rpsbproc.ini` or `-conffile` file is NCBI's report, exit status 2 (line RG probes
// `c_ini_bad`, `c_conf_bad`).
#[test]
fn rpsbproc_starts_without_the_usage_report() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-rpsbproc-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let start = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_start(
            "rpsbproc",
            None,
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
    };
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\n").unwrap();
    let read = start(&[("BLAST_USAGE_REPORT", "abc")]);
    assert_eq!(read.usage_report, Ok(String::new()));
    assert_eq!(read.config, Ok(ConfigLoad::Loaded(String::new())));
    assert_eq!(read.settings.unwrap().blastdb, Ok(Some("/db".to_string())));
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
    let malformed = start(&[]);
    assert_eq!(malformed.usage_report, Ok(String::new()));
    let post = "Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n    T0 \"/export/home/tomcat/TeamCity/Agent2/work/b895474bf74772ad/struct_tree/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = 3)\n\n";
    assert_eq!(
        malformed.config,
        Ok(ConfigLoad::Loaded(
            crate::wording::losat_text(post).into_owned()
        ))
    );
    assert_eq!(
        malformed.settings.unwrap().blastdb,
        Ok(Some("/db".to_string()))
    );
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::write(dir.join("rpsbproc.ini"), "garbage\n").unwrap();
    let report = "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n";
    assert_eq!(
        start(&[]).config,
        Ok(ConfigLoad::Failed(
            crate::wording::losat_text(report).into_owned()
        ))
    );
    let _ = std::fs::remove_dir_all(&dir);
}

// dustmasker has no usage report (dust_mask_app.hpp:47-55, cited at `has_usage_report`):
// `BLAST_USAGE_REPORT` is not read, and `.ncbirc` is read once, by the configuration load,
// whose posts go to the diagnostics as blastdbcmd's do.
#[test]
fn dustmasker_starts_without_the_usage_report() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-dustmasker-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let start = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_start(
            "dustmasker",
            None,
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
    };
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\n").unwrap();
    let read = start(&[("BLAST_USAGE_REPORT", "abc")]);
    assert_eq!(read.usage_report, Ok(String::new()));
    assert_eq!(read.config, Ok(ConfigLoad::Loaded(String::new())));
    assert_eq!(read.settings.unwrap().blastdb, Ok(Some("/db".to_string())));
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
    let malformed = start(&[]);
    assert_eq!(malformed.usage_report, Ok(String::new()));
    assert_eq!(malformed.config, Ok(ConfigLoad::Loaded(critical(3))));
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:195-211,
// c++/src/corelib/ncbiapp.cpp:1247-1279 and c++/src/corelib/metareg.cpp:200-209 (cited at
// `application_start`). NCBI BLAST+ 2.17.0 (g_reg probes of 2026-11-11: run.sh `partial`
// default and BLAST_USAGE_REPORT=false, x_dontuse, x_inidont, x_ini, x_inibothbad, b1full,
// b1ini, inib1r): what each read writes.
#[test]
fn the_start_reads_the_files_as_ncbi_reads_them() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-environment-start-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
    let ini = dir.join("blastn.ini");
    let start = |pairs: &[(&str, &str)]| {
        let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
        all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
        let env = |name: &str| {
            all.iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_start(
            "blastn",
            None,
            all.iter()
                .map(|(n, v)| (OsString::from(n), OsString::from(v))),
            &env,
            &NO_PROCESS,
        )
    };
    let loaded = |text: &str| Ok(ConfigLoad::Loaded(text.to_string()));
    let empty = start(&[]);
    assert_eq!(empty.usage_report, Ok(String::new()));
    assert_eq!(empty.config, loaded(""));
    // A syntax error: both reads report it, and the entries before it are used.
    std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
    let partial = start(&[]);
    assert_eq!(partial.usage_report, Ok(critical(3)));
    assert_eq!(partial.config, loaded(&critical(3)));
    assert_eq!(
        partial.settings.as_ref().unwrap().blastdb,
        Ok(Some("/db".to_string()))
    );
    let late = start(&[("BLAST_USAGE_REPORT", "false")]);
    assert_eq!(late.usage_report, Ok(String::new()));
    assert_eq!(late.config, loaded(&critical(3)));
    for variable in ["NCBI_DONT_USE_NCBIRC", "NCBI_CONFIG__NCBI__DONT_USE_NCBIRC"] {
        let off = start(&[(variable, "")]);
        assert_eq!(off.usage_report, Ok(String::new()));
        assert_eq!(off.config, loaded(""));
        assert_eq!(off.settings.unwrap().blastdb, Ok(None));
    }
    // `[NCBI] DONT_USE_NCBIRC` in the program's file: the usage report's read alone.
    std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC=1\n").unwrap();
    let ini_off = start(&[]);
    assert_eq!(ini_off.usage_report, Ok(critical(3)));
    assert_eq!(ini_off.config, loaded(""));
    // A malformed program's file: exit status 2 with NCBI's report, `.ncbirc` not read.
    std::fs::write(&ini, "garbage\n").unwrap();
    let report = "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n";
    let failed = start(&[]);
    assert_eq!(failed.usage_report, Ok(critical(3)));
    assert_eq!(
        failed.config,
        Ok(ConfigLoad::Failed(
            crate::wording::losat_text(report).into_owned()
        ))
    );
    std::fs::remove_file(&ncbirc).unwrap();
    let failed = start(&[]);
    assert_eq!(failed.usage_report, Ok(String::new()));
    assert_eq!(
        failed.config,
        Ok(ConfigLoad::Failed(
            crate::wording::losat_text(report).into_owned()
        ))
    );
    assert!(application_settings(
        "blastn",
        [(
            OsString::from("NCBI_CONFIG_PATH"),
            dir.as_os_str().to_os_string()
        )]
        .into_iter(),
        &|name: &str| (name == "NCBI_CONFIG_PATH").then(|| dir.as_os_str().to_os_string()),
        &NO_PROCESS,
    )
    .unwrap_err()
    .contains("exit status 2"));
    std::fs::remove_file(&ini).unwrap();
    // A `.ncbirc` of one byte that starts a byte order mark: reported by each read that is
    // not kept from the usage report's.
    std::fs::write(&ncbirc, b"\xEF").unwrap();
    let one = start(&[]);
    assert_eq!(one.usage_report, Ok(error_reading_post().to_string()));
    assert_eq!(one.config, loaded(""));
    assert_eq!(
        start(&[("BLAST_USAGE_REPORT", "0")]).config,
        loaded(error_reading_post())
    );
    std::fs::write(&ini, "[BLAST]\n").unwrap();
    assert_eq!(start(&[]).config, loaded(error_reading_post()));
    // A program's file of that byte is read as empty, without a message.
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::write(&ini, b"\xFF").unwrap();
    let silent = start(&[]);
    assert_eq!(silent.config, loaded(""));
    std::fs::remove_file(&ini).unwrap();
    // LOSAT's refusals: a variable it refuses when there is a message to write; a UTF-16
    // file where it is read; a usage report's entry that NCBI aborts on.
    std::fs::write(&ncbirc, "garbage\n").unwrap();
    let refused = start(&[("DIAG_POST_LEVEL", "Error")]);
    assert!(refused
        .usage_report
        .unwrap_err()
        .contains("DIAG_POST_LEVEL"));
    let refused = start(&[("DIAG_POST_LEVEL", "Error"), ("BLAST_USAGE_REPORT", "off")]);
    assert_eq!(refused.usage_report, Ok(String::new()));
    assert!(refused.config.unwrap_err().contains("DIAG_POST_LEVEL"));
    std::fs::remove_file(&ncbirc).unwrap();
    let quiet = start(&[("DIAG_POST_LEVEL", "Error")]);
    assert_eq!(quiet.usage_report, Ok(String::new()));
    assert_eq!(quiet.config, loaded(""));
    assert!(quiet.settings.unwrap_err().contains("DIAG_POST_LEVEL"));
    std::fs::write(&ncbirc, b"\xFF\xFE[\0B\0]\0").unwrap();
    assert!(start(&[]).usage_report.unwrap_err().contains("UTF-16"));
    assert!(start(&[("BLAST_USAGE_REPORT", "no")])
        .config
        .unwrap_err()
        .contains("UTF-16"));
    assert_eq!(start(&[("NCBI_DONT_USE_NCBIRC", "1")]).config, loaded(""));
    std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
    assert!(start(&[])
        .usage_report
        .unwrap_err()
        .contains("[BLAST] BLAST_USAGE_REPORT"));
    std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = off\n").unwrap();
    assert_eq!(start(&[]).usage_report, Ok(String::new()));
    // The usage report alone reads `[BLAST] BLAST_USAGE_REPORT`: not with
    // `BLAST_USAGE_REPORT` false, nor from the program's file (probes `l1.*`).
    std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
    let off = start(&[("BLAST_USAGE_REPORT", "false")]);
    assert_eq!(off.usage_report, Ok(String::new()));
    assert!(off.settings.is_ok());
    assert!(start(&[("NCBI_DONT_USE_NCBIRC", "")]).settings.is_ok());
    assert!(start(&[("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "yes")])
        .settings
        .is_ok());
    std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
    assert!(start(&[]).usage_report.is_err());
    std::fs::remove_file(&ncbirc).unwrap();
    std::fs::write(&ini, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
    let from_ini = start(&[]);
    assert_eq!(from_ini.usage_report, Ok(String::new()));
    assert!(from_ini.settings.is_ok());
    std::fs::remove_file(&ini).unwrap();
    // The variables: `BLAST_USAGE_REPORT` before any file (probes `m2.*`), the environment's
    // layer of the entry where the report makes its registry.
    for value in ["abc", ""] {
        let aborts = start(&[("BLAST_USAGE_REPORT", value)]);
        assert!(aborts
            .usage_report
            .unwrap_err()
            .contains("the environment variable BLAST_USAGE_REPORT"));
    }
    let layer = [("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "maybe")];
    assert!(start(&layer)
        .usage_report
        .unwrap_err()
        .contains("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT"));
    assert!(start(&[layer[0], ("NCBI_DONT_USE_NCBIRC", "1")])
        .usage_report
        .is_err());
    let off = start(&[layer[0], ("BLAST_USAGE_REPORT", "0")]);
    assert_eq!(off.usage_report, Ok(String::new()));
    assert!(off.settings.is_ok());
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1015-1022 and 1271-1281, and
// c++/src/corelib/metareg.cpp `x_FindRegistry` (cited at `conf_file_path`).
#[test]
fn the_conffile_is_found_as_ncbi_finds_it() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-conffile-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("in_path")).unwrap();
    std::fs::write(dir.join("c.ini"), "[BLAST]\n").unwrap();
    std::fs::write(dir.join("in_path").join("c.ini"), "[BLAST]\n").unwrap();
    let named = dir.join("c.ini");
    let search = [dir.join("in_path")];
    // A name with a directory is that file; the search path is not read.
    assert_eq!(
        conf_file_path(named.as_os_str(), &search, &FileSystem),
        Ok(Some(named.clone()))
    );
    // A directory, a missing file and an empty name are NCBI's post.
    for name in [dir.as_os_str(), OsStr::new(""), OsStr::new("no/such.ini")] {
        assert_eq!(
            conf_file_path(name, &search, &FileSystem),
            Err(ConfFile::NotFound(format!(
                "Critical: (106.23) Registry file \"{}\" not found, aborting.\n",
                name.to_string_lossy()
            ))),
            "{name:?}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1015 (cited at `preparse`).
// NCBI BLAST+ 2.17.0 (env parity review G, F4): `blastn -logfile -version` logs to a file
// named `-version` and prints no version.
#[test]
fn preparse_takes_ncbi_s_words() {
    let words = |line: &[&str]| line.iter().map(|word| word.to_string()).collect::<Vec<_>>();
    let kept = |line: &[&str]| match preparse(&words(line)) {
        Preparsed::Words {
            words,
            conffile,
            dryrun,
        } => (words, conffile, dryrun),
        Preparsed::Version(request) => panic!("version {request:?}"),
    };
    use crate::ncbi_version::VersionRequest;
    assert_eq!(
        preparse(&words(&["-db", "-version"])),
        Preparsed::Version(VersionRequest::Version)
    );
    // The value of `-logfile` or `-conffile` is not read as `-version`; a last `-logfile`
    // is dropped.
    assert_eq!(
        kept(&["-info", "-logfile", "-version"]).0,
        words(&["-info", "-logfile", "-version"])
    );
    assert_eq!(kept(&["-conffile", "-version"]).1, Some("-version".into()));
    assert_eq!(kept(&["-info", "-logfile"]).0, words(&["-info"]));
    assert_eq!(kept(&["--", "-version"]).0, words(&["--", "-version"]));
    assert_eq!(kept(&["-info", "-dryrun"]), (words(&["-info"]), None, true));
    assert_eq!(kept(&["-dryrun", "-conffile", "x"]).1, Some("x".into()));
    // The last `-conffile` is the one NCBI checks and loads.
    assert_eq!(kept(&["-conffile", "a", "-conffile=b"]).1, Some("b".into()));
    assert_eq!(kept(&["-conffile=x"]).1, Some("x".into()));
    assert_eq!(kept(&["-conffilex"]).1, None);
    assert_eq!(kept(&["-info", "-conffile"]).1, None);
    // The first version word is the one printed (line RG item E3).
    assert_eq!(
        preparse(&words(&["-version-full", "-version"])),
        Preparsed::Version(VersionRequest::Full)
    );
    assert_eq!(
        preparse(&words(&["-out", "-version-full-json", "-version-full-xml"])),
        Preparsed::Version(VersionRequest::FullJson)
    );
    assert_eq!(kept(&["-version-full=1"]).0, words(&["-version-full=1"]));
    assert_eq!(
        kept(&["--", "-version-full-xml"]).0,
        words(&["--", "-version-full-xml"])
    );
    assert_eq!(
        kept(&["-logfile", "-version-full"]).0,
        words(&["-logfile", "-version-full"])
    );
    // `-version` comes before a `-conffile` NCBI would check after its scan.
    assert_eq!(
        preparse(&words(&["-conffile", "x", "-version"])),
        Preparsed::Version(VersionRequest::Version)
    );
    // The words of the command line's entry.
    let os = [OsString::from("-logfile"), OsString::from("-version")];
    assert!(matches!(preparse(&os), Preparsed::Words { .. }));
    assert_eq!(log_file(&os), Some(OsStr::new("-version")));
}

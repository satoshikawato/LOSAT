use super::*;
use std::sync::Mutex;

use crate::blastinput::ncbi_environment::TextPlace as Place;

#[derive(Default)]
struct Kept {
    files: Mutex<Vec<(Vec<u8>, Vec<u8>)>>,
    registered: Vec<Vec<u8>>,
    read: Vec<Vec<u8>>,
    out_read: Vec<Vec<u8>>,
}

impl Kept {
    fn with_registered(names: &[&str]) -> Self {
        Self {
            files: Mutex::new(Vec::new()),
            registered: names.iter().map(|name| name.as_bytes().to_vec()).collect(),
            read: Vec::new(),
            out_read: Vec::new(),
        }
    }

    fn with_read_files(read: &[&str], out_read: &[&str]) -> Self {
        Self {
            files: Mutex::new(Vec::new()),
            registered: Vec::new(),
            read: read.iter().map(|name| name.as_bytes().to_vec()).collect(),
            out_read: out_read
                .iter()
                .map(|name| name.as_bytes().to_vec())
                .collect(),
        }
    }
}

impl FileReceiver for Kept {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        self.files
            .lock()
            .unwrap()
            .push((name.to_vec(), bytes.to_vec()));
        Ok(())
    }

    fn registered_files(&self) -> &[Vec<u8>] {
        &self.registered
    }

    fn read_files(&self) -> &[Vec<u8>] {
        &self.read
    }

    fn output_read_files(&self) -> &[Vec<u8>] {
        &self.out_read
    }
}

/// The names of `-logfile` of a receiver: NCBI's special names, the names of directories and
/// of the standard error, the standard output, the null device and the file systems without free
/// space, and a file in memory under the command line's file name (`host_file_name`), against the
/// run's real files (line WQA's sweep, D11 and D19).
#[test]
fn a_receiver_reads_a_log_name_as_the_command_line_on_linux() {
    // Registered folders standing for the folders that exist in the command line's directory for
    // the D11 compatibility matrix (plus `in/q.fa` and `s.fa` for D19's existing-file and folder
    // cases).
    let kept = Arc::new(Kept::with_registered(&[
        "logs/keep",
        "a/b/keep",
        "dev/keep",
        " logs/keep",
        "x/keep",
        "in/q.fa",
        "s.fa",
    ]));
    let files = RunFiles::Receiver(kept.clone());
    let target = |name: Option<&str>| files.log_target(name.map(OsStr::new));
    // Linux's limits: a component of 255 bytes, a path of 4095 bytes given to the system.
    let component_255 = "a".repeat(255);
    let long_component = "a".repeat(256);
    let long_file = format!("logs/{}", "b".repeat(256));
    let long_climbed = format!("{}/../run.log", "a".repeat(256));
    // Round 7: the system follows 40 links in a lookup and refuses the 41st (`ELOOP`).
    let links_40 = format!("{}/tmp/run.log", "/var/lock/../..".repeat(40));
    let links_41 = format!("{}/tmp/run.log", "/var/lock/../..".repeat(41));
    let root_links_40 = format!("{}/tmp/run.log", "/proc/self/root".repeat(20));
    let root_links_42 = format!("{}/tmp/run.log", "/proc/self/root".repeat(21));
    let stdin_links_39 = format!("{}/dev/stdin", "/proc/self/root".repeat(18));
    let stdin_links_41 = format!("{}/dev/stdin", "/proc/self/root".repeat(19));
    let directory_4094 = format!("{}run.log", "x/../".repeat(819));
    let directory_4099 = format!("{}run.log", "x/../".repeat(820));
    let directory_4095 = format!("x{}/run.log", "/.".repeat(2047));
    let directory_4096 = format!("xx{}/run.log", "/.".repeat(2047));
    assert!(matches!(target(None), Ok((DiagTarget::Stderr, None))));
    assert!(matches!(target(Some("-")), Ok((DiagTarget::Stderr, None))));
    assert!(matches!(target(Some("")), Ok((DiagTarget::Discard, None))));
    assert!(matches!(
        target(Some("/dev/null")),
        Ok((DiagTarget::Discard, None))
    ));
    for name in [
        ".",
        "..",
        "../..",
        "/",
        "logs/.",
        "a/b/..",
        "a/./",
        "/tmp/",
        "/tmp",
        "./",
        "//",
        "  ",
        // Case 0 & 1: missing folder, file as folder, and existing folder inside working dir.
        "nodir/run.log",
        "nodir/../run.log",
        "s.fa/run.log",
        "in/q.fa/run.log",
        "in",
        "in/",
        "logs",
        "logs/",
        "logs//",
        "a/b",
        "a/b/",
        "src",
        "src/",
        "/proc/self/cwd/nodir/run.log",
        "/dev/fd",
        "/proc/self/fd/",
        "/home/../usr",
        "/dev/fd/3",
        "/dev/fd/63",
        "/proc/self/fd/3",
        "/proc/self/status",
        "/proc/1/comm",
        "/proc/self/losat",
        "/proc/self/cwd",
        "/proc/self/cwd/../stderr",
        "/dev/null/x",
        "/dev/stderr/x",
        "/dev/fd/2/x",
        "/dev/stdout/../stderr",
        "/dev/fd/1/../2",
        "/dev/pts/0",
        "/dev/pts/ptmx",
        "/dev/mqueue/a/b",
        "/sys/kernel/x",
        "/dev/tty",
        // Round 3: a directory of '/'s only is none, and a user makes no file in a
        // directory of root's and opens none of its other devices.
        "//run.log",
        "///run.log",
        "/run.log",
        "/etc/run.log",
        "/var/log/run.log",
        "/usr/local/run.log",
        "/home/run.log",
        "/root/run.log",
        "/run/run.log",
        "/dev/run.log",
        "/dev/net/run.log",
        "/dev/kmsg",
        "/dev/console",
        "/dev/ptmx/x",
        "/proc/self/root/etc/run.log",
        // Round 4: one rule of owners (only root makes entries in `/`), the links of Linux,
        // the directories of root's that are there, and the system's limits of a path.
        "/data/run.log",
        "/data/logs/run.log",
        "/lost+found/run.log",
        "/snap/run.log",
        "/opt/run.log",
        "/dev/mqueue/a/run.log",
        "/bin/../tmp/run.log",
        "/var/lock/../tmp/run.log",
        "/dev/fd/../../tmp/run.log",
        "/data/../tmp/run.log",
        "/var/lock",
        "/bin",
        long_component.as_str(),
        long_file.as_str(),
        long_climbed.as_str(),
        directory_4096.as_str(),
        directory_4099.as_str(),
        // Round 5: a name ends in a directory that the model knows, and passes none that is
        // not there (a file, a name at the root that is not there).
        "/lib/../tmp/run.log",
        "/var/run/lock/../run.log",
        "/etc/default/../tmp/run.log",
        "/dev/null/../../tmp/run.log",
        "/dev/fd/2/../../../tmp/run.log",
        "/root/x/../../tmp/run.log",
        "/data/x/../../tmp/run.log",
        // Round 6: a name looks nothing up in a directory that the user cannot search (as a
        // user's command line, `EACCES`; measured on the model's host as its first user).
        "/root/../tmp/run.log",
        "/root/./run.log",
        "/var/cache/ldconfig/../../../tmp/run.log",
        "/etc/ssl/private/../../../tmp/run.log",
        "/var/spool/cron/crontabs/../../../../tmp/run.log",
        "/proc/1/fd/../../../tmp/run.log",
        "/proc/1/root/../../../tmp/run.log",
        "/proc/1/task/1/cwd/../../../../../tmp/run.log",
        // Round 7: the thread's directory is a link to a task of the process (`..` is in its
        // `task`), the 41st link is refused, another process's `map_files` and `exe` are
        // closed, and the model's host's Linux 5.15 has no `autogroup` (review r6, all run).
        "/proc/thread-self/../../../tmp/run.log",
        "/proc/thread-self/../../tmp/run.log",
        "/proc/thread-self/fd/../../../tmp/run.log",
        links_41.as_str(),
        root_links_42.as_str(),
        stdin_links_41.as_str(),
        "/proc/1/map_files/../../../tmp/run.log",
        "/proc/1/exe/../../../tmp/run.log",
        "/proc/self/autogroup",
    ] {
        assert!(
            matches!(target(Some(name)), Ok((DiagTarget::Stderr, None))),
            "{name}"
        );
    }
    // Case 2: a file of the run's directory (registered or in `/work`) is refused.
    for existing in [
        "s.fa",
        "./s.fa",
        " s.fa ",
        "s.fa/",
        "in/../s.fa",
        "in/q.fa",
        "/proc/self/cwd/s.fa",
        "Cargo.toml",
    ] {
        let reason = target(Some(existing)).err().expect("existing file refused");
        assert_eq!(
            reason,
            format!(
                "a -logfile file '{existing}' that is a file of the run's directory (LOSAT Web \
                 reads the run's files and writes none of them)"
            ),
            "{existing}"
        );
    }
    for name in [
        "/dev/stderr",
        "//dev/stderr",
        "/dev/./stderr",
        "/tmp/../dev/stderr",
        "/dev/stderr/",
        "/proc/self/root/dev/stderr",
    ] {
        assert!(
            matches!(target(Some(name)), Ok((DiagTarget::StderrFile, None))),
            "{name}"
        );
    }
    for name in [
        "//dev/null",
        "/dev/./null",
        "/tmp/../dev/null",
        "/dev/null/",
        "/dev/fd/2",
        "/dev/fd/1",
        "/dev/fd/0",
        "/dev/stdin",
        "/proc/self/fd/2",
        "/proc/thread-self/fd/1",
        "/dev/fd/../stderr",
        "/dev/fd/../stdout",
        "/dev/fd/../fd/2",
        "/dev/mqueue/q",
        "/dev/full",
        "/dev/zero",
        "/proc/self/comm",
        "/proc/thread-self/comm",
        "/proc/self/oom_score_adj",
        "/proc/self/attr/current",
        "/dev/random",
        "/dev/urandom",
        "/dev/ptmx",
        "/dev/fuse",
        "/dev/net/tun",
        "//dev/ptmx",
        stdin_links_39.as_str(),
    ] {
        assert!(
            matches!(target(Some(name)), Ok((DiagTarget::Discard, None))),
            "{name}"
        );
    }
    // Case 3 & 4: files made, with the command line's host name.
    let names = [
        ("run.log", "run.log"),
        ("logs/run.log", "logs/run.log"),
        ("in/run.log", "in/run.log"),
        ("src/run.log", "src/run.log"),
        ("./run.log", "run.log"),
        ("a/../run.log", "run.log"),
        ("/tmp/run.log", "/tmp/run.log"),
        ("/dev/stdout", "/dev/stdout"),
        ("//dev/stdout", "/dev/stdout"),
        ("/dev/stdout/", "/dev/stdout"),
        ("dev/stderr", "dev/stderr"),
        ("../run.log", "../run.log"),
        ("a/../../run.log", "../run.log"),
        ("../proc/self/fd/2", "../proc/self/fd/2"),
        ("newlog/", "newlog"),
        ("newlog//", "newlog"),
        ("run.log/", "run.log"),
        ("a/newb/", "a/newb"),
        ("/tmp/abs/run.log/", "/tmp/abs/run.log"),
        (" run.log", "run.log"),
        ("run.log ", "run.log"),
        (" logs/run.log", "logs/run.log"),
        ("/proc/self/cwd/run.log", "run.log"),
        ("/proc/self/root/tmp/run.log", "/proc/self/root/tmp/run.log"),
        ("/home/me/run.log", "/home/me/run.log"),
        ("/home/me/logs/run.log", "/home/me/logs/run.log"),
        ("/var/tmp/run.log", "/var/tmp/run.log"),
        ("/dev/shm/run.log", "/dev/shm/run.log"),
        ("/run/lock/run.log", "/run/lock/run.log"),
        // Round 4.
        ("/var/lock/run.log", "/var/lock/run.log"),
        ("/var/run/lock/run.log", "/var/run/lock/run.log"),
        ("/var/crash/run.log", "/var/crash/run.log"),
        ("/var/run/../tmp/run.log", "/var/tmp/run.log"),
        ("/usr/../tmp/run.log", "/tmp/run.log"),
        ("/tmp/a/../run.log", "/tmp/run.log"),
        ("/dev/shm/a/run.log", "/dev/shm/a/run.log"),
        (component_255.as_str(), component_255.as_str()),
        (directory_4094.as_str(), "run.log"),
        (directory_4095.as_str(), "x/run.log"),
        // Round 5: a name climbs out of directories of root's that the model does not list,
        // which a user's command line passes (NCBI BLAST+ 2.17.0 writes the file; review 3).
        ("/etc/default/../../tmp/run.log", "/tmp/run.log"),
        ("/usr/include/../../tmp/run.log", "/tmp/run.log"),
        ("/usr/share/doc/../../../tmp/run.log", "/tmp/run.log"),
        (
            "/usr/lib/x86_64-linux-gnu/../../../tmp/run.log",
            "/tmp/run.log",
        ),
        ("/var/mail/../tmp/run.log", "/var/tmp/run.log"),
        ("/var/backups/../tmp/run.log", "/var/tmp/run.log"),
        ("/dev/disk/../../tmp/run.log", "/tmp/run.log"),
        ("/var/log/../../tmp/run.log", "/tmp/run.log"),
        ("/var/lock/../../tmp/run.log", "/tmp/run.log"),
        // Round 6: the directories around them, and a process's own descriptors, open to the
        // user.
        ("/var/cache/../../tmp/run.log", "/tmp/run.log"),
        ("/etc/ssl/../../tmp/run.log", "/tmp/run.log"),
        ("/proc/1/../../tmp/run.log", "/tmp/run.log"),
        ("/proc/self/fd/../../../tmp/run.log", "/tmp/run.log"),
        // Round 7: the model's host and its first user (review r6): the thread's `..` climbs
        // out of procfs, 40 links, `fdinfo` of Linux 5.15, a directory of the user's group
        // `adm`, a directory at the host's root, and places where the user writes.
        ("/proc/thread-self/../../../../tmp/run.log", "/tmp/run.log"),
        (links_40.as_str(), "/tmp/run.log"),
        (root_links_40.as_str(), root_links_40.as_str()),
        ("/proc/1/fdinfo/../../../tmp/run.log", "/tmp/run.log"),
        (
            "/var/log/unattended-upgrades/../../../tmp/run.log",
            "/tmp/run.log",
        ),
        ("/snap/../tmp/run.log", "/tmp/run.log"),
        ("/run/user/1000/run.log", "/run/user/1000/run.log"),
        ("/var/run/user/1000/run.log", "/var/run/user/1000/run.log"),
        ("/run/shm/run.log", "/run/shm/run.log"),
        ("/mnt/wsl/run.log", "/mnt/wsl/run.log"),
    ];
    for (name, expected_host_name) in names {
        let (log, host_name) = target(Some(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(matches!(log, DiagTarget::Memory(_)), "{name}");
        assert_eq!(
            host_name.as_deref(),
            Some(expected_host_name.as_bytes()),
            "{name}"
        );
        log.write_text(name.as_bytes(), &mut io::sink());
        files.deliver_log(host_name.as_deref(), log).unwrap();
    }
    // A target without a file hands nothing over.
    files.deliver_log(Some(b"-"), DiagTarget::Stderr).unwrap();
    files
        .deliver_log(Some(b"/dev/stderr"), DiagTarget::StderrFile)
        .unwrap();
    files
        .deliver_log(Some(b"/dev/fd/2"), DiagTarget::Discard)
        .unwrap();
    assert_eq!(
        *kept.files.lock().unwrap(),
        names
            .iter()
            .map(|(name, host_name)| (host_name.as_bytes().to_vec(), name.as_bytes().to_vec()))
            .collect::<Vec<_>>()
    );
}

/// A name of 4 KB is read in time linear in its length (review 3 M1, line WQA's D11 round 5,
/// 2026-11-21: resolving each of its directories again from `/` took about a minute in LOSAT
/// Web's serial engine): the walk takes each component a few times, not once per directory
/// before it. The 2047 bytes that NCBI keeps of the first name made absolute are a directory of
/// it, so its messages stay on `cerr`; the second's `..` resolve before the cut, and it names a
/// file. Both as on the command line.
#[test]
fn a_long_log_name_is_read_in_linear_time() {
    let all_dirs: EntryLookup<'_> = &|path, cwd| inside_cwd(path, cwd).map(|_| Entry::Directory);
    let empty: EntryLookup<'_> = &|path, cwd| inside_cwd(path, cwd).map(|_| Entry::Missing);
    let deep = format!("{}run.log", "a/".repeat(2040));
    let up = format!("{}run.log", "../".repeat(1360));
    for (name, entry, is_file) in [
        (deep.clone(), all_dirs, false),
        (deep, empty, false),
        (up.clone(), empty, true),
    ] {
        assert!(name.len() < 4096);
        WALK_STEPS.with(|steps| steps.set(0));
        let started = std::time::Instant::now();
        let read = log_name(name.as_bytes(), Some(entry));
        assert_eq!(matches!(read, LogName::File { .. }), is_file);
        if is_file {
            assert!(matches!(read, LogName::File { name, .. } if name == up.as_bytes()));
        }
        let steps = WALK_STEPS.with(std::cell::Cell::get);
        assert!(steps < 6 * name.len(), "{steps} steps");
        assert!(started.elapsed() < std::time::Duration::from_secs(5));
    }
}

/// The engine carries no list of the system's files (line WQA's D11, Q-D11 = B, 2026-11-21): a
/// name that climbs with `..` out of a file or out of a missing name in root's trees names a
/// file in LOSAT Web, where the command line on Linux refuses the path (`ENOTDIR`, `ENOENT`)
/// and writes the messages to its standard error (`/proc/self/status` too, a file of procfs
/// that a process does not write). docs/web/abi_v2.md and the manual state this difference,
/// so a change to it changes this test and those texts together. A file that the model knows
/// refuses the path as the command line does (`/dev/null`, `/proc/self/comm`).
#[test]
fn a_name_through_a_file_of_roots_names_a_file_as_documented() {
    let empty = |path: &[u8], cwd: &[u8]| inside_cwd(path, cwd).map(|_| Entry::Missing);
    for (name, expected_file) in [
        ("/etc/passwd/../../tmp/x", Some(b"/tmp/x".as_slice())),
        ("/etc/zzz/../../tmp/x", Some(b"/tmp/x".as_slice())),
        ("/usr/bin/ls/../../../tmp/x", Some(b"/tmp/x".as_slice())),
        (
            "/proc/self/status/../../../tmp/x",
            Some(b"/tmp/x".as_slice()),
        ),
        ("/dev/null/../../tmp/x", None),
        ("/proc/self/comm/../../../tmp/x", None),
    ] {
        match (log_name(name.as_bytes(), Some(&empty)), expected_file) {
            (LogName::File { name: got, .. }, Some(want)) => assert_eq!(got, want, "{name}"),
            (LogName::Unopened, None) => {}
            (other, _) => panic!("{name}: unexpected {other:?}"),
        }
    }
}

/// The file that is the standard error writes where `cerr` does, and gets the messages
/// collected before it was opened, which `cerr` has already; a file in memory keeps them, and
/// `cerr` gets none (line WQA's sweep, D11, 2026-11-20).
#[test]
fn the_standard_error_file_writes_to_the_standard_error() {
    use crate::blastinput::ncbi_environment::TextPlace;
    let mut stderr = Vec::new();
    let file = DiagTarget::StderrFile;
    assert_eq!(file.place(), TextPlace::StandardError);
    file.write_collected(b"collected\n", &mut stderr);
    file.write_text(b"posted\n", &mut stderr);
    file.write_to_file(b"error\n");
    assert_eq!(stderr, b"collected\nposted\n");
    let mut stderr = Vec::new();
    DiagTarget::Stderr.write_collected(b"collected\n", &mut stderr);
    assert!(stderr.is_empty());
    let memory = DiagTarget::Memory(Default::default());
    memory.write_collected(b"collected\n", &mut stderr);
    assert!(stderr.is_empty());
    let DiagTarget::Memory(bytes) = memory else {
        unreachable!()
    };
    assert_eq!(bytes.into_inner().unwrap(), b"collected\n");
}

/// The names of `-export_search_strategy` of a receiver (LOSAT Web, line RG's F1w), read as the
/// command line on Linux opens them for writing, the name as typed (`strategy_name`), on the host
/// and as the user of the `-logfile` model: a file of the run where the user makes files, NCBI's
/// error where the command line cannot open the name (a directory, a name ending in `/`, a name in
/// a directory of root's, a component of 256 bytes), the standard output's refusal, no file for the
/// null device and the standard input, and LOSAT Web's refusal of the standard error, of the
/// devices whose writes fail or wait for a reader, and of procfs's and mqueue's files.
#[test]
fn a_strategy_name_is_read_as_the_command_line_opens_it_on_linux() {
    use StrategyTarget as T;
    let files = RunFiles::Receiver(Arc::new(Kept::default()));
    let target = |name: &str| files.strategy_target(OsStr::new(name), &[]);
    let long_component = "a".repeat(256);
    let component_255 = "a".repeat(255);
    for (name, read) in [
        ("s.asn", T::File),
        ("out/s.asn", T::File),
        ("./s.asn", T::File),
        (" s.asn", T::File),
        ("../s.asn", T::File),
        ("a/b/../s.asn", T::File),
        ("/tmp/s.asn", T::File),
        ("/dev/shm/s.asn", T::File),
        ("/var/lock/s.asn", T::File),
        ("/home/me/s.asn", T::File),
        ("/proc/self/cwd/s.asn", T::File),
        ("/proc/thread-self/cwd/s.asn", T::File),
        (component_255.as_str(), T::File),
        ("", T::Inaccessible),
        (".", T::Inaccessible),
        ("..", T::Inaccessible),
        ("out/", T::Inaccessible),
        ("out/.", T::Inaccessible),
        ("/tmp", T::Inaccessible),
        ("/tmp/", T::Inaccessible),
        ("/s.asn", T::Inaccessible),
        ("//s.asn", T::Inaccessible),
        ("/etc/s.asn", T::Inaccessible),
        ("/data/s.asn", T::Inaccessible),
        ("/root/s.asn", T::Inaccessible),
        ("/root/../tmp/s.asn", T::Inaccessible),
        ("/dev/s.asn", T::Inaccessible),
        ("/dev/kmsg", T::Inaccessible),
        ("/dev/null/", T::Inaccessible),
        ("/dev/null/s.asn", T::Inaccessible),
        ("/proc/self/status", T::Inaccessible),
        ("/proc/self/cwd", T::Inaccessible),
        ("/sys/s.asn", T::Inaccessible),
        (long_component.as_str(), T::Inaccessible),
        ("/dev/stdout", T::StandardOutput),
        ("/dev/fd/1", T::StandardOutput),
        ("/proc/self/fd/1", T::StandardOutput),
        ("/proc/thread-self/fd/1", T::StandardOutput),
        ("/dev/null", T::Nowhere),
        ("//dev/null", T::Nowhere),
        ("/dev/zero", T::Nowhere),
        ("/dev/random", T::Nowhere),
        ("/dev/urandom", T::Nowhere),
        ("/dev/stdin", T::Nowhere),
        ("/proc/self/fd/0", T::Nowhere),
        ("/dev/stderr", T::System),
        ("/dev/fd/2", T::System),
        ("/dev/full", T::System),
        ("/dev/ptmx", T::System),
        ("/dev/fuse", T::System),
        ("/dev/net/tun", T::System),
        ("/proc/self/comm", T::System),
        ("/dev/mqueue/q", T::System),
    ] {
        assert_eq!(target(name), read, "{name:?}");
    }
}

/// A strategy name that is the file of another option of the run, as the system resolves both
/// (`-query`, `-subject` and `-window_masker_db` as typed, `-logfile` as NCBI reads it), with
/// either name climbing more: LOSAT Web's refusal names the option. The standard input (`-`) and
/// NCBI's special log names are no file of the run.
#[test]
fn a_strategy_name_of_another_file_of_the_run_names_its_option() {
    use StrategyTarget as T;
    let files = RunFiles::Receiver(Arc::new(Kept::default()));
    let target = |name: &str, others: &[(&'static str, &str)]| {
        let others: Vec<RunFileName<'_>> = others
            .iter()
            .map(|&(option, name)| RunFileName {
                option,
                name: OsStr::new(name),
            })
            .collect();
        files.strategy_target(OsStr::new(name), &others)
    };
    for (name, others, read) in [
        ("q.fa", vec![("query", "q.fa")], T::OtherFile("query")),
        ("./q.fa", vec![("query", "q.fa")], T::OtherFile("query")),
        (
            "/proc/self/cwd/q.fa",
            vec![("query", "q.fa")],
            T::OtherFile("query"),
        ),
        ("q.fa", vec![("query", "../x/../w/q.fa")], T::File),
        (
            "../w/q.fa",
            vec![("query", "../x/../w/q.fa")],
            T::OtherFile("query"),
        ),
        (
            "../../a/q.fa",
            vec![("query", "../../a/q.fa")],
            T::OtherFile("query"),
        ),
        ("../a/q.fa", vec![("query", "../../a/q.fa")], T::File),
        ("q.fa.asn", vec![("query", "q.fa")], T::File),
        (
            "/tmp/s.fa",
            vec![("query", "q.fa"), ("subject", "/tmp/s.fa")],
            T::OtherFile("subject"),
        ),
        ("./-", vec![("query", "-")], T::File),
        (
            "run.log",
            vec![("logfile", "run.log")],
            T::OtherFile("logfile"),
        ),
        (
            "run.log",
            vec![("logfile", " run.log")],
            T::OtherFile("logfile"),
        ),
        ("logs", vec![("logfile", "logs/")], T::OtherFile("logfile")),
        (
            "../run.log",
            vec![("logfile", "../run.log")],
            T::OtherFile("logfile"),
        ),
        ("run.log", vec![("logfile", "../run.log")], T::File),
        ("./-", vec![("logfile", "-")], T::File),
        (
            "wm.obinary",
            vec![("window_masker_db", "wm.obinary")],
            T::OtherFile("window_masker_db"),
        ),
    ] {
        assert_eq!(target(name, &others), read, "{name:?} {others:?}");
    }
}

/// The output file names of a receiver (`RunFiles::create`, `output_name`; line NP WW3a,
/// 2026-11-22): a file in the run's working directory, also in a folder that the run's directory
/// has; no file in a folder that it does not have or in a file of it, for a directory, a name in a
/// directory of root's or one past a directory that a user cannot search; the standard output's
/// names; the lost writes of the null device, of `/dev/zero`, of the random devices and of the
/// standard input (the null device of LOSAT Web's run), and the failed ones of `/dev/full`; the
/// names that LOSAT Web does not write (the standard error, the other devices that a user opens, the
/// files of procfs and the queues of `/dev/mqueue`) and an existing file of the run's directory,
/// which the command line empties; a file of that name outside the working directory where a user
/// makes files.
#[test]
fn a_receiver_reads_an_output_name_as_the_command_line_on_linux() {
    let entry = |path: &[u8], cwd: &[u8]| {
        inside_cwd(path, cwd).map(|rel| match rel {
            b"sub" | b"sub/deep" => Entry::Directory,
            b"in.fa" | b"sub/f" => Entry::File,
            _ => Entry::Missing,
        })
    };
    let name = |name: &str| output_name(name.as_bytes(), &entry);
    let long = "a".repeat(256);
    for (given, expected_host_name) in [
        ("c.counts", "c.counts"),
        (" c.counts", " c.counts"),
        ("./c.counts", "c.counts"),
        ("sub/c.counts", "sub/c.counts"),
        ("sub/deep/c.counts", "sub/deep/c.counts"),
        ("sub/../c.counts", "c.counts"),
        ("./sub/./c.counts", "sub/c.counts"),
        ("../c.counts", "../c.counts"),
        ("sub/../../c.counts", "../c.counts"),
        ("/proc/self/cwd/c.counts", "c.counts"),
        ("/tmp/c.counts", "/tmp/c.counts"),
        ("/home/me/c.counts", "/home/me/c.counts"),
        ("/etc/default/../../tmp/c.counts", "/tmp/c.counts"),
        ("/var/lock/c.counts", "/var/lock/c.counts"),
    ] {
        match name(given) {
            OutputName::File { name: got, .. } => {
                assert_eq!(got, expected_host_name.as_bytes(), "{given:?}");
            }
            other => panic!("{given:?}: expected File, got {other:?}"),
        }
    }
    for (given, expected) in [
        (
            "nodir/c.counts",
            OutputName::Unopened(OpenFailure::NotFound),
        ),
        (
            "nodir/../c.counts",
            OutputName::Unopened(OpenFailure::NotFound),
        ),
        (
            "in.fa/c.counts",
            OutputName::Unopened(OpenFailure::NotADirectory),
        ),
        (
            "/proc/self/cwd/nodir/c.counts",
            OutputName::Unopened(OpenFailure::NotFound),
        ),
        ("sub", OutputName::Unopened(OpenFailure::IsADirectory)),
        ("sub/", OutputName::Unopened(OpenFailure::IsADirectory)),
        (".", OutputName::Unopened(OpenFailure::IsADirectory)),
        ("..", OutputName::Unopened(OpenFailure::IsADirectory)),
        ("sub/..", OutputName::Unopened(OpenFailure::IsADirectory)),
        (
            "/proc/self/cwd",
            OutputName::Unopened(OpenFailure::IsADirectory),
        ),
        ("/tmp", OutputName::Unopened(OpenFailure::IsADirectory)),
        (
            "/c.counts",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        (
            "/etc/c.counts",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        (
            "/data/c.counts",
            OutputName::Unopened(OpenFailure::NotFound),
        ),
        (
            "/root/../tmp/c.counts",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        (
            "/proc/self/status",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        (
            "/dev/fd/3",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        (
            "/dev/kmsg",
            OutputName::Unopened(OpenFailure::PermissionDenied),
        ),
        ("", OutputName::Unopened(OpenFailure::NotFound)),
        (
            long.as_str(),
            OutputName::Unopened(OpenFailure::NameTooLong),
        ),
        ("in.fa", OutputName::Existing),
        ("sub/f", OutputName::Existing),
        ("/proc/self/cwd/in.fa", OutputName::Existing),
        ("/dev/stdout", OutputName::Stdout),
        ("/dev/fd/1", OutputName::Stdout),
        ("/proc/self/fd/1", OutputName::Stdout),
        ("/proc/thread-self/fd/1", OutputName::Stdout),
        ("/dev/stdin", OutputName::Discard),
        ("/dev/null", OutputName::Discard),
        ("//dev/null", OutputName::Discard),
        ("/dev/zero", OutputName::Discard),
        ("/dev/urandom", OutputName::Discard),
        ("/dev/full", OutputName::Full),
        ("/dev/stderr", OutputName::Special),
        ("/dev/fd/2", OutputName::Special),
        ("/dev/ptmx", OutputName::Special),
        ("/proc/self/comm", OutputName::Special),
        ("/dev/mqueue/q", OutputName::Special),
    ] {
        assert_eq!(name(given), expected, "{given:?}");
    }
}

/// `log_name` and `output_name` read the run's directory through one view (`run_entries`,
/// `run_entry`), covering both `/work` entries and registered search names (case 0-4, D19).
#[test]
fn log_and_output_names_share_the_runs_files_view() {
    let entry = |path: &[u8], cwd: &[u8]| {
        inside_cwd(path, cwd).map(|rel| match rel {
            b"db" | b"in" | b"a" => Entry::Directory,
            b"q.fa" | b"in/q.fa" | b"db/keep" | b".ncbirc" => Entry::File,
            _ => Entry::Missing,
        })
    };
    let check_file = |raw: &str, want: &str| match log_name(raw.as_bytes(), Some(&entry)) {
        LogName::File { name, .. } => assert_eq!(name, want.as_bytes(), "{raw:?}"),
        other => panic!("{raw:?}: expected File({want:?}), got {other:?}"),
    };
    // Case 0: folder absent or file as folder -> Unopened; folder present -> File.
    for raw in [
        "logs/run.log",
        "a/../logs/run.log",
        "q.fa/run.log",
        "in/q.fa/run.log",
        "nodir/../run.log",
        "/proc/self/cwd/logs/run.log",
    ] {
        assert_eq!(
            log_name(raw.as_bytes(), Some(&entry)),
            LogName::Unopened,
            "{raw:?}"
        );
    }
    check_file("db/run.log", "db/run.log");
    check_file("in/run.log", "in/run.log");
    check_file("a/../run.log", "run.log");
    // Case 1: folder of the run -> Unopened.
    for raw in ["db", "db/", "in", "in/", "a", "/proc/self/cwd/db"] {
        assert_eq!(
            log_name(raw.as_bytes(), Some(&entry)),
            LogName::Unopened,
            "{raw:?}"
        );
    }
    // Case 2: file of the run -> Existing.
    for raw in [
        "q.fa",
        "./q.fa",
        " q.fa",
        "q.fa/",
        "a/../q.fa",
        "in/q.fa",
        "db/keep",
        ".ncbirc",
        "/proc/self/cwd/q.fa",
    ] {
        assert_eq!(
            log_name(raw.as_bytes(), Some(&entry)),
            LogName::Existing,
            "{raw:?}"
        );
    }
    // Case 3: /dev/stdout -> Stdout; procfs fd/1 stays Nowhere.
    for raw in [
        "/dev/stdout",
        "//dev/stdout",
        "/dev/./stdout",
        "/dev/stdout/",
        "/tmp/../dev/stdout",
        "/proc/self/root/dev/stdout",
    ] {
        assert_eq!(
            log_name(raw.as_bytes(), Some(&entry)),
            LogName::Stdout,
            "{raw:?}"
        );
    }
    for raw in [
        "/dev/fd/1",
        "/proc/self/fd/1",
        "/proc/thread-self/fd/1",
        "/dev/fd/../stdout",
    ] {
        assert_eq!(
            log_name(raw.as_bytes(), Some(&entry)),
            LogName::Nowhere,
            "{raw:?}"
        );
    }
    // `RunFiles::log_host_name` (used by `losat_web2_start` without a run directory).
    for (raw, want) in [
        ("run.log", Some("run.log")),
        ("./run.log", Some("run.log")),
        (" run.log", Some("run.log")),
        ("logs/", Some("logs")),
        ("logs/run.log", Some("logs/run.log")),
        ("/proc/self/cwd/run.log", Some("run.log")),
        ("/etc/default/../../tmp/run.log", Some("/tmp/run.log")),
        ("/dev/stdout", Some("/dev/stdout")),
        ("//dev/stdout", Some("/dev/stdout")),
        ("-", None),
        ("", None),
        ("/dev/null", None),
        (".", None),
        ("/dev/stderr", None),
        ("/dev/fd/2", None),
        ("/etc/run.log", None),
    ] {
        assert_eq!(
            RunFiles::log_host_name(OsStr::new(raw)).as_deref(),
            want.map(str::as_bytes),
            "{raw:?}"
        );
    }
    // Case 4: normalized relative names inside the working tree, NCBI's open path outside.
    check_file(" run.log", "run.log");
    check_file("logs/", "logs");
    check_file("./run.log", "run.log");
    check_file("/proc/self/cwd/run.log", "run.log");
    check_file("../x.log", "../x.log");
    check_file("a/../../x.log", "../x.log");
    check_file("/tmp/run.log", "/tmp/run.log");
    check_file("/etc/default/../../tmp/run.log", "/tmp/run.log");
    // Registered names resolved with `Walk` feed both `log_name` and `output_name`, inside and
    // outside the working directory.
    let receiver = Kept::with_registered(&[
        "./q.fa",
        "sub/../s.fa",
        "in/a.fa",
        "/abs/outside.fa",
        "../up.fa",
    ]);
    let resolved = run_entries(&receiver);
    let cwd1 = working_directory(b"q.fa");
    let cwd2 = working_directory(b"../up.fa");
    assert_eq!(run_entry(&resolved, &cwd1, &cwd1), Some(Entry::Directory));
    assert_eq!(
        run_entry(&resolved, &[cwd1.as_slice(), b"/q.fa"].concat(), &cwd1),
        Some(Entry::File)
    );
    assert_eq!(
        run_entry(&resolved, &[cwd1.as_slice(), b"/s.fa"].concat(), &cwd1),
        Some(Entry::File)
    );
    assert_eq!(
        run_entry(&resolved, &[cwd1.as_slice(), b"/in"].concat(), &cwd1),
        Some(Entry::Directory)
    );
    assert_eq!(
        run_entry(&resolved, &[cwd1.as_slice(), b"/in/a.fa"].concat(), &cwd1),
        Some(Entry::File)
    );
    assert_eq!(
        run_entry(&resolved, &[cwd1.as_slice(), b"/sub"].concat(), &cwd1),
        Some(Entry::Missing)
    );
    assert_eq!(
        run_entry(
            &resolved,
            &[cwd1.as_slice(), b"/outside.fa"].concat(),
            &cwd1
        ),
        Some(Entry::Missing)
    );
    assert_eq!(
        run_entry(&resolved, b"/abs/outside.fa", &cwd1),
        Some(Entry::File)
    );
    assert_eq!(run_entry(&resolved, b"/abs", &cwd1), Some(Entry::Directory));
    assert_eq!(run_entry(&resolved, b"/abs/other.log", &cwd1), None);
    assert_eq!(run_entry(&resolved, b"/\0/up.fa", &cwd2), Some(Entry::File));
    // Direct `log_name` and `output_name` checks with an entry holding outside keys.
    let outside_entry = |path: &[u8], cwd: &[u8]| run_entry(&resolved, path, cwd);
    assert_eq!(
        log_name(b"/abs/outside.fa", Some(&outside_entry)),
        LogName::Existing
    );
    assert_eq!(
        output_name(b"/abs/outside.fa", &outside_entry),
        OutputName::Existing
    );
    assert_eq!(log_name(b"/abs", Some(&outside_entry)), LogName::Unopened);
    assert_eq!(
        log_name(b"/abs/run.log", Some(&outside_entry)),
        LogName::Unopened
    );
    assert_eq!(
        log_name(b"../up.fa", Some(&outside_entry)),
        LogName::Existing
    );
    assert_eq!(
        output_name(b"../up.fa", &outside_entry),
        OutputName::Existing
    );
    assert!(matches!(
        log_name(b"../other.log", Some(&outside_entry)),
        LogName::File { ref name, .. } if name == b"../other.log"
    ));
    // Registered names also reach `OpenFiles::create` through the same view.
    let files = OpenFiles::new(&RunFiles::Receiver(Arc::new(receiver)));
    assert!(matches!(files.create("q.fa"), Ok(Created::Refused(_))));
    assert!(matches!(files.create("in/a.fa"), Ok(Created::Refused(_))));
    assert_eq!(files.create("in").err(), Some(Unopened));
    assert!(matches!(files.create("in/out.txt"), Ok(Created::File(_))));
    assert!(matches!(
        files.log_target(Some(OsStr::new("/dev/stdout"))),
        Ok(ToolLog::Stdout)
    ));
    assert!(matches!(
        files.log_target(Some(OsStr::new("/dev/fd/1"))),
        Ok(ToolLog::Target(DiagTarget::Discard))
    ));
    assert!(files.log_target(Some(OsStr::new("q.fa"))).is_err());
}

/// A receiver's output files (`OpenFiles::create`) in the working directory of this test (the
/// crate's): a new name is a file in memory, which goes to the receiver under its whole name when
/// the run's files end; a file of the directory is refused and a folder of it does not open, nor
/// does a name in a folder that it does not have; the standard output's name is the run's, the
/// null device's writes are lost and `/dev/full`'s fail. `Disk` makes the file on the file system.
#[test]
fn a_receiver_s_output_file_goes_to_the_receiver_when_the_run_ends() {
    let kept = Arc::new(Kept::default());
    let files = OpenFiles::new(&RunFiles::Receiver(kept.clone()));
    let name = "ww3a-receiver-test.counts";
    assert!(!std::path::Path::new(name).exists());
    let Ok(Created::File(mut file)) = files.create(name) else {
        panic!("a file in memory");
    };
    file.write_all(b"counts ").unwrap();
    file.write_all(b"bytes").unwrap();
    file.flush().unwrap();
    drop(file);
    assert!(
        kept.files.lock().unwrap().is_empty(),
        "nothing before the run's files end"
    );
    assert!(matches!(
        files.create("Cargo.toml"),
        Ok(Created::Refused(_))
    ));
    assert_eq!(files.create("src").err(), Some(Unopened));
    assert_eq!(files.create("ww3a-no-such-folder/c").err(), Some(Unopened));
    assert!(matches!(files.create("/dev/stdout"), Ok(Created::Stdout)));
    let Ok(Created::File(mut null)) = files.create("/dev/null") else {
        panic!("the null device");
    };
    assert_eq!(null.write(b"lost").unwrap(), 4);
    let Ok(Created::File(mut full)) = files.create("/dev/full") else {
        panic!("/dev/full");
    };
    assert_eq!(
        full.write_all(b"x").unwrap_err().kind(),
        io::ErrorKind::StorageFull
    );
    drop((null, full, files));
    assert_eq!(
        *kept.files.lock().unwrap(),
        vec![(name.as_bytes().to_vec(), b"counts bytes".to_vec())],
        "one file, none of the devices"
    );
    assert!(!std::path::Path::new(name).exists(), "no file on the disk");
    // The command line's and the library's files are the file system's.
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "ww3a-disk-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("c.counts");
    let disk = OpenFiles::new(&RunFiles::Disk);
    let Ok(Created::File(mut file)) = disk.create(path.to_str().unwrap()) else {
        panic!("a file of the file system");
    };
    file.write_all(b"on disk").unwrap();
    drop(file);
    assert_eq!(std::fs::read(&path).unwrap(), b"on disk");
    assert_eq!(
        disk.create(directory.join("no/c.counts").to_str().unwrap())
            .err(),
        Some(Unopened)
    );
    std::fs::remove_dir_all(&directory).unwrap();
}

/// A receiver that takes its files as they are written (`FileReceiver::open`): the events of each
/// file in order, `open`, each write at its offset, each change of length, `end` (line NP WW3a,
/// D-WW3-6 and D-WW3-9).
#[derive(Default)]
struct Streamed(Arc<Mutex<Vec<String>>>);

struct StreamedFile(String, Arc<Mutex<Vec<String>>>);

impl FileStream for StreamedFile {
    fn write_at(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()> {
        let text = String::from_utf8_lossy(bytes);
        self.1
            .lock()
            .unwrap()
            .push(format!("{} @{offset} {text}", self.0));
        Ok(())
    }

    fn set_len(&mut self, length: u64) -> io::Result<()> {
        self.1
            .lock()
            .unwrap()
            .push(format!("{} length {length}", self.0));
        Ok(())
    }
}

impl Drop for StreamedFile {
    fn drop(&mut self) {
        self.1.lock().unwrap().push(format!("{} end", self.0));
    }
}

impl FileReceiver for Streamed {
    fn receive(&self, name: &[u8], _bytes: &[u8]) -> io::Result<()> {
        panic!("a whole file: {}", String::from_utf8_lossy(name));
    }

    fn open(&self, name: &[u8]) -> Option<Box<dyn FileStream>> {
        let name = String::from_utf8_lossy(name).into_owned();
        self.0.lock().unwrap().push(format!("{name} opened"));
        Some(Box::new(StreamedFile(name, self.0.clone())))
    }
}

/// A receiver that opens its files takes the bytes of an output file and of the file of
/// `-logfile` as they are written, each at its offset, and the files end when the run's files end
/// (the run's `OpenFiles` and every opening dropped), in the order they opened; the names that
/// open no file stay as they are. A receiver that does not open them gets each whole, and the
/// search's `RunFiles::log_target` keeps its file in memory (the files of a search, unchanged).
#[test]
fn a_receiver_that_opens_its_files_takes_them_as_they_are_written() {
    let streamed = Streamed::default();
    let events = streamed.0.clone();
    let files = OpenFiles::new(&RunFiles::Receiver(Arc::new(streamed)));
    let Ok(Created::File(mut file)) = files.create("ww3a-streamed.counts") else {
        panic!("a streamed file");
    };
    file.write_all(b"one").unwrap();
    let Ok(ToolLog::Target(log)) = files.log_target(Some(OsStr::new("ww3a-streamed.log"))) else {
        panic!("a streamed log file");
    };
    assert!(matches!(log, DiagTarget::Stream(_)));
    assert_eq!(log.place(), Place::LogFile);
    log.write_text(b"message", &mut io::sink());
    file.write_all(b"two").unwrap();
    drop(file);
    assert!(matches!(files.create("/dev/stdout"), Ok(Created::Stdout)));
    assert_eq!(files.create("src").err(), Some(Unopened));
    assert!(matches!(
        files.log_target(Some(OsStr::new("-"))),
        Ok(ToolLog::Target(DiagTarget::Stderr))
    ));
    assert_eq!(events.lock().unwrap().len(), 5, "no end before the run's");
    drop((log, files));
    assert_eq!(
        *events.lock().unwrap(),
        [
            "ww3a-streamed.counts opened",
            "ww3a-streamed.counts @0 one",
            "ww3a-streamed.log opened",
            "ww3a-streamed.log @0 message",
            "ww3a-streamed.counts @3 two",
            "ww3a-streamed.counts end",
            "ww3a-streamed.log end",
        ]
    );
    assert!(!std::path::Path::new("ww3a-streamed.counts").exists());
    // The default: whole files, from memory; a search's log in memory.
    let kept = RunFiles::Receiver(Arc::new(Kept::default()));
    assert!(matches!(
        kept.log_target(Some(OsStr::new("ww3a-kept.log"))),
        Ok((DiagTarget::Memory(_), Some(ref name))) if name == b"ww3a-kept.log"
    ));
}

/// A run that opens one file twice (`-logfile x` and then `-out x`, or `./x`) has one file, which
/// the system shares between the openings (POSIX; NCBI BLAST+ 2.17.0 and LOSAT's command line on
/// Linux, line NP WW3a's review M2, D-WW3-9, 2026-11-23): the output's opening empties it and
/// writes from its own position, from 0, and the log's writes go to the file's end, so a message
/// posted before the output's bytes is written over and one posted after them follows them. The
/// receiver that streams gets the change of length and each write at its offset; the one that
/// takes files whole gets one file with the same bytes.
#[test]
fn a_file_opened_twice_is_one_file_with_the_system_s_offsets() {
    let streamed = Streamed::default();
    let events = streamed.0.clone();
    let files = OpenFiles::new(&RunFiles::Receiver(Arc::new(streamed)));
    let Ok(ToolLog::Target(log)) = files.log_target(Some(OsStr::new("ww3a-twice.out"))) else {
        panic!("a streamed log file");
    };
    log.write_text(b"pass 1\n", &mut io::sink());
    let Ok(Created::File(mut out)) = files.create("./ww3a-twice.out") else {
        panic!("the same file");
    };
    log.write_text(b"pass 2\n", &mut io::sink());
    out.write_all(b"counts:").unwrap();
    out.write_all(b"0123456789").unwrap();
    log.write_text(b"done\n", &mut io::sink());
    drop((out, log, files));
    assert_eq!(
        *events.lock().unwrap(),
        [
            "ww3a-twice.out opened",
            "ww3a-twice.out @0 pass 1\n",
            "ww3a-twice.out length 0",
            "ww3a-twice.out @0 pass 2\n",
            "ww3a-twice.out @0 counts:",
            "ww3a-twice.out @7 0123456789",
            "ww3a-twice.out @17 done\n",
            "ww3a-twice.out end",
        ]
    );
    // Whole: the same bytes, once.
    let kept = Arc::new(Kept::default());
    let files = OpenFiles::new(&RunFiles::Receiver(kept.clone()));
    let Ok(ToolLog::Target(log)) = files.log_target(Some(OsStr::new("ww3a-twice.out"))) else {
        panic!("a log file");
    };
    log.write_text(b"pass 1\n", &mut io::sink());
    let Ok(Created::File(mut out)) = files.create("ww3a-twice.out") else {
        panic!("the same file");
    };
    log.write_text(b"a long message\n", &mut io::sink());
    out.write_all(b"short").unwrap();
    drop((out, log, files));
    assert_eq!(
        *kept.files.lock().unwrap(),
        vec![(b"ww3a-twice.out".to_vec(), b"shortg message\n".to_vec())]
    );
}

/// The path that names one file of a run (`file_key`): every name of the same path in the
/// working directory, by the output's reading and by `-logfile`'s (which NCBI trims), and the
/// same place above the working directory whatever the name's depth; other places differ.
#[test]
fn the_names_of_one_path_are_one_file() {
    let out = |name: &str| output_key(name.as_bytes()).unwrap();
    let log = |name: &str| log_key(name.as_bytes()).unwrap();
    for same in ["./x", "sub/../x", "a/b/../../x", "/proc/self/cwd/x", ".//x"] {
        assert_eq!(out(same), out("x"), "{same}");
    }
    assert_eq!(log("x"), out("x"));
    assert_eq!(log(" ./x"), out("x"));
    assert_eq!(out("../x"), out("a/../../x"));
    assert_eq!(out("../../x"), out("../a/../../x"));
    assert_ne!(out("../x"), out("x"));
    assert_ne!(out("../x"), out("../../x"));
    assert_ne!(out("sub/x"), out("x"));
    assert_eq!(out("/tmp/x"), log("/tmp/x"));
    assert_eq!(out("/var/lock/x"), out("/run/lock/x"));
    assert_ne!(out("/tmp/x"), out("x"));
}

/// Registered search inputs outside the working directory (`run_entries`, `run_entry`; line WQA
/// D19 round R1): `log_target` and `OpenFiles::create` see them by `file_key` inside and outside
/// the working directory, treat their parent folders as present even where `directory_at` returns
/// `Directory::None` (`/data` for `/data/q.fa`), keep sibling paths under root directories
/// unopened (`/data/run.log`) and writable user directories open (`/tmp/run.log`, `../run.log`),
/// and ignore non-file directory placeholders (`.`, `..`, `/tmp`, `/`).
#[test]
fn outside_cwd_registered_files_are_seen_by_log_and_output_names() {
    let receiver = Arc::new(Kept::with_registered(&[
        "/tmp/q.fa",
        "../q.fa",
        "/home/me/q.fa",
        "/tmp/s.fa",
        "/data/q.fa",
        "/data/sub/a.fa",
        "../parent_dir/q.fa",
        "/var/lock/l.fa",
        ".",
        "..",
        "/tmp",
        "/",
    ]));
    let run_files = RunFiles::Receiver(receiver.clone());
    let open_files = OpenFiles::new(&run_files);
    // Case 2 outside cwd: refused by both `log_target` and `OpenFiles::create`.
    for existing in [
        "/tmp/q.fa",
        "/tmp/x/../q.fa",
        "../q.fa",
        "./../q.fa",
        "/home/me/q.fa",
        "/tmp/s.fa",
        "/data/q.fa",
        "/data/sub/a.fa",
        "/run/lock/l.fa",
        "/var/lock/l.fa",
    ] {
        assert!(
            run_files.log_target(Some(OsStr::new(existing))).is_err(),
            "log_target should refuse {existing:?}"
        );
        assert!(
            matches!(open_files.create(existing), Ok(Created::Refused(_))),
            "create should refuse {existing:?}"
        );
    }
    // Trailing slash on a registered outside file is stripped by NCBI's `log_file_name`, so
    // `log_target` refuses it too.
    assert!(run_files
        .log_target(Some(OsStr::new("/tmp/q.fa/")))
        .is_err());
    // Ancestor folders and file-as-folder paths outside cwd stay unopened (stderr / Unopened).
    for unopened in [
        "/data",
        "/data/",
        "/data/sub",
        "../parent_dir",
        "/tmp/q.fa/run.log",
        "/data/q.fa/run.log",
    ] {
        assert!(
            matches!(
                run_files.log_target(Some(OsStr::new(unopened))),
                Ok((DiagTarget::Stderr, None))
            ),
            "log_target should keep {unopened:?} on stderr"
        );
        assert_eq!(
            open_files.create(unopened).err(),
            Some(Unopened),
            "create should not open {unopened:?}"
        );
    }
    // Controls: sibling files in writable outside dirs and cwd open; `/data/run.log` stays
    // unopened because `/data` is outside `user_makes`.
    for (raw, want) in [
        ("/tmp/run.log", "/tmp/run.log"),
        ("../run.log", "../run.log"),
        ("run.log", "run.log"),
    ] {
        assert!(
            matches!(
                run_files.log_target(Some(OsStr::new(raw))),
                Ok((DiagTarget::Memory(_), Some(ref name))) if name == want.as_bytes()
            ),
            "log_target control {raw:?}"
        );
        assert!(
            matches!(open_files.create(raw), Ok(Created::File(_))),
            "create control {raw:?}"
        );
    }
    assert!(matches!(
        run_files.log_target(Some(OsStr::new("/data/run.log"))),
        Ok((DiagTarget::Stderr, None))
    ));
    assert_eq!(open_files.create("/data/run.log").err(), Some(Unopened));
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6514-6515
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:745-750
    // `CDirEntry::ConvertToOSPath` calls `NormalizePath` without following symlinks, turning
    // `/proc/self/cwd/../q.fa` into `/proc/self/q.fa`, which fails to open in procfs so NCBI
    // keeps messages on stderr and never opens `../q.fa`.
    assert!(matches!(
        run_files.log_target(Some(OsStr::new("/proc/self/cwd/../q.fa"))),
        Ok((DiagTarget::Stderr, None))
    ));
}

/// Missing input files listed in `FileReceiver::read_files` / `output_read_files` (line WQA D19
/// round R2) refuse `-logfile` / `-out` when `log_name` / `output_name` would otherwise create
/// that file, without making missing parent folders present.
#[test]
fn read_files_refuse_matching_writes_without_creating_folders() {
    let receiver = Arc::new(Kept::with_read_files(
        &["missing.fa", "missing_sub/x.fa", "ids.txt"],
        &["missing.fa"],
    ));
    let run_files = RunFiles::Receiver(receiver);
    let open_files = OpenFiles::new(&run_files);

    // `-logfile` naming a missing read input (also via `./`) is refused.
    assert!(run_files
        .log_target(Some(OsStr::new("missing.fa")))
        .is_err());
    assert!(run_files
        .log_target(Some(OsStr::new("./missing.fa")))
        .is_err());
    assert!(run_files.log_target(Some(OsStr::new("ids.txt"))).is_err());

    // Missing parent folder of a missing input does NOT become present.
    assert!(matches!(
        run_files.log_target(Some(OsStr::new("missing_sub/run.log"))),
        Ok((DiagTarget::Stderr, None))
    ));
    assert!(matches!(
        run_files.log_target(Some(OsStr::new("missing_sub/x.fa"))),
        Ok((DiagTarget::Stderr, None))
    ));

    // `-out` refuses only `output_read_files` (`missing.fa`), not `ids.txt`.
    assert!(matches!(
        open_files.create("missing.fa"),
        Ok(Created::Refused(_))
    ));
    assert!(matches!(
        open_files.create("./missing.fa"),
        Ok(Created::Refused(_))
    ));
    assert!(open_files
        .refused_output(OsStr::new("missing.fa"))
        .is_some());
    assert!(matches!(open_files.create("ids.txt"), Ok(Created::File(_))));
    assert!(open_files.refused_output(OsStr::new("ids.txt")).is_none());
}

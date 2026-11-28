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
/// command line on Linux opens them for writing, the name as typed (`RunFiles::strategy_target`),
/// on the host and as the user of the `-logfile` model, against the run's files (line RG item
/// F1w-r2: the view of `-logfile` and a tool's output, line WQA D19; here registered names, and
/// this test's working directory, the crate's, for `src` and `Cargo.toml`): a file of the run where
/// the user makes files, also in a folder of the run's files; NCBI's error where the command line
/// cannot open the name (a folder that the run's files do not have or that is a file of them, a
/// folder of them, a name ending in `/`, a name in a directory of root's, a component of 256
/// bytes); LOSAT Web's refusal of a file of the run's files; the standard output's refusal, no
/// file for the null device and the standard input, and LOSAT Web's refusal of the standard error,
/// of the devices whose writes fail or wait for a reader, and of procfs's and mqueue's files.
#[test]
fn a_strategy_name_is_read_as_the_command_line_opens_it_on_linux() {
    use StrategyTarget as T;
    // Registered names stand for the folders and files of the command line's directory (as in
    // `a_receiver_reads_a_log_name_as_the_command_line_on_linux`).
    let files = RunFiles::Receiver(Arc::new(Kept::with_registered(&[
        "out/keep", "a/b/keep", "in/q.fa",
    ])));
    let target = |name: &str| files.strategy_target(OsStr::new(name), &[]);
    let file = |name: &str| T::File(name.as_bytes().to_vec());
    let long_component = "a".repeat(256);
    let component_255 = "a".repeat(255);
    for (name, read) in [
        ("s.asn", file("s.asn")),
        ("out/s.asn", file("out/s.asn")),
        ("./s.asn", file("s.asn")),
        (" s.asn", file(" s.asn")),
        ("../s.asn", file("../s.asn")),
        ("a/b/../s.asn", file("a/s.asn")),
        ("a/s.asn", file("a/s.asn")),
        ("in/s.asn", file("in/s.asn")),
        ("in/../s.asn", file("s.asn")),
        ("src/s.asn", file("src/s.asn")),
        ("/tmp/s.asn", file("/tmp/s.asn")),
        ("/dev/shm/s.asn", file("/dev/shm/s.asn")),
        ("/var/lock/s.asn", file("/var/lock/s.asn")),
        // A `..` after a link steps back from the link's target (NCBI 2.17.0 blastp makes
        // `/tmp/s.asn` for `/var/run/../tmp/s.asn`; `/var/run` is `/run`).
        ("/var/run/../tmp/s.asn", file("/tmp/s.asn")),
        ("/var/lock/../shm/s.asn", file("/run/shm/s.asn")),
        ("/proc/self/root/../tmp/s.asn", file("/tmp/s.asn")),
        (
            "/proc/self/root/tmp/s.asn",
            file("/proc/self/root/tmp/s.asn"),
        ),
        ("/home/me/s.asn", file("/home/me/s.asn")),
        ("/proc/self/cwd/s.asn", file("s.asn")),
        ("/proc/self/cwd/out/s.asn", file("out/s.asn")),
        ("/proc/thread-self/cwd/s.asn", file("s.asn")),
        (component_255.as_str(), file(&component_255)),
        ("", T::Inaccessible),
        (".", T::Inaccessible),
        ("..", T::Inaccessible),
        ("out", T::Inaccessible),
        ("out/", T::Inaccessible),
        ("out/.", T::Inaccessible),
        ("in", T::Inaccessible),
        ("a/b", T::Inaccessible),
        ("src", T::Inaccessible),
        ("nodir/s.asn", T::Inaccessible),
        ("nodir/../s.asn", T::Inaccessible),
        ("a/c/../s.asn", T::Inaccessible),
        ("in/q.fa/s.asn", T::Inaccessible),
        ("in/q.fa/", T::Inaccessible),
        ("Cargo.toml/s.asn", T::Inaccessible),
        ("/proc/self/cwd/nodir/s.asn", T::Inaccessible),
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
        ("in/q.fa", T::Existing),
        ("./in/q.fa", T::Existing),
        ("a/../in/q.fa", T::Existing),
        ("/proc/self/cwd/in/q.fa", T::Existing),
        ("out/keep", T::Existing),
        ("Cargo.toml", T::Existing),
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
        // The *start*'s name of the file (`losat_web2_start`, before the run's files) is the run's.
        if let T::File(made) = &read {
            assert_eq!(
                RunFiles::output_host_name(OsStr::new(name)).as_ref(),
                Some(made),
                "{name:?}"
            );
        }
    }
    // Without those folders, as in a new directory that holds only the run's other files, the
    // command line cannot open a name in them (NCBI's `File is not accessible`); the *start* still
    // names the file that a run with the folders makes.
    let bare = RunFiles::Receiver(Arc::new(Kept::default()));
    for (name, made) in [
        ("out/s.asn", "out/s.asn"),
        ("a/b/../s.asn", "a/s.asn"),
        ("in/s.asn", "in/s.asn"),
        ("/proc/self/cwd/out/s.asn", "out/s.asn"),
    ] {
        assert_eq!(
            bare.strategy_target(OsStr::new(name), &[]),
            T::Inaccessible,
            "{name:?}"
        );
        assert_eq!(
            RunFiles::output_host_name(OsStr::new(name)),
            Some(made.as_bytes().to_vec()),
            "{name:?}"
        );
    }
    // A name that makes no file in any run's files has no name in the *start*.
    for name in [
        "",
        ".",
        "..",
        "s.asn/",
        "/tmp",
        "/etc/s.asn",
        "//s.asn",
        "/dev/null",
        "/dev/stdout",
        "/dev/stderr",
        "/dev/full",
        "/proc/self/comm",
        "/dev/mqueue/q",
    ] {
        assert_eq!(
            RunFiles::output_host_name(OsStr::new(name)),
            None,
            "{name:?}"
        );
    }
}

/// A strategy name against the run's files outside the working directory and the files that the
/// search reads (line RG item F1w-r2, as `outside_cwd_registered_files_are_seen_by_log_and_output_names`
/// and `read_files_refuse_matching_writes_without_creating_folders` read `-logfile` and `-out`): a
/// registered file outside the working directory is a file of the run (LOSAT Web's refusal, or the
/// refusal that names its option) and its folder is there; a file that the search reads, present or
/// missing, is LOSAT Web's refusal, without making its folder present.
#[test]
fn a_strategy_name_sees_registered_and_read_files_outside_the_view() {
    use StrategyTarget as T;
    let files = RunFiles::Receiver(Arc::new(Kept {
        registered: vec![b"/data/q.fa".to_vec(), b"../up.fa".to_vec()],
        read: vec![b"tfp.pal".to_vec(), b"missing_sub/x.fa".to_vec()],
        ..Kept::default()
    }));
    let query = [RunFileName {
        option: "query",
        name: OsStr::new("/data/q.fa"),
    }];
    for (name, others, read) in [
        ("/data/q.fa", &[][..], T::Existing),
        ("/data/q.fa", &query[..], T::OtherFile("query")),
        ("/data/./q.fa", &query[..], T::OtherFile("query")),
        ("/data", &[], T::Inaccessible),
        ("/data/", &[], T::Inaccessible),
        ("/data/s.asn", &[], T::Inaccessible),
        ("/data/q.fa/s.asn", &[], T::Inaccessible),
        // `/data` is there, the folder of a registered file: a name may pass it on its way out.
        ("/data/../tmp/s.asn", &[], T::File(b"/tmp/s.asn".to_vec())),
        ("../up.fa", &[], T::Existing),
        ("../s.asn", &[], T::File(b"../s.asn".to_vec())),
        ("tfp.pal", &[], T::Existing),
        ("./tfp.pal", &[], T::Existing),
        ("tfp.pin", &[], T::File(b"tfp.pin".to_vec())),
        ("missing_sub/x.fa", &[], T::Inaccessible),
        ("missing_sub/s.asn", &[], T::Inaccessible),
    ] {
        assert_eq!(
            files.strategy_target(OsStr::new(name), others),
            read,
            "{name:?}"
        );
        // The *start* names each file that the run makes as the run does (`/tmp/s.asn` too,
        // though `/data` is no directory of the model's).
        if let T::File(made) = &read {
            assert_eq!(
                RunFiles::output_host_name(OsStr::new(name)).as_ref(),
                Some(made),
                "{name:?}"
            );
        }
    }
}

/// `-logfile` and `-export_search_strategy` of one name with a `..` after a link are two files, as
/// on the command line (NCBI 2.17.0 blastp with both `/var/run/../tmp/X`: the log `/var/tmp/X`, as
/// NCBI reads `-logfile` lexically, and the strategy file `/tmp/X`, which the system opens as
/// typed, `/var/run` being `/run`): neither is the other's, each under its own name in the run and
/// in the *start*. Without the link, one name is one file (the strategy file is the log's).
#[test]
fn a_log_and_a_strategy_of_one_name_after_a_link_are_two_files() {
    let files = RunFiles::Receiver(Arc::new(Kept::default()));
    for (name, log_name, strategy) in [
        (
            "/var/run/../tmp/X",
            "/var/tmp/X",
            StrategyTarget::File(b"/tmp/X".to_vec()),
        ),
        ("/tmp/X", "/tmp/X", StrategyTarget::OtherFile("logfile")),
    ] {
        let log = [RunFileName {
            option: "logfile",
            name: OsStr::new(name),
        }];
        assert_eq!(
            files.strategy_target(OsStr::new(name), &log),
            strategy,
            "{name:?}"
        );
        let Ok((DiagTarget::Memory(_), Some(logged))) = files.log_target(Some(OsStr::new(name)))
        else {
            panic!("{name:?}: a -logfile file");
        };
        assert_eq!(logged, log_name.as_bytes(), "{name:?}");
        assert_eq!(
            RunFiles::log_host_name(OsStr::new(name)),
            Some(logged),
            "{name:?}"
        );
    }
    assert_eq!(
        RunFiles::output_host_name(OsStr::new("/var/run/../tmp/X")),
        Some(b"/tmp/X".to_vec())
    );
}

/// One view of the run's files for every name that the run writes (line RG item F1w-r2): for a
/// name that NCBI's reading of `-logfile` takes as it is typed (no leading space, no final `/`, no
/// `..` after a link, none of NCBI's special names), the strategy file, the `-logfile` file and a
/// tool's output are each a new file of the run, no file (NCBI's error, the messages on the
/// standard error, a tool's open that fails), or LOSAT Web's refusal of a file of the run, in the
/// same run's files.
#[test]
fn a_strategy_name_a_log_name_and_an_output_name_are_one_reading() {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Read {
        File,
        Unopened,
        Refused,
    }
    let kept = Arc::new(Kept {
        registered: ["out/keep", "a/b/keep", "in/q.fa", "/data/q.fa", "../up.fa"]
            .iter()
            .map(|name| name.as_bytes().to_vec())
            .collect(),
        read: vec![b"tfp.pal".to_vec()],
        out_read: vec![b"tfp.pal".to_vec()],
        ..Kept::default()
    });
    let run_files = RunFiles::Receiver(kept);
    let open_files = OpenFiles::new(&run_files);
    for name in [
        "s.asn",
        "./s.asn",
        "out/s.asn",
        "a/b/../s.asn",
        "a/../s.asn",
        "in/s.asn",
        "src/s.asn",
        "../s.asn",
        "/proc/self/cwd/out/s.asn",
        "/tmp/s.asn",
        "nodir/s.asn",
        "nodir/../s.asn",
        "a/c/../s.asn",
        "in/q.fa/s.asn",
        "Cargo.toml/s.asn",
        "/proc/self/cwd/nodir/s.asn",
        "out",
        "in",
        "a/b",
        "src",
        "/data",
        "/data/s.asn",
        "/etc/s.asn",
        "in/q.fa",
        "./in/q.fa",
        "a/../in/q.fa",
        "/proc/self/cwd/in/q.fa",
        "out/keep",
        "Cargo.toml",
        "/data/q.fa",
        "../up.fa",
        "tfp.pal",
        "./tfp.pal",
    ] {
        let strategy_target = run_files.strategy_target(OsStr::new(name), &[]);
        let strategy = match &strategy_target {
            StrategyTarget::File(_) => Read::File,
            StrategyTarget::Inaccessible => Read::Unopened,
            StrategyTarget::Existing => Read::Refused,
            other => panic!("{name:?}: strategy {other:?}"),
        };
        let log_target = run_files.log_target(Some(OsStr::new(name)));
        let log = match &log_target {
            Ok((DiagTarget::Memory(_), Some(_))) => Read::File,
            Ok((DiagTarget::Stderr, None)) => Read::Unopened,
            Err(_) => Read::Refused,
            Ok(_) => panic!("{name:?}: log target"),
        };
        // A file of both is the command line's one file, under one host name.
        if let (StrategyTarget::File(made), Ok((_, Some(logged)))) = (&strategy_target, &log_target)
        {
            assert_eq!(made, logged, "{name:?}");
        }
        let output = match open_files.create(name) {
            Ok(Created::File(_)) => Read::File,
            Err(Unopened) => Read::Unopened,
            Ok(Created::Refused(_)) => Read::Refused,
            Ok(Created::Stdout) => panic!("{name:?}: output Stdout"),
        };
        assert_eq!((strategy, log), (output, output), "{name:?}");
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
    let file = |name: &str| T::File(name.as_bytes().to_vec());
    for (name, others, read) in [
        ("q.fa", vec![("query", "q.fa")], T::OtherFile("query")),
        ("./q.fa", vec![("query", "q.fa")], T::OtherFile("query")),
        (
            "/proc/self/cwd/q.fa",
            vec![("query", "q.fa")],
            T::OtherFile("query"),
        ),
        ("q.fa", vec![("query", "../x/../w/q.fa")], file("q.fa")),
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
        (
            "../a/q.fa",
            vec![("query", "../../a/q.fa")],
            file("../a/q.fa"),
        ),
        ("q.fa.asn", vec![("query", "q.fa")], file("q.fa.asn")),
        (
            "/tmp/s.fa",
            vec![("query", "q.fa"), ("subject", "/tmp/s.fa")],
            T::OtherFile("subject"),
        ),
        ("./-", vec![("query", "-")], file("-")),
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
        ("run.log", vec![("logfile", "../run.log")], file("run.log")),
        ("./-", vec![("logfile", "-")], file("-")),
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
        ("/var/run/../tmp/c.counts", "/tmp/c.counts"),
        ("/var/lock/../shm/c.counts", "/run/shm/c.counts"),
        ("/proc/self/root/../tmp/c.counts", "/tmp/c.counts"),
        (
            "/proc/self/root/tmp/c.counts",
            "/proc/self/root/tmp/c.counts",
        ),
    ] {
        match name(given) {
            OutputName::File { name: got, .. } => {
                assert_eq!(got, expected_host_name.as_bytes(), "{given:?}");
            }
            other => panic!("{given:?}: expected File, got {other:?}"),
        }
    }
    for (given, expected) in [
        ("nodir/c.counts", OutputName::Unopened),
        ("nodir/../c.counts", OutputName::Unopened),
        ("in.fa/c.counts", OutputName::Unopened),
        ("/proc/self/cwd/nodir/c.counts", OutputName::Unopened),
        ("sub", OutputName::Unopened),
        ("sub/", OutputName::Unopened),
        (".", OutputName::Unopened),
        ("..", OutputName::Unopened),
        ("sub/..", OutputName::Unopened),
        ("/proc/self/cwd", OutputName::Unopened),
        ("/tmp", OutputName::Unopened),
        ("/c.counts", OutputName::Unopened),
        ("/etc/c.counts", OutputName::Unopened),
        ("/data/c.counts", OutputName::Unopened),
        ("/root/../tmp/c.counts", OutputName::Unopened),
        ("/proc/self/status", OutputName::Unopened),
        ("/dev/fd/3", OutputName::Unopened),
        ("/dev/kmsg", OutputName::Unopened),
        ("", OutputName::Unopened),
        (long.as_str(), OutputName::Unopened),
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
        ("/dev/stderr", OutputName::StandardError),
        ("/dev/fd/2", OutputName::StandardError),
        ("/proc/self/fd/2", OutputName::StandardError),
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
        // NCBI reads `-logfile` lexically (a `..` after a link is not the link's parent).
        ("/var/run/../tmp/run.log", Some("/var/tmp/run.log")),
        // A file that a run makes beside a registered `/data/q.fa` (`/data` is no directory of the
        // model's).
        ("/data/../tmp/run.log", Some("/tmp/run.log")),
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
    // A `..` after a link steps back from the link's target: the file is `/tmp/...`.
    let Ok(Created::File(mut linked)) = files.create("/var/run/../tmp/ww3a-linked.counts") else {
        panic!("a file in memory");
    };
    linked.write_all(b"linked").unwrap();
    drop((null, full, linked, files));
    assert_eq!(
        *kept.files.lock().unwrap(),
        vec![
            (name.as_bytes().to_vec(), b"counts bytes".to_vec()),
            (b"/tmp/ww3a-linked.counts".to_vec(), b"linked".to_vec()),
        ],
        "two files, none of the devices"
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
    // Controls: sibling files in writable outside dirs and cwd open, also past the folder of a
    // registered file (`/data/../tmp/run.log`), each under the *start*'s name of it
    // (`log_host_name`); `/data/run.log` stays unopened because `/data` is outside `user_makes`.
    for (raw, want) in [
        ("/tmp/run.log", "/tmp/run.log"),
        ("../run.log", "../run.log"),
        ("run.log", "run.log"),
        ("/data/../tmp/run.log", "/tmp/run.log"),
    ] {
        assert!(
            matches!(
                run_files.log_target(Some(OsStr::new(raw))),
                Ok((DiagTarget::Memory(_), Some(ref name))) if name == want.as_bytes()
            ),
            "log_target control {raw:?}"
        );
        assert_eq!(
            RunFiles::log_host_name(OsStr::new(raw)).as_deref(),
            Some(want.as_bytes()),
            "log_host_name control {raw:?}"
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

    // The check's refusal (`refused_output`) is the run's (`create`) for every name that `create`
    // refuses: the standard error, a device or a file of procfs too (line NP WW4c, 2026-11-28).
    for name in ["missing.fa", "/dev/stderr", "/dev/fd/2", "/proc/self/comm"] {
        let Ok(Created::Refused(reason)) = open_files.create(name) else {
            panic!("{name}: create refuses it");
        };
        assert_eq!(
            open_files.refused_output(OsStr::new(name)),
            Some(reason),
            "{name}"
        );
    }
    for name in [
        "/dev/stdout",
        "/dev/null",
        "/dev/full",
        "new.out",
        "missing_sub/x.out",
    ] {
        assert!(
            open_files.refused_output(OsStr::new(name)).is_none(),
            "{name}"
        );
    }
}

/// A receiver whose run writes its working directory through the file system (makeblastdb of
/// LOSAT Web, `FileReceiver::writes_working_directory`), taking its other files as they are
/// written (`Streamed`).
#[derive(Default)]
struct WritesDirectory(Streamed);

impl FileReceiver for WritesDirectory {
    fn receive(&self, name: &[u8], _bytes: &[u8]) -> io::Result<()> {
        panic!("a whole file: {}", String::from_utf8_lossy(name));
    }

    fn open(&self, name: &[u8]) -> Option<Box<dyn FileStream>> {
        self.0.open(name)
    }

    fn writes_working_directory(&self) -> bool {
        true
    }
}

/// The `-logfile` of a run that writes its working directory (makeblastdb of LOSAT Web; line NP
/// WW4b1 R2, finding 1, 2026-11-28) follows the device model of docs/web/abi_v2.md §4, as
/// windowmasker's and blast_formatter's: the diagnostics' target and the log's own opening
/// (`OpenFiles::log_output`, makeblastdb's `AsOutputFile`) read the name through one owner
/// (`log_name`, `output_name`). The null device loses both; `/dev/stdout` and `/dev/stderr` are
/// the run's standard streams; `/dev/full` fails the log's writes for want of space; a file
/// outside the working directory (`/tmp/mk.log`) is a file of the run that the receiver takes as
/// it is written, the two openings appending to one file; a name in a folder that the directory
/// does not have and a folder leave the diagnostics on the standard error and do not open the
/// log; a device that LOSAT Web does not write is refused, in the run and in a check of the
/// words, before anything opens.
#[test]
fn a_run_that_writes_its_directory_opens_its_log_by_the_device_model() {
    let receiver = WritesDirectory::default();
    let events = receiver.0 .0.clone();
    let run_files = RunFiles::Receiver(Arc::new(receiver));
    let files = OpenFiles::new(&run_files);
    let target = |name: &str| match files.log_target(Some(OsStr::new(name))) {
        Ok(ToolLog::Target(target)) => Some(target),
        Ok(ToolLog::Stdout) => None,
        Err(reason) => panic!("{name}: {reason}"),
    };
    let output = |name: &str| files.log_output(OsStr::new(name));
    // The null device and its spellings: nothing is written anywhere.
    for name in ["/dev/null", "//dev/null", "/dev/zero"] {
        assert!(
            matches!(output(name), Ok(LogOutput::File(OpenedFile::Discard))),
            "{name}"
        );
    }
    assert!(matches!(target("/dev/zero"), Some(DiagTarget::Discard)));
    // The standard streams.
    assert!(target("/dev/stdout").is_none());
    assert!(matches!(output("/dev/stdout"), Ok(LogOutput::Stdout)));
    assert!(matches!(output("/dev/fd/1"), Ok(LogOutput::Stdout)));
    assert!(matches!(
        target("/dev/stderr"),
        Some(DiagTarget::StderrFile)
    ));
    for name in ["/dev/stderr", "/dev/fd/2", "/proc/self/fd/2"] {
        assert!(matches!(output(name), Ok(LogOutput::Stderr)), "{name}");
    }
    // `/dev/full`: the log opens, and its writes fail as the system's.
    let Ok(LogOutput::File(mut full)) = output("/dev/full") else {
        panic!("/dev/full opens");
    };
    assert!(full.write_all(b"line\n").is_err());
    // A file outside the working directory: one file of the run, the diagnostics' opening and
    // the log's appending to it in the order they write.
    let outside = "/tmp/ww4b1-r2-model-mk.log";
    let diagnostics = target(outside).expect("a file");
    assert!(matches!(diagnostics, DiagTarget::Stream(_)));
    let Ok(LogOutput::File(mut log)) = output(outside) else {
        panic!("{outside} opens");
    };
    diagnostics.write_text(b"diagnostic\n", &mut io::sink());
    log.write_all(b"log\n").unwrap();
    drop((diagnostics, log));
    assert!(!std::path::Path::new(outside).exists(), "not on the disk");
    // Names that do not open: the diagnostics stay on the standard error and the log does not
    // open (`src` is a folder of the crate's directory, the tests' working directory).
    for name in ["ww4b1-r2-nodir/mk.log", "src", "src/"] {
        assert!(matches!(target(name), Some(DiagTarget::Stderr)), "{name}");
        assert_eq!(output(name).err(), Some(Unopened), "{name}");
    }
    assert!(!std::path::Path::new("ww4b1-r2-nodir").exists());
    // A device that LOSAT Web does not write: refused by the run and by the check alike.
    for name in ["/dev/ptmx", "/proc/self/comm"] {
        assert!(matches!(output(name), Ok(LogOutput::Refused(_))), "{name}");
        assert!(
            run_files.refused_log_output(OsStr::new(name)).is_some(),
            "{name}"
        );
    }
    for name in [
        "/dev/null",
        "/dev/stdout",
        "/dev/stderr",
        "/dev/full",
        outside,
        "src",
    ] {
        assert_eq!(
            run_files.refused_log_output(OsStr::new(name)),
            None,
            "{name}"
        );
    }
    assert_eq!(
        RunFiles::Disk.refused_log_output(OsStr::new("/dev/ptmx")),
        None
    );
    drop(files);
    assert_eq!(
        *events.lock().unwrap(),
        [
            "/tmp/ww4b1-r2-model-mk.log opened",
            "/tmp/ww4b1-r2-model-mk.log @0 diagnostic\n",
            "/tmp/ww4b1-r2-model-mk.log @11 log\n",
            "/tmp/ww4b1-r2-model-mk.log end",
        ]
    );
}

/// A check of the words of a run that writes its working directory (`RunFiles::log_target`, the
/// `losat_web2_validate` of makeblastdb; line NP WW4b1 R2, finding 2, 2026-11-28) chooses the
/// diagnostics' target as the run does (`OpenFiles::log_target`), from the same reading of the
/// name, and makes nothing: a file of the directory that the run opens is in memory (not made), a
/// name in a missing folder and a folder leave the messages on the standard error, as the run's
/// opening does.
#[test]
fn a_check_of_a_run_that_writes_its_directory_chooses_the_run_s_target() {
    let run_files = RunFiles::Receiver(Arc::new(WritesDirectory::default()));
    let checked = |name: &str| run_files.log_target(Some(OsStr::new(name))).unwrap().0;
    assert!(matches!(
        checked("ww4b1-r2-check.log"),
        DiagTarget::Memory(_)
    ));
    assert!(
        !std::path::Path::new("ww4b1-r2-check.log").exists(),
        "not made"
    );
    for name in ["ww4b1-r2-nodir/mk.log", "src", "src/"] {
        assert!(matches!(checked(name), DiagTarget::Stderr), "{name}");
    }
    assert!(matches!(checked("/dev/stderr"), DiagTarget::StderrFile));
    assert!(matches!(checked("/dev/null"), DiagTarget::Discard));
    assert!(matches!(
        checked("/tmp/ww4b1-r2-check.log"),
        DiagTarget::Memory(_)
    ));
}

/// A receiver whose run writes its working directory, whose host gives a folder for the files
/// that the run makes outside it (makeblastdb of LOSAT Web, `FileReceiver::outside_directory`),
/// taking its other files as they are written (`Streamed`).
struct PlacesOutside(Streamed, std::path::PathBuf);

impl FileReceiver for PlacesOutside {
    fn receive(&self, name: &[u8], _bytes: &[u8]) -> io::Result<()> {
        panic!("a whole file: {}", String::from_utf8_lossy(name));
    }

    fn open(&self, name: &[u8]) -> Option<Box<dyn FileStream>> {
        self.0.open(name)
    }

    fn writes_working_directory(&self) -> bool {
        true
    }

    fn outside_directory(&self) -> Option<std::path::PathBuf> {
        Some(self.1.clone())
    }
}

/// The folder of a database (makeblastdb of LOSAT Web; line NP WW4i, 2026-11-28): a receiver's
/// run answers the questions of `CreateDirectories` on the model of the other tools' outputs
/// (`OpenFiles::folder_is_directory`, `make_folder`, `folder_writable`; the tests' working
/// directory is LOSAT's, with `src` and `Cargo.toml`): the folders of the run's directory are the
/// file system's, those outside it the model's, `..` read as the system reads it (`nodir/..` and
/// `Cargo.toml/..` are no directory). Where the host gives a folder for the files outside the
/// directory, the files of a name whose folder is outside it where the user makes entries are
/// made in a folder of the host's, one per path (`../` and `src/../../` are one), and handed back
/// under their names, the host's folder left empty; a name of the directory, of root's folders,
/// of a folder that the system does not reach, and every name without the host's folder or on
/// `Disk` are not placed.
#[test]
fn a_run_answers_for_its_database_s_folder_and_places_its_files_outside() {
    let directory = OpenFiles::new(&RunFiles::Receiver(Arc::new(WritesDirectory::default())));
    for (d, is_directory, made, writable) in [
        ("..", true, MadeFolder::Outside(true), true),
        ("../nosub", true, MadeFolder::Outside(true), true),
        (
            "src",
            true,
            MadeFolder::InRunDirectory(b"src".to_vec()),
            true,
        ),
        ("src/..", true, MadeFolder::Outside(true), true),
        (
            "nosub",
            false,
            MadeFolder::InRunDirectory(b"nosub".to_vec()),
            false,
        ),
        ("nosub/..", false, MadeFolder::Outside(true), false),
        ("Cargo.toml/..", false, MadeFolder::Outside(true), false),
        ("/tmp/ww4i-x", true, MadeFolder::Outside(true), true),
        ("/proc", true, MadeFolder::Outside(true), false),
        ("/etc/ww4i-x", false, MadeFolder::Outside(false), false),
        ("/work", false, MadeFolder::Outside(false), false),
    ] {
        let d = d.as_bytes();
        assert_eq!(
            directory.folder_is_directory(d),
            Some(is_directory),
            "{d:?}"
        );
        assert_eq!(directory.make_folder(d), Some(made), "{d:?}");
        assert_eq!(directory.folder_writable(d), Some(writable), "{d:?}");
        assert_eq!(OpenFiles::disk().folder_is_directory(d), None);
        assert_eq!(OpenFiles::disk().make_folder(d), None);
        assert_eq!(OpenFiles::disk().folder_writable(d), None);
    }
    assert_eq!(directory.place(b"../db").unwrap(), None, "no host's folder");
    assert_eq!(OpenFiles::disk().place(b"../db").unwrap(), None);

    let tmp = crate::utils::test_dir::TempDir::unique("run-files-outside");
    let host = tmp.to_path_buf();
    let receiver = PlacesOutside(Streamed::default(), host.clone());
    let events = receiver.0 .0.clone();
    let files = OpenFiles::new(&RunFiles::Receiver(Arc::new(receiver)));
    let placed = |name: &str| {
        files
            .place(name.as_bytes())
            .unwrap()
            .map(|path| String::from_utf8(path).unwrap())
    };
    let at = |index: usize, name: &str| Some(format!("{}/{index}/{name}", host.display()));
    for name in [
        "db",
        "src/db",
        "/proc/db",
        "nodir/../db",
        "/etc/ww4i-x/db",
        "../",
    ] {
        assert_eq!(placed(name), None, "{name}");
    }
    assert_eq!(placed("../db"), at(0, "db"));
    assert_eq!(placed("src/../../m"), at(0, "m"));
    assert_eq!(placed("/tmp/ww4i-x/db"), at(1, "db"));
    assert_eq!(placed("../nosub/db"), at(2, "db"));
    std::fs::write(host.join("0/db.pin"), b"p").unwrap();
    std::fs::write(host.join("0/m.gnd"), b"g").unwrap();
    std::fs::write(host.join("1/db.nal"), b"a").unwrap();
    std::fs::create_dir_all(host.join("1/v")).unwrap();
    std::fs::write(host.join("1/v/x"), b"x").unwrap();
    files.hand_back().unwrap();
    assert_eq!(
        std::fs::read_dir(&host).unwrap().count(),
        0,
        "the host's folder is empty"
    );
    drop(files);
    let written: Vec<String> = events
        .lock()
        .unwrap()
        .iter()
        .filter(|event| event.contains(" @"))
        .cloned()
        .collect();
    assert_eq!(
        written,
        [
            "../db.pin @0 p",
            "../m.gnd @0 g",
            "/tmp/ww4i-x/db.nal @0 a",
            "/tmp/ww4i-x/v/x @0 x"
        ]
    );
}

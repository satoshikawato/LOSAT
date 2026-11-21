use super::*;
use std::sync::Mutex;

#[derive(Default)]
struct Kept(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);

impl FileReceiver for Kept {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
        Ok(())
    }
}

/// The names of `-logfile` of a receiver: NCBI's special names, the names of directories and
/// of the standard error, the null device and the file systems without free space, and a file
/// in memory, handed over under its whole name, also with a directory (line WQA's sweep, D11,
/// 2026-11-20: refused before) and as NCBI reads it (its round 2: `logs/` and ` run.log` are
/// files, `/dev/fd/2` and `/proc/self/fd/2` give no message, `/dev/fd/3` and
/// `/proc/self/status` do not open, `/proc/self/cwd/run.log` is `run.log`; its round 3:
/// `//run.log` has no directory, a user's command line opens no new name in a directory of
/// root's (`/etc/run.log`, `/dev/run.log`) and none of its devices but those whose writes are
/// lost (`/dev/ptmx`), and makes files in `/tmp`, `/var/tmp`, `/dev/shm`, `/run/lock`, its
/// home and its own directories; its round 4: the directories of root's are all but the
/// user's (`/data`), the links of `/var/lock`, `/var/run` and merged `/usr` are followed, a
/// directory of root's that is not there is none (`/bin/../tmp`), and the system refuses a
/// component of 256 bytes and a path of 4096; NCBI BLAST+ 2.17.0 and LOSAT's command line on
/// Ubuntu 24.04 as a user, `fix-logfile-dir/r4/cli1.jsonl` and review 2's `cli-names.jsonl`).
#[test]
fn a_receiver_reads_a_log_name_as_the_command_line_on_linux() {
    let kept = Arc::new(Kept::default());
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
    assert!(matches!(target(None), DiagTarget::Stderr));
    assert!(matches!(target(Some("-")), DiagTarget::Stderr));
    assert!(matches!(target(Some("")), DiagTarget::Discard));
    assert!(matches!(target(Some("/dev/null")), DiagTarget::Discard));
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
        assert!(matches!(target(Some(name)), DiagTarget::Stderr), "{name}");
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
            matches!(target(Some(name)), DiagTarget::StderrFile),
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
        assert!(matches!(target(Some(name)), DiagTarget::Discard), "{name}");
    }
    let names = [
        "run.log",
        "logs/run.log",
        "./run.log",
        "a/../run.log",
        "/tmp/run.log",
        "/dev/stdout",
        "dev/stderr",
        "../run.log",
        "../proc/self/fd/2",
        "logs/",
        "logs//",
        "run.log/",
        "a/b/",
        "/tmp/abs/run.log/",
        " run.log",
        "run.log ",
        " logs/run.log",
        "/proc/self/cwd/run.log",
        "/proc/self/root/tmp/run.log",
        "/home/me/run.log",
        "/home/me/logs/run.log",
        "/var/tmp/run.log",
        "/dev/shm/run.log",
        "/run/lock/run.log",
        // Round 4.
        "/var/lock/run.log",
        "/var/run/lock/run.log",
        "/var/crash/run.log",
        "/var/run/../tmp/run.log",
        "/usr/../tmp/run.log",
        "/tmp/a/../run.log",
        "/dev/shm/a/run.log",
        component_255.as_str(),
        directory_4094.as_str(),
        directory_4095.as_str(),
        // Round 5: a name climbs out of directories of root's that the model does not list,
        // which a user's command line passes (NCBI BLAST+ 2.17.0 writes the file; review 3).
        "/etc/default/../../tmp/run.log",
        "/usr/include/../../tmp/run.log",
        "/usr/share/doc/../../../tmp/run.log",
        "/usr/lib/x86_64-linux-gnu/../../../tmp/run.log",
        "/var/mail/../tmp/run.log",
        "/var/backups/../tmp/run.log",
        "/dev/disk/../../tmp/run.log",
        "/var/log/../../tmp/run.log",
        "/var/lock/../../tmp/run.log",
        // Round 6: the directories around them, and a process's own descriptors, open to the
        // user.
        "/var/cache/../../tmp/run.log",
        "/etc/ssl/../../tmp/run.log",
        "/proc/1/../../tmp/run.log",
        "/proc/self/fd/../../../tmp/run.log",
        // Round 7: the model's host and its first user (review r6): the thread's `..` climbs
        // out of procfs, 40 links, `fdinfo` of Linux 5.15, a directory of the user's group
        // `adm`, a directory at the host's root, and places where the user writes.
        "/proc/thread-self/../../../../tmp/run.log",
        links_40.as_str(),
        root_links_40.as_str(),
        "/proc/1/fdinfo/../../../tmp/run.log",
        "/var/log/unattended-upgrades/../../../tmp/run.log",
        "/snap/../tmp/run.log",
        "/run/user/1000/run.log",
        "/var/run/user/1000/run.log",
        "/run/shm/run.log",
        "/mnt/wsl/run.log",
    ];
    for name in names {
        let log = target(Some(name));
        assert!(matches!(log, DiagTarget::Memory(_)), "{name}");
        log.write_text(name.as_bytes(), &mut io::sink());
        files.deliver_log(Some(OsStr::new(name)), log).unwrap();
    }
    // A target without a file hands nothing over.
    files
        .deliver_log(Some(OsStr::new("-")), DiagTarget::Stderr)
        .unwrap();
    files
        .deliver_log(Some(OsStr::new("/dev/stderr")), DiagTarget::StderrFile)
        .unwrap();
    files
        .deliver_log(Some(OsStr::new("/dev/fd/2")), DiagTarget::Discard)
        .unwrap();
    assert_eq!(
        *kept.0.lock().unwrap(),
        names
            .iter()
            .map(|name| (name.as_bytes().to_vec(), name.as_bytes().to_vec()))
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
    for (name, read) in [
        (format!("{}run.log", "a/".repeat(2040)), LogName::Unopened),
        (format!("{}run.log", "../".repeat(1360)), LogName::File),
    ] {
        assert!(name.len() < 4096);
        WALK_STEPS.with(|steps| steps.set(0));
        let started = std::time::Instant::now();
        assert_eq!(log_name(name.as_bytes()), read);
        let steps = WALK_STEPS.with(std::cell::Cell::get);
        assert!(steps < 4 * name.len(), "{steps} steps");
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
    for (name, read) in [
        ("/etc/passwd/../../tmp/x", LogName::File),
        ("/etc/zzz/../../tmp/x", LogName::File),
        ("/usr/bin/ls/../../../tmp/x", LogName::File),
        ("/proc/self/status/../../../tmp/x", LogName::File),
        ("/dev/null/../../tmp/x", LogName::Unopened),
        ("/proc/self/comm/../../../tmp/x", LogName::Unopened),
    ] {
        assert_eq!(log_name(name.as_bytes()), read, "{name}");
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

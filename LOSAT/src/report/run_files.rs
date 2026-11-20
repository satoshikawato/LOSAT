//! The files that a run writes besides its reports: the file of `-logfile` (line RG item E1w-3a).
//! The command line and the library write them to the file system (`RunFiles::Disk`); LOSAT Web,
//! which has none, hands each to its host (`RunFiles::Receiver`), as it hands over the side files
//! of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`). No NCBI counterpart: NCBI writes
//! files.

use std::ffi::OsStr;
use std::io;
use std::sync::Arc;

use crate::blastinput::ncbi_environment::{diag_target, special_log_target, DiagTarget};

/// The receiver of the files of a run that a caller without files takes in place of the files
/// (LOSAT Web's host): the side files of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`) and
/// the files of `RunFiles::Receiver`.
pub trait FileReceiver: Send + Sync {
    /// One file, whole: `name` is the path that the command line writes it to, `bytes` its
    /// contents.
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()>;
}

/// Where the files of a run go.
#[derive(Clone)]
pub enum RunFiles {
    /// The file system, at the names that the arguments give (the command line, the library).
    Disk,
    /// The caller's receiver, which gets each file when the run ends (LOSAT Web).
    Receiver(Arc<dyn FileReceiver>),
}

impl RunFiles {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608 (quoted at
    // `ncbi_environment::diag_target`): no file for `""` and `/dev/null`, `cerr` for `-`, else the
    // file, or `cerr` when its directory does not exist or it cannot be opened (a directory).
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), and any other name as the command line on Linux reads it in a
    /// working directory of its own where the name's directories exist (`log_name`), as LOSAT Web
    /// reads a `-conffile` name as a path of the command line's working directory (line WQA's
    /// sweep, D11, 2026-11-20; it refused a name with a directory before): `cerr` for a directory,
    /// which NCBI cannot open as a file, and for a name that a user's command line cannot open
    /// (`//run.log`, `/etc/run.log`, `/dev/kmsg`); the file that is the standard error
    /// (`DiagTarget::StderrFile`) for `/dev/stderr`; no output for another spelling of `/dev/null`,
    /// for a device whose writes are lost (`/dev/ptmx`) and for a name in a file system that
    /// reports no free space (`/dev/fd/2`); and for any other
    /// name a file in memory under the whole name, which `deliver_log` hands over (`logs/run.log`,
    /// `logs/`, `/tmp/run.log`, and `/dev/stdout`, whose bytes the command line writes into its
    /// standard output among the report's, where their place depends on the report's buffering).
    pub fn log_target(&self, name: Option<&OsStr>) -> DiagTarget {
        match self {
            RunFiles::Disk => diag_target(name),
            RunFiles::Receiver(_) => match (special_log_target(name), name) {
                (Some(target), _) => target,
                (None, None) => DiagTarget::Stderr,
                (None, Some(name)) => match log_name(name.as_encoded_bytes()) {
                    LogName::Unopened => DiagTarget::Stderr,
                    LogName::StandardError => DiagTarget::StderrFile,
                    LogName::Nowhere => DiagTarget::Discard,
                    LogName::File => DiagTarget::Memory(Default::default()),
                },
            },
        }
    }

    /// Hands the file of the `-logfile` value `name`, which `target` holds (`log_target`), to the
    /// receiver. Nothing on `Disk` (the file is written) and for a target without a file.
    pub fn deliver_log(&self, name: Option<&OsStr>, target: DiagTarget) -> io::Result<()> {
        match (self, name, target) {
            (RunFiles::Receiver(receiver), Some(name), DiagTarget::Memory(bytes)) => receiver
                .receive(
                    name.as_encoded_bytes(),
                    &bytes
                        .into_inner()
                        .unwrap_or_else(std::sync::PoisonError::into_inner),
                ),
            _ => Ok(()),
        }
    }
}

/// What a `-logfile` name that is not one of NCBI's special names gives, for a `Receiver`
/// (`RunFiles::log_target`).
#[derive(Debug, PartialEq, Eq)]
enum LogName {
    /// A name that NCBI cannot open, a directory: the messages stay on `cerr`.
    Unopened,
    /// The process's standard error, opened as a file.
    StandardError,
    /// No message: the null device written otherwise than NCBI's special `/dev/null`, or a file
    /// system that reports no free space.
    Nowhere,
    /// A file of that name.
    File,
}

/// The directories of LOSAT Web's Linux (`log_name`) besides the working directory, its
/// ancestors and the name's own: those of the Filesystem Hierarchy Standard at the root and the
/// directories of the file systems that report no free space (`NO_FREE_SPACE`), each with its
/// ancestors.
const DIRECTORIES: [&[u8]; 23] = [
    b"/bin",
    b"/boot",
    b"/dev/fd",
    b"/dev/mqueue",
    b"/dev/pts",
    b"/dev/shm",
    b"/etc",
    b"/home",
    b"/lib",
    b"/lib64",
    b"/media",
    b"/mnt",
    b"/opt",
    b"/proc/self/fd",
    b"/proc/thread-self/fd",
    b"/root",
    b"/run",
    b"/sbin",
    b"/srv",
    b"/sys",
    b"/tmp",
    b"/usr",
    b"/var",
];

/// The directories of LOSAT Web's Linux that are root's (mode 0755 and under it), where a user
/// makes no file and writes none of the system's: `/` itself, `/dev` and those of `DIRECTORIES`
/// at the root that are not file systems of their own, with every directory in them, except
/// `SHARED_DIRECTORIES` and a user's home (`/home/<user>`).
const SYSTEM_DIRECTORIES: [&[u8]; 16] = [
    b"/bin", b"/boot", b"/dev", b"/etc", b"/home", b"/lib", b"/lib64", b"/media", b"/mnt", b"/opt",
    b"/root", b"/run", b"/sbin", b"/srv", b"/usr", b"/var",
];

/// The directories of LOSAT Web's Linux where every user makes files (mode 1777): `/tmp`,
/// `/var/tmp`, `/dev/shm` and `/run/lock` (line WQA's D11 round 3, 2026-11-20: NCBI BLAST+ 2.17.0
/// and LOSAT's command line write their log there as a user); `/dev/mqueue` too, whose names are
/// message queues (`opens`).
const SHARED_DIRECTORIES: [&[u8]; 5] = [
    b"/tmp",
    b"/var/tmp",
    b"/dev/shm",
    b"/dev/mqueue",
    b"/run/lock",
];

/// The file systems of Linux that report no free space (`statfs`'s `f_bavail` is 0): procfs,
/// sysfs, devpts and mqueue, where they are mounted.
const NO_FREE_SPACE: [&[u8]; 4] = [b"/proc", b"/sys", b"/dev/pts", b"/dev/mqueue"];

/// The devices of LOSAT Web's Linux that a user opens for writing (mode 0666), whose writes are
/// lost: the null device written otherwise than NCBI's special `/dev/null` (`//dev/null`),
/// `/dev/full` (whose writes fail, which NCBI does not check), `/dev/zero`, `/dev/random`,
/// `/dev/urandom`, and a new terminal's master (`/dev/ptmx`), FUSE's device (`/dev/fuse`) and the
/// TUN device (`/dev/net/tun`), whose writes fail or no one reads (line WQA's D11 round 3,
/// 2026-11-20: NCBI BLAST+ 2.17.0 and LOSAT's command line write no message there). `/dev/tty`
/// is a device too, which a run without a terminal does not open.
const LOST_DEVICES: [&[u8]; 8] = [
    b"/dev/null",
    b"/dev/full",
    b"/dev/zero",
    b"/dev/random",
    b"/dev/urandom",
    b"/dev/ptmx",
    b"/dev/fuse",
    b"/dev/net/tun",
];

/// The files of procfs that a process may open for writing as a user, of the process
/// (`/proc/self`) and of its thread (`/proc/thread-self`): those of Linux 5.15, which
/// `open(O_WRONLY | O_APPEND | O_CREAT)` opens (line WQA's sweep, D11 round 2, 2026-11-20). The
/// file descriptors 0 to 2 (`fd/0`...) are the standard streams.
const PROCESS_FILES: [&[u8]; 18] = [
    b"attr/current",
    b"attr/exec",
    b"attr/fscreate",
    b"attr/keycreate",
    b"attr/sockcreate",
    b"clear_refs",
    b"comm",
    b"coredump_filter",
    b"gid_map",
    b"loginuid",
    b"mem",
    b"oom_adj",
    b"oom_score_adj",
    b"projid_map",
    b"sched",
    b"timens_offsets",
    b"timerslack_ns",
    b"uid_map",
];
const THREAD_FILES: [&[u8]; 15] = [
    b"attr/current",
    b"attr/exec",
    b"attr/fscreate",
    b"attr/keycreate",
    b"attr/sockcreate",
    b"clear_refs",
    b"comm",
    b"gid_map",
    b"loginuid",
    b"mem",
    b"oom_adj",
    b"oom_score_adj",
    b"projid_map",
    b"sched",
    b"uid_map",
];

/// Whether `path` is `directory` or a directory in it, component by component.
fn within(path: &[u8], directory: &[u8]) -> bool {
    path.strip_prefix(directory)
        .is_some_and(|rest| rest.is_empty() || rest[0] == b'/' || directory.ends_with(b"/"))
}

/// Whether a user's system opens `path` for writing, a path without links whose directories exist:
/// a file that a name may make or write, also a device (No NCBI counterpart). In procfs, only the
/// standard streams and the files of `PROCESS_FILES` and `THREAD_FILES` (a name that procfs does
/// not have is not made there); sysfs and devpts open none of a user's names; an mqueue file is a
/// name in `/dev/mqueue`. In a directory of root's (`SYSTEM_DIRECTORIES`, the root), only the
/// devices of `LOST_DEVICES`: a new name there is not made, and another device or a file of the
/// system's (`/dev/kmsg`, `/dev/console`, `/dev/tty` in a run without a terminal, `/etc/passwd`)
/// is not open to a user; a name in a user's directory is (`/tmp/run.log`, `/home/me/run.log`,
/// `/data/run.log`, whose directory the name gives).
fn opens(path: &[u8]) -> bool {
    if let Some(file) = path.strip_prefix(b"/proc/self/") {
        return matches!(file, b"fd/0" | b"fd/1" | b"fd/2") || PROCESS_FILES.contains(&file);
    }
    if let Some(file) = path.strip_prefix(b"/proc/thread-self/") {
        return THREAD_FILES.contains(&file);
    }
    if let Some(queue) = path.strip_prefix(b"/dev/mqueue/") {
        return !queue.contains(&b'/');
    }
    if [b"/proc".as_slice(), b"/sys", b"/dev/pts"]
        .iter()
        .any(|system| within(path, system))
    {
        return false;
    }
    if LOST_DEVICES.contains(&path) {
        return true;
    }
    let directory = &path[..path.iter().rposition(|byte| *byte == b'/').unwrap_or(0)];
    let home = directory
        .strip_prefix(b"/home/")
        .is_some_and(|user| !user.is_empty());
    let shared = SHARED_DIRECTORIES
        .iter()
        .any(|shared| within(directory, shared));
    home || shared
        || !(directory.is_empty()
            || SYSTEM_DIRECTORIES
                .iter()
                .any(|system| within(directory, system)))
}

/// Whether `path`, a path without links, is a file of LOSAT Web's Linux, which is no directory
/// (`/dev/null/x` and `/dev/stdout/../stderr` name no file): an entry of `/proc/self/fd`, a file
/// of `PROCESS_FILES` or `THREAD_FILES`, a device of `LOST_DEVICES` or `/dev/tty`.
fn known_file(path: &[u8]) -> bool {
    path.strip_prefix(b"/proc/self/fd/")
        .is_some_and(|fd| !fd.is_empty())
        || path
            .strip_prefix(b"/proc/self/")
            .is_some_and(|file| PROCESS_FILES.contains(&file))
        || path
            .strip_prefix(b"/proc/thread-self/")
            .is_some_and(|file| THREAD_FILES.contains(&file))
        || LOST_DEVICES.contains(&path)
        || path == b"/dev/tty"
}

/// What the `-logfile` name `name` gives on the command line on Linux (No NCBI counterpart: LOSAT
/// Web has no file system; line WQA's sweep, D11, 2026-11-20 and its rounds 2 and 3). The name is
/// read as NCBI reads it on Unix (`ncbi_environment::log_file_name`, the command line's reading
/// too) in a working directory of its own, where the directories that the name gives exist (those
/// that NCBI checks and those of the file it opens) on a file system with free space: `logs/run.log`,
/// ` run.log` (`run.log`), `logs/` (the file `logs`) and `/tmp/run.log` are files. A directory of
/// '/'s only (`//run.log`) is none, as `CDirEntry::Reset` leaves it empty. The system then follows
/// the links of Linux in the path that NCBI opens and in the directory whose free space it reads:
/// `/dev/fd` (`/proc/self/fd`), `/dev/stdin`, `/dev/stdout` and `/dev/stderr` (`/proc/self/fd/0`
/// to `2`), `/proc/self/cwd` (the working directory) and `/proc/self/root` (`/`), also of
/// `/proc/thread-self`. The opened path is a directory, which NCBI cannot open (`cerr`), where it
/// is the working directory or one of its ancestors, a directory that the name gives, one of
/// `DIRECTORIES` or an ancestor of any of them (`.`, `logs/.`, `a/b/..`, `/tmp/`,
/// `/proc/self/cwd`); NCBI cannot open a path that the system does not open to a user either
/// (`opens`: `/dev/fd/3`, `/proc/self/status`, `/sys/...`, `/dev/pts/0`, `/dev/tty`,
/// `/run.log`, `/etc/run.log`, `/var/log/run.log`, `/dev/run.log`, `/dev/kmsg`), nor one whose
/// directory is a file (`known_file`: `/dev/null/x`). A file opened in a file system that reports
/// no free space (`NO_FREE_SPACE`) gets no message, as NCBI checks the space before it writes
/// (`/dev/fd/2`, `/proc/self/fd/1`, `/proc/self/comm`, `/dev/mqueue/q`, `/dev/fd/../stderr`,
/// whose directory is `/proc/self`). Otherwise `/dev/stderr` is the standard error, also written
/// otherwise (`//dev/stderr`); the standard input is the null device (`/dev/stdin`), as for a
/// command line run without one; the devices whose writes are lost (`LOST_DEVICES`: `/dev/ptmx`)
/// give no message; `/dev/stdout` is a file of that name, whose bytes the command line writes into
/// its standard output among the report's, and `/proc/self/cwd/run.log` the file `run.log` of the
/// working directory, under the name as typed. The model's command line is a user's (not root),
/// with no terminal, its standard output and error being pipes, as a shell's `|` gives.
fn log_name(name: &[u8]) -> LogName {
    // A working directory deeper than the name's ".." reach, of components no name has (a
    // command line's words have no NUL byte): its ancestors are never the root's directories.
    let climbs = crate::ncbi_str::truncate_spaces(name)
        .split(|byte| *byte == b'/')
        .filter(|part| *part == b"..")
        .count();
    let cwd = b"/\0".repeat(climbs + 1);
    let file = crate::blastinput::ncbi_environment::log_file_name(
        crate::blastinput::ncbi_environment::PathRules::Unix,
        name,
        Some(&cwd),
    )
    .expect("NCBI's Unix reading reads every name");
    let link = |path: &[u8]| -> Option<Vec<u8>> {
        Some(match path {
            b"/dev/fd" | b"/proc/thread-self/fd" => b"/proc/self/fd".to_vec(),
            b"/dev/stdin" => b"/proc/self/fd/0".to_vec(),
            b"/dev/stdout" => b"/proc/self/fd/1".to_vec(),
            b"/dev/stderr" => b"/proc/self/fd/2".to_vec(),
            b"/proc/self/cwd" | b"/proc/thread-self/cwd" => cwd.clone(),
            b"/proc/self/root" | b"/proc/thread-self/root" => b"/".to_vec(),
            _ => return None,
        })
    };
    let resolve = |path: &[u8]| {
        crate::blastinput::ncbi_environment::normalize_path(path, Some(&link))
            .expect("the links of LOSAT Web's Linux end")
    };
    let mut directories: Vec<Vec<u8>> = Vec::new();
    if let Some(directory) = &file.directory {
        if directory.is_empty() {
            return LogName::Unopened;
        }
        let given = if directory.first() == Some(&b'/') {
            directory.clone()
        } else {
            [cwd.as_slice(), b"/", directory].concat()
        };
        directories.extend(
            given
                .iter()
                .enumerate()
                .filter(|(_, byte)| **byte == b'/')
                .map(|(at, _)| resolve(&given[..=at])),
        );
        directories.push(resolve(&given));
    }
    directories.push(cwd.clone());
    directories.extend(DIRECTORIES.iter().map(|directory| directory.to_vec()));
    let path = resolve(&file.path);
    if file.path.is_empty()
        || directories
            .iter()
            .any(|directory| within(directory, &path) || known_file(directory))
        || !opens(&path)
    {
        return LogName::Unopened;
    }
    let space = resolve(&file.space);
    if NO_FREE_SPACE.iter().any(|system| within(&space, system)) {
        return LogName::Nowhere;
    }
    match path.as_slice() {
        b"/proc/self/fd/2" => LogName::StandardError,
        b"/proc/self/fd/0" => LogName::Nowhere,
        path if LOST_DEVICES.contains(&path) => LogName::Nowhere,
        _ => LogName::File,
    }
}

#[cfg(test)]
mod tests {
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
    /// home and its own directories).
    #[test]
    fn a_receiver_reads_a_log_name_as_the_command_line_on_linux() {
        let kept = Arc::new(Kept::default());
        let files = RunFiles::Receiver(kept.clone());
        let target = |name: Option<&str>| files.log_target(name.map(OsStr::new));
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
            "/data/run.log",
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
}

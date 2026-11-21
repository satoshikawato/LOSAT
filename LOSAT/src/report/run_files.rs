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
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608 (cited at
    // `ncbi_environment::diag_target`): no file for `""` and `/dev/null`, `cerr` for `-`, else the
    // file, or `cerr` when its directory does not exist or it cannot be opened (a directory).
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), and any other name as the command line of a user on Linux reads it
    /// in a working directory of the user's own (`log_name`), as LOSAT Web reads a `-conffile`
    /// name as a path of the command line's working directory (line WQA's sweep, D11, 2026-11-20;
    /// it refused a name with a directory before): `cerr` for a directory, which NCBI cannot open
    /// as a file, and for a name that a user's command line cannot open (`//run.log`,
    /// `/etc/run.log`, `/data/run.log`, `/dev/kmsg`, a component of 256 bytes); the file that is
    /// the standard error (`DiagTarget::StderrFile`) for `/dev/stderr`; no output for another
    /// spelling of `/dev/null`, for a device whose writes are lost (`/dev/ptmx`) and for a name in
    /// a file system that reports no free space (`/dev/fd/2`); and for any other name, one in the
    /// working directory, in the user's home or in a directory where every user makes files, a
    /// file in memory under the whole name, which `deliver_log` hands over (`logs/run.log`,
    /// `logs/`, `/tmp/run.log`, `/var/lock/run.log`, and `/dev/stdout`, whose bytes the command
    /// line writes into its standard output among the report's, where their place depends on the
    /// report's buffering).
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

/// The working directory's tree on LOSAT Web's Linux (`log_name`): the working directory is a
/// chain of `"\0"` components in it (`/\0/\0/...`), deeper than a name's `..` reach, so that no
/// name climbs out of it, as a user's command line runs in a directory of the user's own
/// (`~/work`). A command line's words have no NUL byte, so no name gives another path in it.
const WORKING_TREE: &[u8] = b"/\0";

/// The model's host, whose system LOSAT Web's `-logfile` names are read on (`log_name`), and the
/// reference of every measured list here: Ubuntu 24.04.3, kernel 5.15 (WSL2), and its first user
/// (uid 1000, in that user's groups), not root (line WQA's D11 round 7, 2026-11-21: one host and
/// one user, Owner-delegated).
///
/// The directories where the model's user makes files, besides the working directory's tree and a
/// home (`user_makes`): those that the user writes on the model's host, measured as the user
/// (`find -type d -writable`, line WQA's D11 round 7): where every user makes files (mode 1777,
/// `/var/crash` 3777: `/tmp`, `/var/tmp`, `/dev/shm` and `/run/shm`, `/run/lock`, also as
/// `/var/lock`, `/var/crash`, WSL's `/mnt/wsl` and `/mnt/wslg`), the user's own `/run/user/1000`,
/// and WSL's drives of Windows (`/mnt/c`, `/mnt/d`, mode 0777). `/dev/mqueue` (1777 too) holds
/// message queues only (`opens`).
const USER_DIRECTORIES: [&[u8]; 11] = [
    b"/tmp",
    b"/var/tmp",
    b"/dev/shm",
    b"/run/shm",
    b"/run/lock",
    b"/var/crash",
    b"/run/user/1000",
    b"/mnt/wsl",
    b"/mnt/wslg",
    b"/mnt/c",
    b"/mnt/d",
];

/// Whether the user of LOSAT Web's Linux makes entries in `directory`, a path without links: in
/// the working directory's tree (`WORKING_TREE`), in a user's home (`/home/<user>`, which LOSAT
/// Web takes as the user's own, not knowing the user's name) and in `USER_DIRECTORIES`, each with
/// the directories that the user makes in it. Every other directory is root's, as only root makes
/// entries in `/` (line WQA's D11 round 4, 2026-11-20, Owner-delegated: one rule of owners in
/// place of lists of root's directories): `/`, `/data`, `/home`, `/etc`, `/var/log`.
fn user_makes(directory: &[u8]) -> bool {
    within(directory, WORKING_TREE)
        || directory
            .strip_prefix(b"/home/")
            .is_some_and(|user| !user.is_empty())
        || USER_DIRECTORIES
            .iter()
            .any(|place| within(directory, place))
}

/// The directories at the root of LOSAT Web's Linux (`directory_at`): the model's host's (`ls /`,
/// line WQA's D11 round 7), an Ubuntu of merged `/usr` (`/bin`, `/lib`, `/lib32`, `/lib64` and
/// `/sbin` are links, `LINKS`) with its marks `*.usr-is-merged`, Snap's `/snap` and the file
/// system's `/lost+found`. WSL's `/init` is a file; another name at the root is not there
/// (`/data`).
const ROOT_DIRECTORIES: [&[u8]; 20] = [
    b"/bin.usr-is-merged",
    b"/boot",
    b"/dev",
    b"/etc",
    b"/home",
    b"/lib.usr-is-merged",
    b"/lost+found",
    b"/media",
    b"/mnt",
    b"/opt",
    b"/proc",
    b"/root",
    b"/run",
    b"/sbin.usr-is-merged",
    b"/snap",
    b"/srv",
    b"/sys",
    b"/tmp",
    b"/usr",
    b"/var",
];

/// LOSAT Web's process in procfs: `/proc/self` is a link to it (`LINKS`), as to a process id, here
/// a NUL byte, which no name gives (`WORKING_TREE`).
const SELF: &[u8] = b"/proc/\0";

/// LOSAT Web's thread in procfs: `/proc/thread-self` is a link to it (`LINKS`), a task of `SELF`
/// whose id is a NUL byte.
const THREAD: &[u8] = b"/proc/\0/task/\0";

/// The links of LOSAT Web's Linux and their targets (`Walk`), besides `cwd` of `SELF` and of
/// `THREAD`, the working directory: merged `/usr`'s, the Filesystem Hierarchy Standard's
/// `/var/lock` and `/var/run`, the standard streams' (`/dev/fd`, `/dev/stdin`, `/dev/stdout`,
/// `/dev/stderr`) and procfs's (`/proc/self`, `/proc/thread-self` and their `root`).
const LINKS: [(&[u8], &[u8]); 15] = [
    (b"/bin", b"/usr/bin"),
    (b"/lib", b"/usr/lib"),
    (b"/lib32", b"/usr/lib32"),
    (b"/lib64", b"/usr/lib64"),
    (b"/sbin", b"/usr/sbin"),
    (b"/var/lock", b"/run/lock"),
    (b"/var/run", b"/run"),
    (b"/dev/fd", b"/proc/self/fd"),
    (b"/dev/stdin", b"/proc/self/fd/0"),
    (b"/dev/stdout", b"/proc/self/fd/1"),
    (b"/dev/stderr", b"/proc/self/fd/2"),
    (b"/proc/self", SELF),
    (b"/proc/thread-self", THREAD),
    (b"/proc/\0/root", b"/"),
    (b"/proc/\0/task/\0/root", b"/"),
];

/// The links that the system follows in one lookup (`Walk`), each link of the path and of the
/// links' targets: a name through 40 writes its file on the model's host's command line, and one
/// through 41 does not (`ELOOP`; observed in line WQA's D11 review r6, 2026-11-21:
/// `/var/lock/../..` 40 and 41 times, `/proc/self/root`, two links, 20 and 21 times).
const MAX_LINKS: usize = 40;

/// The directories of LOSAT Web's Linux in which its user looks up no name (`searchable`): those
/// that the model's user cannot search on the model's host, measured as that user at every depth
/// but in the user's places and other processes (`find / -type d ! -executable`, line WQA's D11
/// round 7, 2026-11-21). Some are the host's state then: Snap's FUSE connections and the
/// directories of systemd's units in `/run/systemd/propagate`.
const CLOSED_DIRECTORIES: [&[u8]; 39] = [
    b"/etc/credstore",
    b"/etc/credstore.encrypted",
    b"/etc/polkit-1/rules.d",
    b"/etc/ssl/private",
    b"/lost+found",
    b"/proc/tty/driver",
    b"/root",
    b"/run/docker",
    b"/run/systemd/inaccessible/dir",
    b"/run/systemd/propagate/systemd-logind.service",
    b"/run/systemd/propagate/systemd-resolved.service",
    b"/run/systemd/propagate/systemd-timesyncd.service",
    b"/run/systemd/propagate/systemd-udevd.service",
    b"/sys/fs/fuse/connections/86",
    b"/sys/fs/fuse/connections/87",
    b"/sys/fs/fuse/connections/88",
    b"/sys/fs/fuse/connections/89",
    b"/sys/fs/fuse/connections/90",
    b"/sys/fs/fuse/connections/91",
    b"/sys/fs/fuse/connections/92",
    b"/sys/fs/fuse/connections/93",
    b"/sys/fs/fuse/connections/94",
    b"/sys/fs/fuse/connections/95",
    b"/sys/kernel/debug",
    b"/sys/kernel/tracing",
    b"/var/cache/apt/archives/partial",
    b"/var/cache/ldconfig",
    b"/var/cache/private",
    b"/var/lib/apt/lists/partial",
    b"/var/lib/docker",
    b"/var/lib/polkit-1",
    b"/var/lib/private",
    b"/var/lib/snapd/cache",
    b"/var/lib/snapd/cookie",
    b"/var/lib/ubuntu-advantage/apt-esm/var/lib/apt/lists/partial",
    b"/var/lib/update-notifier/package-data-downloads/partial",
    b"/var/log/private",
    b"/var/spool/cron/crontabs",
    b"/var/spool/rsyslog",
];

/// The entries of another process (`/proc/<pid>` and its threads' `/proc/<pid>/task/<tid>`) in
/// which the model's user looks up no name (`searchable`): its descriptors (`fd`, `map_files`) and
/// its links to its executable and its working and root directories (`exe`, `cwd`, `root`), which
/// only the process's user may follow. LOSAT Web's own process is `SELF`; a process of
/// `/proc/<pid>` is root's, as `/proc/1` (measured as the user on the model's host: `stat` of
/// `<entry>/.` is refused there, and allowed in `fdinfo`, `ns`, `attr` and `net`; line WQA's D11
/// round 7, 2026-11-21).
const CLOSED_PROCESS_ENTRIES: [&[u8]; 5] = [b"cwd", b"exe", b"fd", b"map_files", b"root"];

/// Whether a user looks up names in `directory`, a path without links (`Walk`): not in
/// `CLOSED_DIRECTORIES` nor in another process's `CLOSED_PROCESS_ENTRIES`, where the system
/// refuses (`EACCES`), whatever the name (`.` and `..` too).
fn searchable(directory: &[u8]) -> bool {
    let number = |part: &[u8]| !part.is_empty() && part.iter().all(u8::is_ascii_digit);
    let closed_entry = directory.strip_prefix(b"/proc/").is_some_and(|entry| {
        match entry
            .split(|byte| *byte == b'/')
            .collect::<Vec<_>>()
            .as_slice()
        {
            [pid, closed] => number(pid) && CLOSED_PROCESS_ENTRIES.contains(closed),
            [pid, task, tid, closed] => {
                number(pid)
                    && *task == b"task"
                    && number(tid)
                    && CLOSED_PROCESS_ENTRIES.contains(closed)
            }
            _ => false,
        }
    });
    !closed_entry && !CLOSED_DIRECTORIES.contains(&directory)
}

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

/// The files of procfs that a process may open for writing as a user, of the process (`SELF`)
/// and of its thread (`THREAD`): those of the model's host's Linux 5.15 that
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

/// Linux's limits of the paths given to the system (Linux 6.8: include/uapi/linux/limits.h:12-13,
/// `NAME_MAX` and `PATH_MAX`): a component of at most 255 bytes, a path of fewer than 4096 (with
/// its closing NUL). The system refuses a longer one (`ENAMETOOLONG`), which NCBI then cannot open.
const NAME_MAX: usize = 255;
const PATH_MAX: usize = 4096;

/// Whether the system takes `path`, a path given to it as it is (`NAME_MAX`, `PATH_MAX`).
fn fits(path: &[u8]) -> bool {
    path.len() < PATH_MAX
        && path
            .split(|byte| *byte == b'/')
            .all(|part| part.len() <= NAME_MAX)
}

/// Whether `path` is `directory` or a directory in it, component by component.
fn within(path: &[u8], directory: &[u8]) -> bool {
    path.strip_prefix(directory)
        .is_some_and(|rest| rest.is_empty() || rest[0] == b'/' || directory.ends_with(b"/"))
}

/// Whether `path`, a path without links, is a file of procfs that the model knows in `SELF` or
/// `THREAD`: a standard stream (`fd/0` to `fd/2`) or a file of `PROCESS_FILES` or `THREAD_FILES`.
fn is_procfs_file(path: &[u8]) -> bool {
    let stream = |file: &[u8]| matches!(file, b"fd/0" | b"fd/1" | b"fd/2");
    let file = |directory: &[u8]| path.strip_prefix(directory)?.strip_prefix(b"/");
    match (file(THREAD), file(SELF)) {
        (Some(file), _) => stream(file) || THREAD_FILES.contains(&file),
        (None, Some(file)) => stream(file) || PROCESS_FILES.contains(&file),
        (None, None) => false,
    }
}

/// Whether `path`, a path without links, is one of the files that the model knows (`opens`): a
/// device of `LOST_DEVICES` or a file of procfs (`is_procfs_file`). Such a file is no directory
/// (`/dev/null/x`, `/dev/stdout/../stderr`, `/dev/fd/2/..`).
fn is_known_file(path: &[u8]) -> bool {
    LOST_DEVICES.contains(&path) || is_procfs_file(path)
}

/// Whether `directory` holds a file that the model knows (`is_known_file`): `/dev`, `/dev/net`,
/// `SELF` and its `fd` and `attr`, its `task`, `THREAD` and its `fd` and `attr`.
fn holds_known_file(directory: &[u8]) -> bool {
    let holds = |path: &[u8]| path.len() > directory.len() && within(path, directory);
    let in_procfs = |process: &[u8], files: &[&[u8]]| {
        [b"fd/0".as_slice()]
            .iter()
            .chain(files)
            .any(|file| holds(&[process, b"/".as_slice(), *file].concat()))
    };
    LOST_DEVICES.iter().any(|device| holds(device))
        || in_procfs(SELF, PROCESS_FILES.as_slice())
        || in_procfs(THREAD, THREAD_FILES.as_slice())
}

/// What LOSAT Web's Linux has at `directory`, a path without links (`log_name`).
#[derive(Debug, PartialEq, Eq)]
enum Directory {
    /// A directory that the model knows: `/`, one where the user makes entries (`user_makes`),
    /// which the model takes as made where a name gives it (`logs`, `/tmp/a`; the working
    /// directory and its ancestors), a directory at the root (`ROOT_DIRECTORIES`), a mount of
    /// `NO_FREE_SPACE` (`/dev/pts`) and a directory of a file that the model knows (`/dev/net`,
    /// `SELF`'s `fd`).
    Known,
    /// Any other path in root's directories, which may be a directory there or not
    /// (`/etc/default`, `/usr/tmp`): the user cannot change them, and the model does not list
    /// them.
    Unknown,
    /// No directory: a file that the model knows (`/dev/null`, `SELF`'s `fd/1`) or a name at the
    /// root that is not there (`/data`).
    None,
}

/// What the model knows of `directory`, a path without links (`Directory`).
fn directory_at(directory: &[u8]) -> Directory {
    if directory == b"/"
        || user_makes(directory)
        || ROOT_DIRECTORIES.contains(&directory)
        || NO_FREE_SPACE.contains(&directory)
        || holds_known_file(directory)
    {
        Directory::Known
    } else if is_known_file(directory)
        || !ROOT_DIRECTORIES.iter().any(|root| within(directory, root))
    {
        Directory::None
    } else {
        Directory::Unknown
    }
}

/// Whether a user's system opens `path` for writing, a path without links in a directory that
/// exists: a file that a name may make or write, also a device (No NCBI counterpart). In procfs,
/// only the standard streams and the files of `PROCESS_FILES` and `THREAD_FILES` (a name that
/// procfs does not have is not made there); sysfs and devpts open none of a user's names; an
/// mqueue file is a name in `/dev/mqueue`; the devices of `LOST_DEVICES` are open to every user.
/// Any other name opens where the user makes entries (`user_makes`): a new name in a directory of
/// root's is not made, and another device or a file of the system's (`/dev/kmsg`, `/dev/console`,
/// `/dev/tty` in a run without a terminal, `/etc/passwd`) is not open to a user.
fn opens(path: &[u8]) -> bool {
    if within(path, SELF) {
        return is_procfs_file(path);
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
    LOST_DEVICES.contains(&path)
        || user_makes(&path[..path.iter().rposition(|byte| *byte == b'/').unwrap_or(0)])
}

/// A path of LOSAT Web's Linux walked component by component as the system walks it (No NCBI
/// counterpart: the system's resolution, `log_name`): `.` and empty components stay, `..` goes to
/// the parent of the path walked so far (`/` stays), and a link (`LINKS`; `cwd` of `SELF` and
/// `THREAD`, the working directory `cwd`) is replaced by its target, which is walked from `/`. The
/// system refuses the path where a component is looked up in a directory that the user cannot
/// search (`searchable`: `/root/..`, `/root/.`), and at the link past `MAX_LINKS` (also a
/// standard stream of procfs, `fd/0` to `fd/2`, which is a link there). A name is walked once,
/// each component once and each link's target once more; a step's work is at most the length of
/// the path walked so far, under `PATH_MAX` (line WQA's D11 round 5, 2026-11-21: resolving each
/// prefix again was cubic, a minute for a name of 4 KB in LOSAT Web's serial engine).
struct Walk<'a> {
    /// The path walked so far, without links; empty for `/`.
    resolved: Vec<u8>,
    /// The working directory, the target of `cwd` in `SELF`.
    cwd: &'a [u8],
    /// The links that the walk followed.
    links: usize,
    /// Whether the system refused the path.
    refused: bool,
}

impl<'a> Walk<'a> {
    /// The walk of `path` before its first component: at `/` for an absolute path, else at the
    /// working directory `cwd`.
    fn at(path: &[u8], cwd: &'a [u8]) -> Self {
        let mut walk = Walk {
            resolved: Vec::new(),
            cwd,
            links: 0,
            refused: false,
        };
        if path.first() != Some(&b'/') {
            walk.walk(cwd);
        }
        walk
    }

    /// `path` with its links followed and its `.` and `..` taken; `None` where the system refuses
    /// it (`refused`).
    fn resolve(path: &[u8], cwd: &'a [u8]) -> Option<Vec<u8>> {
        let mut walk = Walk::at(path, cwd);
        walk.walk(path);
        (!walk.refused).then(|| walk.path().to_vec())
    }

    /// Takes every component of `path`.
    fn walk(&mut self, path: &[u8]) {
        for part in path.split(|byte| *byte == b'/') {
            self.step(part);
        }
    }

    /// Takes the component `part`.
    fn step(&mut self, part: &[u8]) {
        #[cfg(test)]
        WALK_STEPS.with(|steps| steps.set(steps.get() + 1));
        if part.is_empty() || self.refused {
            return;
        }
        if !searchable(self.path()) {
            self.refused = true;
            return;
        }
        match part {
            b"." => {}
            b".." => {
                let parent = self.resolved.iter().rposition(|byte| *byte == b'/');
                self.resolved.truncate(parent.unwrap_or(0));
            }
            _ => {
                self.resolved.push(b'/');
                self.resolved.extend_from_slice(part);
                let path = self.resolved.as_slice();
                let link = |of: &[u8], name: &[u8]| path.strip_prefix(of) == Some(name);
                let target = if link(SELF, b"/cwd".as_slice()) || link(THREAD, b"/cwd".as_slice()) {
                    Some(self.cwd)
                } else {
                    LINKS
                        .iter()
                        .find(|(link, _)| *link == path)
                        .map(|(_, target)| *target)
                };
                let stream = [SELF, THREAD].iter().any(|of| {
                    [b"/fd/0".as_slice(), b"/fd/1", b"/fd/2"]
                        .iter()
                        .any(|fd| link(of, fd))
                });
                if target.is_some() || stream {
                    self.links += 1;
                    if self.links > MAX_LINKS {
                        self.refused = true;
                        return;
                    }
                }
                if let Some(target) = target {
                    self.resolved.clear();
                    self.walk(target);
                }
            }
        }
    }

    /// The path walked so far (`/` at the root).
    fn path(&self) -> &[u8] {
        if self.resolved.is_empty() {
            b"/"
        } else {
            &self.resolved
        }
    }
}

// The components that `Walk::step` takes on this thread (the test of `log_name`'s cost).
#[cfg(test)]
thread_local! {
    static WALK_STEPS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// What the `-logfile` name `name` gives on the command line on Linux (No NCBI counterpart: LOSAT
/// Web has no file system; line WQA's sweep, D11, 2026-11-20 and its rounds 2 to 7), on the
/// model's host as its first user (`USER_DIRECTORIES`). The name is
/// read as NCBI reads it on Unix (`ncbi_environment::log_file_name`, the command line's reading
/// too) in a working directory of the user's own (`WORKING_TREE`) on a file system with free
/// space: `logs/run.log`, ` run.log` (`run.log`), `logs/` (the file `logs`) and `/tmp/run.log`
/// are files. A directory of '/'s only (`//run.log`) is none, as `CDirEntry::Reset` leaves it
/// empty. The system then follows the links of Linux (`LINKS`, at most `MAX_LINKS`;
/// `/proc/self/cwd` is the working directory) in the directory that NCBI checks, in the path that
/// it opens and in the directory whose free space it reads (`Walk`, once over the name), and
/// refuses a component of more than
/// 255 bytes and a path of 4096 or more (`fits`: `aaa...`, a directory `x/../x/../...` of 4096
/// bytes). The directory that NCBI checks (the name's as written) is walked as the system walks
/// it (`directory_at`): each directory that it passes may be one of root's that the model does
/// not know, which a name passes on its way back out (`/etc/default/../../tmp`, as on the command
/// line), but none that is not there (`/data`, a file: `/dev/null/x`, `/dev/stdout/..`) nor one
/// that the user cannot search (`searchable`: `/root/..`, `/var/cache/ldconfig/..`); the
/// directory where it ends must be one that the model knows, root's or one where the user makes
/// entries, which the user made (`logs`, `/tmp/a`, `/home/me/logs`); not `/bin/..` (`/usr`, then
/// `/usr/tmp`) or `/var/lock/..` (`/run`, then `/run/tmp`). The opened path is a directory, which
/// NCBI cannot open (`cerr`), where it is the working directory or one of its ancestors, or a
/// directory that the walk passed or one of its ancestors (`.`, `logs/.`, `a/b/..`, `/tmp/`,
/// `/proc/self/cwd`); a directory of root's opens to no user (`opens`: `/tmp`, `/proc/self/fd`);
/// NCBI cannot open a path that the system does not open to a user either (`opens`: `/dev/fd/3`,
/// `/proc/self/status`, `/sys/...`, `/dev/pts/0`, `/dev/tty`, a name in a directory of root's:
/// `/run.log`, `/etc/run.log`, `/var/log/run.log`, `/dev/run.log`, `/dev/kmsg`). A file opened in
/// a file system that reports no free space (`NO_FREE_SPACE`) gets no message, as NCBI checks
/// the space before it writes (`/dev/fd/2`, `/proc/self/fd/1`, `/proc/self/comm`,
/// `/dev/mqueue/q`, `/dev/fd/../stderr`, whose directory is `/proc/self`). Otherwise
/// `/dev/stderr` is the standard error, also written otherwise (`//dev/stderr`); the standard
/// input is the null device (`/dev/stdin`), as for a command line run without one; the devices
/// whose writes are lost (`LOST_DEVICES`: `/dev/ptmx`) give no message; `/dev/stdout` is a file
/// of that name, whose bytes the command line writes into its standard output among the
/// report's, `/proc/self/cwd/run.log` the file `run.log` of the working directory, and
/// `/var/lock/run.log` a file in `/run/lock`, each under the name as typed. The model's command
/// line is a user's (not root), with no terminal, its standard output and error being pipes, as a
/// shell's `|` gives.
fn log_name(name: &[u8]) -> LogName {
    // A working directory deeper than the name's ".." reach (`WORKING_TREE`).
    let climbs = crate::ncbi_str::truncate_spaces(name)
        .split(|byte| *byte == b'/')
        .filter(|part| *part == b"..")
        .count();
    let cwd = WORKING_TREE.repeat(climbs + 1);
    let file = crate::blastinput::ncbi_environment::log_file_name(
        crate::blastinput::ncbi_environment::PathRules::Unix,
        name,
        Some(&cwd),
    )
    .expect("NCBI's Unix reading reads every name");
    // The directories that the walk of the checked directory passes, each a directory of the
    // system: the opened path is none of them nor an ancestor of one.
    let mut directories: Vec<Vec<u8>> = Vec::new();
    if let Some(directory) = &file.directory {
        if directory.is_empty() || !fits(directory) {
            return LogName::Unopened;
        }
        let mut walk = Walk::at(directory, &cwd);
        for part in directory.split(|byte| *byte == b'/') {
            walk.step(part);
            if walk.refused || directory_at(walk.path()) == Directory::None {
                return LogName::Unopened;
            }
            if directories.last().map(Vec::as_slice) != Some(walk.path()) {
                directories.push(walk.path().to_vec());
            }
        }
        if directory_at(walk.path()) != Directory::Known {
            return LogName::Unopened;
        }
    }
    directories.push(cwd.clone());
    let Some(path) = Walk::resolve(&file.path, &cwd) else {
        return LogName::Unopened;
    };
    if file.path.is_empty()
        || !fits(&file.path)
        || directories.iter().any(|directory| within(directory, &path))
        || !opens(&path)
    {
        return LogName::Unopened;
    }
    // A directory whose free space the system does not read keeps the file (NCBI ignores the
    // error); its walk is a part of the checked directory's, which the system took.
    let space = Walk::resolve(&file.space, &cwd).unwrap_or_default();
    if NO_FREE_SPACE.iter().any(|system| within(&space, system)) {
        return LogName::Nowhere;
    }
    let stream = |descriptor: &[u8]| path.strip_prefix(SELF) == Some(descriptor);
    if stream(b"/fd/2".as_slice()) {
        LogName::StandardError
    } else if stream(b"/fd/0".as_slice()) || LOST_DEVICES.contains(&path.as_slice()) {
        LogName::Nowhere
    } else {
        LogName::File
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
}

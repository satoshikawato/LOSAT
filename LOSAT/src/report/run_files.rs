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

/// The directories of LOSAT Web's Linux where every user makes files (mode 1777, `/var/crash`
/// 3777): `/tmp` and `/var/tmp` of the Filesystem Hierarchy Standard, `/dev/shm`, Ubuntu's
/// `/run/lock` (also `/var/lock`, `LINKS`) and Ubuntu's `/var/crash` (line WQA's D11 rounds 3 and
/// 4, 2026-11-20: NCBI BLAST+ 2.17.0 and LOSAT's command line write their log there as a user on
/// Ubuntu 24.04). `/dev/mqueue` (1777 too) holds message queues only (`opens`).
const SHARED_DIRECTORIES: [&[u8]; 5] = [
    b"/tmp",
    b"/var/tmp",
    b"/dev/shm",
    b"/run/lock",
    b"/var/crash",
];

/// Whether the user of LOSAT Web's Linux makes entries in `directory`, a path without links: in
/// the working directory's tree (`WORKING_TREE`), in a user's home (`/home/<user>`, which LOSAT
/// Web takes as the user's own, not knowing the user's name) and in the directories where every
/// user makes files (`SHARED_DIRECTORIES`), each with the directories that the user makes in it.
/// Every other directory is root's, as only root makes entries in `/` (line WQA's D11 round 4,
/// 2026-11-20, Owner-delegated: one rule of owners in place of lists of root's directories):
/// `/`, `/data`, `/home`, `/etc`, `/var/log`, `/run/user/<uid>` (a user's under systemd only).
fn user_makes(directory: &[u8]) -> bool {
    within(directory, WORKING_TREE)
        || directory
            .strip_prefix(b"/home/")
            .is_some_and(|user| !user.is_empty())
        || SHARED_DIRECTORIES
            .iter()
            .any(|shared| within(directory, shared))
}

/// The directories at the root of LOSAT Web's Linux, an Ubuntu of merged `/usr` (`directory_at`):
/// those that the Filesystem Hierarchy Standard 3.0 requires there and its Linux annex (`/proc`,
/// `/sys`), with `/home` and `/root` (`/bin`, `/lib`, `/lib64` and `/sbin` are links, `LINKS`).
/// Another name at the root is not there (`/data`, `/snap`, `/lost+found`, which only root makes).
const ROOT_DIRECTORIES: [&[u8]; 15] = [
    b"/boot", b"/dev", b"/etc", b"/home", b"/media", b"/mnt", b"/opt", b"/proc", b"/root", b"/run",
    b"/srv", b"/sys", b"/tmp", b"/usr", b"/var",
];

/// The links of LOSAT Web's Linux and their targets (`log_name`), besides `/proc/self/cwd` and
/// `/proc/thread-self/cwd`, the working directory: merged `/usr`'s, the Filesystem Hierarchy
/// Standard's `/var/lock` and `/var/run`, the standard streams' (`/dev/fd`, `/dev/stdin`,
/// `/dev/stdout`, `/dev/stderr`) and procfs's (a thread's descriptors are the process's).
const LINKS: [(&[u8], &[u8]); 13] = [
    (b"/bin", b"/usr/bin"),
    (b"/lib", b"/usr/lib"),
    (b"/lib64", b"/usr/lib64"),
    (b"/sbin", b"/usr/sbin"),
    (b"/var/lock", b"/run/lock"),
    (b"/var/run", b"/run"),
    (b"/dev/fd", b"/proc/self/fd"),
    (b"/dev/stdin", b"/proc/self/fd/0"),
    (b"/dev/stdout", b"/proc/self/fd/1"),
    (b"/dev/stderr", b"/proc/self/fd/2"),
    (b"/proc/thread-self/fd", b"/proc/self/fd"),
    (b"/proc/self/root", b"/"),
    (b"/proc/thread-self/root", b"/"),
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
/// `open(O_WRONLY | O_APPEND | O_CREAT)` opens (line WQA's sweep, D11 round 2, 2026-11-20), and
/// the process's `autogroup` (mode 0644), which the kernels of Ubuntu 24.04 have (automatic group
/// scheduling; Linux 6.8: fs/proc/base.c, `tgid_base_stuff`) and the measured one did not (line
/// WQA's D11 round 5, 2026-11-21). The file descriptors 0 to 2 (`fd/0`...) are the standard
/// streams.
const PROCESS_FILES: [&[u8]; 19] = [
    b"attr/current",
    b"attr/exec",
    b"attr/fscreate",
    b"attr/keycreate",
    b"attr/sockcreate",
    b"autogroup",
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

/// Whether `path`, a path without links, is one of the files that the model knows (`opens`): a
/// device of `LOST_DEVICES`, a standard stream (`/proc/self/fd/0` to `2`) or a file of procfs's
/// `PROCESS_FILES` and `THREAD_FILES`. Such a file is no directory (`/dev/null/x`,
/// `/dev/stdout/../stderr`, `/dev/fd/2/..`).
fn is_known_file(path: &[u8]) -> bool {
    LOST_DEVICES.contains(&path)
        || path.strip_prefix(b"/proc/self/").is_some_and(|file| {
            matches!(file, b"fd/0" | b"fd/1" | b"fd/2") || PROCESS_FILES.contains(&file)
        })
        || path
            .strip_prefix(b"/proc/thread-self/")
            .is_some_and(|file| THREAD_FILES.contains(&file))
}

/// Whether `directory` holds a file that the model knows (`is_known_file`): `/dev`, `/dev/net`,
/// `/proc/self` and its `fd` and `attr`, `/proc/thread-self` and its `attr`.
fn holds_known_file(directory: &[u8]) -> bool {
    let holds = |path: &[u8]| path.len() > directory.len() && within(path, directory);
    LOST_DEVICES.iter().any(|device| holds(device))
        || [b"fd/0".as_slice()]
            .iter()
            .chain(&PROCESS_FILES)
            .any(|file| holds(&[b"/proc/self/".as_slice(), file].concat()))
        || THREAD_FILES
            .iter()
            .any(|file| holds(&[b"/proc/thread-self/".as_slice(), file].concat()))
}

/// What LOSAT Web's Linux has at `directory`, a path without links (`log_name`).
#[derive(Debug, PartialEq, Eq)]
enum Directory {
    /// A directory that the model knows: `/`, one where the user makes entries (`user_makes`),
    /// which the model takes as made where a name gives it (`logs`, `/tmp/a`; the working
    /// directory and its ancestors), a directory at the root (`ROOT_DIRECTORIES`), a mount of
    /// `NO_FREE_SPACE` (`/dev/pts`) and a directory of a file that the model knows (`/dev/net`,
    /// `/proc/self/fd`).
    Known,
    /// Any other path in root's directories, which may be a directory there or not
    /// (`/etc/default`, `/usr/tmp`): the user cannot change them, and the model does not list
    /// them.
    Unknown,
    /// No directory: a file that the model knows (`/dev/null`, `/proc/self/fd/1`), a name at the
    /// root that is not there (`/data`, `/snap`), or a path in `/root`, whose entries a user
    /// cannot reach (mode 0700).
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
        || within(directory, b"/root")
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
    LOST_DEVICES.contains(&path)
        || user_makes(&path[..path.iter().rposition(|byte| *byte == b'/').unwrap_or(0)])
}

/// A path of LOSAT Web's Linux walked component by component as the system walks it (No NCBI
/// counterpart: the system's resolution, `log_name`): `.` and empty components stay, `..` goes to
/// the parent of the path walked so far (`/` stays), and a link (`LINKS`; `/proc/self/cwd` and
/// `/proc/thread-self/cwd`, the working directory `cwd`) is replaced by its target, which is
/// walked from `/`. Each step costs the component's length, so a name is read in time linear in
/// its length (line WQA's D11 round 5, 2026-11-21: resolving each prefix again was cubic, a minute
/// for a name of 4 KB in LOSAT Web's serial engine).
struct Walk<'a> {
    /// The path walked so far, without links; empty for `/`.
    resolved: Vec<u8>,
    /// The working directory, `/proc/self/cwd`'s target.
    cwd: &'a [u8],
}

impl<'a> Walk<'a> {
    /// The walk of `path` before its first component: at `/` for an absolute path, else at the
    /// working directory `cwd`.
    fn at(path: &[u8], cwd: &'a [u8]) -> Self {
        let mut walk = Walk {
            resolved: Vec::new(),
            cwd,
        };
        if path.first() != Some(&b'/') {
            walk.walk(cwd);
        }
        walk
    }

    /// `path` with its links followed and its `.` and `..` taken.
    fn resolve(path: &[u8], cwd: &'a [u8]) -> Vec<u8> {
        let mut walk = Walk::at(path, cwd);
        walk.walk(path);
        walk.path().to_vec()
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
        match part {
            b"" | b"." => {}
            b".." => {
                let parent = self.resolved.iter().rposition(|byte| *byte == b'/');
                self.resolved.truncate(parent.unwrap_or(0));
            }
            _ => {
                self.resolved.push(b'/');
                self.resolved.extend_from_slice(part);
                let target = match self.resolved.as_slice() {
                    b"/proc/self/cwd" | b"/proc/thread-self/cwd" => Some(self.cwd),
                    path => LINKS
                        .iter()
                        .find(|(link, _)| *link == path)
                        .map(|(_, target)| *target),
                };
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
/// Web has no file system; line WQA's sweep, D11, 2026-11-20 and its rounds 2 to 5). The name is
/// read as NCBI reads it on Unix (`ncbi_environment::log_file_name`, the command line's reading
/// too) in a working directory of the user's own (`WORKING_TREE`) on a file system with free
/// space: `logs/run.log`, ` run.log` (`run.log`), `logs/` (the file `logs`) and `/tmp/run.log`
/// are files. A directory of '/'s only (`//run.log`) is none, as `CDirEntry::Reset` leaves it
/// empty. The system then follows the links of Linux (`LINKS`; `/proc/self/cwd` is the working
/// directory) in the directory that NCBI checks, in the path that it opens and in the directory
/// whose free space it reads (`Walk`, once over the name), and refuses a component of more than
/// 255 bytes and a path of 4096 or more (`fits`: `aaa...`, a directory `x/../x/../...` of 4096
/// bytes). The directory that NCBI checks (the name's as written) is walked as the system walks
/// it (`directory_at`): each directory that it passes may be one of root's that the model does
/// not know, which a name passes on its way back out (`/etc/default/../../tmp`, as on the command
/// line), but none that is not there (`/data`, a file: `/dev/null/x`, `/dev/stdout/..`); the
/// directory where it ends must be one that the model knows, root's or one where the user makes
/// entries, which the user made (`logs`, `/tmp/a`, `/home/me/logs`); not `/bin/..` (`/usr`, then
/// `/usr/tmp`) or `/var/lock/..` (`/run`, then `/run/tmp`). The opened path is a directory, which
/// NCBI cannot open (`cerr`), where it is the working directory or one of its ancestors, or a
/// directory that the walk passed or one of its ancestors (`.`, `logs/.`, `a/b/..`, `/tmp/`,
/// `/proc/self/cwd`); a directory of root's opens to no user (`opens`: `/tmp`, `/proc/self/fd`); NCBI
/// cannot open a path that the system does not open to a user either (`opens`: `/dev/fd/3`,
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
            if directory_at(walk.path()) == Directory::None {
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
    let path = Walk::resolve(&file.path, &cwd);
    if file.path.is_empty()
        || !fits(&file.path)
        || directories.iter().any(|directory| within(directory, &path))
        || !opens(&path)
    {
        return LogName::Unopened;
    }
    let space = Walk::resolve(&file.space, &cwd);
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
            "/run/user/1000/run.log",
            "/var/run/user/1000/run.log",
            "/run/shm/run.log",
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
            // not there (a file, a name at the root that is not there, an entry of `/root`).
            "/lib/../tmp/run.log",
            "/var/run/lock/../run.log",
            "/etc/default/../tmp/run.log",
            "/dev/null/../../tmp/run.log",
            "/dev/fd/2/../../../tmp/run.log",
            "/root/x/../../tmp/run.log",
            "/data/x/../../tmp/run.log",
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
            // Round 5: Ubuntu 24.04's kernels have `autogroup` (automatic group scheduling).
            "/proc/self/autogroup",
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

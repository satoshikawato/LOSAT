//! makeblastdb's local time on Windows against Microsoft's C runtime (UCRT), which NCBI BLAST+'s
//! Windows build reads it from: `CTime(CTime::eCurrent)` calls `localtime`, which takes the zone
//! from `TZ` in UCRT's form or from the system (`LOSAT/src/writedb/windows_zone.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1722-1734 (`x_SetTime`: `localtime`),
//! c++/src/objtools/blast/seqdb_writer/build_db.cpp:1104-1105 (the log's current time),
//! c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:88-97 (the database's date).
//!
//! The oracle is the UCRT of the machine that runs the test: this test binary starts itself
//! (`ucrt_child`) with `TZ` in its environment, as makeblastdb is started, and prints what
//! `_localtime64_s` gives. Each process sets UCRT's zone once, at its first `localtime`.
//!
//! The instants of the vectors lie inside `localtime`'s first branch, more than three days after
//! 1970-01-01 00:00:00 and more than three days before its last second, 3001-01-19 07:59:59 UTC
//! (ucrt/time/localtime.cpp:58; corecrt_internal_time.h:98). That is where LOSAT's port is
//! defined. The branch for the three days at either end is not ported: the current time is
//! never there, and `windows_zone.rs` gives UTC instead of one UCRT version's edge behaviour.
//! The `ucrt/time` sources of the Windows SDKs 10.0.22621.0 and 10.0.26100.0 are identical.

#![cfg(windows)]

use std::ffi::OsString;
use std::os::windows::ffi::OsStringExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// UCRT's zones and local times with `TZ` set, from `ucrtbase.dll` (scripts/ucrt_tz_vectors.py);
/// the unit tests of `windows_zone.rs` compare LOSAT's port with the same file.
const VECTORS: &str = include_str!("unit/ucrt_local_time_vectors.tsv");

/// The variable that makes `ucrt_child` print (the instants, separated by spaces).
const CHILD_INSTANTS: &str = "LOSAT_TEST_UCRT_INSTANTS";

#[repr(C)]
#[derive(Default)]
struct Tm {
    tm_sec: i32,
    tm_min: i32,
    tm_hour: i32,
    tm_mday: i32,
    tm_mon: i32,
    tm_year: i32,
    tm_wday: i32,
    tm_yday: i32,
    tm_isdst: i32,
}

// The C runtime that Rust programs link on Windows (`ucrtbase.dll`); `long` is 32 bits.
extern "C" {
    fn _localtime64_s(tm: *mut Tm, time: *const i64) -> i32;
    fn _get_timezone(seconds: *mut i32) -> i32;
    fn _get_daylight(hours: *mut i32) -> i32;
    fn _get_dstbias(seconds: *mut i32) -> i32;
}

/// Days from 1970-01-01 to the first day of `year` (Gregorian).
fn days_before_year(year: i64) -> i64 {
    let y = year - 1;
    365 * (year - 1970) + (y / 4 - y / 100 + y / 400) - (1969 / 4 - 1969 / 100 + 1969 / 400)
}

/// UCRT's local time of `t` in this process: the offset from UTC and the calendar fields.
fn ucrt_local(t: i64) -> Option<(i64, Tm)> {
    let mut tm = Tm::default();
    // SAFETY: `tm` and `t` live for the call, which writes `tm` only.
    if unsafe { _localtime64_s(&mut tm, &t) } != 0 {
        return None;
    }
    let local = (days_before_year(1900 + i64::from(tm.tm_year)) + i64::from(tm.tm_yday)) * 86_400
        + i64::from(tm.tm_hour) * 3600
        + i64::from(tm.tm_min) * 60
        + i64::from(tm.tm_sec);
    Some((local - t, tm))
}

const MONTHS: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

/// The two forms NCBI writes (`CTime::AsString()` in the log, "b d, Y  H:m P" for the date).
fn ncbi_forms(tm: &Tm) -> (String, String) {
    let log = format!(
        "{:02}/{:02}/{:04} {:02}:{:02}:{:02}",
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_year + 1900,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec
    );
    let date = format!(
        "{} {}, {:04}  {}:{:02} {}",
        MONTHS[tm.tm_mon as usize],
        tm.tm_mday,
        tm.tm_year + 1900,
        (tm.tm_hour + 11) % 12 + 1,
        tm.tm_min,
        if tm.tm_hour < 12 { "AM" } else { "PM" }
    );
    (log, date)
}

/// Started by `in_child` only: prints UCRT's zone and the local times of the instants of
/// `LOSAT_TEST_UCRT_INSTANTS` for this process's `TZ`.
#[test]
fn ucrt_child() {
    let Some(instants) = std::env::var_os(CHILD_INSTANTS) else {
        return;
    };
    // With one test thread libtest writes "test windows_local_time::ucrt_child ... " before the
    // test runs and ends that line only after it (library/test/src/formatters/pretty.rs,
    // `write_test_start`): the records start on a line of their own.
    println!();
    for t in instants.to_str().unwrap().split(' ') {
        let t: i64 = t.parse().unwrap();
        match ucrt_local(t) {
            Some((offset, tm)) => {
                let (log, date) = ncbi_forms(&tm);
                println!("UCRT\ttime\t{t}\t{offset}\t{log}\t{date}");
            }
            None => println!("UCRT\ttime\t{t}\tnone"),
        }
    }
    let (mut timezone, mut daylight, mut dstbias) = (0, 0, 0);
    // SAFETY: each call writes the one value it is given.
    unsafe {
        _get_timezone(&mut timezone);
        _get_daylight(&mut daylight);
        _get_dstbias(&mut dstbias);
    }
    println!("UCRT\tzone\t{timezone}\t{daylight}\t{dstbias}");
}

/// The output of `ucrt_child` in a process of this binary.
struct Child {
    /// The records (the lines that start with `UCRT\t`, without it).
    lines: Vec<String>,
    /// The whole stdout and stderr, for failure messages.
    raw: String,
}

/// `ucrt_child` in a process of this binary with `TZ` (`None`: unset) and `instants`.
fn in_child(tz: Option<&[u16]>, instants: &[i64]) -> Child {
    let mut command = Command::new(std::env::current_exe().unwrap());
    command.args([
        "windows_local_time::ucrt_child",
        "--exact",
        "--nocapture",
        "--test-threads",
        "1",
    ]);
    match tz {
        Some(units) => command.env("TZ", OsString::from_wide(units)),
        None => command.env_remove("TZ"),
    };
    let list: Vec<String> = instants.iter().map(i64::to_string).collect();
    command.env(CHILD_INSTANTS, list.join(" "));
    let out = command.output().unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    let raw = format!(
        "child {:?} with TZ {tz:?}, {CHILD_INSTANTS}={:?}, {}\n--- stdout\n{stdout}--- stderr\n{}",
        command.get_program(),
        list.join(" "),
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out.status.success(), "{raw}");
    let lines = stdout
        .lines()
        .filter_map(|l| l.strip_prefix("UCRT\t").map(str::to_string))
        .collect();
    Child { lines, raw }
}

/// A `TZ` value of the vector file as UTF-16 units: `\\` and `\u{XXXX}` escapes.
fn vector_tz(text: &str) -> Vec<u16> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(c) = rest.chars().next() {
        if let Some(r) = rest.strip_prefix("\\\\") {
            out.push(u16::from(b'\\'));
            rest = r;
        } else if let Some(r) = rest.strip_prefix("\\u{") {
            let end = r.find('}').unwrap();
            out.push(u16::from_str_radix(&r[..end], 16).unwrap());
            rest = &r[end + 1..];
        } else {
            out.push(u16::try_from(u32::from(c)).unwrap());
            rest = &rest[c.len_utf8()..];
        }
    }
    out
}

// The UCRT of this machine gives the zones and times of the vectors, which LOSAT's port gives
// (`windows_zone.rs`, `ucrt_tz_values_match_the_c_runtime`).
#[test]
fn ucrt_of_this_machine_gives_the_vectors() {
    let mut zones: Vec<(String, String, Vec<(i64, String)>)> = Vec::new();
    for line in VECTORS.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "zone" => zones.push((f[1].to_string(), f[2..5].join("\t"), Vec::new())),
            "time" => {
                let zone = zones.last_mut().unwrap();
                assert_eq!(zone.0, f[1], "{line}");
                zone.2.push((f[2].parse().unwrap(), f[3].to_string()));
            }
            _ => {}
        }
    }
    assert!(zones.len() >= 70, "{} zones", zones.len());
    let mut bad = Vec::new();
    // The child's whole output for the first zone that differs.
    let mut first_raw = None;
    for (tz, globals, times) in &zones {
        let instants: Vec<i64> = times.iter().map(|(t, _)| *t).collect();
        let child = in_child(Some(&vector_tz(tz)), &instants);
        let bad_before = bad.len();
        let want_zone = format!("zone\t{globals}");
        if !child.lines.contains(&want_zone) {
            bad.push(format!("{tz}: {want_zone:?}, UCRT {:?}", child.lines));
        }
        for (t, want) in times {
            let got = child
                .lines
                .iter()
                .find_map(|l| l.strip_prefix(&format!("time\t{t}\t")))
                .map(|rest| rest.split('\t').next().unwrap().to_string());
            if got.as_deref() != Some(want.as_str()) {
                bad.push(format!("{tz} at {t}: {want}, UCRT {got:?}"));
            }
        }
        if bad.len() > bad_before && first_raw.is_none() {
            first_raw = Some(child.raw);
        }
    }
    assert!(
        bad.is_empty(),
        "\n{}\n{}",
        bad.join("\n"),
        first_raw.unwrap_or_default()
    );
}

fn work_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-windows-local-time-{tag}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn now() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs() as i64
}

/// Runs makeblastdb and blastdbcmd -info in `dir` with `TZ`; returns the log's current time,
/// the database's date and the seconds the run took place in.
fn build(dir: &Path, tz: Option<&[u16]>) -> (String, String, Vec<i64>) {
    // Test protein: the first residues of E. coli AsnA (RefSeq NP_418200.1).
    std::fs::write(dir.join("in.fa"), ">p1\nMKTAYIAKQRQISFVKSHFSRQ\n").unwrap();
    let losat = env!("CARGO_BIN_EXE_losat");
    let mut command = Command::new(losat);
    command.current_dir(dir).args([
        "makeblastdb",
        "-in",
        "in.fa",
        "-dbtype",
        "prot",
        "-out",
        "db",
    ]);
    match tz {
        Some(units) => command.env("TZ", OsString::from_wide(units)),
        None => command.env_remove("TZ"),
    };
    let start = now();
    let out = command.output().unwrap();
    let end = now();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let log = stdout
        .lines()
        .find_map(|l| l.strip_prefix("Building a new DB, current time: "))
        .unwrap_or_else(|| panic!("no current time in {stdout:?}"))
        .to_string();
    let info = Command::new(losat)
        .current_dir(dir)
        .args(["blastdbcmd", "-db", "db", "-dbtype", "prot", "-info"])
        .output()
        .unwrap();
    assert!(
        info.status.success(),
        "{}",
        String::from_utf8_lossy(&info.stderr)
    );
    let info = String::from_utf8(info.stdout).unwrap();
    let date = info
        .lines()
        .find_map(|l| l.strip_prefix("Date: "))
        .and_then(|l| l.split('\t').next())
        .unwrap_or_else(|| panic!("no date in {info:?}"))
        .to_string();
    (log, date, (start..=end).collect())
}

// makeblastdb's log time and the database's date are UCRT's local time of a second of the run,
// with `TZ` in UCRT's forms (an offset, a daylight-time name, a white space that `wcstol` skips,
// an offset of 68 years) and unset (the system zone; UTC on GitHub's runners).
#[test]
fn makeblastdb_writes_the_c_runtimes_local_time() {
    let values: [Option<&str>; 6] = [
        Some("JST-9"),
        Some("PST8PDT"),
        Some("AAA-12:34:56BBB"),
        Some("AAA\u{3000}8"),
        Some("AAA-596523BBB"),
        None,
    ];
    for (i, tz) in values.into_iter().enumerate() {
        let units: Option<Vec<u16>> = tz.map(|v| v.encode_utf16().collect());
        let dir = work_dir(&i.to_string());
        let (log, date, seconds) = build(&dir, units.as_deref());
        let child = in_child(units.as_deref(), &seconds);
        let forms: Vec<(&str, &str)> = child
            .lines
            .iter()
            .filter_map(|l| l.strip_prefix("time\t"))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                assert_eq!(f.len(), 4, "{tz:?}: no local time for {l}\n{}", child.raw);
                (f[2], f[3])
            })
            .collect();
        assert_eq!(forms.len(), seconds.len(), "{tz:?}\n{}", child.raw);
        assert!(
            forms.iter().any(|(l, _)| *l == log),
            "{tz:?}: log {log:?}, UCRT {forms:?}\n{}",
            child.raw
        );
        assert!(
            forms.iter().any(|(_, d)| *d == date),
            "{tz:?}: date {date:?}, UCRT {forms:?}\n{}",
            child.raw
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}

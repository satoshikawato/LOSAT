//! `losat makeblastdb` on a full disk: the index file of a volume (and every other output file
//! that passes the limit) is cut short, against the outcomes of makeblastdb 2.17.0 (the official
//! Linux x64 release) under the same limit, frozen in `fixtures/makeblastdb_failure/fsize_rows.txt`.
//!
//! makeblastdb opens the database it built (`CSeqDB`) before it writes the metadata file, and
//! NCBI's reader takes the index file as it finds it: a file that ends early reads as zeros to the
//! end of its last page, so a file of 1 or 2 bytes is "Not a valid version 4 or 5 database.", one
//! of 7 bytes holds the wrong sequence type (a protein database), one of 8 or 9 bytes is a database
//! with no title, no date and no sequences, and one that ends inside the date is a date that does
//! not read: makeblastdb ends with the exception of `CTime`, exit status 255, and writes no
//! metadata file. An index file of no bytes is mapped as a null pointer and NCBI's process is
//! killed by SIGSEGV (rc -11 in the fixture), which LOSAT refuses ("NCBI crashes on it"). A header
//! field that reaches past the file's last page reads whatever memory follows, so what NCBI does
//! depends on its memory (it need not crash): LOSAT refuses that too, in other words, which no row
//! of the fixture covers (it needs a field of 4 KiB; `seqdb::damaged_files_tests` has it). An alias
//! file that is cut short names no database ("No database names were found in alias file [...]",
//! without the file's path).
//!
//! Version 5 databases are not in the fixture: under `prlimit --fsize` the LMDB file fails first
//! (`mdb_env_open: File too large`) in both programs, before the index file is read. They were
//! compared with a write limiter in front of the index file (the probes of line LIB-O, item O1e).
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1296-1318,
//! c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:134-236 (`CSeqDBIdxFile`).
//!
//! The index files and the alias file hold the build's date ("Nov 21, 2026  5:07 AM"), and a cut at a
//! fixed number of bytes lands on another byte of the file when the date has one byte more or less.
//! The date's length depends on the digits of the hour and of the day of the month (the month's
//! name has 3 letters in NCBI's table and the year 4 digits, so they do not change it). The hour is
//! pinned by the zone the test gives LOSAT (`pinned_zone`: the local hour is 5 or 17, one digit on
//! the 12-hour clock); the day cannot be pinned (a zone moves the local date by less than a day), so
//! the fixture holds the outcomes of NCBI for both classes of day (`DayClass`), frozen under a fixed
//! clock (`faketime` shim, time zone UTC) at 2026-11-21 05:07 (two digits) and 2026-11-05 05:07
//! (one digit): 15 of the 18 cases are the same in both, the three that differ (m6 limits 95, 100
//! and 111, the cut lands inside the date of the alias file or of a volume's index) are in the
//! fixture twice, the rows of the one-digit day tagged `d1`. The test picks the rows by the length
//! of the date that LOSAT writes (`date_class`, unit-tested on the edge dates below).
//!
//! The limit is `prlimit --fsize` with SIGXFSZ ignored (the child of `sh` keeps the ignored signal);
//! the test is skipped where `prlimit` is missing. The cases that cut the date string of the index
//! inside the digits of the day or the time (the outcome of which depends on the day and the hour)
//! are in the date table of `seqdb::date`'s tests.

#![cfg(target_os = "linux")]

use std::collections::BTreeMap;
use std::process::Command;

/// The input of the kind n4: 6 nucleotide sequences of 300 bases (the input of the frozen run).
const NUC: &str = include_str!("fixtures/makeblastdb_failure/n4_input.fa");

const ROWS: &str = include_str!("fixtures/makeblastdb_failure/fsize_rows.txt");

const RICH_P: &str = ">gi|1001|ref|NP_000001.1| alpha\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n>gnl|ti|12345 trace one\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA\n>sp|P55555|CCC_HUMAN Gamma\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARR\n";

/// 130 proteins of 60 residues.
fn many() -> String {
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut out = String::new();
    for i in 0..130usize {
        out.push_str(&format!(">sp|P{i:05}|S{i} desc {i}\n"));
        out.extend((0..60usize).map(|j| AA[(i * 7 + j * 3 + j * j % 11) % 20] as char));
        out.push('\n');
    }
    out
}

/// A row of the fixture: what NCBI's makeblastdb did.
#[derive(Clone, Default)]
struct Row {
    kind: String,
    limit: u64,
    /// The row was frozen at a one-digit day of the month (tag `d1`); an untagged row holds for
    /// every day unless a `d1` row of the same case exists.
    one_digit_day: bool,
    rc: i32,
    err: String,
    out: String,
    files: String,
    pjs: String,
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

fn rows() -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    for line in ROWS.lines().filter(|line| !line.starts_with('#')) {
        let (field, value) = line.split_once('\t').unwrap();
        if field == "case" {
            let mut fields = value.split('\t');
            let kind = fields.next().unwrap();
            let limit = fields.next().unwrap();
            let one_digit_day = match fields.next() {
                None => false,
                Some("d1") => true,
                Some(other) => panic!("unknown tag {other}"),
            };
            rows.push(Row {
                kind: kind.to_string(),
                limit: limit.parse().unwrap(),
                one_digit_day,
                ..Row::default()
            });
            continue;
        }
        let row = rows.last_mut().unwrap();
        match field {
            "rc" => row.rc = value.parse().unwrap(),
            "err" => row.err = unescape(value),
            "out" => row.out = unescape(value),
            "files" => row.files = value.to_string(),
            "pjs" => row.pjs = unescape(value),
            other => panic!("unknown field {other}"),
        }
    }
    rows
}

/// The rows to compare for a date of `class`: the untagged rows, each replaced by the row of the
/// same case tagged `d1` where the class is `OneDigit` and such a row exists.
fn rows_for(class: DayClass) -> Vec<Row> {
    let all = rows();
    let mut chosen: Vec<Row> = Vec::new();
    for row in all.iter().filter(|row| !row.one_digit_day) {
        let tagged = all.iter().find(|other| {
            other.one_digit_day && other.kind == row.kind && other.limit == row.limit
        });
        chosen.push(match (class, tagged) {
            (DayClass::OneDigit, Some(tagged)) => tagged.clone(),
            _ => row.clone(),
        });
    }
    chosen
}

fn has_prlimit() -> bool {
    Command::new("prlimit")
        .arg("--version")
        .output()
        .is_ok_and(|output| output.status.success())
}

/// The time and the duration of the log and the build's date of the metadata masked.
fn masked(text: &str, d: &str) -> String {
    let mut out = text.replace(d, "<D>");
    for (start, end, fill) in [
        ("current time: ", "\n", "<T>"),
        (" sequences in ", " seconds.", "<S>"),
        ("\"last-updated\": \"", "\"", "<ISO>"),
    ] {
        let mut at = 0;
        while let Some(i) = out[at..].find(start) {
            let from = at + i + start.len();
            let Some(j) = out[from..].find(end) else {
                break;
            };
            out.replace_range(from..from + j, fill);
            at = from + fill.len();
        }
    }
    out
}

/// The class of the date of the build by the number of digits of its day of the month: the cut of
/// a file at a fixed byte lands on another byte when the date is one byte longer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DayClass {
    /// Days 1 to 9 ("Nov 5, 2026").
    OneDigit,
    /// Days 10 to 31 ("Nov 21, 2026").
    TwoDigits,
}

/// How many hours ahead of UTC (`UTC-k` in POSIX `TZ`) the zone is whose local hour at `now` (Unix
/// seconds) is 5 or 17: one digit on the 12-hour clock of the date ("5:07 AM"). At that hour a
/// run of the test does not reach midnight, so the day of the date does not change during it.
fn pinned_hours_ahead(now: u64) -> u64 {
    let utc_hour = (now / 3600) % 24;
    (5 + 12 - utc_hour % 12) % 12
}

fn pinned_zone_at(now: u64) -> String {
    format!("UTC-{}", pinned_hours_ahead(now))
}

/// The day of the month (1 to 31) of the Gregorian date `seconds` after 1970-01-01 00:00:00.
fn day_of_month(seconds: u64) -> u32 {
    // Days since 0000-03-01, with the leap day last in the year (the usual civil-date algorithm).
    let z = (seconds / 86_400) as i64 + 719_468;
    let day_of_era = z.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_from_march = (5 * day_of_year + 2) / 153;
    (day_of_year - (153 * month_from_march + 2) / 5 + 1) as u32
}

/// The class of the date that LOSAT writes at `now` (Unix seconds) in the zone `pinned_zone_at(now)`.
fn date_class(now: u64) -> DayClass {
    let local = now + 3600 * pinned_hours_ahead(now);
    if day_of_month(local) < 10 {
        DayClass::OneDigit
    } else {
        DayClass::TwoDigits
    }
}

fn unix_now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_secs())
}

/// The date string that did not parse, in the message of `CTime` ("Time string '...' is too
/// short"), masked: it holds the month's name of the day the fixture was frozen.
fn mask_date(text: &str) -> String {
    let start = "Time string '";
    match text.find(start) {
        Some(i) => {
            let from = i + start.len();
            match text[from..].find('\'') {
                Some(j) => format!("{}<DATE>{}", &text[..from], &text[from + j..]),
                None => text.to_string(),
            }
        }
        None => text.to_string(),
    }
}

fn dated(name: &str) -> bool {
    [".pin", ".nin", ".pjs", ".njs", ".pal", ".nal"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// The number of files of the output directory, then `name:size` of its index, alias and metadata
/// files (`*` for the size of a file that holds the build's date), as the fixture holds them.
fn listing(dir: &std::path::Path) -> String {
    let mut entries = BTreeMap::new();
    for entry in std::fs::read_dir(dir).unwrap().map_while(Result::ok) {
        let name = entry.file_name().to_string_lossy().into_owned();
        let size = entry.metadata().unwrap().len();
        entries.insert(
            name.clone(),
            if dated(&name) {
                "*".to_string()
            } else {
                size.to_string()
            },
        );
    }
    let chosen: Vec<String> = entries
        .iter()
        .filter(|(name, _)| {
            [".pin", ".nin", ".pal", ".nal", ".pjs", ".njs"]
                .iter()
                .any(|ext| name.ends_with(ext))
        })
        .map(|(n, s)| format!("{n}:{s}"))
        .collect();
    format!("{} {}", entries.len(), chosen.join(" "))
}

#[test]
fn a_full_disk_leaves_the_index_as_ncbi_reads_it() {
    if !has_prlimit() {
        eprintln!("prlimit is not installed: skipped");
        return;
    }
    let mut compared = 0;
    let now = unix_now();
    let zone = pinned_zone_at(now);
    for row in rows_for(date_class(now)) {
        let tmp = std::fs::canonicalize(std::env::temp_dir())
            .unwrap()
            .join(format!(
                "losat-mkdb-failed-index-{}-{}-{}",
                row.kind,
                row.limit,
                std::process::id()
            ));
        let _ = std::fs::remove_dir_all(&tmp);
        let out_dir = tmp.join("o");
        std::fs::create_dir_all(&out_dir).unwrap();
        let (input, dbtype, index, args): (String, &str, &str, Vec<&str>) = match row.kind.as_str()
        {
            "p4" => (
                RICH_P.to_string(),
                "prot",
                "pin",
                vec!["-hash_index", "-blastdb_version", "4"],
            ),
            "n4" => (
                NUC.to_string(),
                "nucl",
                "nin",
                vec!["-blastdb_version", "4"],
            ),
            "m6" => (
                many(),
                "prot",
                "pin",
                vec!["-max_file_sz", "2000B", "-blastdb_version", "4"],
            ),
            other => panic!("unknown kind {other}"),
        };
        std::fs::write(tmp.join("in.fa"), input).unwrap();
        let mut command = Command::new("sh");
        command
            .current_dir(&tmp)
            .env("HOME", &tmp)
            .env("TZ", &zone)
            .env_remove("BLASTDB")
            .arg("-c")
            .arg(format!(
                "trap '' XFSZ; exec prlimit --fsize={} -- \"$@\"",
                row.limit
            ))
            .arg("sh")
            .arg(env!("CARGO_BIN_EXE_losat"))
            .args([
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                dbtype,
                "-parse_seqids",
            ])
            .args(&args)
            .args(["-out", "o/db"]);
        let output = command.output().expect("run losat makeblastdb");
        let d = tmp.to_str().unwrap();
        let status = output.status.code().unwrap_or(-1);
        let stderr = masked(&String::from_utf8_lossy(&output.stderr), d);
        let stdout = masked(&String::from_utf8_lossy(&output.stdout), d);
        let files = listing(&out_dir);
        let name = format!("{} limit {}", row.kind, row.limit);
        if row.rc == -11 {
            // NCBI's process is killed by SIGSEGV where it maps the index file of no bytes: LOSAT
            // refuses (explicit rejection, exit status 1) and leaves the files NCBI leaves.
            assert_eq!(status, 1, "{name}");
            assert!(
                stderr.starts_with(&format!("losat makeblastdb: <D>/o/db.{index}: "))
                    && stderr.ends_with("(NCBI crashes on it)\n"),
                "{name}: {stderr}"
            );
        } else {
            assert_eq!(status, row.rc, "{name}");
            let expected_err =
                String::from_utf8(losat::wording::losat_bytes(row.err.as_bytes()).into_owned())
                    .unwrap();
            assert_eq!(mask_date(&stderr), mask_date(&expected_err), "{name}");
        }
        assert_eq!(stdout, row.out, "{name}");
        let pjs: String = ["db.pjs", "db.njs"]
            .iter()
            .filter_map(|n| std::fs::read_to_string(out_dir.join(n)).ok())
            .collect();
        assert_eq!(masked(&pjs, d), row.pjs, "{name}");
        assert_eq!(files, row.files, "{name}");
        let _ = std::fs::remove_dir_all(&tmp);
        compared += 1;
    }
    assert_eq!(compared, 18);
}

/// The edge dates of the class choice. The epoch values are `calendar.timegm` of the UTC date.
#[test]
fn the_date_class_follows_the_local_day_at_the_edges() {
    use DayClass::{OneDigit, TwoDigits};
    // (UTC time, zone, class): `UTC-k` puts the local hour at 5 or 17.
    for (now, zone, class) in [
        // 2026-11-09 04:59:59 UTC: 1 h ahead is 05:59 on day 9.
        (1_794_182_400 + 4 * 3600 + 3599, "UTC-1", OneDigit),
        // 2026-11-09 00:00:00 UTC: 5 h ahead is 05:00 on day 9.
        (1_794_182_400, "UTC-5", OneDigit),
        // 2026-11-09 23:59:59 UTC: the zone is 6 h ahead, local 05:59 on day 10.
        (1_794_268_799, "UTC-6", TwoDigits),
        // 2026-11-10 00:00:00 UTC: local 05:00 on day 10.
        (1_794_268_800, "UTC-5", TwoDigits),
        // 2026-10-31 23:59:59 UTC: local 05:59 on 1 November (6 h ahead).
        (1_793_491_199, "UTC-6", OneDigit),
        // 2026-11-01 00:00:00 UTC.
        (1_793_491_200, "UTC-5", OneDigit),
        // 2026-12-31 23:59:59 UTC: local 05:59 on 1 January 2027.
        (1_798_761_599, "UTC-6", OneDigit),
        // 2027-01-01 00:00:00 UTC.
        (1_798_761_600, "UTC-5", OneDigit),
        // 2024-02-28 23:59:59 UTC: local 05:59 on 29 February 2024 (a leap day).
        (1_709_164_799, "UTC-6", TwoDigits),
        // 2024-02-29 12:00:00 UTC: local 17:00 (UTC+5) on day 29.
        (1_709_208_000, "UTC-5", TwoDigits),
        // 2024-03-01 00:00:00 UTC.
        (1_709_251_200, "UTC-5", OneDigit),
        // 2026-09-09 19:00:00 UTC: 10 h ahead is 05:00 on day 10.
        (1_788_980_400, "UTC-10", TwoDigits),
    ] {
        assert_eq!(pinned_zone_at(now), zone, "zone at {now}");
        assert_eq!(date_class(now), class, "class at {now}");
    }
}

/// Every hour of the years 2024 to 2030 (leap years included): the local hour of the pinned zone is
/// 5 or 17, and the class is the number of digits of the day of a calendar walked day by day.
#[test]
fn the_pinned_zone_and_the_class_hold_on_every_hour_of_every_day() {
    fn month_length(year: u64, month: usize) -> u64 {
        let leap =
            year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400));
        let february = if leap { 29 } else { 28 };
        [31, february, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month] as u64
    }
    // The day of the month of each day from 2024-01-01 (day number 19,723 since the epoch) to the
    // end of 2030, plus one more month.
    let first_day = 19_723u64;
    let mut walked: Vec<u64> = Vec::new();
    let (mut year, mut month, mut day) = (2024u64, 0usize, 1u64);
    while year < 2031 {
        walked.push(day);
        day += 1;
        if day > month_length(year, month) {
            day = 1;
            month += 1;
            if month == 12 {
                month = 0;
                year += 1;
            }
        }
    }
    walked.extend(1..=31);
    let mut checked = 0u64;
    for utc_index in 0..walked.len() - 31 {
        for hour in 0..24u64 {
            let now = (first_day + utc_index as u64) * 86_400 + hour * 3600 + 1234;
            let ahead = pinned_hours_ahead(now);
            let local_hour = (hour + ahead) % 24;
            assert!(
                local_hour == 5 || local_hour == 17,
                "hour {local_hour} at {now}"
            );
            // The local time is at most 11 hours after the UTC time: the same day or the next.
            let local_day = walked[utc_index + ((hour + ahead) / 24) as usize];
            let expected = if local_day < 10 {
                DayClass::OneDigit
            } else {
                DayClass::TwoDigits
            };
            assert_eq!(date_class(now), expected, "at {now}");
            checked += 1;
        }
    }
    assert_eq!(checked, 24 * (365 * 5 + 366 * 2));
}

/// The fixture's `d1` rows each replace an untagged row, and the two classes differ only there.
#[test]
fn the_one_digit_day_rows_replace_untagged_rows() {
    let all = rows();
    let tagged: Vec<&Row> = all.iter().filter(|row| row.one_digit_day).collect();
    assert_eq!(tagged.len(), 3);
    for row in &tagged {
        let base = all
            .iter()
            .filter(|other| {
                !other.one_digit_day && other.kind == row.kind && other.limit == row.limit
            })
            .count();
        assert_eq!(base, 1, "{} {}", row.kind, row.limit);
    }
    assert_eq!(rows_for(DayClass::TwoDigits).len(), 18);
    assert_eq!(rows_for(DayClass::OneDigit).len(), 18);
}

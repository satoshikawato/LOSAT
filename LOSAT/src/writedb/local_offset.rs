//! The local time zone of the clock: how far the local time of an instant is from UTC, as the C
//! library's `localtime` gives it to NCBI's builds (`super::time::LocalTime::now`, the dates that
//! makeblastdb writes into a database and its log).
//!
//! On Windows it is UCRT's `localtime` (`super::windows_zone`). On every other target LOSAT reads
//! the zone that glibc reads for NCBI's Linux build, and the crate jiff reads the zone file (TZif,
//! RFC 8536) or the POSIX rule (Owner 2026-11-20 16:27 JST 「TZ=A」). Where glibc finds the zone
//! was observed with each form of `TZ` (Python's `time.localtime`, the fixture
//! `LOSAT/tests/unit/local_time_glibc_vectors.tsv`):
//!
//! - `TZ` unset: the zone file `/etc/localtime`.
//! - `TZ` set: one leading `:` is dropped, and what is left empty is UTC. Otherwise it names a
//!   zone file: itself when it starts with `/`, else the file of that name under the directory
//!   `TZDIR` (`/usr/share/zoneinfo` when `TZDIR` is unset or empty). When that is no zone file,
//!   it is a POSIX rule (`JST-9`, `EST5EDT,M3.2.0,M11.1.0`).
//! - A zone file with leap seconds (tzdata's `right/` zones) gives the time less its leap seconds,
//!   as glibc's does: LOSAT reads the file's table of leap seconds itself (RFC 8536 §3.2), which
//!   jiff does not use.
//! - A zone that jiff cannot read is UTC, as is an instant outside jiff's range. The forms of `TZ`
//!   that glibc reads and jiff does not are in `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md`.

/// The offset in seconds of the local time from UTC at `unix_seconds`, or `None` for UTC
/// (`super::time::LocalTime::now` then gives UTC).
pub fn local_utc_offset(unix_seconds: i64) -> Option<i64> {
    #[cfg(windows)]
    {
        super::windows_zone::windows_local_utc_offset(unix_seconds)
    }
    #[cfg(not(windows))]
    {
        zone_of(
            std::env::var_os("TZ").as_deref(),
            std::env::var_os("TZDIR").as_deref(),
        )?
        .offset_at(unix_seconds)
    }
}

#[cfg(not(windows))]
use jiff::tz::TimeZone;

/// A zone that `TZ` names: jiff's reading of its zone file or POSIX rule, and the leap seconds of
/// the zone file.
#[cfg(not(windows))]
struct Zone {
    zone: TimeZone,
    /// (instant, correction) of each leap second (RFC 8536 §3.2, "leap-second records"): after
    /// the instant, the local time is `correction` seconds less. At the instant itself glibc gives
    /// the second 60 of the minute, which LOSAT's `LocalTime` writes as the next minute (tzdata
    /// 2026e lists no leap second after 2016).
    leap_seconds: Vec<(i64, i64)>,
}

#[cfg(not(windows))]
impl Zone {
    /// The offset of the local time from UTC at `unix_seconds`, or `None` outside jiff's range of
    /// instants.
    fn offset_at(&self, unix_seconds: i64) -> Option<i64> {
        let instant = jiff::Timestamp::from_second(unix_seconds).ok()?;
        let leap = self.leap_seconds.partition_point(|&(at, _)| at < unix_seconds);
        let correction = leap.checked_sub(1).map_or(0, |i| self.leap_seconds[i].1);
        Some(i64::from(self.zone.to_offset(instant).seconds()) - correction)
    }
}

/// The zone of `TZ` and `TZDIR` (the module's rules), or `None` for UTC.
#[cfg(not(windows))]
fn zone_of(tz: Option<&std::ffi::OsStr>, tzdir: Option<&std::ffi::OsStr>) -> Option<Zone> {
    let Some(tz) = tz else {
        return zone_file("/etc/localtime".as_ref());
    };
    let tz = tz.as_encoded_bytes();
    let name = tz.strip_prefix(b":").unwrap_or(tz);
    if name.is_empty() {
        return None;
    }
    let mut path = Vec::new();
    if !name.starts_with(b"/") {
        match tzdir.map(std::ffi::OsStr::as_encoded_bytes) {
            Some(dir) if !dir.is_empty() => path.extend_from_slice(dir),
            _ => path.extend_from_slice(b"/usr/share/zoneinfo"),
        }
        path.push(b'/');
    }
    path.extend_from_slice(name);
    // SAFETY: `path` joins whole `OsStr`s (`tz` after an ASCII `:`, `TZDIR`) and ASCII.
    let path = unsafe { std::ffi::OsStr::from_encoded_bytes_unchecked(&path) };
    zone_file(path.as_ref()).or_else(|| {
        let zone = TimeZone::posix(std::str::from_utf8(name).ok()?).ok()?;
        Some(Zone { zone, leap_seconds: Vec::new() })
    })
}

/// The zone of the zone file at `path`, or `None` where it is not one that jiff reads.
#[cfg(not(windows))]
fn zone_file(path: &std::path::Path) -> Option<Zone> {
    use std::io::Read;
    // A zone file is a few kilobytes; the bound keeps a device (`TZ=/dev/zero`) from being read
    // without end.
    const MAX_ZONE_FILE: u64 = 1 << 20;
    let mut data = Vec::new();
    std::fs::File::open(path)
        .ok()?
        .take(MAX_ZONE_FILE + 1)
        .read_to_end(&mut data)
        .ok()?;
    if data.len() as u64 > MAX_ZONE_FILE {
        return None;
    }
    let zone = TimeZone::tzif("localtime", &data).ok()?;
    Some(Zone { zone, leap_seconds: leap_seconds(&data)? })
}

/// The leap-second records of TZif data that jiff has read (RFC 8536 §3): those of the data
/// block of version 2 and later (64-bit instants) when there is one, else of the version 1 block.
#[cfg(not(windows))]
fn leap_seconds(data: &[u8]) -> Option<Vec<(i64, i64)>> {
    // The header: "TZif", the version, 15 bytes, and the counts isutcnt, isstdcnt, leapcnt,
    // timecnt, typecnt and charcnt (32-bit big-endian).
    fn counts(block: &[u8]) -> Option<[usize; 6]> {
        let mut counts = [0; 6];
        for (i, count) in counts.iter_mut().enumerate() {
            let bytes = block.get(20 + 4 * i..24 + 4 * i)?;
            *count = u32::from_be_bytes(bytes.try_into().ok()?) as usize;
        }
        Some(counts)
    }
    // A data block: transition times (`time_size` bytes each), their type indices (1 byte), the
    // local time types (6 bytes), the designations, then the leap-second records (an instant of
    // `time_size` bytes and a 32-bit correction).
    fn block_len(counts: [usize; 6], time_size: usize) -> usize {
        let [isut, isstd, leap, time, types, chars] = counts;
        44 + time * (time_size + 1) + types * 6 + chars + leap * (time_size + 4) + isstd + isut
    }
    let v1 = counts(data)?;
    let (block, time_size) = if *data.get(4)? >= b'2' {
        (data.get(block_len(v1, 4)..)?, 8)
    } else {
        (data, 4)
    };
    let counts = counts(block)?;
    let [_, _, leap, time, types, chars] = counts;
    let start = 44 + time * (time_size + 1) + types * 6 + chars;
    let records = block.get(start..start + leap * (time_size + 4))?;
    Some(
        records
            .chunks_exact(time_size + 4)
            .map(|record| {
                let (at, correction) = record.split_at(time_size);
                let at = match time_size {
                    8 => i64::from_be_bytes(at.try_into().unwrap()),
                    _ => i64::from(i32::from_be_bytes(at.try_into().unwrap())),
                };
                (at, i64::from(i32::from_be_bytes(correction.try_into().unwrap())))
            })
            .collect(),
    )
}

/// The offset of the POSIX rule `spec` at `unix_seconds` (`None` where jiff does not read it),
/// which the tests of `super::windows_zone` compare UCRT's system zones with.
#[cfg(test)]
pub(super) fn offset_from_posix_tz(spec: &str, unix_seconds: i64) -> Option<i32> {
    let zone = jiff::tz::TimeZone::posix(spec).ok()?;
    let instant = jiff::Timestamp::from_second(unix_seconds).ok()?;
    Some(zone.to_offset(instant).seconds())
}

#[cfg(all(test, not(windows)))]
mod tests {
    use super::*;
    use std::ffi::OsStr;

    /// glibc's local time of each form of `TZ` and `TZDIR` (scripts/local_time_glibc_vectors.py).
    const VECTORS: &str = include_str!("../../tests/unit/local_time_glibc_vectors.tsv");
    /// The zone files of the forms (`{dir}` in the vectors).
    const ZONES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/unit/local_time");

    #[derive(Clone, Copy, PartialEq, Debug)]
    enum Differs {
        /// LOSAT's time is UTC at every instant.
        Utc,
        /// jiff reads the rule otherwise than glibc.
        Jiff,
    }

    /// The forms that glibc reads otherwise than LOSAT (docs/product_decisions/
    /// PD-LOSAT-NCBI-DEFECTS.md, makeblastdb's dates): every other form gives glibc's time.
    const DIFFERING: &[(&str, Differs)] = &[
        // Numbers outside POSIX's ranges, which glibc brings into them.
        ("XXX+25", Differs::Utc),
        ("XXX5:70", Differs::Utc),
        ("XXX99999", Differs::Utc),
        ("XXX-168", Differs::Utc),
        ("EST5EDT,M13.2.0,M11.1.0", Differs::Utc),
        ("EST5EDT,M3.6.0,M11.1.0", Differs::Utc),
        ("EST5EDT,M3.2.7,M11.1.0", Differs::Utc),
        ("EST5EDT,J0,J300", Differs::Utc),
        ("EST5EDT,366,300", Differs::Utc),
        ("EST5EDT,M3.2.0/168,M11.1.0", Differs::Utc),
        ("XXX-24YYY,0/0,J365/25", Differs::Jiff),
        // Two signs.
        ("AAA--8", Differs::Utc),
        ("JST--9", Differs::Utc),
        ("JST+-9", Differs::Utc),
        // A rule cut short.
        ("EST5EDT,M3.2.0", Differs::Utc),
        ("EST5EDT,", Differs::Utc),
        ("EST5EDT,M3,M11.1.0", Differs::Utc),
        ("EST5EDT,M3.2,M11.1.0", Differs::Utc),
        ("EST5EDT,,M11.1.0", Differs::Utc),
        ("EST5EDT,m3.2.0,m11.1.0", Differs::Utc),
        // Daylight time without a rule or a zone file of the name (glibc: the United States').
        ("EST5EDT4", Differs::Utc),
        ("XXX5YYY", Differs::Utc),
        ("XXX5YYY4", Differs::Utc),
        ("XXX-9YYY", Differs::Utc),
        ("JST-9junk", Differs::Utc),
        // Minutes or seconds of one digit.
        ("XXX1:2", Differs::Utc),
        ("XXX-5:3", Differs::Utc),
        ("XXX5YYY4:3,M3.2.0,M11.1.0", Differs::Utc),
        // Text after the rule.
        ("JST-9 ", Differs::Utc),
        ("JST-9:", Differs::Utc),
        ("JST-9,", Differs::Utc),
        ("EST5EDT,M3.2.0,M11.1.0 ", Differs::Utc),
        ("EST5EDT,M3.2.0,M11.1.0,M12.1.0", Differs::Utc),
        // A zone file whose footer disagrees with its last transition (zic -b slim, tzdata 2026e).
        ("slim/America/Ojinaga", Differs::Utc),
    ];

    fn losat(tz: Option<&str>, tzdir: Option<&str>, unix_seconds: i64) -> Option<i64> {
        zone_of(tz.map(OsStr::new), tzdir.map(OsStr::new))?.offset_at(unix_seconds)
    }

    /// (instant, glibc's offset) of a row of the vectors or of the server.
    fn points(row: &str) -> Vec<(i64, i64)> {
        row.split(' ')
            .map(|point| {
                let (t, offset) = point.split_once(':').unwrap();
                (t.parse().unwrap(), offset.parse().unwrap())
            })
            .collect()
    }

    /// The forms of the vectors: (TZ, TZDIR, points), `{dir}` made ZONES.
    fn vectors() -> Vec<(String, Option<String>, Vec<(i64, i64)>)> {
        VECTORS
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                assert_eq!((fields.len(), fields[0]), (4, "form"), "{line}");
                let tzdir = (fields[2] != "-").then(|| fields[2].replace("{dir}", ZONES));
                (fields[1].replace("{dir}", ZONES), tzdir, points(fields[3]))
            })
            .collect()
    }

    /// The instants where LOSAT's time is not glibc's, unless `DIFFERING` lists `form` (where
    /// LOSAT's time must differ as listed).
    fn compare(
        form: &str,
        tz: Option<&str>,
        tzdir: Option<&str>,
        glibc: &[(i64, i64)],
    ) -> Option<String> {
        let losat: Vec<Option<i64>> = glibc.iter().map(|&(t, _)| losat(tz, tzdir, t)).collect();
        let differing = DIFFERING.iter().find(|(f, _)| *f == form).map(|&(_, d)| d);
        let unequal: Vec<_> = glibc
            .iter()
            .zip(&losat)
            .filter(|&(&(_, offset), losat)| losat.unwrap_or(0) != offset)
            .map(|(&(t, offset), losat)| (t, offset, *losat))
            .collect();
        match differing {
            None if unequal.is_empty() => None,
            None => Some(format!(
                "TZ={tz:?} TZDIR={tzdir:?}: {} of {} instants differ from glibc, e.g. \
                 (instant, glibc, LOSAT) {:?}",
                unequal.len(),
                glibc.len(),
                &unequal[..unequal.len().min(3)]
            )),
            Some(_) if unequal.is_empty() => Some(format!("{form:?} is glibc's time: drop it from DIFFERING")),
            Some(Differs::Utc) if losat.iter().all(Option::is_none) => None,
            Some(Differs::Jiff) if losat.iter().all(Option::is_some) => None,
            Some(differs) => Some(format!("{form:?} is not {differs:?}: LOSAT {losat:?}")),
        }
    }

    #[test]
    fn local_time_is_glibcs() {
        let vectors = vectors();
        let failures: Vec<String> = vectors
            .iter()
            .filter_map(|(tz, tzdir, glibc)| compare(tz, Some(tz.as_str()), tzdir.as_deref(), glibc))
            .collect();
        assert!(failures.is_empty(), "{}", failures.join("\n"));
        for (form, _) in DIFFERING {
            assert!(
                vectors.iter().any(|(tz, _, _)| tz == form),
                "{form:?} of DIFFERING is not a form of the vectors"
            );
        }
    }

    /// The host's glibc against LOSAT, on Linux: the forms of the vectors, the forms of the host
    /// (`TZ` unset, names under /usr/share/zoneinfo) and each zone file there (Python's
    /// `time.localtime` through scripts/local_time_glibc_vectors.py --serve).
    #[cfg(target_os = "linux")]
    #[test]
    fn local_time_is_the_hosts_glibcs() {
        use std::io::Write;
        use std::process::{Command, Stdio};
        let script = concat!(env!("CARGO_MANIFEST_DIR"), "/../scripts/local_time_glibc_vectors.py");
        let libc = Command::new("python3")
            .args(["-I", "-c", "import platform; print(platform.libc_ver()[0])"])
            .output()
            .expect("python3 runs");
        if String::from_utf8_lossy(&libc.stdout).trim() != "glibc" {
            eprintln!("skipped: the C library of this Linux is not glibc");
            return;
        }
        let mut jobs: Vec<(Option<String>, Option<String>)> = vectors()
            .into_iter()
            .map(|(tz, tzdir, _)| (Some(tz), tzdir))
            .collect();
        for tz in [
            None,
            Some("Asia/Tokyo"),
            Some(":Asia/Tokyo"),
            Some("asia/tokyo"),
            Some("/usr/share/zoneinfo/Asia/Tokyo"),
            Some("../zoneinfo/Asia/Tokyo"),
            Some("localtime"),
            Some("posixrules"),
            Some("/etc/localtime"),
            Some("/dev/null"),
            Some("/dev/zero"),
        ] {
            jobs.push((tz.map(str::to_owned), None));
        }
        jobs.push((Some("Asia/Tokyo".to_owned()), Some(String::new())));
        let root = std::path::Path::new("/usr/share/zoneinfo");
        let mut dirs = vec![root.to_path_buf()];
        while let Some(dir) = dirs.pop() {
            for entry in std::fs::read_dir(&dir).into_iter().flatten().flatten() {
                let path = entry.path();
                // Not through a link to a directory (Debian's `posix` -> `.`).
                if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
                    dirs.push(path);
                } else if std::fs::read(&path).is_ok_and(|data| data.starts_with(b"TZif")) {
                    let name = path.strip_prefix(root).unwrap().to_str().unwrap().to_owned();
                    jobs.push((Some(name), None));
                }
            }
        }
        let input: String = jobs
            .iter()
            .map(|(tz, tzdir)| {
                let field = |v: &Option<String>| v.clone().unwrap_or_else(|| "\0".to_owned());
                format!("{}\t{}\t*\t*\n", field(tz), field(tzdir))
            })
            .collect();
        let mut server = Command::new("python3")
            .args(["-I", script, "--serve"])
            .env_remove("TZ")
            .env_remove("TZDIR")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("python3 runs");
        let mut stdin = server.stdin.take().unwrap();
        let writer = std::thread::spawn(move || stdin.write_all(input.as_bytes()));
        let output = server.wait_with_output().unwrap();
        writer.join().unwrap().unwrap();
        assert!(output.status.success());
        let rows = String::from_utf8(output.stdout).unwrap();
        let rows: Vec<&str> = rows.lines().collect();
        assert_eq!(rows.len(), jobs.len());
        let failures: Vec<String> = jobs
            .iter()
            .zip(rows)
            .filter_map(|((tz, tzdir), row)| {
                let form = tz.as_deref().unwrap_or("<unset>");
                compare(form, tz.as_deref(), tzdir.as_deref(), &points(row))
            })
            .collect();
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}

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
//! - A zone that jiff cannot read is UTC, as is an instant outside jiff's range. Where glibc's
//!   time is another (forms of `TZ` that glibc reads and jiff does not, rules whose change falls
//!   in another year in UTC, zone files that break RFC 8536) is in
//!   `docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md`.

/// How far the local time of an instant is from UTC (`local_utc_offset`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LocalOffset {
    /// The local calendar time is that of the instant plus `seconds`.
    pub seconds: i64,
    /// The instant is a leap second of the zone file (tzdata's `right/` zones), which the C
    /// library writes as second 60: the calendar time of `seconds` is second 59 of that minute.
    pub leap_second: bool,
}

/// The local time of `unix_seconds` as an offset from UTC, or `None` for UTC
/// (`super::time::LocalTime::now` then gives UTC).
pub fn local_utc_offset(unix_seconds: i64) -> Option<LocalOffset> {
    #[cfg(windows)]
    {
        super::windows_zone::windows_local_utc_offset(unix_seconds).map(|seconds| LocalOffset {
            seconds,
            leap_second: false,
        })
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
    /// (instant, correction) of each leap second (RFC 8536 §3.2, "leap-second records"): from
    /// the instant on, the local time is `correction` seconds less than the zone's offset gives.
    leap_seconds: Vec<(i64, i64)>,
}

#[cfg(not(windows))]
impl Zone {
    /// The local time at `unix_seconds`, or `None` outside jiff's range of instants. As glibc
    /// gives it (observed): the correction of the last record at or before the instant, and at
    /// the instant of a record that adds a second, second 60 (a record that takes one away skips
    /// a second instead).
    fn offset_at(&self, unix_seconds: i64) -> Option<LocalOffset> {
        let instant = jiff::Timestamp::from_second(unix_seconds).ok()?;
        let in_force = self
            .leap_seconds
            .partition_point(|&(at, _)| at <= unix_seconds);
        let correction =
            |records: usize| records.checked_sub(1).map_or(0, |i| self.leap_seconds[i].1);
        let leap_second = in_force > 0
            && self.leap_seconds[in_force - 1].0 == unix_seconds
            && correction(in_force) > correction(in_force - 1);
        Some(LocalOffset {
            seconds: i64::from(self.zone.to_offset(instant).seconds()) - correction(in_force),
            leap_second,
        })
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
        Some(Zone {
            zone,
            leap_seconds: Vec::new(),
        })
    })
}

/// The zone of the zone file at `path`, or `None` where it is not one that jiff reads.
#[cfg(not(windows))]
fn zone_file(path: &std::path::Path) -> Option<Zone> {
    let data = read_tzif(std::fs::File::open(path).ok()?)?;
    let zone = TimeZone::tzif("localtime", &data).ok()?;
    Some(Zone {
        zone,
        leap_seconds: leap_seconds(&data)?,
    })
}

/// The length of a TZif header: "TZif", the version, 15 bytes, and the counts isutcnt, isstdcnt,
/// leapcnt, timecnt, typecnt and charcnt (32-bit big-endian; RFC 8536 §3.1).
#[cfg(not(windows))]
const HEADER: usize = 44;

/// The counts of the TZif header at the start of `header`.
#[cfg(not(windows))]
fn counts(header: &[u8]) -> Option<[usize; 6]> {
    let mut counts = [0; 6];
    for (i, count) in counts.iter_mut().enumerate() {
        let bytes = header.get(20 + 4 * i..24 + 4 * i)?;
        *count = usize::try_from(u32::from_be_bytes(bytes.try_into().ok()?)).ok()?;
    }
    Some(counts)
}

/// The length of a header and its data block (RFC 8536 §3.2): transition times (`time_size`
/// bytes each), their type indices (1 byte), the local time types (6 bytes), the designations,
/// the leap-second records (an instant of `time_size` bytes and a 32-bit correction) and the
/// standard/wall and UT/local indicators (1 byte each).
#[cfg(not(windows))]
fn block_len(counts: [usize; 6], time_size: usize) -> Option<usize> {
    let [isut, isstd, leap, time, types, chars] = counts;
    [
        time.checked_mul(time_size + 1)?,
        types.checked_mul(6)?,
        chars,
        leap.checked_mul(time_size + 4)?,
        isstd,
        isut,
    ]
    .into_iter()
    .try_fold(HEADER, usize::checked_add)
}

/// The TZif data at the start of `file`, as far as its headers say (RFC 8536 §3): the first header
/// and block, and for every version but 0 (as glibc and jiff read them) the second header and
/// block and the footer. A file that does not start with "TZif" is read no further
/// (`TZ=/dev/zero`), a block only as far as the file goes, and the footer, a `TZ` string between
/// two newlines, up to 64 KiB (jiff reads past the footer no further).
#[cfg(not(windows))]
fn read_tzif(mut file: std::fs::File) -> Option<Vec<u8>> {
    use std::io::Read;
    const FOOTER: u64 = 1 << 16;
    let mut read = |data: &mut Vec<u8>, len: usize| -> Option<()> {
        let got = (&mut file)
            .take(u64::try_from(len).ok()?)
            .read_to_end(data)
            .ok()?;
        (got == len).then_some(())
    };
    let mut data = Vec::new();
    read(&mut data, HEADER)?;
    if !data.starts_with(b"TZif") {
        return None;
    }
    let len = block_len(counts(&data)?, 4)?;
    read(&mut data, len - HEADER)?;
    if data[4] == 0 {
        return Some(data);
    }
    let second = data.len();
    read(&mut data, HEADER)?;
    let len = block_len(counts(&data[second..])?, 8)?;
    read(&mut data, len - HEADER)?;
    file.take(FOOTER).read_to_end(&mut data).ok()?;
    Some(data)
}

/// The leap-second records of TZif data that jiff has read (RFC 8536 §3): those of the second
/// data block (64-bit instants) for every version but 0, else of the first block.
#[cfg(not(windows))]
fn leap_seconds(data: &[u8]) -> Option<Vec<(i64, i64)>> {
    let (block, time_size) = if *data.get(4)? != 0 {
        (data.get(block_len(counts(data)?, 4)?..)?, 8)
    } else {
        (data, 4)
    };
    let [_, _, leap, time, types, chars] = counts(block)?;
    let start = block_len([0, 0, 0, time, types, chars], time_size)?;
    let records = block.get(start..start.checked_add(leap.checked_mul(time_size + 4)?)?)?;
    Some(
        records
            .chunks_exact(time_size + 4)
            .map(|record| {
                let (at, correction) = record.split_at(time_size);
                let at = match time_size {
                    8 => i64::from_be_bytes(at.try_into().unwrap()),
                    _ => i64::from(i32::from_be_bytes(at.try_into().unwrap())),
                };
                (
                    at,
                    i64::from(i32::from_be_bytes(correction.try_into().unwrap())),
                )
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
        /// jiff reads the rule otherwise than glibc: at `instants` instants of the vectors LOSAT's
        /// offset is `offset` and not glibc's.
        Jiff { offset: i64, instants: usize },
    }

    /// The forms that glibc reads otherwise than LOSAT (docs/product_decisions/
    /// PD-LOSAT-DATABASE-SEARCH.md, the creation date): every other form gives glibc's time.
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
        // A change that falls in another year in UTC than in local time, which jiff keeps in its
        // year in UTC (glibc puts it where the rule says).
        (
            "XXX-24YYY,0/0,J365/25",
            Differs::Jiff {
                offset: 90000,
                instants: 13,
            },
        ),
        (
            "AAA-9BBB,0,365",
            Differs::Jiff {
                offset: 36000,
                instants: 2,
            },
        ),
        (
            "AAA-9BBB,0/0,J365/25",
            Differs::Jiff {
                offset: 36000,
                instants: 2,
            },
        ),
        (
            "AAA-9BBB,J60,365",
            Differs::Jiff {
                offset: 36000,
                instants: 14,
            },
        ),
        (
            "AAA1BBB,J60,365/-3",
            Differs::Jiff {
                offset: 0,
                instants: 14,
            },
        ),
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
        // Daylight time without a rule (glibc: the rules of the zone file `posixrules` under
        // `TZDIR`, or without it those of the United States since 2007).
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

    /// A local time as the vectors write glibc's: the calendar time less the instant (second 60
    /// counted as the next minute's second 0), and whether it is second 60.
    type Local = (i64, bool);
    /// (instant, local time) of a form.
    type Points = Vec<(i64, Local)>;

    /// LOSAT's local time at `unix_seconds`.
    fn losat(tz: &str, tzdir: Option<&str>, unix_seconds: i64) -> Option<Local> {
        let offset =
            zone_of(Some(OsStr::new(tz)), tzdir.map(OsStr::new))?.offset_at(unix_seconds)?;
        Some((
            offset.seconds + i64::from(offset.leap_second),
            offset.leap_second,
        ))
    }

    /// (instant, (glibc's offset, second 60)) of a row of the vectors: `instant:offset`, or
    /// `instant:offset:60` where glibc writes the second as 60.
    fn points(row: &str) -> Points {
        row.split(' ')
            .map(|point| {
                let mut fields = point.split(':');
                let t = fields.next().unwrap().parse().unwrap();
                let offset = fields.next().unwrap().parse().unwrap();
                let leap_second = match fields.next() {
                    None => false,
                    Some("60") => true,
                    Some(other) => panic!("{point}: {other}"),
                };
                (t, (offset, leap_second))
            })
            .collect()
    }

    /// The forms of the vectors: (TZ, TZDIR, points), `{dir}` made ZONES.
    fn vectors() -> Vec<(String, Option<String>, Points)> {
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

    /// The instants where LOSAT's time is not glibc's, unless `DIFFERING` lists `tz` (where
    /// LOSAT's time must differ as listed).
    fn compare(tz: &str, tzdir: Option<&str>, glibc: &[(i64, Local)]) -> Option<String> {
        let losat: Vec<Option<Local>> = glibc.iter().map(|&(t, _)| losat(tz, tzdir, t)).collect();
        let differing = DIFFERING.iter().find(|(f, _)| *f == tz).map(|&(_, d)| d);
        let unequal: Vec<_> = glibc
            .iter()
            .zip(&losat)
            .filter(|&(&(_, glibc), losat)| losat.unwrap_or_default() != glibc)
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
            Some(_) if unequal.is_empty() => {
                Some(format!("{tz:?} is glibc's time: drop it from DIFFERING"))
            }
            Some(Differs::Utc) if losat.iter().all(Option::is_none) => None,
            Some(Differs::Jiff { offset, instants })
                if unequal.len() == instants
                    && unequal.iter().all(|&(_, _, l)| l == Some((offset, false))) =>
            {
                None
            }
            Some(differs) => Some(format!("{tz:?} is not {differs:?}: LOSAT {losat:?}")),
        }
    }

    #[test]
    fn local_time_is_glibcs() {
        let vectors = vectors();
        let failures: Vec<String> = vectors
            .iter()
            .filter_map(|(tz, tzdir, glibc)| compare(tz, tzdir.as_deref(), glibc))
            .collect();
        assert!(failures.is_empty(), "{}", failures.join("\n"));
        for (form, _) in DIFFERING {
            assert!(
                vectors.iter().any(|(tz, _, _)| tz == form),
                "{form:?} of DIFFERING is not a form of the vectors"
            );
        }
    }

    /// The files that glibc reads by default, which the vectors cannot name: `/etc/localtime`
    /// with `TZ` unset, and a zone name under `/usr/share/zoneinfo` with `TZDIR` unset or empty
    /// (where the host has these files; else both sides are UTC).
    #[test]
    fn the_default_files_are_glibcs() {
        let offsets = |zone: Option<Zone>| -> Vec<Option<LocalOffset>> {
            [0, 1_767_225_600, 1_785_000_000]
                .iter()
                .map(|&t| zone.as_ref().and_then(|zone| zone.offset_at(t)))
                .collect()
        };
        assert_eq!(
            offsets(zone_of(None, Some(OsStr::new("/nonexistent")))),
            offsets(zone_file("/etc/localtime".as_ref()))
        );
        let tokyo = offsets(zone_file("/usr/share/zoneinfo/Asia/Tokyo".as_ref()));
        for tzdir in [None, Some(""), Some("/usr/share/zoneinfo")] {
            let zone = zone_of(Some(OsStr::new("Asia/Tokyo")), tzdir.map(OsStr::new));
            assert_eq!(offsets(zone), tokyo, "TZDIR={tzdir:?}");
        }
    }
}

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
#[path = "local_offset_tests.rs"]
mod tests;

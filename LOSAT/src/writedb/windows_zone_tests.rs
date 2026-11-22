use super::*;
use crate::writedb::local_offset::offset_from_posix_tz;

/// UCRT's `_tzset` and `_localtime64_s` of `ucrtbase.dll` with `TZ` set
/// (scripts/ucrt_tz_vectors.py; LOSAT/tests/windows_local_time.rs checks them again against
/// the UCRT of the Windows runner).
const UCRT_VECTORS: &str = include_str!("../../tests/unit/ucrt_local_time_vectors.tsv");

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

/// The characters above U+00FF of the vector file that `GetStringTypeW` calls white space.
fn vector_spaces() -> Vec<(u16, bool)> {
    UCRT_VECTORS
        .lines()
        .filter_map(|l| l.strip_prefix("space\t"))
        .map(|l| {
            let (unit, space) = l.split_once('\t').unwrap();
            (u16::from_str_radix(unit, 16).unwrap(), space == "1")
        })
        .collect()
}

fn utf16(s: &str) -> Vec<u16> {
    s.encode_utf16().collect()
}

/// `iswspace` above U+00FF for the tests: the characters the vector file names.
fn test_space(c: u16) -> bool {
    vector_spaces().contains(&(c, true))
}

fn tz_zone(tz: &str) -> UcrtZone {
    ucrt_zone_from_tz(&utf16(tz), &test_space)
}

#[test]
fn ucrt_tz_values_match_the_c_runtime() {
    let (mut zones, mut times, mut bad) = (0, 0, Vec::new());
    for line in UCRT_VECTORS.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "space" => {}
            "zone" => {
                zones += 1;
                let z = ucrt_zone_from_tz(&vector_tz(f[1]), &test_space);
                let want = (
                    f[2].parse::<i64>().unwrap(),
                    f[3] == "1",
                    f[4].parse::<i64>().unwrap(),
                );
                if (z.timezone, z.daylight, z.dstbias) != want || z.changes.is_some() {
                    bad.push(format!("{line} -> {z:?}"));
                }
            }
            "time" => {
                times += 1;
                let z = ucrt_zone_from_tz(&vector_tz(f[1]), &test_space);
                let t: i64 = f[2].parse().unwrap();
                let want = (f[3] != "none").then(|| f[3].parse::<i64>().unwrap());
                if z.offset_at(t) != want {
                    bad.push(format!("{line} -> {:?}", z.offset_at(t)));
                }
            }
            kind => panic!("unknown row {kind}"),
        }
    }
    assert!(
        zones >= 70 && times >= 1_800,
        "{zones} zones, {times} times"
    );
    assert!(
        bad.is_empty(),
        "{} of {times} rows differ:\n{}",
        bad.len(),
        bad[..bad.len().min(30)].join("\n")
    );
}

#[cfg(windows)]
#[test]
fn windows_white_space_is_the_c_runtimes() {
    for (unit, space) in vector_spaces() {
        assert_eq!(windows_char_is_space(unit), space, "U+{unit:04X}");
    }
}

// The forms of UCRT's documentation (`tzn[+|-]hh[:mm[:ss]][dzn]`, west of UTC positive) and
// the readings its code gives the forms of other C libraries.
#[test]
fn ucrt_tz_forms() {
    let zone = |tz: &str| {
        let z = tz_zone(tz);
        (z.timezone, z.daylight)
    };
    assert_eq!(zone("JST-9"), (-32_400, false));
    assert_eq!(zone("PST8PDT"), (28_800, true));
    assert_eq!(zone("IST-5:30"), (-19_800, false));
    assert_eq!(zone("AAA12:34:56BBB"), (45_296, true));
    assert_eq!(zone("AAA+8"), (28_800, false));
    // The name is three characters, whatever they are.
    assert_eq!(zone("AKST9AKDT"), (0, true));
    assert_eq!(zone("123-9"), (-32_400, false));
    assert_eq!(zone("AB"), (0, false));
    assert_eq!(zone(""), (0, false));
    // `wcstol` takes white space and a second sign; the scan after it stops there, and what
    // is left makes a daylight-time zone.
    assert_eq!(zone("AAA -8"), (-28_800, true));
    assert_eq!(zone("AAA--8"), (28_800, true));
    assert_eq!(zone("AAA++8"), (0, false));
    assert_eq!(zone("AAA\u{3000}8"), (28_800, true));
    assert_eq!(zone("AAA\u{200B}8"), (0, true));
    // Digits of other scripts count for `wcstol`, not for the scan.
    assert_eq!(zone("AAA\u{FF18}"), (28_800, true));
    assert_eq!(zone("AAA8:\u{FF13}\u{FF10}"), (30_600, true));
    // POSIX rules and Olson names are not read: any rest means the United States rules.
    assert_eq!(zone("AAA8,M3.2.0,M11.1.0"), (28_800, true));
    assert_eq!(zone(":Asia/Tokyo"), (0, true));
    assert_eq!(zone("EST5EDT4"), (18_000, true));
    // A NUL ends the value.
    assert_eq!(
        ucrt_zone_from_tz(&[u16::from(b'J'), 0, u16::from(b'X')], &test_space).timezone,
        0
    );
}

#[test]
fn ucrt_wcstol_reads_as_the_c_runtime() {
    let n = |s: &str| ucrt_wcstol(&utf16(s), &test_space);
    assert_eq!(n("42x"), 42);
    assert_eq!(n(" \t\u{A0}\u{85}-7"), -7);
    assert_eq!(n("+"), 0);
    assert_eq!(n("-"), 0);
    assert_eq!(n("x1"), 0);
    assert_eq!(n("\u{0669}\u{06F1}\u{0967}\u{1819}"), 9119);
    // The ranges end at the tenth digit.
    assert_eq!(n("\u{066A}"), 0);
    assert_eq!(n("2147483647"), i32::MAX);
    assert_eq!(n("2147483648"), i32::MAX);
    assert_eq!(n("-2147483648"), i32::MIN);
    assert_eq!(n("-2147483649"), i32::MIN);
    assert_eq!(n("99999999999999999999999"), i32::MAX);
    assert_eq!(n("-99999999999999999999999"), i32::MIN);
    // U+FFFF is `WEOF`: never white space, even where the system would say so.
    assert!(!ucrt_iswspace(0xFFFF, &|_| true));
    assert!(ucrt_iswspace(0x3000, &|_| true));
    assert!(!ucrt_iswspace(0x00A1, &|_| true));
}

/// The United States rules of `_isindst` with `TZ` set: 2007 and later from the second
/// Sunday of March to the first Sunday of November, before 2007 from the first Sunday of
/// April to the last Sunday of October, at 02:00 local time.
#[test]
fn ucrt_tz_dst_follows_the_us_rules() {
    let pst = tz_zone("PST8PDT");
    for (change, before, after) in [
        (1_772_964_000, -28_800, -25_200), // 2026-03-08 10:00 UTC
        (1_793_523_600, -25_200, -28_800), // 2026-11-01 09:00 UTC
        (1_143_972_000, -28_800, -25_200), // 2006-04-02 10:00 UTC
        (1_162_112_400, -25_200, -28_800), // 2006-10-29 09:00 UTC
    ] {
        assert_eq!(pst.offset_at(change - 1), Some(before), "{change}");
        assert_eq!(pst.offset_at(change), Some(after), "{change}");
    }
    // Without a daylight name there is no change.
    assert_eq!(tz_zone("MST7").offset_at(1_782_864_000), Some(-25_200));
}

/// `TIME_ZONE_INFORMATION` of the registry of Windows 11 (`TZI`, 2026-11-18).
fn win_zone(
    bias: i32,
    standard: (u16, u16, u16, u16, u16),
    daylight: (u16, u16, u16, u16, u16),
    daylight_bias: i32,
) -> WinZone {
    let change = |(month, day, day_of_week, hour, minute): (u16, u16, u16, u16, u16)| WinChange {
        month,
        day,
        day_of_week,
        hour,
        minute,
        ..WinChange::default()
    };
    WinZone {
        bias,
        standard_date: change(standard),
        standard_bias: 0,
        daylight_date: change(daylight),
        daylight_bias,
    }
}

/// Windows zones and the POSIX rule strings of the same changes, as jiff reads them
/// (`offset_from_posix_tz`; jiff reads these strings as glibc does,
/// LOSAT/tests/unit/local_time_glibc_vectors.tsv).
fn system_zones() -> Vec<(&'static str, WinZone, &'static str)> {
    vec![
        (
            "Pacific",
            win_zone(480, (11, 1, 0, 2, 0), (3, 2, 0, 2, 0), -60),
            "PST8PDT,M3.2.0,M11.1.0",
        ),
        (
            "AUS Eastern",
            win_zone(-600, (4, 1, 0, 3, 0), (10, 1, 0, 2, 0), -60),
            "AEST-10AEDT,M10.1.0,M4.1.0/3",
        ),
        (
            "GMT",
            win_zone(0, (10, 5, 0, 2, 0), (3, 5, 0, 1, 0), -60),
            "GMT0BST,M3.5.0/1,M10.5.0",
        ),
        (
            "Newfoundland",
            win_zone(210, (11, 1, 0, 2, 0), (3, 2, 0, 2, 0), -60),
            "NST3:30NDT,M3.2.0,M11.1.0",
        ),
        (
            "Lord Howe",
            win_zone(-630, (4, 1, 0, 2, 0), (10, 1, 0, 2, 0), -30),
            "LHST-10:30LHDT-11,M10.1.0,M4.1.0",
        ),
        (
            "Chatham Islands",
            win_zone(-765, (4, 1, 0, 3, 45), (9, 5, 0, 2, 45), -60),
            "CHAST-12:45CHADT,M9.5.0/2:45,M4.1.0/3:45",
        ),
        (
            "Morocco",
            win_zone(0, (9, 3, 0, 2, 0), (3, 4, 0, 2, 0), -60),
            "XXX0YYY,M3.4.0,M9.3.0",
        ),
        (
            "Tokyo",
            win_zone(-540, (0, 0, 0, 0, 0), (0, 0, 0, 0, 0), -60),
            "JST-9",
        ),
        (
            "UTC",
            win_zone(0, (0, 0, 0, 0, 0), (0, 0, 0, 0, 0), 0),
            "UTC0",
        ),
    ]
}

// UCRT places the changes of the year of the local standard time, glibc (and jiff) those of
// the UTC year: they meet for zones whose changes are not at the turn of the year.
#[test]
fn ucrt_system_zones_change_as_their_rules_say() {
    for (name, info, posix) in system_zones() {
        let zone = ucrt_zone_from_system(Some(&info));
        let mut t = 1_735_689_600; // 2025-01-01
        while t < 1_830_297_600 {
            // 2028-01-01
            for u in [t - 1, t, t + 1] {
                assert_eq!(
                    zone.offset_at(u),
                    offset_from_posix_tz(posix, u).map(i64::from),
                    "{name} at {u}"
                );
            }
            t += 900;
        }
    }
}

#[test]
fn ucrt_system_zone_values() {
    let zones = system_zones();
    let zone = |name: &str| {
        let info = zones.iter().find(|z| z.0 == name).unwrap().1;
        ucrt_zone_from_system(Some(&info))
    };
    let p = zone("Pacific");
    assert_eq!((p.timezone, p.daylight, p.dstbias), (28_800, true, -3_600));
    assert_eq!(p.offset_at(1_772_964_000 - 1), Some(-28_800));
    assert_eq!(p.offset_at(1_772_964_000), Some(-25_200));
    let a = zone("AUS Eastern");
    assert_eq!(a.offset_at(1_775_318_400 - 1), Some(39_600));
    assert_eq!(a.offset_at(1_775_318_400), Some(36_000));
    assert_eq!(a.offset_at(1_791_043_200 - 1), Some(36_000));
    assert_eq!(a.offset_at(1_791_043_200), Some(39_600));
    // A daylight bias without a daylight date is no daylight time.
    let j = zone("Tokyo");
    assert_eq!((j.timezone, j.daylight, j.dstbias), (-32_400, false, 0));
    let l = zone("Lord Howe");
    assert_eq!((l.timezone, l.daylight, l.dstbias), (-37_800, true, -1_800));
    // A failed GetTimeZoneInformation leaves the start values: Pacific time, US rules.
    let f = ucrt_zone_from_system(None);
    assert_eq!(f, UCRT_START);
    assert_eq!(f.offset_at(1_772_964_000), Some(-25_200));
    // The standard bias counts only with a standard date, and the daylight bias is then
    // relative to it.
    let mut info = win_zone(300, (11, 1, 0, 2, 0), (3, 2, 0, 2, 0), -60);
    info.standard_bias = 15;
    let z = ucrt_zone_from_system(Some(&info));
    assert_eq!((z.timezone, z.dstbias), (18_900, -4_500));
    info.standard_date.month = 0;
    assert_eq!(ucrt_zone_from_system(Some(&info)).timezone, 18_000);
}

#[test]
fn ucrt_cvtdate_forms() {
    let at = |year, month, day_of_week, day, hour| WinChange {
        year,
        month,
        day_of_week,
        day,
        hour,
        ..WinChange::default()
    };
    // 2026 (126): the second Sunday of March is March 8 (day 66), the last Sunday of
    // October is October 25 (day 297), the first Sunday of November is November 1 (304).
    assert_eq!(ucrt_cvtdate(&at(0, 3, 0, 2, 2), 126), Some((66, 7_200_000)));
    assert_eq!(
        ucrt_cvtdate(&at(0, 10, 0, 5, 2), 126),
        Some((297, 7_200_000))
    );
    assert_eq!(
        ucrt_cvtdate(&at(0, 11, 0, 1, 2), 126),
        Some((304, 7_200_000))
    );
    // 2024 (124, leap): the last Thursday of February is February 29 (day 59).
    assert_eq!(ucrt_cvtdate(&at(0, 2, 4, 5, 0), 124), Some((59, 0)));
    // An absolute date (`wYear` not 0) is that day of the month in every year.
    assert_eq!(ucrt_cvtdate(&at(2026, 3, 0, 21, 0), 124), Some((80, 0)));
    assert_eq!(ucrt_cvtdate(&at(2026, 3, 0, 21, 0), 126), Some((79, 0)));
    // Month 13 of an absolute date is the last entry of the day table; months beyond the
    // tables are undefined reads.
    assert_eq!(ucrt_cvtdate(&at(1, 13, 0, 1, 0), 126), Some((365, 0)));
    assert_eq!(ucrt_cvtdate(&at(1, 14, 0, 1, 0), 126), None);
    assert_eq!(ucrt_cvtdate(&at(1, 0, 0, 1, 0), 126), None);
    assert_eq!(ucrt_cvtdate(&at(0, 13, 0, 1, 0), 126), None);
    assert_eq!(ucrt_cvtdate(&at(0, 0, 0, 1, 0), 126), None);
}

// `int` and `long` overflow in UCRT's reading of `TZ` beyond 596,523 hours: LOSAT takes
// the mathematical value.
#[test]
fn ucrt_offsets_beyond_32_bits_are_mathematical() {
    let z = tz_zone("AAA-596524");
    assert_eq!(z.timezone, -2_147_486_400);
    assert_eq!(z.offset_at(1_700_000_000), Some(2_147_486_400));
    // Hours beyond `long` are `LONG_MAX` hours: no local time.
    let z = tz_zone("AAA99999999999");
    assert_eq!(z.timezone, i64::from(i32::MAX) * 3600);
    assert_eq!(z.offset_at(1_700_000_000), None);
    assert_eq!(tz_zone("AAA-0:35791395").timezone, -2_147_483_700);
}

#[test]
fn ucrt_localtime_range() {
    let utc = tz_zone("UTC");
    // The first and last three days of the range take a branch that is not ported.
    assert_eq!(utc.offset_at(3 * UCRT_DAY_SEC), None);
    assert_eq!(utc.offset_at(3 * UCRT_DAY_SEC + 1), Some(0));
    assert_eq!(
        utc.offset_at(UCRT_MAX_TIME64 - 3 * UCRT_DAY_SEC - 1),
        Some(0)
    );
    assert_eq!(utc.offset_at(UCRT_MAX_TIME64 - 3 * UCRT_DAY_SEC), None);
    assert_eq!(utc.offset_at(-1), None);
    // `gmtime` of the local time must be within its range.
    let west = tz_zone("AAA13");
    assert_eq!(west.offset_at(3 * UCRT_DAY_SEC + 1), Some(-46_800));
    let east = tz_zone("AAA-596523");
    assert_eq!(east.offset_at(1_700_000_000), Some(2_147_482_800));
    let far_west = tz_zone("AAA596523");
    assert_eq!(far_west.offset_at(1_700_000_000), None);
    assert_eq!(far_west.offset_at(2_147_439_600), Some(-2_147_482_800));
}

// Within its range `_gmtime64_s` is the Gregorian calendar of `LocalTime`.
#[test]
fn ucrt_gmtime_is_the_gregorian_calendar() {
    let mut t = UCRT_MIN_LOCAL_TIME;
    let mut step = 1i64;
    while t <= UCRT_MAX_TIME64 + UCRT_MAX_LOCAL_TIME {
        let tm = ucrt_gmtime(t).unwrap();
        let lt = crate::writedb::time::LocalTime::from_unix_seconds(t);
        let days_before_month = if ucrt_is_leap_year(tm.year) {
            &UCRT_LPDAYS
        } else {
            &UCRT_DAYS
        };
        let yday = i64::from(days_before_month[lt.month as usize - 1]) + i64::from(lt.day);
        assert_eq!(
            (tm.year + 1900, tm.yday, tm.hour, tm.min, tm.sec),
            (
                i64::from(lt.year),
                yday,
                i64::from(lt.hour),
                i64::from(lt.minute),
                i64::from(lt.second)
            ),
            "{t}"
        );
        t += step;
        step = (step * 3 + 7) % 9_999_991;
    }
    assert_eq!(ucrt_gmtime(UCRT_MIN_LOCAL_TIME - 1), None);
    assert_eq!(ucrt_gmtime(UCRT_MAX_TIME64 + UCRT_MAX_LOCAL_TIME + 1), None);
}

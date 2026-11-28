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
    let offset = zone_of(Some(OsStr::new(tz)), tzdir.map(OsStr::new))?.offset_at(unix_seconds)?;
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

/// The rules that LOSAT Web's host puts into `TZ` for the browser's zone in 2026
/// (web/packages/losat/src/infra/reactor/time-zone.ts and its test; line NP WW4b1, D-WW4b-2,
/// 2026-11-27): read as `TZ` (no zone file has the rule's name, so jiff reads the POSIX rule), each
/// gives its zone's offset at every hour of 2026 (where the host has the zone files).
#[test]
fn the_rules_of_the_web_host_are_their_zones() {
    const RULES: &[(&str, &str)] = &[
        ("Asia/Tokyo", "<+09>-9"),
        ("Asia/Kolkata", "<+0530>-5:30"),
        ("UTC", "<+00>0"),
        ("America/Phoenix", "<-07>7"),
        ("America/New_York", "<-05>5<-04>4,M3.2.0/2,M11.1.0/2"),
        ("Europe/London", "<+00>0<+01>-1,M3.5.0/1,M10.5.0/2"),
        ("Europe/Berlin", "<+01>-1<+02>-2,M3.5.0/2,M10.5.0/3"),
        ("Australia/Sydney", "<+10>-10<+11>-11,M10.1.0/2,M4.1.0/3"),
        (
            "Australia/Lord_Howe",
            "<+1030>-10:30<+11>-11,M10.1.0/2,M4.1.0/2",
        ),
    ];
    const START_2026: i64 = 1_767_225_600;
    for &(name, rule) in RULES {
        let Some(file) = zone_file(format!("/usr/share/zoneinfo/{name}").as_ref()) else {
            eprintln!("no zone file for {name}: its rule is not compared");
            continue;
        };
        let zone = zone_of(Some(OsStr::new(rule)), None).expect(rule);
        for hour in 0..366 * 24 {
            let at = START_2026 + hour * 3600;
            assert_eq!(zone.offset_at(at), file.offset_at(at), "{rule} at {at}");
        }
    }
}

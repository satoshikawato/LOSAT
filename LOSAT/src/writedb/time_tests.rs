use super::*;

fn at(month: u32, day: u32, hour: u32, minute: u32) -> LocalTime {
    LocalTime {
        year: 2026,
        month,
        day,
        hour,
        minute,
        second: 7,
    }
}

// The forms makeblastdb 2.17.0 wrote: "Oct 12, 2026  3:22 PM" in an
// index file and "10/12/2026 20:29:48" in its log.
#[test]
fn index_date_has_a_12_hour_clock_without_leading_zero() {
    assert_eq!(at(10, 12, 15, 22).index_date(), "Oct 12, 2026  3:22 PM");
    assert_eq!(at(10, 12, 20, 29).index_date(), "Oct 12, 2026  8:29 PM");
    assert_eq!(at(3, 9, 9, 5).index_date(), "Mar 9, 2026  9:05 AM");
    assert_eq!(at(1, 1, 0, 0).index_date(), "Jan 1, 2026  12:00 AM");
    assert_eq!(at(12, 31, 12, 30).index_date(), "Dec 31, 2026  12:30 PM");
    assert_eq!(at(12, 31, 23, 59).index_date(), "Dec 31, 2026  11:59 PM");
    assert_eq!(at(6, 5, 10, 1).index_date(), "Jun 5, 2026  10:01 AM");
}

#[test]
fn log_time_is_zero_padded() {
    assert_eq!(at(10, 12, 20, 29).log_time(), "10/12/2026 20:29:07");
    assert_eq!(at(3, 9, 9, 5).log_time(), "03/09/2026 09:05:07");
}

#[test]
fn the_reader_parses_the_date_the_writer_formats() {
    let earlier = at(3, 9, 9, 5).index_date();
    let later = at(10, 12, 15, 22).index_date();
    assert_eq!(
        crate::seqdb::date::database_date(&[earlier.as_str(), later.as_str()]).unwrap(),
        later
    );
}

// The calendar dates of Python's `datetime` for the same instants.
#[test]
fn unix_seconds_become_calendar_times() {
    let cases: [(i64, (i32, u32, u32, u32, u32, u32)); 17] = [
        (0, (1970, 1, 1, 0, 0, 0)),
        (-1, (1969, 12, 31, 23, 59, 59)),
        (86399, (1970, 1, 1, 23, 59, 59)),
        (86400, (1970, 1, 2, 0, 0, 0)),
        (951782399, (2000, 2, 28, 23, 59, 59)),
        (951782400, (2000, 2, 29, 0, 0, 0)),
        (951868800, (2000, 3, 1, 0, 0, 0)),
        (1798761599, (2026, 12, 31, 23, 59, 59)),
        (1798761600, (2027, 1, 1, 0, 0, 0)),
        (4102444800, (2100, 1, 1, 0, 0, 0)),
        (4107542400, (2100, 3, 1, 0, 0, 0)),
        (-2208988800, (1900, 1, 1, 0, 0, 0)),
        (253402300799, (9999, 12, 31, 23, 59, 59)),
        (1791807000, (2026, 10, 12, 12, 10, 0)),
        (68169599, (1972, 2, 28, 23, 59, 59)),
        (68256000, (1972, 3, 1, 0, 0, 0)),
        (-62135596800, (1, 1, 1, 0, 0, 0)),
    ];
    for (seconds, (year, month, day, hour, minute, second)) in cases {
        assert_eq!(
            LocalTime::from_unix_seconds(seconds),
            LocalTime {
                year,
                month,
                day,
                hour,
                minute,
                second
            },
            "{seconds}"
        );
    }
}

#[test]
fn the_clock_gives_a_plausible_time() {
    let now = LocalTime::now().expect("a clock");
    assert!(now.year >= 2026 && (1..=12).contains(&now.month) && (1..=31).contains(&now.day));
    assert!(now.hour < 24 && now.minute < 60 && now.second < 62);
}

use super::*;

#[test]
fn the_metadata_date_is_iso_8601() {
    // probe J185-J188 of line BDC: the dates of copies of a version 4 database.
    assert_eq!(
        iso8601_date_time("Oct 29, 2026  12:13 AM").as_deref(),
        Ok("2026-10-29T00:13:00")
    );
    assert_eq!(
        iso8601_date_time("Oct 29, 2026  12:13 PM").as_deref(),
        Ok("2026-10-29T12:13:00")
    );
    assert_eq!(
        iso8601_date_time("Jan 1, 2026  1:05 AM").as_deref(),
        Ok("2026-01-01T01:05:00")
    );
    assert_eq!(
        iso8601_date_time("Oct 29, 2026  7:13 PM").as_deref(),
        Ok("2026-10-29T19:13:00")
    );
    // J189, J190: NCBI's constructor throws (exit status 255).
    assert!(iso8601_date_time("Feb 30, 2026  1:05 PM").is_err());
    assert!(iso8601_date_time("not a date").is_err());
}

// Date strings as makeblastdb 2.17.0 writes them.
const PM_322: &str = "Oct 12, 2026  3:22 PM";
const AM_101: &str = "Oct 12, 2026  1:01 AM";

fn time(year: u32, month: u32, day: u32, hour: u32, min: u32) -> Time {
    Time {
        year,
        month,
        day,
        hour,
        min,
    }
}

fn parsed(text: &str) -> Time {
    parse_blast_db_date(text).unwrap()
}

fn time_exception(text: &str) -> String {
    match parse_blast_db_date(text) {
        Err(ParseFailure::TimeException(exception)) => {
            String::from_utf8_lossy(exception.get_msg()).into_owned()
        }
        other => panic!("{text:?}: expected a CTimeException, got {other:?}"),
    }
}

fn too_short(text: &str) -> String {
    format!("Time string '{text}' is too short for time format 'b d, Y  H:m P'")
}

fn too_long(text: &str) -> String {
    format!("Time string '{text}' is too long for time format 'b d, Y  H:m P'")
}

fn invalid(text: &str) -> String {
    format!("Unable to convert string '{text}' to CTime")
}

fn date_of(volume_dates: &[&str]) -> String {
    database_date(volume_dates).unwrap()
}

// ---- database_date (CSeqDBImpl::GetDate) ----

#[test]
fn no_volume_gives_empty_string() {
    // The loop does not run; `date` is still the empty string.
    assert_eq!(date_of(&[]), "");
}

#[test]
fn single_volume_is_returned_without_parsing() {
    // `date.empty()` is true on the only iteration, so no CTime is built.
    assert_eq!(date_of(&[PM_322]), PM_322);
    assert_eq!(date_of(&["not a date"]), "not a date");
}

#[test]
fn equal_dates_are_not_compared() {
    assert_eq!(date_of(&[PM_322, PM_322, PM_322]), PM_322);
    // `d != date` is false, so even a string that CTime rejects is not parsed.
    assert_eq!(date_of(&["not a date", "not a date"]), "not a date");
}

#[test]
fn later_second_date_replaces_first() {
    // Different minute.
    assert_eq!(
        date_of(&["Oct 12, 2026  3:22 PM", "Oct 12, 2026  3:23 PM"]),
        "Oct 12, 2026  3:23 PM"
    );
    // AM against PM: 1:01 AM is hour 1, 3:22 PM is hour 15.
    assert_eq!(date_of(&[AM_101, PM_322]), PM_322);
    // PM beats a larger AM hour number: 11:59 AM is hour 11, 1:00 PM is hour 13.
    assert_eq!(
        date_of(&["Oct 12, 2026 11:59 AM", "Oct 12, 2026  1:00 PM"]),
        "Oct 12, 2026  1:00 PM"
    );
    // Different day.
    assert_eq!(
        date_of(&["Oct 12, 2026  3:22 PM", "Oct 13, 2026  1:01 AM"]),
        "Oct 13, 2026  1:01 AM"
    );
    // Different month: months compare by number (Sep = 9, Oct = 10), not by name.
    assert_eq!(
        date_of(&["Sep 30, 2026 11:59 PM", "Oct 1, 2026 12:00 AM"]),
        "Oct 1, 2026 12:00 AM"
    );
    // Different year.
    assert_eq!(
        date_of(&["Dec 31, 2025 11:59 PM", "Jan 1, 2026 12:00 AM"]),
        "Jan 1, 2026 12:00 AM"
    );
}

#[test]
fn earlier_second_date_keeps_first() {
    assert_eq!(date_of(&[PM_322, AM_101]), PM_322);
    assert_eq!(
        date_of(&["Oct 13, 2026  1:01 AM", "Oct 12, 2026  3:22 PM"]),
        "Oct 13, 2026  1:01 AM"
    );
    assert_eq!(
        date_of(&["Jan 1, 2026 12:00 AM", "Dec 31, 2025 11:59 PM"]),
        "Jan 1, 2026 12:00 AM"
    );
}

#[test]
fn same_time_written_differently_keeps_first() {
    // The strings differ, both parse to the same time, and `t2 > t1` is false.
    assert_eq!(
        date_of(&["Oct 12, 2026  3:22 PM", "oct 12, 2026 03:22 pm"]),
        "Oct 12, 2026  3:22 PM"
    );
}

#[test]
fn latest_of_several_volumes_wins() {
    assert_eq!(
        date_of(&[
            AM_101,
            "Oct 14, 2026  9:00 AM",
            PM_322,
            "Oct 14, 2026  9:00 AM"
        ]),
        "Oct 14, 2026  9:00 AM"
    );
}

#[test]
fn twelve_am_is_hour_zero() {
    // `SET_HOUR` stores the hour modulo 12 (so 12 becomes 0), and AM adds nothing.
    assert_eq!(parsed("Oct 12, 2026 12:30 AM"), time(2026, 10, 12, 0, 30));
    assert_eq!(
        date_of(&["Oct 12, 2026 12:30 AM", "Oct 12, 2026  1:00 AM"]),
        "Oct 12, 2026  1:00 AM"
    );
    assert_eq!(
        date_of(&["Oct 12, 2026  1:00 AM", "Oct 12, 2026 12:30 AM"]),
        "Oct 12, 2026  1:00 AM"
    );
}

#[test]
fn twelve_pm_is_hour_twelve() {
    // `SET_HOUR` stores the hour modulo 12 (so 12 becomes 0), and PM adds 12.
    assert_eq!(parsed("Oct 12, 2026 12:30 PM"), time(2026, 10, 12, 12, 30));
    assert_eq!(
        date_of(&["Oct 12, 2026 11:59 AM", "Oct 12, 2026 12:30 PM"]),
        "Oct 12, 2026 12:30 PM"
    );
    assert_eq!(
        date_of(&["Oct 12, 2026 12:30 PM", "Oct 12, 2026  1:00 PM"]),
        "Oct 12, 2026  1:00 PM"
    );
    assert_eq!(
        date_of(&["Oct 12, 2026 12:30 AM", "Oct 12, 2026 12:30 PM"]),
        "Oct 12, 2026 12:30 PM"
    );
}

#[test]
fn hour_zero_is_accepted_like_twelve() {
    // CHECK_RANGE2_HOUR12 allows 0..=12, and 0 % 12 == 12 % 12.
    assert_eq!(parsed("Oct 12, 2026  0:30 AM"), time(2026, 10, 12, 0, 30));
    assert_eq!(parsed("Oct 12, 2026  0:30 PM"), time(2026, 10, 12, 12, 30));
}

#[test]
fn empty_first_date_is_replaced_without_comparison() {
    // `date.empty()` is tested on every iteration.
    assert_eq!(date_of(&["", PM_322]), PM_322);
    assert_eq!(date_of(&["", "", AM_101, PM_322]), PM_322);
    // The replacement is not parsed at that point.
    assert_eq!(date_of(&["", "not a date"]), "not a date");
}

#[test]
fn empty_later_date_is_never_later() {
    // An empty string gives an empty CTime (all fields 0) without any exception, and
    // year 0 is not greater than the year of t1.
    assert_eq!(date_of(&[PM_322, ""]), PM_322);
    assert_eq!(date_of(&[PM_322, "", AM_101]), PM_322);
}

#[test]
fn number_without_digit_is_caught_and_keeps_current_date() {
    // t2: `NStr::StringToLong` on an empty string throws `CStringException`, caught by `GetDate`.
    for bad in [
        "Oct xx, 2026  3:22 PM", // day
        "Oct 12, year  3:22 PM", // year
        "Oct 12, 2026  x:22 PM", // hour
        "Oct 12, 2026  3:xx PM", // minute
        "Oct 12, 2026  3:",      // minute at the end of the string
        "Oct ",                  // day at the end of the string
    ] {
        assert!(
            matches!(
                parse_blast_db_date(bad),
                Err(ParseFailure::StringException(_))
            ),
            "{bad:?}"
        );
        assert_eq!(date_of(&[PM_322, bad]), PM_322, "{bad:?}");
        // The next volume is still compared.
        assert_eq!(
            date_of(&[PM_322, bad, "Oct 13, 2026  1:01 AM"]),
            "Oct 13, 2026  1:01 AM",
            "{bad:?}"
        );
    }
}

#[test]
fn current_date_without_digit_is_caught_and_stays() {
    // t1 throws CStringException on every later comparison, so the first date stays
    // although it is not a valid date and a valid one follows.
    let bad = "Oct xx, 2026  3:22 PM";
    assert_eq!(date_of(&[bad, PM_322]), bad);
    assert_eq!(date_of(&[bad, PM_322, AM_101]), bad);
}

#[test]
fn time_exception_is_an_error() {
    // t2: CTimeException is not a CStringException, so it leaves GetDate.
    let err = database_date(&[PM_322, "not a date"]).unwrap_err();
    assert_eq!(err.date, "not a date");
    assert_eq!(err.reason, too_short("not a date"));
    assert_eq!(err.exception.get_msg(), too_short("not a date").as_bytes());
    assert_eq!(
        err.to_string(),
        "BLAST database volume date 'not a date' is rejected: \
             Time string 'not a date' is too short for time format 'b d, Y  H:m P'"
    );

    // t1: the current date is rejected.
    let err = database_date(&["not a date", PM_322]).unwrap_err();
    assert_eq!(err.date, "not a date");
    assert_eq!(err.reason, too_short("not a date"));

    // t1 is built even when d is empty.
    let err = database_date(&["not a date", ""]).unwrap_err();
    assert_eq!(err.date, "not a date");

    // An unparsed first date fails at the first differing volume, wherever it is.
    let err = database_date(&["", "Feb 30, 2026  1:00 AM", PM_322]).unwrap_err();
    assert_eq!(err.date, "Feb 30, 2026  1:00 AM");
    assert_eq!(err.reason, invalid("Feb 30, 2026  1:00 AM"));
}

#[test]
fn current_date_is_parsed_before_new_date() {
    let no_digit = "Oct xx, 2026  3:22 PM";
    // t1 throws CTimeException before t2 is built.
    let err = database_date(&["not a date", no_digit]).unwrap_err();
    assert_eq!(err.date, "not a date");
    // t1 throws CStringException, which is caught, so t2 is never built.
    assert_eq!(date_of(&[no_digit, "not a date"]), no_digit);
}

// ---- parse_blast_db_date (CTime::x_Init with "b d, Y  H:m P") ----

#[test]
fn format_has_only_ported_characters() {
    // The `_` arm of the format loop is the literal branch; it must see only `,` and `:`.
    for &f in BLAST_DB_DATE_FORMAT.as_bytes() {
        assert!(b"bdYHmP,: ".contains(&f), "{:?}", f as char);
    }
}

#[test]
fn makeblastdb_dates_parse() {
    assert_eq!(parsed(PM_322), time(2026, 10, 12, 15, 22));
    assert_eq!(parsed(AM_101), time(2026, 10, 12, 1, 1));
    assert_eq!(parsed("Jan 1, 2026 11:05 AM"), time(2026, 1, 1, 11, 5));
}

#[test]
fn empty_string_is_the_empty_time() {
    // x_Init returns before reading the format.
    assert_eq!(parse_blast_db_date(""), Ok(time(0, 0, 0, 0, 0)));
}

#[test]
fn month_names_in_any_case() {
    let names = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    for (i, name) in names.iter().enumerate() {
        let month = i as u32 + 1;
        for name in [name.to_string(), name.to_uppercase(), name.to_lowercase()] {
            let text = format!("{name} 1, 2026  1:00 AM");
            assert_eq!(parsed(&text), time(2026, month, 1, 1, 0), "{text:?}");
        }
    }
}

#[test]
fn unknown_or_full_month_name_is_too_short() {
    // No name matches: nothing is consumed and `d` then finds no white space.
    for text in [
        "Okt 12, 2026  3:22 PM",
        "10 12, 2026  3:22 PM",
        "Oc",
        // "Oct" matches, and "ober" is not white space.
        "October 12, 2026  3:22 PM",
    ] {
        assert_eq!(time_exception(text), too_short(text));
    }
}

#[test]
fn am_pm_in_any_case() {
    for (word, hour) in [("AM", 3), ("am", 3), ("Am", 3), ("PM", 15), ("pm", 15)] {
        let text = format!("Oct 12, 2026  3:22 {word}");
        assert_eq!(parsed(&text), time(2026, 10, 12, hour, 22), "{text:?}");
    }
}

#[test]
fn missing_am_pm_after_white_space_leaves_twelve_hour_value() {
    // The white space before `P` is skipped, `P` matches neither word and consumes
    // nothing, the string is at its end, and hour_format stays e24.
    assert_eq!(parsed("Oct 12, 2026  3:22 "), time(2026, 10, 12, 3, 22));
    assert_eq!(parsed("Oct 12, 2026 12:22 "), time(2026, 10, 12, 0, 22));
}

#[test]
fn white_space_follows_the_format() {
    // A format space needs one or more white-space characters of any kind.
    assert_eq!(
        parsed("Oct 12, 2026 3:22 PM"),
        time(2026, 10, 12, 15, 22),
        "one space after the year"
    );
    assert_eq!(
        parsed("Oct \t12,\n 2026 \x0B\x0C\r 3:22   PM"),
        time(2026, 10, 12, 15, 22)
    );
    // Missing white space where the format has it.
    for text in [
        "Oct12, 2026  3:22 PM",
        "Oct 12,2026  3:22 PM",
        "Oct 12, 20263:22 PM",
        "Oct 12, 2026  3:22PM",
        "Oct 12, 2026  3:22",
    ] {
        assert_eq!(time_exception(text), too_short(text));
    }
    // White space where the format has none.
    for text in [
        " Oct 12, 2026  3:22 PM",
        "Oct 12 , 2026  3:22 PM",
        "Oct 12, 2026  3 :22 PM",
        "Oct 12, 2026  3: 22 PM",
    ] {
        assert_eq!(time_exception(text), too_short(text));
    }
}

#[test]
fn literals_must_match() {
    for text in [
        "Oct 12 2026  3:22 PM",
        "Oct 12; 2026  3:22 PM",
        "Oct 12, 2026  3.22 PM",
        "Oct 12, 2026  3",
    ] {
        assert_eq!(time_exception(text), too_short(text));
    }
}

#[test]
fn digit_counts() {
    // One digit is enough for day, hour and minute.
    assert_eq!(parsed("Oct 2, 2026  3:5 PM"), time(2026, 10, 2, 15, 5));
    // Leading zeros count as digits.
    assert_eq!(parsed("Oct 02, 2026 03:05 PM"), time(2026, 10, 2, 15, 5));
    // `d`, `H` and `m` read at most two digits, `Y` at most four; the next digit then
    // meets a literal or a format space.
    for text in [
        "Oct 012, 2026  3:22 PM",
        "Oct 12, 20260  3:22 PM",
        "Oct 12, 2026  003:22 PM",
        "Oct 12, 2026  3:022 PM",
    ] {
        assert_eq!(time_exception(text), too_short(text));
    }
}

#[test]
fn trailing_characters_are_too_long() {
    for text in [
        "Oct 12, 2026  3:22 PM ",
        "Oct 12, 2026  3:22 PMx",
        "Oct 12, 2026  3:22 PM\n",
        // `P` consumes nothing, so "XM" is left over.
        "Oct 12, 2026  3:22 XM",
        "Oct 12, 2026  3:22 P",
    ] {
        assert_eq!(time_exception(text), too_long(text));
    }
}

#[test]
fn values_out_of_range() {
    assert_eq!(
        time_exception("Oct 0, 2026  3:22 PM"),
        "Day value '0' is out of range"
    );
    assert_eq!(
        time_exception("Oct 32, 2026  3:22 PM"),
        "Day value '32' is out of range"
    );
    assert_eq!(
        time_exception("Oct 12, 1582  3:22 PM"),
        "Year value '1582' is out of range"
    );
    assert_eq!(
        time_exception("Oct 12, 26  3:22 PM"),
        "Year value '26' is out of range"
    );
    assert_eq!(
        time_exception("Oct 12, 2026 13:22 PM"),
        "Hour value '13' is out of range"
    );
    assert_eq!(
        time_exception("Oct 12, 2026  3:60 PM"),
        "Minute value '60' is out of range"
    );
    // The range check comes before the rest of the string is read.
    assert_eq!(
        time_exception("Oct 99, garbage"),
        "Day value '99' is out of range"
    );
    // Without the `:`, `H` reads the two digits "32" and fails its range check before
    // the literal is compared.
    assert_eq!(
        time_exception("Oct 12, 2026  322 PM"),
        "Hour value '32' is out of range"
    );
}

#[test]
fn year_boundaries_and_twelve_bit_field() {
    assert_eq!(parsed("Oct 12, 1583  3:22 PM"), time(1583, 10, 12, 15, 22));
    assert_eq!(parsed("Oct 12, 4095  3:22 PM"), time(4095, 10, 12, 15, 22));
    // SET_YEAR keeps value & 0xFFF: 4096 -> 0, which IsValid rejects (Year() < 1583).
    let text = "Oct 12, 4096  3:22 PM";
    assert_eq!(time_exception(text), invalid(text));
    // 5679 = 4096 + 1583 -> 1583, and 9999 = 8192 + 1807 -> 1807.
    assert_eq!(parsed("Oct 12, 5679  3:22 PM"), time(1583, 10, 12, 15, 22));
    assert_eq!(parsed("Oct 12, 9999  3:22 PM"), time(1807, 10, 12, 15, 22));
    // So such a year compares as the reduced one.
    assert_eq!(
        date_of(&["Oct 12, 2026  3:22 PM", "Oct 12, 9999  3:22 PM"]),
        "Oct 12, 2026  3:22 PM"
    );
}

#[test]
fn day_must_exist_in_month() {
    assert_eq!(parsed("Feb 29, 2024  1:00 AM"), time(2024, 2, 29, 1, 0));
    assert_eq!(parsed("Feb 29, 2000  1:00 AM"), time(2000, 2, 29, 1, 0));
    assert_eq!(parsed("Apr 30, 2026  1:00 AM"), time(2026, 4, 30, 1, 0));
    assert_eq!(parsed("Dec 31, 2026  1:00 AM"), time(2026, 12, 31, 1, 0));
    for text in [
        "Feb 29, 2026  1:00 AM",
        "Feb 29, 1900  1:00 AM",
        "Feb 30, 2024  1:00 AM",
        "Apr 31, 2026  1:00 AM",
        "Nov 31, 2026  1:00 AM",
    ] {
        assert_eq!(time_exception(text), invalid(text));
    }
}

#[test]
fn non_ascii_text_is_rejected_without_panic() {
    // Bytes above 0x7F are neither white space, digits nor letters of the names.
    let text = "Oct 12, 2026  3:22 PM\u{00a0}";
    assert_eq!(time_exception(text), too_long(text));
    let text = "\u{6642}\u{9593} 12, 2026  3:22 PM";
    assert_eq!(time_exception(text), too_short(text));
    // No digit after the white space: CStringException.
    assert!(matches!(
        parse_blast_db_date("Oct \u{ff11}\u{ff12}, 2026  3:22 PM"),
        Err(ParseFailure::StringException(_))
    ));
}

/// The rows of `tests/fixtures/makeblastdb_failure/index_dates.tsv`: NCBI's `blastdbcmd
/// -metadata` on index files whose date is a prefix of the written one (a failed write).
#[test]
fn the_dates_of_a_cut_short_index_are_read_as_ncbi_reads_them() {
    let table = include_str!("../../tests/fixtures/makeblastdb_failure/index_dates.tsv");
    let mut rows = 0;
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        let result = iso8601_date_time_exception(columns[0]);
        match columns[1] {
            "ok" => assert_eq!(result.as_deref(), Ok(columns[2]), "{:?}", columns[0]),
            _ => {
                let exception = result.expect_err(columns[0]);
                let site = exception.site();
                assert_eq!(
                    (
                        site.line.to_string(),
                        site.kind,
                        site.function.unwrap_or(""),
                        exception.get_msg()
                    ),
                    (
                        columns[2].to_string(),
                        columns[3],
                        columns[4],
                        columns[5].as_bytes()
                    ),
                    "{:?}",
                    columns[0]
                );
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 95);
}

//! The local time of a build (`CTime(CTime::eCurrent)`) and its two text
//! forms: the creation date of an index file and the time of the log.

use super::local_offset::local_utc_offset;

/// A local calendar time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalTime {
    pub year: i32,
    /// 1 to 12.
    pub month: u32,
    /// 1 to 31.
    pub day: u32,
    /// 0 to 23.
    pub hour: u32,
    pub minute: u32,
    pub second: u32,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:102-105
// ```c++
// static const char* kMonthAbbr[12] = {
//     "Jan", "Feb", "Mar", "Apr", "May", "Jun",
//     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
// };
// ```
const MONTH_ABBR: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

impl LocalTime {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1722-1734
    // ```c++
    //     // Bind values to internal variables
    //     struct tm *t;
    //
    // #ifdef HAVE_LOCALTIME_R
    //     struct tm temp;
    //     if (GetTimeZone() == eLocal) {
    //         localtime_r(&timer, &temp);
    //     } else {
    //         gmtime_r(&timer, &temp);
    //     }
    //     t = &temp;
    // #else
    //     t = ( GetTimeZone() == eLocal ) ? localtime(&timer) : gmtime(&timer);
    // ```
    // `CTime(CTime::eCurrent)` is the local time of the C library. LOSAT
    // calls no C library: the zone's offset from UTC comes from
    // `local_offset`, which reads the zone as the C library does, and where
    // the zone cannot be found (Windows, WASI, an unreadable zone file) the
    // time is UTC.
    /// The current local time, or `None` where the platform has no clock.
    pub fn now() -> Option<Self> {
        let seconds = unix_seconds()?;
        let offset = local_utc_offset(seconds).unwrap_or(0);
        Some(Self::from_unix_seconds(seconds + i64::from(offset)))
    }

    /// The calendar time `seconds` after 1970-01-01 00:00:00 (Gregorian
    /// calendar, days of 86,400 seconds).
    pub fn from_unix_seconds(seconds: i64) -> Self {
        let days = seconds.div_euclid(86_400);
        let second_of_day = seconds.rem_euclid(86_400);
        // Days since 0000-03-01, in 400-year eras that start on a 1 March.
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let day_of_era = z.rem_euclid(146_097);
        let year_of_era =
            (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
        let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
        // Months counted from March, so that the leap day is the last day.
        let month_from_march = (5 * day_of_year + 2) / 153;
        let day = day_of_year - (153 * month_from_march + 2) / 5 + 1;
        let month = if month_from_march < 10 {
            month_from_march + 3
        } else {
            month_from_march - 9
        };
        let year = year_of_era + era * 400 + i64::from(month <= 2);
        Self {
            year: year as i32,
            month: month as u32,
            day: day as u32,
            hour: (second_of_day / 3_600) as u32,
            minute: (second_of_day % 3_600 / 60) as u32,
            second: (second_of_day % 60) as u32,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:88-97
    // ```c++
    //     CTime now(CTime::eCurrent);
    //
    //     m_Date = now.AsString("b d, Y  ");
    //     string t = now.AsString("H:m P");
    //
    //     if (t[0] == '0') {
    //         t.assign(t, 1, t.size() - 1);
    //     }
    //
    //     m_Date += t;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1592-1617
    // ```c++
    //         case 'Y': s_AddZeroPadInt(str, t->Year(), 4);       break;
    //         case 'M': s_AddZeroPadInt2(str, t->Month());        break;
    //         case 'b': str += kMonthAbbr[t->Month()-1];          break;
    //         case 'B': str += kMonthFull[t->Month()-1];          break;
    //         case 'D': s_AddZeroPadInt2(str, t->Day());          break;
    //         case 'd': s_AddZeroPadInt(str, t->Day(),1);         break;
    //         case 'h': s_AddZeroPadInt2(str, t->Hour());         break;
    //         case 'H': s_AddZeroPadInt2(str, (t->Hour()+11) % 12+1);
    //         ...
    //         case 'm': s_AddZeroPadInt2(str, t->Minute());       break;
    //         case 's': s_AddZeroPadInt2(str, t->Second());       break;
    //         ...
    //         case 'P': str += ( t->Hour() < 12) ? "AM" : "PM" ;  break;
    // ```
    /// The creation date of the index files and of the alias file, for
    /// example `Oct 12, 2026  3:22 PM`.
    pub fn index_date(&self) -> String {
        let month = MONTH_ABBR[(self.month.clamp(1, 12) - 1) as usize];
        let hour12 = (self.hour + 11) % 12 + 1;
        let half = if self.hour < 12 { "AM" } else { "PM" };
        format!(
            "{month} {}, {:04}  {hour12}:{:02} {half}",
            self.day, self.year, self.minute
        )
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:121
    // ```c++
    // static const char* kDefaultFormatTime      = "M/D/Y h:m:s";
    // ```
    /// `CTime::AsString()` with the default format, for example
    /// `10/12/2026 20:29:48`.
    pub fn log_time(&self) -> String {
        format!(
            "{:02}/{:02}/{:04} {:02}:{:02}:{:02}",
            self.month, self.day, self.year, self.hour, self.minute, self.second
        )
    }
}

/// Seconds since 1970-01-01 00:00:00 UTC by the system clock.
#[cfg(any(unix, windows, target_os = "wasi"))]
fn unix_seconds() -> Option<i64> {
    let elapsed = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?;
    i64::try_from(elapsed.as_secs()).ok()
}

#[cfg(not(any(unix, windows, target_os = "wasi")))]
fn unix_seconds() -> Option<i64> {
    None
}

#[cfg(test)]
mod tests {
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
}

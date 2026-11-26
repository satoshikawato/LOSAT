//! The local time of a build (`CTime` with `CTime::eCurrent`) and its two text
//! forms: the creation date of an index file and the time of the log; and the
//! same instant as `time` gives it when called with a null pointer (the seed of makeblastdb's
//! `-randomize`).

use super::local_offset::local_utc_offset;

/// A local calendar time, and the instant it is.
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
    /// What `time` gives with a null pointer: the instant, in seconds since 1970-01-01 00:00:00
    /// UTC.
    pub unix_seconds: i64,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:102-105
const MONTH_ABBR: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

impl LocalTime {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1722-1734
    // `CTime` with `CTime::eCurrent` is the local time of the C library. LOSAT
    // calls no C library: the zone's offset from UTC comes from
    // `local_offset`, which finds the zone where the C library finds it
    // (glibc on Unix and on WASI, from the environment and the directories
    // that the host gives, read by the crate jiff; UCRT on Windows), and
    // where no zone is found or read (no zone file, a WASI host that gives
    // neither `TZ` nor a zone file, a `TZ` that jiff does not read) the time
    // is UTC. A leap second of a zone file is second 60, as glibc writes it.
    /// The current local time, or `None` where the platform has no clock.
    pub fn now() -> Option<Self> {
        let seconds = unix_seconds()?;
        let offset = local_utc_offset(seconds).unwrap_or_default();
        let mut time = Self::from_unix_seconds(seconds + offset.seconds);
        if offset.leap_second {
            time.second = 60;
        }
        time.unix_seconds = seconds;
        Some(time)
    }

    /// The calendar time `seconds` after 1970-01-01 00:00:00 (Gregorian
    /// calendar, days of 86,400 seconds): the time in UTC of the instant `seconds`.
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
            unix_seconds: seconds,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:88-97
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1592-1617
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1595
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1592-1617
    /// `CTime::AsString` in the format `CSeqDB::kBlastDbDateFormat` (`b d, Y  H:m P`): the day
    /// without padding, the hour of 12 with two digits, for example `Nov 5, 2026  07:00 AM`. The
    /// creation date of a binary seqid list (`WriteBlastSeqidlistFile`); the index files' date
    /// drops the hour's zero (`index_date`).
    pub fn blast_db_date(&self) -> String {
        let month = MONTH_ABBR[(self.month.clamp(1, 12) - 1) as usize];
        let hour12 = (self.hour + 11) % 12 + 1;
        let half = if self.hour < 12 { "AM" } else { "PM" };
        format!(
            "{month} {}, {:04}  {hour12:02}:{:02} {half}",
            self.day, self.year, self.minute
        )
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:121
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
#[path = "time_tests.rs"]
mod tests;

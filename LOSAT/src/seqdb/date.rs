//! Date of a BLAST database, chosen from the date strings of its volumes.
//!
//! Port of `CSeqDBImpl::GetDate` of NCBI BLAST+ 2.17.0 and of the part of `CTime` string
//! parsing that it reaches with the format `"b d, Y  H:m P"`. No other format symbol, format
//! flag or time zone is ported.
//!
//! What NCBI does, as read from the source cited below:
//!
//! - The first non-empty volume date becomes the database date. A later volume date that
//!   differs from it as a string replaces it only when `CTime` reads it as a later time.
//! - `b` is one of `Jan` .. `Dec` in any letter case. `d` is up to two digits in 1..=31. `Y` is
//!   up to four digits, at least 1583, kept in a 12-bit field. `H` is up to two digits in
//!   0..=12, kept modulo 12. `m` is up to two digits in 0..=59. `P` is `AM` or `PM` in any
//!   letter case; `PM` adds 12 to the hour, so `12 AM` is hour 0 and `12 PM` is hour 12.
//! - A run of spaces in the format needs at least one white-space character in the string and
//!   skips all of them. Where the format has no space, the string must have none. `,` and `:`
//!   must be present. The string must end where the format ends.
//! - A number with no digit makes `NStr::StringToLong` throw `CStringException`, which
//!   `GetDate` catches: the current date stays. Every other failure is a `CTimeException`,
//!   which `GetDate` does not catch: [`database_date`] returns an error.

use std::fmt;

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::{is_space, string_to_int8};

/// Error of `database_date`: NCBI's CTime rejects one of two differing volume dates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DateError {
    /// The volume date string that `CTime` rejects.
    pub date: String,
    /// The message of the `CTimeException` that NCBI throws for it.
    pub reason: String,
    /// The exception itself, with the place that NCBI's report names.
    pub exception: Box<NcbiException>,
}

impl fmt::Display for DateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "BLAST database volume date '{}' is rejected: {}",
            self.date, self.reason
        )
    }
}

impl std::error::Error for DateError {}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1595
const BLAST_DB_DATE_FORMAT: &str = "b d, Y  H:m P";

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:102-105
const MONTH_ABBR: [&[u8; 3]; 12] = [
    b"Jan", b"Feb", b"Mar", b"Apr", b"May", b"Jun", b"Jul", b"Aug", b"Sep", b"Oct", b"Nov", b"Dec",
];

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:97-99
const DAYS_IN_MONTH: [u32; 12] = [31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

/// How `CTime::CTime(const string&, const CTimeFormat&)` fails for the BLAST database format.
///
/// The two exception classes derive from `CCoreException` separately, so
/// `catch(CStringException &)` in `CSeqDBImpl::GetDate` catches only the first.
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:4499
// NCBI reference (598d8ae6): c++/include/corelib/ncbitime.hpp:2074
#[derive(Debug, Clone, PartialEq, Eq)]
enum ParseFailure {
    /// `NStr::StringToLong` threw this `CStringException`: a number has no digit.
    StringException(NcbiException),
    /// `CTime::x_Init` threw this `CTimeException`.
    TimeException(NcbiException),
}

impl ParseFailure {
    /// `CTime::x_Init` throws `CTimeException` `kind` (`eFormat`, `eArgument` or `eConvert`) at
    /// `line` of `ncbitime.cpp`.
    fn time(kind: &'static str, line: u32, reason: String) -> Self {
        ParseFailure::TimeException(NcbiException::new(
            ThrowSite {
                file: "src/corelib/ncbitime.cpp",
                line,
                kind,
                function: Some("ncbi::CTime::x_Init"),
            },
            reason,
        ))
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:537-540
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HourFormat {
    E24,
    Am,
    Pm,
}

/// The `CTime` fields that the BLAST database format sets.
///
/// Seconds and nanoseconds stay zero, and the time zone is `eLocal` (the default of the
/// constructor) for every time built from this format, which has no `Z`, `z` or `o` symbol.
// NCBI reference (598d8ae6): c++/include/corelib/ncbitime.hpp:1248-1252
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct Time {
    year: u32,
    month: u32,
    day: u32,
    hour: u32,
    min: u32,
}

impl Time {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:2211-2215
    fn is_leap(&self) -> bool {
        (self.year.is_multiple_of(4) && !self.year.is_multiple_of(100))
            || self.year.is_multiple_of(400)
    }

    /// `CTime::IsValid` for a time that the whole format has filled.
    ///
    /// The other tests of `IsValid` cannot fail here: the day is at least 1, so the time is
    /// not empty; the time zone is `eLocal`; the hour is at most 11 + 12; the minute passed
    /// its range check; seconds and nanoseconds are zero. The month test cannot fail either
    /// (an unmatched month name stops the format loop before validation); it is kept because
    /// it guards the table index.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:2023-2033
    fn is_valid(&self) -> bool {
        if self.year < 1583 {
            return false;
        }
        if self.month < 1 || self.month > 12 {
            return false;
        }
        if self.month == 2 {
            if self.day < 1 || self.day > (if self.is_leap() { 29 } else { 28 }) {
                return false;
            }
        } else if self.day < 1 || self.day > DAYS_IN_MONTH[(self.month - 1) as usize] {
            return false;
        }
        true
    }

    /// `*self > *t` of `CTime::operator>`.
    ///
    /// `tmp.ToTime(GetTimeZone())` changes nothing because both times are `eLocal`
    /// (`ToTime` converts only `if (GetTimeZone() != tz)`), and `Second()` and `NanoSecond()`
    /// are zero on both sides. An empty time (all fields zero, from an empty string) is
    /// compared by the same fields and is never later than a valid one.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:2137-2162
    fn is_later_than(&self, t: &Time) -> bool {
        (self.year, self.month, self.day, self.hour, self.min)
            > (t.year, t.month, t.day, t.hour, t.min)
    }
}

/// `*sss` of a NUL-terminated C string: the byte at `pos`, or 0 at the end.
fn byte_at(s: &[u8], pos: usize) -> u8 {
    s.get(pos).copied().unwrap_or(0)
}

/// Reads up to `len` decimal digits at `*sss` and converts them as `NStr::StringToLong` does.
///
/// With no digit, `value_str` is empty, `StringToInt8` finds `pos == pos0`, and it throws
/// `CStringException` because the flags (0) lack `fConvErr_NoThrow`.
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:742-746
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:572-575,860-862
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:658-666
// `StringToInt8` with flags 0 is `ncbi_str::string_to_int8`; 4 digits at most do not overflow
// a 64-bit `long`.
fn read_number(s: &[u8], sss: &mut usize, len: usize) -> Result<i64, ParseFailure> {
    let start = *sss;
    while *sss - start < len && byte_at(s, *sss).is_ascii_digit() {
        *sss += 1;
    }
    string_to_int8(&s[start..*sss], 0).map_err(ParseFailure::StringException)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:149-156
// `line` is the line of the macro's use in `x_Init` (751 for the year, 772, 782, 788:
// ncbitime.cpp).
fn check_range(value: i64, what: &str, min: i64, max: i64, line: u32) -> Result<(), ParseFailure> {
    if value < min || value > max {
        return Err(ParseFailure::time(
            "CTimeException::eArgument",
            line,
            format!("{what} value '{value}' is out of range"),
        ));
    }
    Ok(())
}

/// `CTime t(text, CSeqDB::kBlastDbDateFormat)`: the fields of the time, or the exception.
///
/// The `const char*` format converts to `CTimeFormat` with `fDefault` flags, which `SetFormat`
/// turns into `fFormat_Simple | fMatch_Strict`: every format character can be a format
/// symbol, white space is matched (not ignored), and neither a short string nor a short
/// format is accepted. `x_Init` runs with `err_action = eErr_Throw`, and its own
/// `catch (const CException&)` rethrows in that mode.
// NCBI reference (598d8ae6): c++/include/corelib/ncbitime.hpp:209,472-474,1199-1200
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:361-371
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1017-1029
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:941-947
fn parse_blast_db_date(text: &str) -> Result<Time, ParseFailure> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:509-515
    let mut time = Time::default();
    if text.is_empty() {
        return Ok(time);
    }

    let fmt = BLAST_DB_DATE_FORMAT;
    let s = text.as_bytes();
    let mut sss = 0usize;
    let mut hour_format = HourFormat::E24;
    let mut is_12hour = false;
    let mut is_format_space = false;

    // `format_matched` is false when the loop stops at a `break`, that is, when NCBI leaves
    // its loop with `*fff != '\0'`.
    let format_matched = 'format: {
        for &f in fmt.as_bytes() {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:545,549-578
            // `x_Init` calls the C library's `isspace` (`ncbi_str::is_space`): the NCBI C++
            // Toolkit code that BLAST+ runs does not call `setlocale`, so the "C" locale is in
            // effect.
            if is_space(f) {
                is_format_space = true;
                continue;
            }
            if is_space(byte_at(s, sss)) {
                if is_format_space {
                    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:460-461
                    while is_space(byte_at(s, sss)) {
                        sss += 1;
                    }
                } else {
                    break 'format false;
                }
            } else if is_format_space {
                break 'format false;
            }
            is_format_space = false;

            match f {
                // A name that matches none of the twelve leaves the month at 0 and consumes
                // nothing; the next format symbol then stops the loop on non-matching spaces.
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:617-635
                b'b' => {
                    for (i, name) in MONTH_ABBR.iter().enumerate() {
                        let namelen = name.len();
                        let matches = s
                            .get(sss..sss + namelen)
                            .is_some_and(|w| w.eq_ignore_ascii_case(*name));
                        if matches {
                            sss += namelen;
                            time.month = (i as u32 + 1) & 0x0F;
                            break;
                        }
                    }
                }

                // Text that is neither AM nor PM consumes nothing and leaves `e24`.
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:721-730
                b'P' => {
                    let word = s.get(sss..sss + 2);
                    if word.is_some_and(|w| w.eq_ignore_ascii_case(b"AM")) {
                        hour_format = HourFormat::Am;
                        sss += 2;
                    } else if word.is_some_and(|w| w.eq_ignore_ascii_case(b"PM")) {
                        hour_format = HourFormat::Pm;
                        sss += 2;
                    }
                }

                // `Y` reads up to four digits. The 12-bit field keeps `value & 0xFFF`, so a
                // year of 4096 or more is stored reduced (5679 becomes 1583).
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:735-737,750-754
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:166,177
                b'Y' => {
                    let value = read_number(s, &mut sss, 4)?;
                    check_range(value, "Year", 1583, i64::from(i32::MAX), 751)?;
                    time.year = (value & 0xFFF) as u32;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:770-775
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:168,179
                b'd' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Day", 1, 31, 772)?;
                    time.day = (value & 0x1F) as u32;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:781-786
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:170,180
                b'H' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Hour", 0, 12, 782)?;
                    time.hour = ((value % 12) & 0x1F) as u32;
                    is_12hour = true;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:787-791
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:171,181
                b'm' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Minute", 0, 59, 788)?;
                    time.min = (value & 0x3F) as u32;
                }

                // `,` and `:`, the characters of the format that are neither spaces nor in
                // kFormatSymbolsTime.
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:129,604-611
                _ => {
                    if f == byte_at(s, sss) {
                        sss += 1;
                        continue;
                    }
                    break 'format false;
                }
            }
        }
        true
    };

    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:849-851
    if is_12hour && hour_format == HourFormat::Pm {
        time.hour = (time.hour + 12) & 0x1F;
    }

    // NCBI skips trailing white space here only `if (isspace(*sss) && is_format_space)`.
    // That condition is never true for this format: the format ends in the symbol `P`, so a
    // completed loop leaves `is_format_space` false, and no `break` above is taken with both
    // a white-space character at `sss` and `is_format_space` true. Trailing white space is
    // therefore "too long".
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:857-859,863-872
    if !format_matched {
        return Err(ParseFailure::time(
            "CTimeException::eFormat",
            865,
            format!("Time string '{text}' is too short for time format '{fmt}'"),
        ));
    }
    if byte_at(s, sss) != 0 {
        return Err(ParseFailure::time(
            "CTimeException::eFormat",
            870,
            format!("Time string '{text}' is too long for time format '{fmt}'"),
        ));
    }

    // A completed loop has set year, month, day and time, so `ptcache` is 2222 and no field
    // is taken from the current time.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:880-884
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:932-934
    if !time.is_valid() {
        return Err(ParseFailure::time(
            "CTimeException::eConvert",
            933,
            format!("Unable to convert string '{text}' to CTime"),
        ));
    }
    Ok(time)
}

/// Date of a database from its volumes' date strings, in volume order
/// (CSeqDBImpl::GetDate). An empty slice gives an empty string.
///
/// Each string must already be cut at its first NUL byte (`x_FixString`).
///
/// The first non-empty date is taken as it is. A later date that differs from the current
/// one replaces it when `CTime` reads it as later. When `NStr::StringToLong` throws
/// `CStringException` for either of the two (a number without a digit), NCBI catches it and
/// the current date stays. When `CTime` throws `CTimeException` (any other mismatch with the
/// format, a value out of range, or an invalid date), NCBI does not catch it and `GetDate`
/// fails; this function then returns `Err` with the rejected string and NCBI's message.
/// The current date is parsed before the new one, and the new one is not parsed when the
/// current one already failed.
///
/// NCBI also caches the result in `m_Date`; that is left to the caller.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1010-1035
pub fn database_date(volume_dates: &[&str]) -> Result<String, DateError> {
    let mut date = String::new();

    for &d in volume_dates {
        if date.is_empty() {
            date = d.to_owned();
        } else if d != date {
            let t1 = match parse_blast_db_date(&date) {
                Ok(time) => time,
                Err(ParseFailure::StringException(_)) => continue,
                Err(ParseFailure::TimeException(exception)) => {
                    return Err(DateError {
                        date,
                        reason: String::from_utf8_lossy(exception.get_msg()).into_owned(),
                        exception: Box::new(exception),
                    });
                }
            };
            let t2 = match parse_blast_db_date(d) {
                Ok(time) => time,
                Err(ParseFailure::StringException(_)) => continue,
                Err(ParseFailure::TimeException(exception)) => {
                    return Err(DateError {
                        date: d.to_owned(),
                        reason: String::from_utf8_lossy(exception.get_msg()).into_owned(),
                        exception: Box::new(exception),
                    });
                }
            };

            if t2.is_later_than(&t1) {
                date = d.to_owned();
            }
        }
    }

    Ok(date)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1730-1732
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1512-1520
// The ISO 8601 form is `Y-M-DTh:m:s` (ncbitime.cpp:398, `eISO8601_DateTimeSec`); the
// stored date has no seconds. An empty date is the empty time, written as the empty
// string. The constructor has checked the time a date string gives (`IsValid`).
/// `CTime(date, "b d, Y  H:m P").AsString(eISO8601_DateTimeSec)`, the `last-updated` of
/// `blastdbcmd -metadata`: `YYYY-MM-DDThh:mm:00`. `Err` with NCBI's message where the
/// constructor throws (`CTimeException`, or `CStringException` for a number without digits).
pub fn iso8601_date_time(date: &str) -> Result<String, String> {
    match iso8601(date) {
        Ok(text) => Ok(text),
        Err(ParseFailure::TimeException(exception)) => {
            Err(String::from_utf8_lossy(exception.get_msg()).into_owned())
        }
        Err(ParseFailure::StringException(_)) => {
            Err(format!("a number of the date '{date}' has no digit"))
        }
    }
}

/// [`iso8601_date_time`] with NCBI's exception, as `CSeqDB::GetDBMetaData` lets it leave
/// (`CTime`'s constructor throws `CTimeException`, or `CStringException` for a number without
/// digits): NCBI's applications end with it, exit status 255.
pub fn iso8601_date_time_exception(date: &str) -> Result<String, NcbiException> {
    iso8601(date).map_err(|failure| match failure {
        ParseFailure::TimeException(exception) | ParseFailure::StringException(exception) => {
            exception
        }
    })
}

fn iso8601(date: &str) -> Result<String, ParseFailure> {
    if date.is_empty() {
        return Ok(String::new());
    }
    let time = parse_blast_db_date(date)?;
    Ok(format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:00",
        time.year, time.month, time.day, time.hour, time.min
    ))
}

#[cfg(test)]
#[path = "date_tests.rs"]
mod tests;

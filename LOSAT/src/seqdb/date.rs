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

/// Error of `database_date`: NCBI's CTime rejects one of two differing volume dates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DateError {
    /// The volume date string that `CTime` rejects.
    pub date: String,
    /// The message of the `CTimeException` that NCBI throws for it.
    pub reason: String,
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
// ```c++
// const char* CSeqDB::kBlastDbDateFormat = "b d, Y  H:m P";
// ```
const BLAST_DB_DATE_FORMAT: &str = "b d, Y  H:m P";

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:102-105
// ```c++
// static const char* kMonthAbbr[12] = {
//     "Jan", "Feb", "Mar", "Apr", "May", "Jun",
//     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
// };
// ```
const MONTH_ABBR: [&[u8; 3]; 12] = [
    b"Jan", b"Feb", b"Mar", b"Apr", b"May", b"Jun", b"Jul", b"Aug", b"Sep", b"Oct", b"Nov", b"Dec",
];

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:97-99
// ```c++
// static int s_DaysInMonth[12] = {
//     31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
// };
// ```
const DAYS_IN_MONTH: [u32; 12] = [31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

/// How `CTime::CTime(const string&, const CTimeFormat&)` fails for the BLAST database format.
///
/// The two exception classes derive from `CCoreException` separately, so
/// `catch(CStringException &)` in `CSeqDBImpl::GetDate` catches only the first.
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:4499
// ```c++
// class NCBI_XNCBI_EXPORT CStringException : public CParseTemplException<CCoreException>
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbitime.hpp:2074
// ```c++
// class NCBI_XNCBI_EXPORT CTimeException : public CCoreException
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
enum ParseFailure {
    /// `NStr::StringToLong` threw `CStringException`: a number has no digit.
    StringException,
    /// `CTime::x_Init` threw `CTimeException` with this message.
    TimeException(String),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:537-540
// ```c++
//     enum EHourFormat{
//         e24, eAM, ePM
//     };
//     EHourFormat hour_format = e24;
// ```
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
// ```c++
//         unsigned int  year        NCBI_TIME_BITFIELD(12);  // 4 digits
//         unsigned char month       NCBI_TIME_BITFIELD( 4);  // 0..12
//         unsigned char day         NCBI_TIME_BITFIELD( 5);  // 0..31
//         unsigned char hour        NCBI_TIME_BITFIELD( 5);  // 0..23
//         unsigned char min         NCBI_TIME_BITFIELD( 6);  // 0..59
// ```
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
    // ```c++
    // bool CTime::IsLeap(void) const
    // {
    //     int year = Year();
    //     return (year % 4 == 0  &&  year % 100 != 0)  ||  year % 400 == 0;
    // }
    // ```
    fn is_leap(&self) -> bool {
        (self.year % 4 == 0 && self.year % 100 != 0) || self.year % 400 == 0
    }

    /// `CTime::IsValid` for a time that the whole format has filled.
    ///
    /// The other tests of `IsValid` cannot fail here: the day is at least 1, so the time is
    /// not empty; the time zone is `eLocal`; the hour is at most 11 + 12; the minute passed
    /// its range check; seconds and nanoseconds are zero. The month test cannot fail either
    /// (an unmatched month name stops the format loop before validation); it is kept because
    /// it guards the table index.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:2023-2033
    // ```c++
    //     if (Year() < 1583) // first Gregorian date February 24, 1582
    //         return false;
    //     if (Month()  < 1  ||  Month()  > 12)
    //         return false;
    //     if (Month() == 2) {
    //         if (Day() < 1 ||  Day() > (IsLeap() ? 29 : 28))
    //             return false;
    //     } else {
    //         if (Day() < 1 ||  Day() > s_DaysInMonth[Month() - 1])
    //             return false;
    //     }
    // ```
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
    // ```c++
    // bool CTime::operator> (const CTime& t) const
    // {
    //     CTime tmp(t);
    //     if ( !tmp.IsEmptyDate() ) {
    //         tmp.ToTime(GetTimeZone());
    //     }
    //     if (Year()   > tmp.Year())
    //         return true;
    //     if (Year()   < tmp.Year())
    //         return false;
    //     if (Month()  > tmp.Month())
    //         return true;
    //     if (Month()  < tmp.Month())
    //         return false;
    //     if (Day()    > tmp.Day())
    //         return true;
    //     if (Day()    < tmp.Day())
    //         return false;
    //     if (Hour()   > tmp.Hour())
    //         return true;
    //     if (Hour()   < tmp.Hour())
    //         return false;
    //     if (Minute() > tmp.Minute())
    //         return true;
    //     if (Minute() < tmp.Minute())
    //         return false;
    // ```
    fn is_later_than(&self, t: &Time) -> bool {
        (self.year, self.month, self.day, self.hour, self.min)
            > (t.year, t.month, t.day, t.hour, t.min)
    }
}

/// C `isspace` in the "C" locale: space, `\t`, `\n`, `\v`, `\f`, `\r`.
///
/// `x_Init` calls the C library's `isspace`. The NCBI C++ Toolkit code that BLAST+ runs does
/// not call `setlocale`, so the "C" locale is in effect.
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:460-461
// ```c++
// #define SKIP_SPACES(s) \
//     while (isspace((unsigned char)(*s))) ++s;
// ```
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
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
// ```c++
//         for ( ; len  &&  *sss  &&  isdigit((unsigned char)(*sss));  --len) {
//             *s++ = *sss++;
//         }
//         *s = '\0';
//         long value = NStr::StringToLong(value_str);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:572-575,860-862
// ```c++
//         if ( !(flags & NStr::fConvErr_NoThrow) ) {                    \
//             err_guard.Throw();                                        \
//             NCBI_THROW2(CStringException, eConvert,                   \
//                         err_guard.Message(str, #to_type, msg), pos);  \
// ...
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
// ```
fn read_number(s: &[u8], sss: &mut usize, mut len: usize) -> Result<i64, ParseFailure> {
    let mut value: i64 = 0;
    let mut digits = 0;
    while len > 0 && byte_at(s, *sss).is_ascii_digit() {
        value = value * 10 + i64::from(byte_at(s, *sss) - b'0');
        *sss += 1;
        digits += 1;
        len -= 1;
    }
    if digits == 0 {
        return Err(ParseFailure::StringException);
    }
    Ok(value)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:149-156
// ```c++
// #define CHECK_RANGE2(value, what, min, max, err_action) \
//     if ( value < min  ||  value > max ) {  \
//         if ( err_action == eErr_Throw ) { \
//             NCBI_THROW(CTimeException, eArgument, \
//                        what " value '" + \
//                        NStr::Int8ToString((Int8)value) + "' is out of range"); \
//         } else { return false; } \
//     }
// ```
fn check_range(value: i64, what: &str, min: i64, max: i64) -> Result<(), ParseFailure> {
    if value < min || value > max {
        return Err(ParseFailure::TimeException(format!(
            "{what} value '{value}' is out of range"
        )));
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
// ```c++
//     CTimeFormat(const char* fmt, TFlags flags = fDefault);
// ...
//     explicit CTime(const string& str, const CTimeFormat& fmt = kEmptyStr,
//                    ETimeZone tz = eLocal,
//                    ETimeZonePrecision tzp = eTZPrecisionDefault);
// ...
//     bool x_Init(const string& str, const CTimeFormat& fmt,
//                 EErrAction err_action = eErr_Throw);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:361-371
// ```c++
//     if ((flags & f) == 0) {
//         flags |= fFormat_Simple;  // default
//     }
// ...
//     if ((flags & (fMatch_Strict | fMatch_Weak)) == 0) {
//         flags |= fMatch_Strict;  // default
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:1017-1029
// ```c++
//     memset(&m_Data, 0, sizeof(m_Data));
// ...
//     if ( !m_Data.tz ) {
//         m_Data.tz = tz;
//     }
//     m_Data.tzprec = tzp;
//     if (fmt.IsEmpty()) {
//         x_Init(str, GetFormat());
//     } else {
//         x_Init(str, fmt);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:941-947
// ```c++
//     } catch (const CException&) {
//         if (err_action == eErr_Throw) {
//             throw;
//         } else {
//             return false;
//         }
//     }
// ```
fn parse_blast_db_date(text: &str) -> Result<Time, ParseFailure> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:509-515
    // ```c++
    //     Clear();
    // ...
    //     if ( str.empty() ) {
    //         return true;
    //     }
    // ```
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
            // ```c++
            //     bool is_ignore_spaces = ((format.GetFlags() & CTimeFormat::fMatch_IgnoreSpaces) != 0);
            // ...
            //     for (fff = fmt.c_str();  *fff != '\0';  ++fff) {
            // ...
            //             if (isspace((unsigned char)(*fff))) {
            //                 is_format_space = true;
            //                 continue;
            //             }
            //             if (is_ignore_spaces) {
            // ...
            //             } else {
            //                 // match spaces
            //                 if (isspace((unsigned char)(*sss))) {
            //                     if (is_format_space) {
            //                         // Skip space symbols in time string
            //                         SKIP_SPACES(sss);
            //                     } else {
            //                         break;  // error: non-matching spaces
            //                     }
            //                 } else {
            //                     if (is_format_space) {
            //                         break;  // error: non-matching spaces
            //                     }
            //                 }
            //                 is_format_space = false;
            //             }
            // ```
            if is_c_space(f) {
                is_format_space = true;
                continue;
            }
            if is_c_space(byte_at(s, sss)) {
                if is_format_space {
                    while is_c_space(byte_at(s, sss)) {
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
                // ```c++
                //         if (*fff == 'b'  ||  *fff == 'B') {
                //             const char** name;
                //             if (*fff == 'b') {
                //                 name = kMonthAbbr;
                // ...
                //             for (unsigned char i = 0;  i < 12;  ++i) {
                //                 size_t namelen = strlen(*name);
                //                 if (NStr::strncasecmp(sss, *name, namelen) == 0) {
                //                     sss += namelen;
                //                     SET_MONTH(i + 1);
                //                     break;
                //                 }
                //                 ++name;
                //             }
                //             is_month_present = true;
                //             continue;
                //         }
                // ```
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
                // ```c++
                //         if (*fff == 'p'  ||  *fff == 'P') {
                //             if (NStr::strncasecmp(sss, "AM", 2) == 0) {
                //                 hour_format = eAM;
                //                 sss += 2;
                //             } else if (NStr::strncasecmp(sss, "PM", 2) == 0) {
                //                 hour_format = ePM;
                //                 sss += 2;
                //             }
                //             continue;
                //         }
                // ```
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
                // ```c++
                //         size_t len = 2;
                //         switch (*fff) {
                //             case 'Y': len = 4; break;
                // ...
                //         case 'Y':
                //             CHECK_RANGE2_YEAR(value, err_action);
                //             SET_YEAR(value);
                //             is_year_present = true;
                //             break;
                // ```
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:166,177
                // ```c++
                // #define CHECK_RANGE2_YEAR(value, err)   CHECK_RANGE2(value, "Year", 1583, kMax_Int, err)
                // #define SET_YEAR(value)  m_Data.year  = static_cast<unsigned int>((value) & 0xFFF)
                // ```
                b'Y' => {
                    let value = read_number(s, &mut sss, 4)?;
                    check_range(value, "Year", 1583, i64::from(i32::MAX))?;
                    time.year = (value & 0xFFF) as u32;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:770-775
                // ```c++
                //         case 'D':
                //         case 'd':
                //             CHECK_RANGE2_DAY(value, err_action);
                //             SET_DAY(value);
                //             is_day_present = true;
                //             break;
                // ```
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:168,179
                // ```c++
                // #define CHECK_RANGE2_DAY(value, err)    CHECK_RANGE2(value, "Day", 1, 31, err)
                // #define SET_DAY(value)   m_Data.day   = static_cast<unsigned char>((value) & 0x1F)
                // ```
                b'd' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Day", 1, 31)?;
                    time.day = (value & 0x1F) as u32;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:781-786
                // ```c++
                //         case 'H':
                //             CHECK_RANGE2_HOUR12(value, err_action);
                //             SET_HOUR(value % 12);
                //             is_12hour = true;
                //             is_time_present = true;
                //             break;
                // ```
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:170,180
                // ```c++
                // #define CHECK_RANGE2_HOUR12(value, err) CHECK_RANGE2(value, "Hour", 0, 12, err)
                // #define SET_HOUR(value)  m_Data.hour  = static_cast<unsigned char>((value) & 0x1F)
                // ```
                b'H' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Hour", 0, 12)?;
                    time.hour = ((value % 12) & 0x1F) as u32;
                    is_12hour = true;
                }

                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:787-791
                // ```c++
                //         case 'm':
                //             CHECK_RANGE2_MIN(value, err_action);
                //             SET_MIN(value);
                //             is_time_present = true;
                //             break;
                // ```
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:171,181
                // ```c++
                // #define CHECK_RANGE2_MIN(value, err)    CHECK_RANGE2(value, "Minute", 0, 59, err)
                // #define SET_MIN(value)   m_Data.min   = static_cast<unsigned char>((value) & 0x3F)
                // ```
                b'm' => {
                    let value = read_number(s, &mut sss, 2)?;
                    check_range(value, "Minute", 0, 59)?;
                    time.min = (value & 0x3F) as u32;
                }

                // `,` and `:`, the characters of the format that are neither spaces nor in
                // kFormatSymbolsTime.
                // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:129,604-611
                // ```c++
                // static const char* kFormatSymbolsTime = "yYMbBDdhHmsSlrgGzZwWpPo";
                // ...
                //             if (strchr(kFormatSymbolsTime, *fff) == 0) {
                //                 // Match non-format symbols
                //                 if (*fff == *sss) {
                //                     ++sss;
                //                     continue;
                //                 }
                //                 break;  // error: non-matching symbols
                //             }
                // ```
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
    // ```c++
    //     if (is_12hour  &&  hour_format == ePM) {
    //         SET_HOUR(m_Data.hour + 12);
    //     }
    // ```
    if is_12hour && hour_format == HourFormat::Pm {
        time.hour = (time.hour + 12) & 0x1F;
    }

    // NCBI skips trailing white space here only `if (isspace(*sss) && is_format_space)`.
    // That condition is never true for this format: the format ends in the symbol `P`, so a
    // completed loop leaves `is_format_space` false, and no `break` above is taken with both
    // a white-space character at `sss` and `is_format_space` true. Trailing white space is
    // therefore "too long".
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:857-859,863-872
    // ```c++
    //         if (isspace((unsigned char)(*sss))  &&  is_format_space) {
    //             SKIP_SPACES(sss);
    //         }
    // ...
    //     if (*fff != '\0'  &&
    //         !(format.GetFlags() & CTimeFormat::fMatch_ShortTime)) {
    //         X_INIT_ERROR(eFormat, "Time string '" + str +
    //                               "' is too short for time format '" + fmt + "'");
    //     }
    //     if (*sss != '\0'  &&
    //         !(format.GetFlags() & CTimeFormat::fMatch_ShortFormat)) {
    //         X_INIT_ERROR(eFormat, "Time string '" + str +
    //                               "' is too long for time format '" + fmt + "'");
    //     }
    // ```
    if !format_matched {
        return Err(ParseFailure::TimeException(format!(
            "Time string '{text}' is too short for time format '{fmt}'"
        )));
    }
    if byte_at(s, sss) != 0 {
        return Err(ParseFailure::TimeException(format!(
            "Time string '{text}' is too long for time format '{fmt}'"
        )));
    }

    // A completed loop has set year, month, day and time, so `ptcache` is 2222 and no field
    // is taken from the current time.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:880-884
    // ```c++
    //     int ptcache = 0;
    //     ptcache += (is_year_present  ? 2000 : 1000);
    //     ptcache += (is_month_present ? 200 : 100);
    //     ptcache += (is_day_present   ? 20 : 10);
    //     ptcache += (is_time_present  ? 2 : 1);
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbitime.cpp:932-934
    // ```c++
    //     if ( !IsValid() ) {
    //         X_INIT_ERROR(eConvert, "Unable to convert string '" + str + "' to CTime");
    //     }
    // ```
    if !time.is_valid() {
        return Err(ParseFailure::TimeException(format!(
            "Unable to convert string '{text}' to CTime"
        )));
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
// ```c++
//     string date;
//
//     for(int i = 0; i < m_VolSet.GetNumVols(); i++) {
//         string d = x_FixString( m_VolSet.GetVol(i)->GetDate() );
//
//         if (date.empty()) {
//             date = d;
//         } else if (d != date) {
//             try {
//                 CTime t1(date, CSeqDB::kBlastDbDateFormat);
//                 CTime t2(d, CSeqDB::kBlastDbDateFormat);
//
//                 if (t2 > t1) {
//                     date.swap(d);
//                 }
//             }
//             catch(CStringException &) {
//                 // Here I think it is better to pick any valid date
//                 // than to propagate a string exception.
//             }
//         }
//     }
// ...
//     return date;
// ```
pub fn database_date(volume_dates: &[&str]) -> Result<String, DateError> {
    let mut date = String::new();

    for &d in volume_dates {
        if date.is_empty() {
            date = d.to_owned();
        } else if d != date {
            let t1 = match parse_blast_db_date(&date) {
                Ok(time) => time,
                Err(ParseFailure::StringException) => continue,
                Err(ParseFailure::TimeException(reason)) => {
                    return Err(DateError { date, reason });
                }
            };
            let t2 = match parse_blast_db_date(d) {
                Ok(time) => time,
                Err(ParseFailure::StringException) => continue,
                Err(ParseFailure::TimeException(reason)) => {
                    return Err(DateError {
                        date: d.to_owned(),
                        reason,
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

#[cfg(test)]
mod tests {
    use super::*;

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
            Err(ParseFailure::TimeException(reason)) => reason,
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
        // SET_HOUR(value % 12) stores 0 for 12, and AM adds nothing.
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
        // SET_HOUR(value % 12) stores 0 for 12, and PM adds 12.
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
        // t2: NStr::StringToLong("") throws CStringException, caught by GetDate.
        for bad in [
            "Oct xx, 2026  3:22 PM", // day
            "Oct 12, year  3:22 PM", // year
            "Oct 12, 2026  x:22 PM", // hour
            "Oct 12, 2026  3:xx PM", // minute
            "Oct 12, 2026  3:",      // minute at the end of the string
            "Oct ",                  // day at the end of the string
        ] {
            assert_eq!(
                parse_blast_db_date(bad),
                Err(ParseFailure::StringException),
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
        assert_eq!(
            err,
            DateError {
                date: "not a date".to_owned(),
                reason: too_short("not a date"),
            }
        );
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
        assert_eq!(
            parse_blast_db_date("Oct \u{ff11}\u{ff12}, 2026  3:22 PM"),
            Err(ParseFailure::StringException)
        );
    }
}

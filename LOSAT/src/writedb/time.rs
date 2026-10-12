//! The local time of a build (`CTime(CTime::eCurrent)`) and its two text
//! forms: the creation date of an index file and the time of the log.

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
    /// The current local time, or `None` where the platform has no clock or
    /// no local time (`wasm32-unknown-unknown`).
    pub fn now() -> Option<Self> {
        now_impl()
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
// `CTime(CTime::eCurrent)` is the local time of the C library. The standard
// library of Rust has no local time, so the C library is called directly
// (`localtime_r`; `localtime_s` on Windows).
#[cfg(any(unix, target_os = "wasi"))]
fn now_impl() -> Option<LocalTime> {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_secs() as libc::time_t;
    // SAFETY: `tm` is plain data for which all-zero bytes are a valid value.
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    // SAFETY: both pointers refer to live, properly aligned values owned by
    // this frame; `localtime_r` writes only into `tm`.
    let result = unsafe { libc::localtime_r(&seconds, &mut tm) };
    if result.is_null() {
        return None;
    }
    Some(from_tm(&tm))
}

#[cfg(windows)]
fn now_impl() -> Option<LocalTime> {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_secs() as libc::time_t;
    // SAFETY: `tm` is plain data for which all-zero bytes are a valid value.
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    // SAFETY: both pointers refer to live, properly aligned values owned by
    // this frame; `localtime_s` writes only into `tm`.
    let status = unsafe { libc::localtime_s(&mut tm, &seconds) };
    if status != 0 {
        return None;
    }
    Some(from_tm(&tm))
}

#[cfg(not(any(unix, windows, target_os = "wasi")))]
fn now_impl() -> Option<LocalTime> {
    None
}

#[cfg(any(unix, windows, target_os = "wasi"))]
fn from_tm(tm: &libc::tm) -> LocalTime {
    LocalTime {
        year: tm.tm_year + 1900,
        month: (tm.tm_mon + 1) as u32,
        day: tm.tm_mday as u32,
        hour: tm.tm_hour as u32,
        minute: tm.tm_min as u32,
        second: tm.tm_sec as u32,
    }
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

    #[test]
    fn the_clock_gives_a_plausible_time() {
        let now = LocalTime::now().expect("a clock");
        assert!(now.year >= 2026 && (1..=12).contains(&now.month) && (1..=31).contains(&now.day));
        assert!(now.hour < 24 && now.minute < 60 && now.second < 62);
    }
}

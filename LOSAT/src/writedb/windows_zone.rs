// Code of Microsoft's C runtime (UCRT), from the sources that Microsoft publishes under the
// MIT License in the NuGet package Microsoft.Windows.SDK.CRTSource 10.0.22621.3
// (ucrt/...; also in THIRD_PARTY_NOTICES, section 6 (UCRT)). The package's LICENSE:
//
// MIT License
//
// Copyright (c) Microsoft Corporation.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE

//! Local UTC offset on Windows, as Microsoft's C runtime (UCRT) gives it to NCBI's Windows
//! build (`localtime`, see `super::local_offset`): `TZ` in UCRT's own form
//! (`tzn[+|-]hh[:mm[:ss]][dzn]`, with the United States DST rules), else the system zone of
//! `GetTimeZoneInformation` with its DST rules. The one entry is `windows_local_utc_offset`,
//! which `super::local_offset::local_utc_offset` calls on Windows. Every target compiles this
//! module, and its unit tests run on all of them.

// The code below follows the sources of UCRT that Microsoft publishes under the MIT license
// in the NuGet package Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (`ucrt/...` below).
// NCBI's Windows build calls UCRT's `localtime` (see `super::local_offset::local_utc_offset`).
// The global state that `_tzset` leaves is a `UcrtZone`; `localtime` of one instant is
// `UcrtZone::offset_at`.
//
// UCRT computes the zone in `int` and `long` (32 bits). Where an absurd `TZ` value (more than
// 596,523 hours) or a corrupt system zone would overflow them, the values here are the
// mathematical ones (the C++ overflow is undefined; LOSAT does not reproduce overflows, as for
// NCBI's own), and `localtime`'s range checks apply to them: past its range the offset is
// `None`.

/// The time zone state of UCRT after `_tzset`: `_timezone` (seconds west of UTC of standard
/// time), `_daylight`, `_dstbias` (seconds added to `_timezone` in daylight time) and where
/// `_isindst` takes the DST changes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct UcrtZone {
    timezone: i64,
    daylight: bool,
    dstbias: i64,
    /// `[DaylightDate, StandardDate]` of the system zone; `None` for the United States rules.
    changes: Option<[WinChange; 2]>,
}

/// A `SYSTEMTIME` of `TIME_ZONE_INFORMATION`: a DST change as `_isindst` reads it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct WinChange {
    year: u16,
    month: u16,
    day_of_week: u16,
    day: u16,
    hour: u16,
    minute: u16,
    second: u16,
    milliseconds: u16,
}

/// The fields of `TIME_ZONE_INFORMATION` that `_tzset` reads (biases in minutes west of UTC).
#[derive(Clone, Copy, Debug, Default)]
struct WinZone {
    bias: i32,
    standard_date: WinChange,
    standard_bias: i32,
    daylight_date: WinChange,
    daylight_bias: i32,
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/timeset.cpp:67-69
/// The start values of the time zone globals of UCRT.
const UCRT_START: UcrtZone = UcrtZone {
    timezone: 8 * 3600,
    daylight: true,
    dstbias: -3600,
    changes: None,
};

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:145-200
/// `_tzset` without `TZ`: the zone of `GetTimeZoneInformation` (`None`: the call failed, and
/// the globals keep their start values with the United States rules). The first `_tzset` of
/// a process is the one `localtime` makes, so the values before it are the start values.
#[cfg_attr(not(windows), allow(dead_code))]
fn ucrt_zone_from_system(info: Option<&WinZone>) -> UcrtZone {
    let Some(z) = info else {
        return UCRT_START;
    };
    let mut timezone = i64::from(z.bias) * 60;
    if z.standard_date.month != 0 {
        timezone += i64::from(z.standard_bias) * 60;
    }
    let (daylight, dstbias) = if z.daylight_date.month != 0 && z.daylight_bias != 0 {
        (
            true,
            (i64::from(z.daylight_bias) - i64::from(z.standard_bias)) * 60,
        )
    } else {
        (false, 0)
    };
    UcrtZone {
        timezone,
        daylight,
        dstbias,
        changes: Some([z.daylight_date, z.standard_date]),
    }
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:264-335
/// `_tzset` with a non-empty `TZ` (`tz`, UTF-16 as `_wgetenv_s` gives it): any first three
/// characters are the zone's name, an offset in hours, minutes and seconds west of UTC
/// follows, and anything after it makes the zone use the United States DST rules. Each number
/// is read by `wcstol` (white space, a sign, digits of any script), and the scan that skips it
/// takes ASCII digits only (and `+` after the hours). `_dstbias` is not set here: it keeps its
/// start value. `space_above_latin1` is `iswspace` of the characters above U+00FF.
#[cfg_attr(not(windows), allow(dead_code))]
fn ucrt_zone_from_tz(tz: &[u16], space_above_latin1: &dyn Fn(u16) -> bool) -> UcrtZone {
    // The value is a C string: it ends at its first NUL.
    let tz = &tz[..tz.iter().position(|&c| c == 0).unwrap_or(tz.len())];
    let at = |i: usize| tz.get(i).copied().unwrap_or(0);
    let is_digit = |c: u16| (u16::from(b'0')..=u16::from(b'9')).contains(&c);
    let number = |from: usize| i64::from(ucrt_wcstol(&tz[from..], space_above_latin1));
    let mut p = tz.len().min(3);
    let negative = at(p) == u16::from(b'-');
    p += usize::from(negative);
    let mut timezone = number(p) * 3600;
    while at(p) == u16::from(b'+') || is_digit(at(p)) {
        p += 1;
    }
    if at(p) == u16::from(b':') {
        p += 1;
        timezone += number(p) * 60;
        while is_digit(at(p)) {
            p += 1;
        }
        if at(p) == u16::from(b':') {
            p += 1;
            timezone += number(p);
            while is_digit(at(p)) {
                p += 1;
            }
        }
    }
    if negative {
        timezone = -timezone;
    }
    UcrtZone {
        timezone,
        daylight: at(p) != 0,
        dstbias: UCRT_START.dstbias,
        changes: None,
    }
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/inc/corecrt_internal_strtox.h:248-345
/// `wcstol(s, &end, 10)` of UCRT (`ucrt/convert/strtox.cpp:237-243`): leading white space, an
/// optional sign, then decimal digits; a value outside `long` (32 bits) is `LONG_MIN` or
/// `LONG_MAX` by its sign; no digit is 0.
fn ucrt_wcstol(s: &[u16], space_above_latin1: &dyn Fn(u16) -> bool) -> i32 {
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut i = 0;
    while ucrt_iswspace(at(i), space_above_latin1) {
        i += 1;
    }
    let negative = at(i) == u16::from(b'-');
    if negative || at(i) == u16::from(b'+') {
        i += 1;
    }
    let mut number: Option<u64> = None;
    while let Some(digit) = ucrt_wide_digit(at(i)) {
        // Past 2^33 every value is out of range; the bound keeps the sum from growing.
        number = Some((number.unwrap_or(0) * 10 + u64::from(digit)).min(1 << 33));
        i += 1;
    }
    let Some(number) = number else {
        return 0;
    };
    let value = if negative {
        -(number as i64)
    } else {
        number as i64
    };
    value.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/inc/corecrt_internal_strtox.h:88-158
/// The decimal digit of `c` for `wcstol` in base 10 (letters are digits 10 and up, and so end
/// a decimal number).
fn ucrt_wide_digit(c: u16) -> Option<u32> {
    const ZEROS: [u16; 16] = [
        0x0660, 0x06F0, 0x0966, 0x09E6, 0x0A66, 0x0AE6, 0x0B66, 0x0C66, 0x0CE6, 0x0D66, 0x0E50,
        0x0ED0, 0x0F20, 0x1040, 0x17E0, 0x1810,
    ];
    let digit = |zero: u16| (zero..zero + 10).contains(&c).then(|| u32::from(c - zero));
    if c < 0x0030 {
        return None;
    }
    if c < 0x0030 + 10 {
        return digit(0x0030);
    }
    if c < 0xFF10 {
        // The ranges are in increasing order: the first zero above `c` ends the search.
        for zero in ZEROS {
            if c < zero {
                return None;
            }
            if let Some(d) = digit(zero) {
                return Some(d);
            }
        }
        return None;
    }
    digit(0xFF10)
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/convert/iswctype.cpp:42-54
// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/locale/ctype.cpp:53-57, 76, 177, 204
// (`_wctype`, the entries with `_SPACE`)
/// `iswspace(c)` of UCRT (`is_space` of `wcstol`, `corecrt_internal_strtox.h:185-188`):
/// UCRT's table up to U+00FF, `space_above_latin1` (the system's character type) above it,
/// and never `WEOF` (U+FFFF).
fn ucrt_iswspace(c: u16, space_above_latin1: &dyn Fn(u16) -> bool) -> bool {
    match c {
        0x09..=0x0D | 0x20 | 0x85 | 0xA0 => true,
        0x00..=0xFF | 0xFFFF => false,
        _ => space_above_latin1(c),
    }
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/days.cpp:14-24
const UCRT_DAYS: [i32; 13] = [-1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364];
const UCRT_LPDAYS: [i32; 13] = [-1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/inc/corecrt_internal_time.h:153-175
/// Years count from 1900 (`tm_year`); division truncates toward zero, as in C.
fn ucrt_is_leap_year(yr: i64) -> bool {
    (yr % 4 == 0 && yr % 100 != 0) || (yr + 1900) % 400 == 0
}

fn ucrt_elapsed_leap_years(yr: i64) -> i64 {
    (yr - 1) / 4 - (yr - 1) / 100 + (yr + 299) / 400 - 17
}

const UCRT_DAY_SEC: i64 = 24 * 60 * 60;
const UCRT_YEAR_SEC: i64 = 365 * UCRT_DAY_SEC;
const UCRT_MS_PER_DAY: i64 = UCRT_DAY_SEC * 1000;

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/inc/corecrt_internal_time.h:37-40, 98
const UCRT_MAX_LOCAL_TIME: i64 = 14 * 60 * 60;
const UCRT_MIN_LOCAL_TIME: i64 = -12 * 60 * 60;
const UCRT_MAX_TIME64: i64 = 0x7_9358_2aff;

/// The fields of a `tm` that `_isindst` reads: `tm_year` (from 1900), `tm_yday` (from 0) and
/// the time of day.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct UcrtTm {
    year: i64,
    yday: i64,
    hour: i64,
    min: i64,
    sec: i64,
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/gmtime.cpp:54-77, 97-134
/// `_gmtime64_s` of `caltim`; `None` outside its range. Within it the calendar is the
/// Gregorian one of `LocalTime::from_unix_seconds`.
fn ucrt_gmtime(caltim: i64) -> Option<UcrtTm> {
    if !(UCRT_MIN_LOCAL_TIME..=UCRT_MAX_TIME64 + UCRT_MAX_LOCAL_TIME).contains(&caltim) {
        return None;
    }
    let mut c = caltim;
    let mut year = c / UCRT_YEAR_SEC + 70;
    c -= (year - 70) * UCRT_YEAR_SEC;
    c -= ucrt_elapsed_leap_years(year) * UCRT_DAY_SEC;
    if c < 0 {
        c += UCRT_YEAR_SEC;
        year -= 1;
        if ucrt_is_leap_year(year) {
            c += UCRT_DAY_SEC;
        }
    }
    let yday = c / UCRT_DAY_SEC;
    c -= yday * UCRT_DAY_SEC;
    let hour = c / 3600;
    c -= hour * 3600;
    let min = c / 60;
    Some(UcrtTm {
        year,
        yday,
        hour,
        min,
        sec: c - min * 60,
    })
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:446-496
// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:546-559
// (the choice of the form)
/// `cvtdate`: the day of the year (from 0) and the millisecond of the day of the change `c`
/// in year `yr` (from 1900): with `wYear` 0, the `wDay`-th (5: the last) weekday `wDayOfWeek`
/// of month `wMonth`; otherwise day `wDay` of month `wMonth`. `None` for a month whose day
/// table entry lies outside UCRT's tables (the C++ read is undefined).
fn ucrt_cvtdate(c: &WinChange, yr: i64) -> Option<(i64, i64)> {
    let days = if ucrt_is_leap_year(yr) {
        &UCRT_LPDAYS
    } else {
        &UCRT_DAYS
    };
    let month = usize::from(c.month);
    let yearday = if c.year == 0 {
        if !(1..=12).contains(&month) {
            return None;
        }
        let (week, dayofweek) = (i64::from(c.day), i64::from(c.day_of_week));
        let mut yearday = 1 + i64::from(days[month - 1]);
        // `_BASE_DOW`: 1970-01-01 was a Thursday (`corecrt_internal_time.h:34`).
        let monthdow = (yearday + (yr - 70) * 365 + ucrt_elapsed_leap_years(yr) + 4) % 7;
        yearday += if monthdow <= dayofweek {
            (dayofweek - monthdow) + (week - 1) * 7
        } else {
            (dayofweek - monthdow) + week * 7
        };
        if week == 5 && yearday > i64::from(days[month]) {
            yearday -= 7;
        }
        yearday
    } else {
        if !(1..=13).contains(&month) {
            return None;
        }
        i64::from(days[month - 1]) + i64::from(c.day)
    };
    let ms = i64::from(c.milliseconds)
        + 1000 * (i64::from(c.second) + 60 * (i64::from(c.minute) + 60 * i64::from(c.hour)));
    Some((yearday, ms))
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:611-651
/// The United States DST changes of year `yr` (from 1900) in the form of `GetTimeZoneInformation`.
fn ucrt_us_changes(yr: i64) -> [WinChange; 2] {
    let ((start_month, start_week), (end_month, end_week)) = if 107 > yr {
        ((4, 1), (10, 5))
    } else {
        ((3, 2), (11, 1))
    };
    let change = |month, week| WinChange {
        month,
        day: week,
        hour: 2,
        ..WinChange::default()
    };
    [change(start_month, start_week), change(end_month, end_week)]
}

impl UcrtZone {
    // Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:498-512
    // (cvtdate, end of DST)
    // Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:655-685
    /// `_isindst` of the local standard time `tm`. The changes are those of `tm`'s year.
    fn is_in_dst(&self, tm: &UcrtTm) -> Option<bool> {
        let [start, end] = self.changes.unwrap_or_else(|| ucrt_us_changes(tm.year));
        let (start_yd, start_ms) = ucrt_cvtdate(&start, tm.year)?;
        let (mut end_yd, mut end_ms) = ucrt_cvtdate(&end, tm.year)?;
        end_ms += self.dstbias * 1000;
        if end_ms < 0 {
            end_ms += UCRT_MS_PER_DAY;
            end_yd -= 1;
        } else if end_ms >= UCRT_MS_PER_DAY {
            end_ms -= UCRT_MS_PER_DAY;
            end_yd += 1;
        }
        if start_yd < end_yd {
            if tm.yday < start_yd || tm.yday > end_yd {
                return Some(false);
            }
            if tm.yday > start_yd && tm.yday < end_yd {
                return Some(true);
            }
        } else {
            if tm.yday < end_yd || tm.yday > start_yd {
                return Some(true);
            }
            if tm.yday > end_yd && tm.yday < start_yd {
                return Some(false);
            }
        }
        let ms = 1000 * (tm.sec + 60 * tm.min + 3600 * tm.hour);
        Some(if tm.yday == start_yd {
            ms >= start_ms
        } else {
            ms < end_ms
        })
    }

    // Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/localtime.cpp:46-80
    /// `_localtime64_s` of `t` as an offset: the local time is the calendar time of `t` plus
    /// the offset. `None` where `localtime` fails (NCBI's `CTime` then throws). The first and
    /// last three days of `localtime`'s range (to 1970-01-04 00:00:00 and from 3001-01-16
    /// 07:59:59 UTC on) take another branch of `localtime.cpp`, which is not ported: the current
    /// time is never there, and `None` (UTC) claims no UCRT version's behaviour there.
    #[cfg_attr(not(windows), allow(dead_code))]
    fn offset_at(&self, t: i64) -> Option<i64> {
        if !(t > 3 * UCRT_DAY_SEC && t < UCRT_MAX_TIME64 - 3 * UCRT_DAY_SEC) {
            return None;
        }
        let mut ltime = t - self.timezone;
        let tm = ucrt_gmtime(ltime)?;
        if self.daylight && self.is_in_dst(&tm)? {
            ltime -= self.dstbias;
            ucrt_gmtime(ltime)?;
        }
        Some(ltime - t)
    }
}

// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/tzset.cpp:347-363
// `get_tz_environment_variable` reads `TZ` with `_wgetenv_s` (`tzset.cpp:76-103`), from the
// CRT's wide environment, which UCRT makes from the process environment block
// (`ucrt/env/environment_initialization.cpp:166-188`); the name is compared without case
// (`ucrt/env/getenv.cpp:44`). `std::env::var_os` reads the same block, without case, and keeps
// its UTF-16 units. `localtime` runs `_tzset` once per process, at its first call (`__tzset`,
// `tzset.cpp:391-416`); LOSAT asks at each time it prints, with the same environment and zone.
/// The offset of UCRT's `localtime` at `unix_seconds` from `TZ` and the system zone.
#[cfg(windows)]
pub(super) fn windows_local_utc_offset(unix_seconds: i64) -> Option<i64> {
    use std::os::windows::ffi::OsStrExt;
    let tz: Vec<u16> = std::env::var_os("TZ")
        .map(|v| v.encode_wide().collect())
        .unwrap_or_default();
    let zone = if tz.first().is_some_and(|&c| c != 0) {
        ucrt_zone_from_tz(&tz, &windows_char_is_space)
    } else {
        ucrt_zone_from_system(windows_time_zone().as_ref())
    };
    zone.offset_at(unix_seconds)
}

/// `GetTimeZoneInformation`, the call of `tzset_from_system_nolock`; `None` when it fails.
#[cfg(windows)]
fn windows_time_zone() -> Option<WinZone> {
    use windows_sys::Win32::Foundation::SYSTEMTIME;
    use windows_sys::Win32::System::Time::{
        GetTimeZoneInformation, TIME_ZONE_ID_INVALID, TIME_ZONE_INFORMATION,
    };
    let mut info = TIME_ZONE_INFORMATION::default();
    // SAFETY: `info` is a TIME_ZONE_INFORMATION that lives for the call, which fills it.
    if unsafe { GetTimeZoneInformation(&mut info) } == TIME_ZONE_ID_INVALID {
        return None;
    }
    let change = |s: SYSTEMTIME| WinChange {
        year: s.wYear,
        month: s.wMonth,
        day_of_week: s.wDayOfWeek,
        day: s.wDay,
        hour: s.wHour,
        minute: s.wMinute,
        second: s.wSecond,
        milliseconds: s.wMilliseconds,
    };
    Some(WinZone {
        bias: info.Bias,
        standard_date: change(info.StandardDate),
        standard_bias: info.StandardBias,
        daylight_date: change(info.DaylightDate),
        daylight_bias: info.DaylightBias,
    })
}

// `__acrt_GetStringTypeW` is declared with `GetStringTypeW`'s parameters
// (`ucrt/inc/corecrt_internal.h:570-575`); its body is not among the shipped sources, and
// `GetStringTypeW` is the system call of that signature.
/// `iswspace` above U+00FF: `CT_CTYPE1` of the character has `C1_SPACE` (`iswctype.cpp:48-54`).
#[cfg(windows)]
fn windows_char_is_space(c: u16) -> bool {
    use windows_sys::Win32::Globalization::{GetStringTypeW, C1_SPACE, CT_CTYPE1};
    let mut char_type: u16 = 0;
    // SAFETY: the source is the one UTF-16 unit `c` (`cchSrc` 1), which lives for the call,
    // and `char_type` takes the one word the call writes for it.
    let ok = unsafe { GetStringTypeW(CT_CTYPE1, &c, 1, &mut char_type) };
    ok != 0 && u32::from(char_type) & C1_SPACE != 0
}

#[cfg(test)]
mod tests {
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
        let change =
            |(month, day, day_of_week, hour, minute): (u16, u16, u16, u16, u16)| WinChange {
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

    /// Windows zones and the POSIX rule strings of the same changes, which the glibc reader
    /// above takes (it matches glibc, `glibc_tsv`).
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

    // UCRT places the changes of the year of the local standard time, glibc those of the UTC
    // year: they meet for zones whose changes are not at the turn of the year.
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
}

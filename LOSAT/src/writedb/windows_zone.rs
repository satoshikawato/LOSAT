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
/// `wcstol` of `s`, an end pointer and base 10, of UCRT (`ucrt/convert/strtox.cpp:237-243`):
/// leading white space, an optional sign, then decimal digits; a value outside `long` (32 bits) is
/// `LONG_MIN` or `LONG_MAX` by its sign; no digit is 0.
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
/// `iswspace` of `c`, of UCRT (`is_space` of `wcstol`, `corecrt_internal_strtox.h:185-188`):
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
#[path = "windows_zone_tests.rs"]
mod tests;

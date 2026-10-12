//! Local UTC offset of the system time zone, found with the standard library only.
//!
//! The zone is chosen the way the C library does it: the `TZ` environment variable
//! (POSIX, "Other Environment Variables"), else `/etc/localtime`. A zone is either a
//! TZif file (RFC 9636) or a POSIX TZ rule string. Leap-second records of a TZif file
//! (the `right/` zones) are ignored: they never change the UTC offset of a type.
//!
//! The rule-string reader reproduces the C library's reader, including how it treats
//! malformed strings and its use of the UTC year of the instant when it places the
//! DST transitions, so that both give the same answer for every input.
#![cfg_attr(not(unix), allow(dead_code))]

use std::ffi::OsStr;
use std::path::Path;

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
/// Seconds east of UTC of the system's local time zone at `unix_seconds`, as the C
/// library's `localtime_r` gives it on Unix. `None` when the zone cannot be determined
/// (the caller then uses UTC) and on platforms that are not Unix.
pub fn local_utc_offset(unix_seconds: i64) -> Option<i32> {
    #[cfg(unix)]
    {
        let read = |p: &Path| std::fs::read(p).ok();
        resolve(
            std::env::var_os("TZ").as_deref(),
            std::env::var_os("TZDIR").as_deref(),
            &read,
            unix_seconds,
        )
    }
    #[cfg(not(unix))]
    {
        let _ = unix_seconds;
        None
    }
}

const DEFAULT_TZDIR: &str = "/usr/share/zoneinfo";

/// The zone the C library picks for `TZ` (`tz`; `None` when unset) and `TZDIR`, then the
/// offset at `unix_seconds`. `read` returns the bytes of a file or `None` when it cannot
/// be read. `None` only when the value is not valid UTF-8 or names an undefined rule.
fn resolve(
    tz: Option<&OsStr>,
    tzdir: Option<&OsStr>,
    read: &dyn Fn(&Path) -> Option<Vec<u8>>,
    unix_seconds: i64,
) -> Option<i32> {
    let dir = match tzdir {
        Some(d) if !d.is_empty() => Path::new(d),
        _ => Path::new(DEFAULT_TZDIR),
    };
    let name = match tz {
        None => "/etc/localtime",
        Some(v) => {
            let v = v.to_str()?;
            // An empty value means UTC; one leading ':' is dropped (implementation-defined syntax).
            let v = v.strip_prefix(':').unwrap_or(v);
            if v.is_empty() {
                return Some(0);
            }
            v
        }
    };
    let path = if name.starts_with('/') {
        Path::new(name).to_path_buf()
    } else {
        dir.join(name)
    };
    // A readable TZif file wins; any other value is read as a rule string, which gives
    // UTC when it is not one (the unset case `/etc/localtime` has no rule-string reading).
    if let Some(z) = read(&path).as_deref().and_then(parse_tzif) {
        return z.offset_at(unix_seconds);
    }
    if tz.is_none() {
        return Some(0);
    }
    let p = parse_posix(name.as_bytes());
    if p.no_rule {
        // A DST name without rules takes the transitions of `posixrules` with the
        // offsets of the string.
        let rules = read(&dir.join("posixrules"));
        if let Some(z) = rules.as_deref().and_then(parse_tzif) {
            if let Some(o) = z.offset_with_rules(&p, unix_seconds) {
                return Some(o);
            }
        }
    }
    p.offset_at(unix_seconds)
}

// ---------------------------------------------------------------- TZif (RFC 9636)

/// A TZif file with its counts checked; slices point into the file bytes.
struct Tzif<'a> {
    /// Transition times, `tsize` bytes each (4 in version 1, 8 otherwise).
    times: &'a [u8],
    tsize: usize,
    /// Local time type index of each transition.
    idx: &'a [u8],
    /// Local time types: utoff (4 bytes), isdst, abbreviation index.
    types: &'a [u8],
    /// The POSIX TZ string of the footer (version 2 and later).
    footer: Option<&'a [u8]>,
}

fn be_u32(d: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_be_bytes(
        d.get(at..at.checked_add(4)?)?.try_into().ok()?,
    ))
}

fn parse_tzif(d: &[u8]) -> Option<Tzif<'_>> {
    // Header: magic, version, 15 reserved bytes, six counts (RFC 9636 section 3.1).
    let header = |at: usize| -> Option<([u64; 6], usize)> {
        if d.get(at..at.checked_add(4)?)? != b"TZif" {
            return None;
        }
        let mut c = [0u64; 6];
        for (i, v) in c.iter_mut().enumerate() {
            *v = u64::from(be_u32(d, at + 20 + 4 * i)?);
        }
        Some((c, at + 44))
    };
    let (mut c, mut at) = header(0)?;
    let mut tsize = 4usize;
    let wide = d[4] >= b'2';
    // [isutcnt, isstdcnt, leapcnt, timecnt, typecnt, charcnt]
    let block = |c: &[u64; 6], tsize: u64| {
        c[3] * (tsize + 1) + c[4] * 6 + c[5] + c[2] * (tsize + 4) + c[1] + c[0]
    };
    if wide {
        at = at.checked_add(usize::try_from(block(&c, 4)).ok()?)?;
        (c, at) = header(at)?;
        tsize = 8;
    }
    let (n, k, chars) = (
        usize::try_from(c[3]).ok()?,
        usize::try_from(c[4]).ok()?,
        usize::try_from(c[5]).ok()?,
    );
    let end = at.checked_add(usize::try_from(block(&c, tsize as u64)).ok()?)?;
    if k == 0 || d.len() < end {
        return None;
    }
    let times = &d[at..at + n * tsize];
    let idx = &d[at + n * tsize..at + n * (tsize + 1)];
    let types = &d[at + n * (tsize + 1)..at + n * (tsize + 1) + k * 6];
    let footer = if wide {
        // Version 2 and later end with "\n" TZ string "\n"; a file without that is malformed.
        // An empty string means that no rule follows the last transition.
        let line = d.get(end..)?.strip_prefix(b"\n")?;
        let spec = &line[..line.iter().position(|&b| b == b'\n')?];
        Some(spec).filter(|s| !s.is_empty())
    } else {
        None
    };
    let z = Tzif {
        times,
        tsize,
        idx,
        types,
        footer,
    };
    let sorted = (1..n).all(|i| z.time(i - 1) <= z.time(i));
    let types_ok = types
        .chunks(6)
        .all(|t| t[4] <= 1 && usize::from(t[5]) < chars && be_u32(t, 0) != Some(0x8000_0000));
    (sorted && types_ok && idx.iter().all(|&i| usize::from(i) < k)).then_some(z)
}

/// The number of leading indices below `n` for which `f` holds (`f` is true, then false).
fn count_while(n: usize, f: impl Fn(usize) -> bool) -> usize {
    let (mut lo, mut hi) = (0, n);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if f(mid) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

impl Tzif<'_> {
    fn n(&self) -> usize {
        self.idx.len()
    }
    fn time(&self, i: usize) -> i64 {
        let b = &self.times[i * self.tsize..(i + 1) * self.tsize];
        if self.tsize == 4 {
            i64::from(i32::from_be_bytes([b[0], b[1], b[2], b[3]]))
        } else {
            i64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
        }
    }
    fn utoff(&self, ty: usize) -> i32 {
        let b = &self.types[ty * 6..];
        i32::from_be_bytes([b[0], b[1], b[2], b[3]])
    }
    fn is_dst(&self, ty: usize) -> bool {
        self.types[ty * 6 + 4] != 0
    }
    /// Type used before the first transition: the first standard-time type, else type 0.
    fn first_type(&self) -> usize {
        (0..self.types.len() / 6)
            .find(|&t| !self.is_dst(t))
            .unwrap_or(0)
    }

    /// The offset at `t`: the type of the last transition at or before `t`; before the
    /// first transition (or with none) the first standard-time type; from the last
    /// transition on, the footer rule when there is one, else the last type.
    fn offset_at(&self, t: i64) -> Option<i32> {
        let n = self.n();
        if n == 0 || t < self.time(0) {
            return Some(self.utoff(self.first_type()));
        }
        if t >= self.time(n - 1) {
            if let Some(f) = self.footer {
                return parse_posix(f).offset_at(t);
            }
        }
        // Number of transitions at or before `t`; at least 1 here.
        let lo = count_while(n, |i| self.time(i) <= t);
        Some(self.utoff(usize::from(self.idx[lo - 1])))
    }

    /// The `posixrules` reading of a rule string without rules. Only the flavor of each
    /// transition counts: standard time takes the offset of `p`, DST the DST offset. The C
    /// library moves the transitions by the offset change of the new flavor: a change to
    /// DST by the standard offset of `p` minus that of the last standard-time transition
    /// of the file, a change to standard time by the DST offset of `p`. Before the first
    /// transition the answer is standard time; from the last on the footer of this file
    /// decides with its own offsets.
    fn offset_with_rules(&self, p: &Posix, t: i64) -> Option<i32> {
        let n = self.n();
        if self.types.len() / 6 < 2 {
            return None;
        }
        let dst_at = |i: usize| self.is_dst(usize::from(self.idx[i]));
        let last_std = (0..n)
            .rev()
            .find(|&i| !dst_at(i))
            .map_or(0, |i| self.utoff(usize::from(self.idx[i])));
        let moved = |i: usize| {
            let shift = if dst_at(i) {
                i64::from(p.off[0]) - i64::from(last_std)
            } else {
                i64::from(p.off[1])
            };
            self.time(i) + shift
        };
        if n == 0 || t < moved(0) {
            return Some(p.off[0]);
        }
        if let (true, Some(f)) = (t >= moved(n - 1), self.footer) {
            return parse_posix(f).offset_at(t);
        }
        let lo = count_while(n, |i| moved(i) <= t);
        Some(p.off[usize::from(dst_at(lo - 1))])
    }
}

/// Offset of the TZif file `data` at `unix_seconds`; `None` when the file is malformed.
#[cfg(test)]
fn offset_from_tzif(data: &[u8], unix_seconds: i64) -> Option<i32> {
    parse_tzif(data)?.offset_at(unix_seconds)
}

// ------------------------------------------------- POSIX TZ strings (POSIX, TZ variable)

/// One DST change rule. `kind`: 0 = `n` (day 0..=365, leap days counted), 1 = `Jn`
/// (day 1..=365, leap day not counted), 2 = `Mm.n.d`.
#[derive(Clone, Copy, Default)]
struct Rule {
    kind: u8,
    m: u16,
    n: u16,
    d: u16,
    secs: i64,
}

/// A parsed rule string in the state the C library keeps it in: a part that does not
/// parse leaves its fields at zero (or at what was read before the failure) instead of
/// failing the whole string.
#[derive(Default)]
struct Posix {
    /// Seconds east of UTC: standard time, DST.
    off: [i32; 2],
    rule: [Rule; 2],
    /// A DST name without rules was given.
    no_rule: bool,
}

/// C `sscanf("%hu")`: skips white space, takes an optional sign, wraps to 16 bits.
fn scan_hu(s: &[u8], p: &mut usize) -> Option<u16> {
    let mut i = *p;
    while i < s.len() && (s[i] == b' ' || (b'\t'..=b'\r').contains(&s[i])) {
        i += 1;
    }
    let neg = s.get(i) == Some(&b'-');
    if neg || s.get(i) == Some(&b'+') {
        i += 1;
    }
    let start = i;
    let mut v: Option<u64> = Some(0);
    while let Some(&c) = s.get(i).filter(|c| c.is_ascii_digit()) {
        v = v.and_then(|v| v.checked_mul(10)?.checked_add(u64::from(c - b'0')));
        i += 1;
    }
    if i == start {
        return None;
    }
    *p = i;
    // An overflowing number reads as the largest value, whatever its sign.
    Some(match v {
        None => u64::MAX as u16,
        Some(v) if neg => v.wrapping_neg() as u16,
        Some(v) => v as u16,
    })
}

/// `hh[:mm[:ss]]` as `sscanf("%hu:%hu:%hu")` reads it; fields not read keep `init`.
fn scan_hms(s: &[u8], p: &mut usize, init: [u16; 3]) -> Option<[u16; 3]> {
    let mut v = init;
    let mut i = *p;
    v[0] = scan_hu(s, &mut i)?;
    *p = i;
    for k in 1..3 {
        if s.get(i) != Some(&b':') {
            break;
        }
        i += 1;
        match scan_hu(s, &mut i) {
            Some(x) => v[k] = x,
            None => break,
        }
        *p = i;
    }
    Some(v)
}

fn parse_name(s: &[u8], p: &mut usize) -> bool {
    let quoted = s.get(*p) == Some(&b'<');
    let start = *p + usize::from(quoted);
    let mut i = start;
    while let Some(&c) = s.get(i) {
        if c.is_ascii_alphabetic() || (quoted && (c.is_ascii_digit() || c == b'+' || c == b'-')) {
            i += 1;
        } else {
            break;
        }
    }
    if (quoted && s.get(i) != Some(&b'>')) || i - start < 3 {
        return false;
    }
    *p = i + usize::from(quoted);
    true
}

/// `[+-]hh[:mm[:ss]]`; POSIX counts west of UTC as positive. Hours are held to 24,
/// minutes and seconds to 59. A missing DST offset is one hour ahead of standard time.
fn parse_offset(s: &[u8], p: &mut usize, which: usize, st: &mut Posix) -> bool {
    let first = s.get(*p).copied();
    if which == 0 && !matches!(first, Some(b'+' | b'-' | b'0'..=b'9')) {
        return false;
    }
    let sign = match first {
        Some(b'-') => 1,
        _ => -1,
    };
    if matches!(first, Some(b'+' | b'-')) {
        *p += 1;
    }
    match scan_hms(s, p, [0; 3]) {
        Some([h, m, sec]) => {
            st.off[which] = sign
                * (i32::from(h.min(24)) * 3600
                    + i32::from(m.min(59)) * 60
                    + i32::from(sec.min(59)));
            true
        }
        None if which == 0 => {
            st.off = [0; 2];
            false
        }
        None => {
            st.off[1] = st.off[0] + 3600;
            true
        }
    }
}

fn parse_rule(s: &[u8], p: &mut usize, which: usize, st: &mut Posix) -> bool {
    let mut i = *p + usize::from(s.get(*p) == Some(&b','));
    let r = &mut st.rule[which];
    match s.get(i) {
        Some(&c) if c == b'J' || c.is_ascii_digit() => {
            r.kind = u8::from(c == b'J');
            i += usize::from(c == b'J');
            if !s.get(i).is_some_and(u8::is_ascii_digit) {
                return false;
            }
            let digits = s[i..].iter().take_while(|c| c.is_ascii_digit()).count();
            let d = s[i..i + digits].iter().try_fold(0u64, |a, &c| {
                a.checked_mul(10)?.checked_add(u64::from(c - b'0'))
            });
            match d {
                Some(d) if d <= 365 && !(r.kind == 1 && d == 0) => r.d = d as u16,
                _ => return false,
            }
            i += digits;
        }
        Some(b'M') => {
            r.kind = 2;
            i += 1;
            let mut ok = true;
            for k in 0..3 {
                if k > 0 {
                    if s.get(i) != Some(&b'.') {
                        ok = false;
                        break;
                    }
                    i += 1;
                }
                let Some(v) = scan_hu(s, &mut i) else {
                    ok = false;
                    break;
                };
                *[&mut r.m, &mut r.n, &mut r.d][k] = v;
            }
            if !ok || !(1..=12).contains(&r.m) || !(1..=5).contains(&r.n) || r.d > 6 {
                return false;
            }
        }
        // The United States rules since 2007: second Sunday of March to first Sunday of November.
        None => {
            r.kind = 2;
            (r.m, r.n, r.d) = if which == 0 { (3, 2, 0) } else { (11, 1, 0) };
        }
        Some(_) => return false,
    }
    match s.get(i) {
        None | Some(b',') => r.secs = 7200,
        Some(b'/') => {
            i += 1;
            if i >= s.len() {
                return false;
            }
            let neg = s[i] == b'-';
            i += usize::from(neg);
            let [h, m, sec] = scan_hms(s, &mut i, [2, 0, 0]).unwrap_or([2, 0, 0]);
            let v = i64::from(h) * 3600 + i64::from(m) * 60 + i64::from(sec);
            r.secs = if neg { -v } else { v };
        }
        Some(_) => return false,
    }
    *p = i;
    true
}

/// `std offset [dst [offset] [,start[/time],end[/time]]]`. Never fails: a string that
/// does not start with a valid standard zone gives UTC.
fn parse_posix(s: &[u8]) -> Posix {
    let mut st = Posix::default();
    let mut p = 0;
    if parse_name(s, &mut p) && parse_offset(s, &mut p, 0, &mut st) {
        if p < s.len() {
            if parse_name(s, &mut p) {
                parse_offset(s, &mut p, 1, &mut st);
                st.no_rule = p == s.len() || s[p..] == [b','];
            }
            if parse_rule(s, &mut p, 0, &mut st) {
                parse_rule(s, &mut p, 1, &mut st);
            }
        } else {
            st.off[1] = st.off[0];
        }
    }
    st
}

fn is_leap(y: i64) -> bool {
    y % 4 == 0 && (y % 100 != 0 || y % 400 == 0)
}

/// The UTC year of `t` (proleptic Gregorian calendar).
fn utc_year(t: i64) -> Option<i64> {
    if t.unsigned_abs() > 1 << 52 {
        return None;
    }
    let z = t.div_euclid(86400) + 719468;
    let doe = z.rem_euclid(146097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let month = (5 * doy + 2) / 153;
    Some(yoe + z.div_euclid(146097) * 400 + i64::from(month >= 10))
}

impl Rule {
    /// The UTC instant of the change in `year`, for a zone `off` seconds east of UTC.
    /// Years up to 1970 count from 1970-01-01. `None` for a month outside 1..=12,
    /// whose reading is undefined.
    fn change(&self, year: i64, off: i32) -> Option<i64> {
        let mut t = 0;
        if year > 1970 {
            t = ((year - 1970) * 365 + (year - 1 - 1968) / 4 - (year - 1 - 1900) / 100
                + (year - 1 - 1600) / 400)
                * 86400;
        }
        let d = i64::from(self.d);
        match self.kind {
            1 => t += (d - 1 + i64::from(d >= 60 && is_leap(year))) * 86400,
            0 => t += d * 86400,
            _ => {
                const CUM: [i64; 13] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];
                let m = usize::from(self.m);
                if !(1..=12).contains(&m) {
                    return None;
                }
                let leap = i64::from(is_leap(year));
                let first = CUM[m - 1] + if m > 2 { leap } else { 0 };
                let len = CUM[m] - CUM[m - 1] + if m == 2 { leap } else { 0 };
                t += first * 86400;
                // Zeller's congruence: weekday of the first day of the month.
                let m1 = (m as i64 + 9) % 12 + 1;
                let y0 = if m <= 2 { year - 1 } else { year };
                let (y1, y2) = (y0 / 100, y0 % 100);
                let dow = (((26 * m1 - 2) / 10 + 1 + y2 + y2 / 4 + y1 / 4 - 2 * y1) % 7 + 7) % 7;
                let mut day = d - dow;
                if day < 0 {
                    day += 7;
                }
                for _ in 1..self.n {
                    if day + 7 >= len {
                        break;
                    }
                    day += 7;
                }
                t += day * 86400;
            }
        }
        Some(t - i64::from(off) + self.secs)
    }
}

impl Posix {
    fn offset_at(&self, t: i64) -> Option<i32> {
        let y = utc_year(t)?;
        let c0 = self.rule[0].change(y, self.off[0])?;
        let c1 = self.rule[1].change(y, self.off[1])?;
        // The DST period wraps the year end when its start falls after its end.
        let dst = if c0 > c1 {
            t < c1 || t >= c0
        } else {
            t >= c0 && t < c1
        };
        Some(self.off[usize::from(dst)])
    }
}

/// Offset of a POSIX TZ string at `unix_seconds`. A DST name without rules uses the
/// United States rules here (no `posixrules` file is consulted).
#[cfg(test)]
fn offset_from_posix_tz(spec: &str, unix_seconds: i64) -> Option<i32> {
    parse_posix(spec.as_bytes()).offset_at(unix_seconds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::PathBuf;

    // Asia/Tokyo and Pacific/Norfolk from the system zone files (fat format, version 2).
    const TOKYO: &[u8] = &[
        84, 90, 105, 102, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4,
        0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 4, 0, 0, 0, 12, 128, 0, 0, 0, 215, 62, 2, 112, 215, 237,
        89, 240, 216, 248, 250, 112, 217, 205, 59, 240, 219, 7, 0, 240, 219, 173, 29, 240, 220,
        230, 226, 240, 221, 140, 255, 240, 3, 1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 131, 3, 0, 0, 0, 0,
        140, 160, 1, 4, 0, 0, 126, 144, 0, 8, 0, 0, 126, 144, 0, 8, 76, 77, 84, 0, 74, 68, 84, 0,
        74, 83, 84, 0, 0, 0, 0, 1, 0, 0, 0, 1, 84, 90, 105, 102, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 4, 0, 0, 0, 12,
        255, 255, 255, 255, 101, 194, 164, 112, 255, 255, 255, 255, 215, 62, 2, 112, 255, 255, 255,
        255, 215, 237, 89, 240, 255, 255, 255, 255, 216, 248, 250, 112, 255, 255, 255, 255, 217,
        205, 59, 240, 255, 255, 255, 255, 219, 7, 0, 240, 255, 255, 255, 255, 219, 173, 29, 240,
        255, 255, 255, 255, 220, 230, 226, 240, 255, 255, 255, 255, 221, 140, 255, 240, 3, 1, 2, 1,
        2, 1, 2, 1, 2, 0, 0, 131, 3, 0, 0, 0, 0, 140, 160, 1, 4, 0, 0, 126, 144, 0, 8, 0, 0, 126,
        144, 0, 8, 76, 77, 84, 0, 74, 68, 84, 0, 74, 83, 84, 0, 0, 0, 0, 1, 0, 0, 0, 1, 10, 74, 83,
        84, 45, 57, 10,
    ];
    const NORFOLK: &[u8] = &[
        84, 90, 105, 102, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7,
        0, 0, 0, 0, 0, 0, 0, 43, 0, 0, 0, 7, 0, 0, 0, 30, 128, 0, 0, 0, 220, 65, 248, 128, 9, 15,
        202, 104, 9, 181, 231, 104, 86, 15, 230, 104, 93, 152, 175, 240, 94, 136, 160, 240, 95,
        120, 145, 240, 96, 104, 130, 240, 97, 88, 115, 240, 98, 72, 100, 240, 99, 56, 85, 240, 100,
        40, 70, 240, 101, 24, 55, 240, 102, 17, 99, 112, 103, 1, 84, 112, 103, 241, 69, 112, 104,
        225, 54, 112, 105, 209, 39, 112, 106, 193, 24, 112, 107, 177, 9, 112, 108, 160, 250, 112,
        109, 144, 235, 112, 110, 128, 220, 112, 111, 112, 205, 112, 112, 105, 248, 240, 113, 89,
        233, 240, 114, 73, 218, 240, 115, 57, 203, 240, 116, 41, 188, 240, 117, 25, 173, 240, 118,
        9, 158, 240, 118, 249, 143, 240, 119, 233, 128, 240, 120, 217, 113, 240, 121, 201, 98, 240,
        122, 185, 83, 240, 123, 178, 127, 112, 124, 162, 112, 112, 125, 146, 97, 112, 126, 130, 82,
        112, 127, 114, 67, 112, 127, 255, 255, 255, 1, 2, 3, 4, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6,
        5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 6, 0, 0, 157,
        120, 0, 0, 0, 0, 157, 128, 0, 4, 0, 0, 161, 184, 0, 10, 0, 0, 175, 200, 1, 16, 0, 0, 161,
        184, 0, 10, 0, 0, 154, 176, 0, 22, 0, 0, 168, 192, 1, 26, 76, 77, 84, 0, 43, 49, 49, 49,
        50, 0, 43, 49, 49, 51, 48, 0, 43, 49, 50, 51, 48, 0, 43, 49, 49, 0, 43, 49, 50, 0, 0, 0, 0,
        1, 1, 1, 1, 84, 90, 105, 102, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 43, 0, 0, 0, 7, 0, 0, 0, 30, 255, 255, 255, 255, 126, 54,
        23, 136, 255, 255, 255, 255, 220, 65, 248, 128, 0, 0, 0, 0, 9, 15, 202, 104, 0, 0, 0, 0, 9,
        181, 231, 104, 0, 0, 0, 0, 86, 15, 230, 104, 0, 0, 0, 0, 93, 152, 175, 240, 0, 0, 0, 0, 94,
        136, 160, 240, 0, 0, 0, 0, 95, 120, 145, 240, 0, 0, 0, 0, 96, 104, 130, 240, 0, 0, 0, 0,
        97, 88, 115, 240, 0, 0, 0, 0, 98, 72, 100, 240, 0, 0, 0, 0, 99, 56, 85, 240, 0, 0, 0, 0,
        100, 40, 70, 240, 0, 0, 0, 0, 101, 24, 55, 240, 0, 0, 0, 0, 102, 17, 99, 112, 0, 0, 0, 0,
        103, 1, 84, 112, 0, 0, 0, 0, 103, 241, 69, 112, 0, 0, 0, 0, 104, 225, 54, 112, 0, 0, 0, 0,
        105, 209, 39, 112, 0, 0, 0, 0, 106, 193, 24, 112, 0, 0, 0, 0, 107, 177, 9, 112, 0, 0, 0, 0,
        108, 160, 250, 112, 0, 0, 0, 0, 109, 144, 235, 112, 0, 0, 0, 0, 110, 128, 220, 112, 0, 0,
        0, 0, 111, 112, 205, 112, 0, 0, 0, 0, 112, 105, 248, 240, 0, 0, 0, 0, 113, 89, 233, 240, 0,
        0, 0, 0, 114, 73, 218, 240, 0, 0, 0, 0, 115, 57, 203, 240, 0, 0, 0, 0, 116, 41, 188, 240,
        0, 0, 0, 0, 117, 25, 173, 240, 0, 0, 0, 0, 118, 9, 158, 240, 0, 0, 0, 0, 118, 249, 143,
        240, 0, 0, 0, 0, 119, 233, 128, 240, 0, 0, 0, 0, 120, 217, 113, 240, 0, 0, 0, 0, 121, 201,
        98, 240, 0, 0, 0, 0, 122, 185, 83, 240, 0, 0, 0, 0, 123, 178, 127, 112, 0, 0, 0, 0, 124,
        162, 112, 112, 0, 0, 0, 0, 125, 146, 97, 112, 0, 0, 0, 0, 126, 130, 82, 112, 0, 0, 0, 0,
        127, 114, 67, 112, 0, 0, 0, 0, 127, 255, 255, 255, 1, 2, 3, 4, 5, 6, 5, 6, 5, 6, 5, 6, 5,
        6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 6,
        0, 0, 157, 120, 0, 0, 0, 0, 157, 128, 0, 4, 0, 0, 161, 184, 0, 10, 0, 0, 175, 200, 1, 16,
        0, 0, 161, 184, 0, 10, 0, 0, 154, 176, 0, 22, 0, 0, 168, 192, 1, 26, 76, 77, 84, 0, 43, 49,
        49, 49, 50, 0, 43, 49, 49, 51, 48, 0, 43, 49, 50, 51, 48, 0, 43, 49, 49, 0, 43, 49, 50, 0,
        0, 0, 0, 1, 1, 1, 1, 10, 60, 43, 49, 49, 62, 45, 49, 49, 60, 43, 49, 50, 62, 44, 77, 49,
        48, 46, 49, 46, 48, 44, 77, 52, 46, 49, 46, 48, 47, 51, 10,
    ];

    /// (instant, offset) pairs taken from the C library (`TZ=:<file> date -d @t +%::z`).
    const TOKYO_GLIBC: &[(i64, i32)] = &[
        (-3000000000, 33539),
        (-2587712401, 33539),
        (-2587712400, 32400),
        (-577962001, 36000),
        (-577962000, 32400),
        (-577961999, 32400),
        (0, 32400),
        (1000000000, 32400),
        (1767225600, 32400),
        (1782864000, 32400),
        (2147483648, 32400),
        (4118083200, 32400),
        (16725225600, 32400),
    ];
    const NORFOLK_GLIBC: &[(i64, i32)] = &[
        (-3000000000, 40312),
        (-2177493113, 40312),
        (-2177493112, 40320),
        (0, 41400),
        (1000000000, 41400),
        (1767225600, 43200),
        (1775314799, 43200),
        (1775314800, 39600),
        (1782864000, 39600),
        (1791039599, 39600),
        (1791039600, 43200),
        (2147483646, 43200),
        (2147483647, 43200),
        (2147483648, 43200),
        (2216818799, 43200),
        (2216818800, 39600),
        (2233148399, 39600),
        (2233148400, 43200),
        (4118083200, 39600),
        (16725225600, 43200),
    ];

    /// (TZ string, [(instant, offset)]) taken from the C library (`TZ=<string> date -d @t +%::z`).
    const POSIX_GLIBC: &[(&str, &[(i64, i32)])] = &[
        (
            "JST-9",
            &[
                (-615513600, 32400),
                (-1, 32400),
                (1767225600, 32400),
                (1782864000, 32400),
                (2210198400, 32400),
                (16740864000, 32400),
            ],
        ),
        (
            "<+0330>-3:30",
            &[
                (-615513600, 12600),
                (-1, 12600),
                (1767225600, 12600),
                (1782864000, 12600),
                (2210198400, 12600),
                (16740864000, 12600),
            ],
        ),
        (
            "<-03>3",
            &[
                (-615513600, -10800),
                (-1, -10800),
                (1767225600, -10800),
                (1782864000, -10800),
                (2210198400, -10800),
                (16740864000, -10800),
            ],
        ),
        (
            "UTC0",
            &[
                (-615513600, 0),
                (-1, 0),
                (1767225600, 0),
                (1782864000, 0),
                (2210198400, 0),
                (16740864000, 0),
            ],
        ),
        (
            "EST5EDT,M3.2.0,M11.1.0",
            &[
                (-615513600, -18000),
                (-1, -18000),
                (1767225600, -18000),
                (1772953199, -18000),
                (1772953200, -14400),
                (1782864000, -14400),
                (1793512799, -14400),
                (1793512800, -18000),
                (2210198400, -18000),
                (16740864000, -14400),
            ],
        ),
        (
            "CET-1CEST,M3.5.0,M10.5.0/3",
            &[
                (-615513600, 3600),
                (-1, 3600),
                (1767225600, 3600),
                (1774745999, 3600),
                (1774746000, 7200),
                (1782864000, 7200),
                (1792889999, 7200),
                (1792890000, 3600),
                (2210198400, 3600),
                (16740864000, 7200),
            ],
        ),
        (
            "AEST-10AEDT,M10.1.0,M4.1.0/3",
            &[
                (-615513600, 39600),
                (-1, 39600),
                (1767225600, 39600),
                (1775318399, 39600),
                (1775318400, 36000),
                (1782864000, 36000),
                (1791043199, 36000),
                (1791043200, 39600),
                (2210198400, 39600),
                (16740864000, 36000),
            ],
        ),
        (
            "<-03>3<-02>,M3.5.0/-2,M10.5.0/-1",
            &[
                (-615513600, -10800),
                (-1, -10800),
                (1767225600, -10800),
                (1774745999, -10800),
                (1774746000, -7200),
                (1782864000, -7200),
                (1792889999, -7200),
                (1792890000, -10800),
                (2210198400, -10800),
                (16740864000, -7200),
            ],
        ),
        (
            "IST-2IDT,M3.4.4/26,M10.5.0",
            &[
                (-615513600, 7200),
                (-1, 7200),
                (1767225600, 7200),
                (1774569599, 7200),
                (1774569600, 10800),
                (1782864000, 10800),
                (1792882799, 10800),
                (1792882800, 7200),
                (2210198400, 7200),
                (16740864000, 10800),
            ],
        ),
        (
            "ABC5DEF4,J60,J300",
            &[
                (-615513600, -18000),
                (-1, -18000),
                (1767225600, -18000),
                (1772348399, -18000),
                (1772348400, -14400),
                (1782864000, -14400),
                (1793080799, -14400),
                (1793080800, -18000),
                (2210198400, -18000),
                (16740864000, -14400),
            ],
        ),
        (
            "ABC5DEF4,59,300",
            &[
                (-615513600, -18000),
                (-1, -18000),
                (1767225600, -18000),
                (1772348399, -18000),
                (1772348400, -14400),
                (1782864000, -14400),
                (1793167199, -14400),
                (1793167200, -18000),
                (2210198400, -18000),
                (16740864000, -14400),
            ],
        ),
        (
            "ABC5DEF4,M3.2.0/167,M11.1.0/-167",
            &[
                (-615513600, -18000),
                (-1, -18000),
                (1767225600, -18000),
                (1773547199, -18000),
                (1773547200, -14400),
                (1782864000, -14400),
                (1792904399, -14400),
                (1792904400, -18000),
                (2210198400, -18000),
                (16740864000, -14400),
            ],
        ),
        (
            "<+13>-13<+14>,M9.5.0/3,M4.1.0/4",
            &[
                (-615513600, 50400),
                (-1, 50400),
                (1767225600, 50400),
                (1775311199, 50400),
                (1775311200, 46800),
                (1782864000, 46800),
                (1790431199, 46800),
                (1790431200, 50400),
                (2210198400, 50400),
                (16740864000, 46800),
            ],
        ),
        (
            "ABC5:",
            &[
                (-615513600, 0),
                (-1, 0),
                (1767225599, 0),
                (1767225600, -18000),
                (1767243599, -18000),
                (1767243600, 0),
                (1782864000, 0),
                (1798761599, 0),
                (1798761600, -18000),
                (1798779599, -18000),
                (1798779600, 0),
                (2210198400, 0),
                (16740864000, 0),
            ],
        ),
        (
            "ABC5DEF4,M3.2.0/x,M11.1.0",
            &[
                (-615513600, -14400),
                (-1, -14400),
                (1767225600, -14400),
                (1767239999, -14400),
                (1767240000, -18000),
                (1772953199, -18000),
                (1772953200, -14400),
                (1782864000, -14400),
                (1798775999, -14400),
                (1798776000, -18000),
                (2210198400, -18000),
                (16740864000, -14400),
            ],
        ),
        (
            "ABC5DEF4,J366,0",
            &[
                (-615513600, -18000),
                (-1, -14400),
                (1767225599, -18000),
                (1767225600, -14400),
                (1767239999, -14400),
                (1767240000, -18000),
                (1782864000, -18000),
                (1798761599, -18000),
                (1798761600, -14400),
                (1798775999, -14400),
                (1798776000, -18000),
                (2210198400, -18000),
                (16740864000, -18000),
            ],
        ),
        (
            "ABC-24:60:61DEF,J1,J365",
            &[
                (-615513600, 89999),
                (-1, 93599),
                (1767052800, 93599),
                (1767052801, 89999),
                (1767225599, 89999),
                (1767225600, 93599),
                (1782864000, 93599),
                (1798588800, 93599),
                (1798588801, 89999),
                (1798761599, 89999),
                (1798761600, 93599),
                (2210198400, 93599),
                (16740864000, 93599),
            ],
        ),
        (
            "FOO",
            &[
                (-615513600, 0),
                (-1, 0),
                (1767225600, 0),
                (1782864000, 0),
                (2210198400, 0),
                (16740864000, 0),
            ],
        ),
        (
            "Nonexistent/Zone",
            &[
                (-615513600, 0),
                (-1, 0),
                (1767225600, 0),
                (1782864000, 0),
                (2210198400, 0),
                (16740864000, 0),
            ],
        ),
        (
            "",
            &[
                (-615513600, 0),
                (-1, 0),
                (1767225600, 0),
                (1782864000, 0),
                (2210198400, 0),
                (16740864000, 0),
            ],
        ),
    ];

    /// Builds a TZif file. Version 0 has only the 32-bit block; the 32-bit block of a later
    /// version holds the transitions that fit it. `footer` of `None` is an empty string.
    fn tzif(
        version: u8,
        types: &[(i32, bool)],
        trans: &[(i64, u8)],
        footer: Option<&str>,
    ) -> Vec<u8> {
        let abbr: Vec<u8> = (0..types.len())
            .flat_map(|i| format!("T{i}\0").into_bytes())
            .collect();
        let block = |wide: bool| {
            let keep: Vec<_> = trans
                .iter()
                .filter(|(t, _)| wide || i32::try_from(*t).is_ok())
                .collect();
            let mut b = b"TZif".to_vec();
            b.push(version);
            b.extend([0u8; 15]);
            for c in [0, 0, 0, keep.len(), types.len(), abbr.len()] {
                b.extend((c as u32).to_be_bytes());
            }
            for (t, _) in &keep {
                if wide {
                    b.extend(t.to_be_bytes());
                } else {
                    b.extend((*t as i32).to_be_bytes());
                }
            }
            b.extend(keep.iter().map(|(_, i)| *i));
            for (i, (off, dst)) in types.iter().enumerate() {
                b.extend(off.to_be_bytes());
                b.extend([u8::from(*dst), 3 * i as u8]);
            }
            b.extend(&abbr);
            b
        };
        let mut out = block(false);
        if version != 0 {
            out.extend(block(true));
            out.extend(format!("\n{}\n", footer.unwrap_or("")).into_bytes());
        }
        out
    }

    fn fake_fs(files: &[(&str, &[u8])]) -> HashMap<PathBuf, Vec<u8>> {
        files
            .iter()
            .map(|(p, d)| (PathBuf::from(p), d.to_vec()))
            .collect()
    }

    fn resolve_in(
        fs: &HashMap<PathBuf, Vec<u8>>,
        tz: Option<&str>,
        tzdir: Option<&str>,
        t: i64,
    ) -> Option<i32> {
        let read = |p: &Path| fs.get(p).cloned();
        resolve(tz.map(OsStr::new), tzdir.map(OsStr::new), &read, t)
    }

    #[test]
    fn posix_strings_match_glibc() {
        for (spec, cases) in POSIX_GLIBC {
            for &(t, want) in *cases {
                assert_eq!(
                    offset_from_posix_tz(spec, t),
                    Some(want),
                    "TZ={spec:?} t={t}"
                );
            }
        }
    }

    #[test]
    fn posix_us_rules_change_at_the_second() {
        // 2026-03-08 02:00 EST = 07:00 UTC; 2026-11-01 02:00 EDT = 06:00 UTC.
        let tz = "EST5EDT,M3.2.0,M11.1.0";
        assert_eq!(offset_from_posix_tz(tz, 1772953199), Some(-18000));
        assert_eq!(offset_from_posix_tz(tz, 1772953200), Some(-14400));
        assert_eq!(offset_from_posix_tz(tz, 1793512799), Some(-14400));
        assert_eq!(offset_from_posix_tz(tz, 1793512800), Some(-18000));
    }

    #[test]
    fn posix_southern_dst_wraps_the_year_end() {
        let tz = "AEST-10AEDT,M10.1.0,M4.1.0/3";
        assert_eq!(offset_from_posix_tz(tz, 1767225600), Some(39600)); // January
        assert_eq!(offset_from_posix_tz(tz, 1782864000), Some(36000)); // July
    }

    #[test]
    fn posix_offsets_are_west_positive_and_clamped() {
        assert_eq!(offset_from_posix_tz("JST-9", 0), Some(32400));
        assert_eq!(offset_from_posix_tz("<-03>3", 0), Some(-10800));
        assert_eq!(offset_from_posix_tz("<+0545>-5:45", 0), Some(20700));
        assert_eq!(offset_from_posix_tz("ABC-99", 0), Some(86400));
        assert_eq!(offset_from_posix_tz("ABC-9:59:60", 0), Some(35999));
        assert_eq!(offset_from_posix_tz("ABC-65537", 0), Some(3600)); // 16-bit wrap
    }

    #[test]
    fn posix_malformed_strings_read_as_the_c_library_does() {
        assert_eq!(offset_from_posix_tz("", 0), Some(0));
        assert_eq!(offset_from_posix_tz("FOO", 0), Some(0));
        assert_eq!(offset_from_posix_tz("Nonexistent/Zone", 0), Some(0));
        assert_eq!(offset_from_posix_tz("<>-9", 0), Some(0));
        // A DST name that cannot be read leaves the DST offset at 0 and both rules at day 0.
        assert_eq!(offset_from_posix_tz("ABC5:", 1767225600), Some(-18000));
        assert_eq!(offset_from_posix_tz("ABC5:", 1782864000), Some(0));
    }

    #[test]
    fn posix_month_outside_the_table_is_undefined() {
        assert_eq!(offset_from_posix_tz("ABC5DEF4,M13.2.0,M11.1.0", 0), None);
        assert_eq!(offset_from_posix_tz("ABC5DEF4,M0.2.0,M11.1.0", 0), None);
    }

    #[test]
    fn posix_dst_name_without_rule_uses_us_rules() {
        let with = "AAA3BBB2,M3.2.0,M11.1.0";
        for t in [
            1700000000, 1772960399, 1772960400, 1782864000, 1793505599, 1793505600,
        ] {
            assert_eq!(
                offset_from_posix_tz("AAA3BBB2", t),
                offset_from_posix_tz(with, t),
                "t={t}"
            );
        }
    }

    #[test]
    fn tzif_tokyo_matches_glibc() {
        for &(t, want) in TOKYO_GLIBC {
            assert_eq!(offset_from_tzif(TOKYO, t), Some(want), "t={t}");
        }
    }

    #[test]
    fn tzif_norfolk_footer_covers_dst_across_the_year_end() {
        for &(t, want) in NORFOLK_GLIBC {
            assert_eq!(offset_from_tzif(NORFOLK, t), Some(want), "t={t}");
        }
    }

    #[test]
    fn tzif_before_first_transition_takes_the_first_standard_type() {
        let z = tzif(
            b'2',
            &[(7200, true), (3600, false), (-3600, false)],
            &[(1000, 2)],
            None,
        );
        assert_eq!(offset_from_tzif(&z, 999), Some(3600));
        assert_eq!(offset_from_tzif(&z, 1000), Some(-3600));
        // Without a footer the last transition's type holds for ever.
        assert_eq!(offset_from_tzif(&z, 5_000_000), Some(-3600));
    }

    #[test]
    fn tzif_without_transitions_ignores_the_footer() {
        let z = tzif(b'2', &[(7200, true), (3600, false)], &[], Some("JST-9"));
        assert_eq!(offset_from_tzif(&z, 1_000_000_000), Some(3600));
    }

    #[test]
    fn tzif_footer_rule_applies_from_the_last_transition() {
        let z = tzif(
            b'2',
            &[(0, false), (3600, true)],
            &[(1000, 1)],
            Some("EST5EDT,M3.2.0,M11.1.0"),
        );
        assert_eq!(offset_from_tzif(&z, 999), Some(0));
        assert_eq!(offset_from_tzif(&z, 1000), Some(-18000));
        assert_eq!(offset_from_tzif(&z, 1_785_000_000), Some(-14400));
        assert_eq!(offset_from_tzif(&z, 1_767_225_600), Some(-18000));
    }

    #[test]
    fn tzif_version_one_reads_the_32_bit_block() {
        let z = tzif(
            0,
            &[(0, false), (3600, true)],
            &[(1000, 1), (2000, 0)],
            None,
        );
        assert_eq!(offset_from_tzif(&z, 999), Some(0));
        assert_eq!(offset_from_tzif(&z, 1500), Some(3600));
        assert_eq!(offset_from_tzif(&z, 2000), Some(0));
        assert_eq!(offset_from_tzif(&z, 1 << 33), Some(0));
    }

    #[test]
    fn tzif_later_versions_read_the_64_bit_block() {
        let z = tzif(b'2', &[(0, false), (3600, true)], &[(1 << 33, 1)], None);
        assert_eq!(offset_from_tzif(&z, (1 << 33) - 1), Some(0));
        assert_eq!(offset_from_tzif(&z, 1 << 33), Some(3600));
    }

    #[test]
    fn tzif_malformed_data_is_none() {
        assert_eq!(offset_from_tzif(b"", 0), None);
        assert_eq!(offset_from_tzif(b"TZif", 0), None);
        assert_eq!(
            offset_from_tzif(
                b"not a zone file at all, but long enough to hold a header........",
                0
            ),
            None
        );
        let ok = |z: &[u8]| offset_from_tzif(z, 5).is_some();
        assert!(ok(&tzif(b'2', &[(0, false)], &[(1, 0)], None)));
        assert!(!ok(&tzif(b'2', &[], &[], None)), "no types");
        assert!(
            !ok(&tzif(b'2', &[(0, false)], &[(1, 1)], None)),
            "type index out of range"
        );
        assert!(
            !ok(&tzif(b'2', &[(0, false)], &[(5, 0), (4, 0)], None)),
            "transitions out of order"
        );
        assert!(
            !ok(&tzif(b'2', &[(i32::MIN, false)], &[], None)),
            "utoff -2^31"
        );
        let mut z = tzif(b'2', &[(0, false)], &[(1, 0)], None);
        let n = z.len();
        z.truncate(n - 1);
        assert!(!ok(&z), "footer without its final newline");
        z.truncate(n - 2);
        assert!(!ok(&z), "no footer");
        let mut bad_flag = tzif(b'2', &[(0, true)], &[], None);
        let flag = bad_flag.windows(2).position(|w| w == [1, 0]).unwrap();
        bad_flag[flag] = 2;
        assert!(!ok(&bad_flag), "isdst 2");
    }

    #[test]
    fn tzif_huge_counts_do_not_allocate_or_loop() {
        let mut h = b"TZif2".to_vec();
        h.extend([0u8; 15]);
        for _ in 0..6 {
            h.extend(u32::MAX.to_be_bytes());
        }
        assert_eq!(offset_from_tzif(&h, 0), None);
        let mut h2 = h.clone();
        h2.extend(vec![0u8; 100]);
        assert_eq!(offset_from_tzif(&h2, 0), None);
    }

    #[test]
    fn tzif_every_prefix_and_every_corruption_is_survived() {
        for z in [TOKYO, NORFOLK] {
            for n in 0..z.len() {
                assert_eq!(offset_from_tzif(&z[..n], 1_700_000_000), None, "prefix {n}");
            }
            let mut c = z.to_vec();
            for i in 0..z.len() {
                for v in [0x00u8, 0x01, 0x7f, 0x80, 0xff, z[i] ^ 0x55] {
                    c[i] = v;
                    for t in [i64::MIN, -3_000_000_000, 0, 1_700_000_000, i64::MAX] {
                        let _ = offset_from_tzif(&c, t);
                    }
                }
                c[i] = z[i];
            }
        }
    }

    #[test]
    fn resolve_unset_reads_etc_localtime() {
        let fs = fake_fs(&[("/etc/localtime", TOKYO)]);
        assert_eq!(resolve_in(&fs, None, None, 0), Some(32400));
        // The TZ value is not looked at when unset, and TZDIR does not move /etc/localtime.
        assert_eq!(resolve_in(&fs, None, Some("/elsewhere"), 0), Some(32400));
        assert_eq!(
            resolve_in(&fake_fs(&[]), None, None, 0),
            Some(0),
            "missing file is UTC"
        );
        assert_eq!(
            resolve_in(&fake_fs(&[("/etc/localtime", b"junk")]), None, None, 0),
            Some(0)
        );
    }

    #[test]
    fn resolve_empty_value_is_utc() {
        let fs = fake_fs(&[
            ("/etc/localtime", TOKYO),
            ("/usr/share/zoneinfo/Asia/Tokyo", TOKYO),
        ]);
        assert_eq!(resolve_in(&fs, Some(""), None, 0), Some(0));
        assert_eq!(resolve_in(&fs, Some(":"), None, 0), Some(0));
    }

    #[test]
    fn resolve_names_are_files_under_tzdir() {
        let fs = fake_fs(&[
            ("/usr/share/zoneinfo/Asia/Tokyo", TOKYO),
            ("/zi/Asia/Tokyo", NORFOLK),
            ("/abs/Tokyo", TOKYO),
        ]);
        for tz in [
            "Asia/Tokyo",
            ":Asia/Tokyo",
            "/usr/share/zoneinfo/Asia/Tokyo",
            ":/abs/Tokyo",
            "/abs/Tokyo",
        ] {
            assert_eq!(resolve_in(&fs, Some(tz), None, 0), Some(32400), "TZ={tz}");
        }
        assert_eq!(
            resolve_in(&fs, Some("Asia/Tokyo"), Some("/zi"), 0),
            Some(41400),
            "TZDIR is used"
        );
        assert_eq!(
            resolve_in(&fs, Some("Asia/Tokyo"), Some(""), 0),
            Some(32400),
            "empty TZDIR is the default"
        );
        assert_eq!(
            resolve_in(&fs, Some("/usr/share/zoneinfo/Asia/Tokyo"), Some("/zi"), 0),
            Some(32400)
        );
    }

    #[test]
    fn resolve_other_values_are_rule_strings_or_utc() {
        let fs = fake_fs(&[
            ("/etc/localtime", TOKYO),
            ("/etc/passwd", b"root:x:0:0:root:/root:/bin/bash\n"),
        ]);
        let at = |tz| resolve_in(&fs, Some(tz), None, 1_000_000_000);
        assert_eq!(at("JST-9"), Some(32400));
        assert_eq!(at("<+0330>-3:30"), Some(12600));
        assert_eq!(at("EST5EDT,M3.2.0,M11.1.0"), Some(-14400));
        assert_eq!(at("Nonexistent/Zone"), Some(0));
        assert_eq!(at(":Nonexistent"), Some(0));
        assert_eq!(at("../../etc/passwd"), Some(0));
        assert_eq!(
            at("/etc/passwd"),
            Some(0),
            "a readable file that is not TZif"
        );
        assert_eq!(at("/etc/localtime"), Some(32400));
    }

    #[test]
    fn resolve_a_zone_file_wins_over_the_rule_string_reading() {
        let fs = fake_fs(&[(
            "/usr/share/zoneinfo/EST5EDT",
            &tzif(b'2', &[(3600, false)], &[], None),
        )]);
        assert_eq!(
            resolve_in(&fs, Some("EST5EDT"), None, 1_000_000_000),
            Some(3600)
        );
        // The same name that is not a TZif file is a rule string.
        let fs = fake_fs(&[("/usr/share/zoneinfo/EST5EDT", b"nonsense")]);
        assert_eq!(
            resolve_in(&fs, Some("EST5EDT"), None, 1_000_000_000),
            Some(-14400)
        );
    }

    #[test]
    fn resolve_dst_name_without_rule_reads_posixrules() {
        // A New York like rule file: the transitions move by the offsets of the string.
        let rules = tzif(
            b'2',
            &[(-18000, false), (-14400, true)],
            &[(1772953200, 1), (1793512800, 0), (1804402800, 1)],
            None,
        );
        let fs = fake_fs(&[("/usr/share/zoneinfo/posixrules", &rules)]);
        for (t, want) in [
            (1700000000, -10800),
            (1772953200, -10800),
            (1782864000, -7200),
            (1793512799, -10800),
            (1900000000, -7200),
        ] {
            assert_eq!(
                resolve_in(&fs, Some("AAA3BBB2"), None, t),
                Some(want),
                "t={t}"
            );
        }
        // Without the file, or with a file of one type, the United States rules apply.
        for fs in [
            fake_fs(&[]),
            fake_fs(&[(
                "/usr/share/zoneinfo/posixrules",
                &tzif(b'2', &[(0, false)], &[], None),
            )]),
        ] {
            assert_eq!(
                resolve_in(&fs, Some("AAA3BBB2"), None, 1772960399),
                Some(-7200)
            );
            assert_eq!(
                resolve_in(&fs, Some("AAA3BBB2"), None, 1793505599),
                Some(-7200)
            );
            assert_eq!(
                resolve_in(&fs, Some("AAA3BBB2"), None, 1793505600),
                Some(-10800)
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn resolve_non_utf8_value_is_unknown() {
        use std::os::unix::ffi::OsStrExt;
        let read = |_: &Path| None;
        assert_eq!(
            resolve(Some(OsStr::from_bytes(b"\xff\xfe")), None, &read, 0),
            None
        );
    }

    #[test]
    fn local_utc_offset_follows_the_platform() {
        assert_eq!(local_utc_offset(1_700_000_000).is_some(), cfg!(unix));
    }

    fn real_read(p: &Path) -> Option<Vec<u8>> {
        std::fs::read(p).ok()
    }

    /// Compares every line of the TSV named by `LOSAT_LOCAL_OFFSET_CHECK` with the code.
    /// Lines: `kind<TAB>value<TAB>unix_seconds<TAB>offset_seconds`; `kind` is `file`
    /// (a TZif path), `tz` (a TZ value, real file system, default TZDIR), `tzdir`
    /// (TZDIR, U+001F, TZ value) or `unset`. See scripts/db2_tz_check.py.
    #[test]
    #[ignore]
    fn glibc_tsv() {
        let Some(path) = std::env::var_os("LOSAT_LOCAL_OFFSET_CHECK") else {
            panic!("LOSAT_LOCAL_OFFSET_CHECK names no TSV file");
        };
        let text = std::fs::read_to_string(path).unwrap();
        let (mut total, mut bad) = (0usize, 0usize);
        let mut undefined = std::collections::BTreeSet::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let (t, want): (i64, i32) = (f[2].parse().unwrap(), f[3].parse().unwrap());
            let got = match f[0] {
                "file" => std::fs::read(f[1])
                    .ok()
                    .and_then(|d| offset_from_tzif(&d, t)),
                "tz" => resolve(Some(OsStr::new(f[1])), None, &real_read, t),
                "unset" => resolve(None, None, &real_read, t),
                "tzdir" => {
                    let (dir, tz) = f[1].split_once('\u{1f}').unwrap();
                    resolve(Some(OsStr::new(tz)), Some(OsStr::new(dir)), &real_read, t)
                }
                k => panic!("unknown kind {k}"),
            };
            total += 1;
            if got.is_none() && f[0] == "tz" {
                // A rule with a month outside 1..=12 reads memory outside the C library's
                // day table; the answer there is not defined, so none is given.
                undefined.insert(f[1].to_string());
            } else if got != Some(want) {
                bad += 1;
                if bad <= 40 {
                    println!("DISAGREE {line} -> {got:?}");
                }
            }
        }
        println!(
            "checked {total} lines, {bad} disagreements, {} TZ values with an undefined rule month",
            undefined.len()
        );
        assert_eq!(bad, 0);
    }
}

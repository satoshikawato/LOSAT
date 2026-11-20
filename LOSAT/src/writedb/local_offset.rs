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
//!
//! The WASI targets read the zone the same way, from what their host gives the module: its
//! environment (`TZ`, `TZDIR`) and the files of its preopened directories. WASI preview 1
//! has no call for the time zone.
//!
//! On Windows the zone is chosen the way Microsoft's C runtime (UCRT) does it: `TZ` in UCRT's
//! own form (`tzn[+|-]hh[:mm[:ss]][dzn]`, with the United States DST rules), else the system
//! zone of `GetTimeZoneInformation` with its DST rules (section "Windows" below). Each target
//! compiles both readers (the unit tests run them all); each uses its own.
#![cfg_attr(windows, allow(dead_code))]

use std::ffi::OsStr;
use std::io::Read;
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
// NCBI reference (598d8ae6): c++/src/build-system/cmake/CMakeChecks.basic-checks.cmake:116
// ```cmake
// check_function_exists(localtime_r HAVE_LOCALTIME_R)
// ```
// The reference platform is Unix (glibc). A WASI module has no time zone call (WASI
// preview 1); it sees the environment and the directories its host gives it, and reading
// `TZ`, `TZDIR` and `/etc/localtime` there as glibc reads them on the host gives the time
// that NCBI's makeblastdb prints on that host. Microsoft's C runtime has no `localtime_r`, so
// NCBI's Windows build takes the `localtime` branch above: UCRT's `_localtime64`, whose zone is
// `TZ` in UCRT's form or the system zone (`ucrt_local_utc_offset`).
/// Seconds east of UTC of the local time zone at `unix_seconds`, as the C library's
/// `localtime_r` gives it on Unix, from the process's environment and file system (on
/// WASI, the host's environment and preopened directories), and as UCRT's `localtime` gives it
/// on Windows. `None` when the zone cannot be determined or the C library gives no local time
/// (the caller then uses UTC).
pub fn local_utc_offset(unix_seconds: i64) -> Option<i64> {
    #[cfg(not(windows))]
    {
        resolve(
            std::env::var_os("TZ").as_deref(),
            std::env::var_os("TZDIR").as_deref(),
            &read_zone_file,
            unix_seconds,
        )
        .map(i64::from)
    }
    #[cfg(windows)]
    {
        ucrt_local_utc_offset(unix_seconds)
    }
}

const DEFAULT_TZDIR: &str = "/usr/share/zoneinfo";

/// The largest zone file that is read. The files of the time zone database are
/// under 10 kB; `TZ` is the user's, and may name anything.
const MAX_ZONE_FILE_BYTES: u64 = 1 << 20;

/// The bytes of `input` when it ends within `cap` bytes; `None` when it is longer or
/// cannot be read. At most `cap + 1` bytes are taken from `input`.
fn read_bounded(input: impl Read, cap: u64) -> Option<Vec<u8>> {
    let mut bytes = Vec::new();
    input
        .take(cap.saturating_add(1))
        .read_to_end(&mut bytes)
        .ok()?;
    (bytes.len() as u64 <= cap).then_some(bytes)
}

/// The bytes of the zone file `path`, or `None` when it cannot be read. Only a
/// regular file of at most `MAX_ZONE_FILE_BYTES` is read: a `TZ` that names a device
/// without an end (`/dev/zero`) or a huge file must not be read into memory. The
/// caller treats an unreadable file as the C library treats a file that is not a
/// zone file.
fn read_zone_file(path: &Path) -> Option<Vec<u8>> {
    if !std::fs::metadata(path).ok()?.is_file() {
        return None;
    }
    read_bounded(std::fs::File::open(path).ok()?, MAX_ZONE_FILE_BYTES)
}

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
            // Both values come from the file: the sum must not overflow.
            self.time(i).saturating_add(shift)
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

// ------------------------------------------------------- Windows: Microsoft's C runtime (UCRT)
//
// The code below follows the sources of UCRT that Microsoft ships with the Windows SDK
// (10.0.26100.0 and 10.0.22621.0 have the same files, `ucrt/...` below). NCBI's Windows build
// calls UCRT's `localtime` (see `local_utc_offset`). The global state that `_tzset` leaves is a
// `UcrtZone`; `localtime` of one instant is `UcrtZone::offset_at`.
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

// UCRT reference (10.0.26100.0): ucrt/time/timeset.cpp:67-69
// ```c++
//     _timezone.initialize(8 * 3600L);
//     _daylight.initialize(1);
//     _dstbias.initialize (-3600);
// ```
/// The start values of the time zone globals of UCRT.
const UCRT_START: UcrtZone = UcrtZone {
    timezone: 8 * 3600,
    daylight: true,
    dstbias: -3600,
    changes: None,
};

// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:145-200
// ```c++
//     if (GetTimeZoneInformation(&tz_info) != 0xFFFFFFFF)
//     {
//         // Record that the API was used:
//         tz_api_used = 1;
//
//         // Derive _timezone value from Bias and StandardBias fields.
//         timezone = tz_info.Bias * 60;
//
//         if (tz_info.StandardDate.wMonth != 0)
//             timezone += tz_info.StandardBias * 60;
//
//         // Check to see if there is a daylight time bias. Since the StandardBias
//         // has been added into _timezone, it must be compensated for in the
//         // value computed for _dstbias:
//         if (tz_info.DaylightDate.wMonth != 0 && tz_info.DaylightBias != 0)
//         {
//             daylight = 1;
//             dstbias = (tz_info.DaylightBias - tz_info.StandardBias) * 60;
//         }
//         else
//         {
//             daylight = 0;
//             ...
//             dstbias = 0;
//         }
//         ...
//     }
//
//     _set_timezone(timezone);
//     _set_daylight(daylight);
//     _set_dstbias(dstbias);
// ```
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

// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:264-335
// ```c++
//     rsize_t const tzname_length = 3;
//     ...
//     // Skip first few characters if present.
//     for (rsize_t i = 0; i < tzname_length; ++i)
//     {
//         if (*tz_env)
//         {
//             ++tz_env;
//         }
//     }
//
//     // The time difference is of the form:
//     //     [+|-]hh[:mm[:ss]]
//     // Check for the minus sign first:
//     bool const is_negative_difference = *tz_env == L'-';
//     if (is_negative_difference)
//     {
//         ++tz_env;
//     }
//     ...
//     // process, then skip over, the hours
//     timezone = wcstol(tz_env, &dummy, decimal_base) * 3600;
//     while (*tz_env == '+' || (*tz_env >= L'0' && *tz_env <= L'9'))
//     {
//         ++tz_env;
//     }
//
//     // Check if minutes were specified:
//     if (*tz_env == L':')
//     {
//         // Process, then skip over, the minutes
//         timezone += wcstol(++tz_env, &dummy, decimal_base) * 60;
//         while (*tz_env >= L'0' && *tz_env <= L'9')
//         {
//             ++tz_env;
//         }
//
//         // Check if seconds were specified:
//         if (*tz_env == L':')
//         {
//             // Process, then skip over, the seconds:
//             timezone += wcstol(++tz_env, &dummy, decimal_base);
//             while (*tz_env >= L'0' && *tz_env <= L'9')
//             {
//                 ++tz_env;
//             }
//         }
//     }
//
//     if (is_negative_difference)
//     {
//         timezone = -timezone;
//     }
//
//     // Finally, check for a DST zone suffix:
//     daylight = *tz_env ? 1 : 0;
//     ...
//     _set_timezone(timezone);
//     _set_daylight(daylight);
// ```
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

// UCRT reference (10.0.26100.0): ucrt/inc/corecrt_internal_strtox.h:248-345
// ```c++
//         while (is_space(c, loc))
//         {
//             c = source.get();
//         }
//     ...
//     // Optional sign (+ or -):
//     if (c == '-')
//     {
//         flags |= FL_NEGATIVE;
//     }
//
//     if (c == '-' || c == '+')
//     {
//         c = source.get();
//     }
//     ...
//     for (;;)
//     {
//         unsigned const digit{parse_digit(c)};
//         if (digit >= static_cast<unsigned>(base))
//         {
//             ...
//             break;
//         }
//         ...
//     }
//     ...
//     // If we failed to read any digits, there's no number to be read:
//     if ((flags & FL_READ_DIGIT) == 0)
//     {
//         source.restore_state(initial_state);
//         return 0;
//     }
//
//     if (is_overflow_condition(flags, number))
//     {
//         ptd.get_errno().set(ERANGE);
//         ...
//         else if (flags & FL_NEGATIVE)
//         {
//             return minimum_signed_value(UnsignedInteger());
//         }
//         else
//         {
//             return maximum_signed_value(UnsignedInteger());
//         }
//     }
// ```
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

// UCRT reference (10.0.26100.0): ucrt/inc/corecrt_internal_strtox.h:88-158
// ```c++
// __forceinline int __cdecl wide_character_to_digit(wchar_t const c) throw()
// {
//     #define DIGIT_RANGE_TEST(zero) \
//         if (c < zero)              \
//             return -1;             \
//                                    \
//         if (c < zero + 10)         \
//             return c - zero;
//
//     DIGIT_RANGE_TEST(0x0030)        // 0030;DIGIT ZERO
//     if (c < 0xFF10)                 // FF10;FULLWIDTH DIGIT ZERO
//     {
//         DIGIT_RANGE_TEST(0x0660)    // 0660;ARABIC-INDIC DIGIT ZERO
//         ... (06F0, 0966, 09E6, 0A66, 0AE6, 0B66, 0C66, 0CE6, 0D66, 0E50, 0ED0, 0F20,
//              1040, 17E0)
//         DIGIT_RANGE_TEST(0x1810)    // 1810;MONGOLIAN DIGIT ZERO
//
//         return -1;
//     }
//
//     #undef DIGIT_RANGE_TEST
//
//     if (c < 0xFF10 + 10)            // FF10;FULLWIDTH DIGIT ZERO
//         return c - 0xFF10;
//
//     return -1;
// }
// ...
// __forceinline unsigned __cdecl parse_digit(wchar_t const c) throw()
// {
//     int const value = wide_character_to_digit(c);
//     if (value != -1)
//         return static_cast<unsigned>(value);
//
//     if (__ascii_iswalpha(c))
//         return static_cast<unsigned>(__ascii_towupper(c) - 'A' + 10);
//
//     return static_cast<unsigned>(-1);
// }
// ```
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

// UCRT reference (10.0.26100.0): ucrt/convert/iswctype.cpp:42-54
// ```c++
//     if (c == WEOF)
//         return 0;
//
//     if (c < 256)
//         return static_cast<int>(_pwctype[c] & mask);
//
//     wchar_t const wide_character = c;
//     wint_t char_type = 0;
//     if (__acrt_GetStringTypeW(CT_CTYPE1, &wide_character, 1, &char_type) == 0)
//         return 0;
//
//     return static_cast<int>(char_type & mask);
// ```
// UCRT reference (10.0.26100.0): ucrt/locale/ctype.cpp:53-57, 76, 177, 204 (`_wctype`, the
// entries with `_SPACE`)
// ```c++
//     _SPACE | _CONTROL ,             // 09 (HT)
//     ... (0A, 0B, 0C)
//     _SPACE | _CONTROL ,             // 0D (CR)
//     _SPACE | _BLANK ,               // 20 SPACE
//     _SPACE | _CONTROL ,             // 85 (NEL / NEXT LINE)
//     _SPACE ,                        // A0 (NBSP / NO-BREAK SPACE)
// ```
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

// UCRT reference (10.0.26100.0): ucrt/time/days.cpp:14-24
// ```c++
// extern "C" int const _days[] =
// {
//     -1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364
// };
//
// extern "C" int const _lpdays[] =
// {
//     -1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
// };
// ```
const UCRT_DAYS: [i32; 13] = [-1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364];
const UCRT_LPDAYS: [i32; 13] = [-1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];

// UCRT reference (10.0.26100.0): ucrt/inc/corecrt_internal_time.h:153-175
// ```c++
// bool __cdecl __crt_time_is_leap_year(TimeType const yr) throw()
// {
//     if (yr % 4 == 0 && yr % 100 != 0)
//         return true;
//
//     if ((yr + 1900) % 400 == 0)
//         return true;
//
//     return false;
// }
// ...
// TimeType __cdecl __crt_time_elapsed_leap_years(TimeType const yr) throw()
// {
//     static TimeType const leap_years_between_1900_and_1970 = 17;
//
//     TimeType const elapsed_leap_years = ((yr - 1) / 4) - ((yr - 1) / 100) + ((yr + 299) / 400);
//
//     return elapsed_leap_years - leap_years_between_1900_and_1970;
// }
// ```
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

// UCRT reference (10.0.26100.0): ucrt/inc/corecrt_internal_time.h:37-40, 98
// ```c++
// #define _MAX_LOCAL_TIME (14 * 60 * 60)
// ...
// #define _MIN_LOCAL_TIME (-12 * 60 * 60)
// ...
//         max_time_t = 0x793582affLL,
// ```
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

// UCRT reference (10.0.26100.0): ucrt/time/gmtime.cpp:54-77, 97-134
// ```c++
//     // Determine the years since 1900. Start by ignoring leap years:
//     int tmptim = static_cast<int>(caltim / _YEAR_SEC) + 70;
//     caltim -= static_cast<__time64_t>(tmptim - 70) * _YEAR_SEC;
//
//     // Correct for elapsed leap years:
//     caltim -= static_cast<__time64_t>(__crt_time_elapsed_leap_years(tmptim)) * _DAY_SEC;
//
//     // If we have underflowed the __time64_t range (i.e., if caltim < 0),
//     // back up one year, adjusting the correction if necessary.
//     if (caltim < 0)
//     {
//         caltim += static_cast<__time64_t>(_YEAR_SEC);
//         tmptim--;
//         if (__crt_time_is_leap_year(tmptim))
//         {
//             caltim += _DAY_SEC;
//             ...
//         }
//     }
//     ...
//     _VALIDATE_RETURN_ERRCODE_NOEXC(caltim >= _MIN_LOCAL_TIME, EINVAL)
//     ...
//     _VALIDATE_RETURN_ERRCODE_NOEXC(caltim <= time_traits::max_time_t + _MAX_LOCAL_TIME, EINVAL)
//     ...
//     ptm->tm_yday = static_cast<int>(caltim / _DAY_SEC);
//     caltim -= static_cast<TimeType>(ptm->tm_yday) * _DAY_SEC;
//     ...
//     ptm->tm_hour = static_cast<int>(caltim / 3600);
//     caltim -= static_cast<TimeType>(ptm->tm_hour) * 3600L;
//
//     ptm->tm_min = static_cast<int>(caltim / 60);
//     ptm->tm_sec = static_cast<int>(caltim - (ptm->tm_min) * 60);
// ```
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

// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:446-496
// ```c++
//     if (datetype == date_type::day_in_month)
//     {
//         // Figure out the year-day of the start of the month:
//         yearday = 1 + (__crt_time_is_leap_year(year)
//             ? _lpdays[month - 1]
//             : _days[month - 1]);
//
//         // Figureo ut the day of the week of the start of the month:
//         monthdow = (yearday + ((year - 70) * 365) +
//                     __crt_time_elapsed_leap_years(year) + _BASE_DOW) % 7;
//
//         // Figure out the year-day of the transition date:
//         if (monthdow <= dayofweek)
//             yearday += (dayofweek - monthdow) + (week - 1) * 7;
//         else
//             yearday += (dayofweek - monthdow) + week * 7;
//
//         // We may have to adjust the calculation above if week == 5 (meaning the
//         // last instance of the day in the month).  Check if the year falls
//         // beyond after month and adjust accordingly:
//         int const days_to_compare = __crt_time_is_leap_year(year)
//             ? _lpdays[month]
//             : _days[month];
//
//         if (week == 5 && yearday > days_to_compare)
//         {
//             yearday -= 7;
//         }
//     }
//     else
//     {
//         yearday = __crt_time_is_leap_year(year)
//             ? _lpdays[month - 1]
//             : _days[month - 1];
//
//         yearday += date;
//     }
//
//     if (trantype == transition_type::start_of_dst)
//     {
//         dststart.yd = yearday;
//         dststart.ms = msec + (1000 * (sec + 60 * (min + 60 * hour)));
//         ...
//     }
//     else // end_of_dst
//     {
//         dstend.yd = yearday;
//         dstend.ms = msec + (1000 * (sec + 60 * (min + 60 * hour)));
// ```
// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:546-559 (the choice of the form)
// ```c++
//             if (tz_info.DaylightDate.wYear == 0)
//             {
//                 cvtdate(
//                     transition_type::start_of_dst,
//                     date_type::day_in_month,
//                     tb->tm_year,
//                     tz_info.DaylightDate.wMonth,
//                     tz_info.DaylightDate.wDay,
//                     tz_info.DaylightDate.wDayOfWeek,
//                     0,
//                     ...
// ```
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

// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:611-651
// ```c++
//             // The GetTimeZoneInformation API was not used, or failed.  We use
//             // the USA Daylight Savings Time rules as a fallback.
//             int startmonth = 3; // March
//             int startweek  = 2; // Second week
//             int endmonth   = 11;// November
//             int endweek    = 1; // First week
//
//             // The rules changed in 2007:
//             if (107 > tb->tm_year)
//             {
//                 startmonth = 4; // April
//                 startweek  = 1; // first week
//                 endmonth   = 10;// October
//                 endweek    = 5; // last week
//             }
//
//             cvtdate(
//                 transition_type::start_of_dst,
//                 date_type::day_in_month,
//                 tb->tm_year,
//                 startmonth,
//                 startweek,
//                 0, // Sunday
//                 0,
//                 2, // 02:00 (2 AM)
//                 ...
// ```
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
    // UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:498-512 (cvtdate, end of DST)
    // ```c++
    //         // The converted date is still a DST date.  We must convert to a standard
    //         // (local) date while being careful the millisecond field does not
    //         // overflow or underflow
    //         _ERRCHECK(_get_dstbias(&dstbias));
    //         dstend.ms += (dstbias * 1000);
    //         if (dstend.ms < 0)
    //         {
    //             dstend.ms += milliseconds_per_day;
    //             dstend.yd--;
    //         }
    //         else if (dstend.ms >= milliseconds_per_day)
    //         {
    //             dstend.ms -= milliseconds_per_day;
    //             dstend.yd++;
    //         }
    // ```
    // UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:655-685
    // ```c++
    //     // Handle simple cases first:
    //     if (dststart.yd < dstend.yd)
    //     {
    //         // Northern hemisphere ordering:
    //         if (tb->tm_yday < dststart.yd || tb->tm_yday > dstend.yd)
    //             return 0;
    //
    //         if (tb->tm_yday > dststart.yd && tb->tm_yday < dstend.yd)
    //             return 1;
    //     }
    //     else
    //     {
    //         // Southern hemisphere ordering:
    //         if (tb->tm_yday < dstend.yd || tb->tm_yday > dststart.yd)
    //             return 1;
    //
    //         if (tb->tm_yday > dstend.yd && tb->tm_yday < dststart.yd)
    //             return 0;
    //     }
    //
    //     long const ms = 1000 * (tb->tm_sec + 60 * tb->tm_min + 3600 * tb->tm_hour);
    //
    //     if (tb->tm_yday == dststart.yd)
    //     {
    //
    //         return ms >= dststart.ms ? 1 : 0;
    //     }
    //     else
    //     {
    //         return ms < dstend.ms ? 1 : 0;
    //     }
    // ```
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

    // UCRT reference (10.0.26100.0): ucrt/time/localtime.cpp:46-80
    // ```c++
    //     _VALIDATE_RETURN_ERRCODE_NOEXC(*ptime >= 0,                       EINVAL);
    //     _VALIDATE_RETURN_ERRCODE_NOEXC(*ptime <= time_traits::max_time_t, EINVAL);
    //
    //     __tzset();
    //     ...
    //     if (*ptime > 3 * _DAY_SEC && *ptime < time_traits::max_time_t - 3 * _DAY_SEC)
    //     {
    //         ...
    //         TimeType ltime = *ptime - timezone;
    //
    //         errno_t status0 = time_traits::gmtime_s(ptm, &ltime);
    //         if (status0 != 0)
    //             return status0;
    //
    //         // Check and adjust for daylight savings time:
    //         if (daylight && _isindst(ptm))
    //         {
    //             ltime -= dstbias;
    //
    //             errno_t const status1 = time_traits::gmtime_s(ptm, &ltime);
    //             if (status1 != 0)
    //                 return status1;
    //
    //             ptm->tm_isdst = 1;
    //         }
    //     }
    // ```
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

// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:347-363
// ```c++
//     // Get the value of the TZ environment variable:
//     wchar_t local_env_buffer[local_env_buffer_size];
//     wchar_t* const tz_env = get_tz_environment_variable(local_env_buffer);
//     ...
//     if (tz_env == nullptr || tz_env[0] == '\0')
//         return tzset_from_system_nolock();
//
//     return tzset_from_environment_nolock(tz_env);
// ```
// `get_tz_environment_variable` reads `TZ` with `_wgetenv_s` (`tzset.cpp:76-103`), from the
// CRT's wide environment, which UCRT makes from the process environment block
// (`ucrt/env/environment_initialization.cpp:166-188`); the name is compared without case
// (`ucrt/env/getenv.cpp:44`). `std::env::var_os` reads the same block, without case, and keeps
// its UTF-16 units. `localtime` runs `_tzset` once per process, at its first call (`__tzset`,
// `tzset.cpp:391-416`); LOSAT asks at each time it prints, with the same environment and zone.
/// The offset of UCRT's `localtime` at `unix_seconds` from `TZ` and the system zone.
#[cfg(windows)]
fn ucrt_local_utc_offset(unix_seconds: i64) -> Option<i64> {
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
// (`ucrt/inc/corecrt_internal.h:577-582`); its body is not among the shipped sources, and
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

    // Unix: `TZ` unset or a zone, else UTC; Windows: `TZ` in UCRT's form or the system zone.
    #[test]
    fn local_utc_offset_is_known_on_every_platform() {
        assert!(local_utc_offset(1_700_000_000).is_some());
    }

    #[test]
    fn a_zone_file_is_read_up_to_the_cap() {
        assert_eq!(read_bounded(&b"abc"[..], 3), Some(b"abc".to_vec()));
        assert_eq!(read_bounded(&b"abc"[..], 4), Some(b"abc".to_vec()));
        assert_eq!(read_bounded(&b"abcd"[..], 3), None);
        assert_eq!(read_bounded(&b""[..], 0), Some(Vec::new()));
        // An input without an end is given up after `cap + 1` bytes.
        assert_eq!(read_bounded(std::io::repeat(0), 16), None);
        assert_eq!(read_bounded(std::io::repeat(0), MAX_ZONE_FILE_BYTES), None);
    }

    #[cfg(unix)]
    #[test]
    fn a_device_is_not_read_as_a_zone_file() {
        assert_eq!(read_zone_file(Path::new("/dev/zero")), None);
        assert_eq!(read_zone_file(Path::new("/dev/null")), None);
        assert_eq!(read_zone_file(Path::new("/")), None);
        // `TZ=:/dev/zero` is then a rule string that is none: UTC.
        for tz in [":/dev/zero", "/dev/zero"] {
            assert_eq!(
                resolve(Some(OsStr::new(tz)), None, &read_zone_file, 0),
                Some(0)
            );
        }
    }

    #[test]
    fn posixrules_transitions_at_the_ends_of_the_range_do_not_overflow() {
        // The standard-time transition moves by the DST offset of the string (-7200) and
        // the DST transition by +7200: both sums leave the range of the type.
        let rules = tzif(
            b'2',
            &[(-18000, false), (-14400, true)],
            &[(i64::MIN, 0), (i64::MAX, 1)],
            None,
        );
        let fs = fake_fs(&[("/usr/share/zoneinfo/posixrules", &rules)]);
        assert_eq!(resolve_in(&fs, Some("AAA3BBB2"), None, 0), Some(-10800));
        assert_eq!(
            resolve_in(&fs, Some("AAA3BBB2"), None, i64::MAX),
            Some(-7200)
        );
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

    // ------------------------------------------------------------ Windows (UCRT)

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
            let rules = parse_posix(posix.as_bytes());
            let mut t = 1_735_689_600; // 2025-01-01
            while t < 1_830_297_600 {
                // 2028-01-01
                for u in [t - 1, t, t + 1] {
                    assert_eq!(
                        zone.offset_at(u),
                        rules.offset_at(u).map(i64::from),
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

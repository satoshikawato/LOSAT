//! The `-t_*_pct` settings of the windowmasker application: the percentages of the `##pct:`
//! metadata lines of a counts file, which `CSeqMaskerIstatFactory::create` reads before the
//! counts and which replace the masker's thresholds (`s_GetCountByPct`). blastn passes none
//! (`seq_masker.hpp:178-181`, -1 each), so its counts are read without this step.

use super::istat::{extract_uint, is_c_space};
use super::{MaskerParams, StatType, WinMaskError};

/// The four percentages of the windowmasker application (`-t_low_pct`, `-t_extend_pct`,
/// `-t_thres_pct`, `-t_high_pct`); `None` for one not given (NCBI's -1.0).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PctSettings {
    pub min_pct: Option<f64>,
    pub extend_pct: Option<f64>,
    pub thres_pct: Option<f64>,
    pub max_pct: Option<f64>,
}

impl PctSettings {
    /// Whether any percentage is given.
    pub fn any(&self) -> bool {
        self.min_pct.is_some()
            || self.extend_pct.is_some()
            || self.thres_pct.is_some()
            || self.max_pct.is_some()
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:206-214
/// (`count_map` is not empty here.)
fn count_by_pct(pct: f64, count_map: &[f64]) -> u32 {
    count_map
        .iter()
        .position(|&value| value >= pct)
        .unwrap_or(count_map.len() - 1) as u32
}

/// Why LOSAT refuses a `##pct:` line. One byte: with a `&'static str` here `PctLine` became the
/// largest variant of `WinMaskError`, whose discriminant then moved into a niche of its `String`,
/// so that every `Ok` of the masker's count lookup (`UnitCounts::at`, for each unit of each window
/// of blastn `-window_masker_db`) was written and tested as a 64-bit constant instead of a tag
/// byte: +0.16 % instructions in blastn with binary counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PctWhy {
    /// Its words are numbers that C++ libraries read differently.
    RuntimeWord,
    /// A binary line whose largest count is `u32::MAX`.
    MaxCountWraps,
    /// A binary line whose table does not fit the memory.
    NoMemory,
}

impl PctWhy {
    /// The words of the refusal that follow "the percentage line '<line>' of the counts file".
    pub(super) fn text(self) -> &'static str {
        match self {
            PctWhy::RuntimeWord => "holds a number that C++ libraries read differently (an infinity, a NaN, a hexadecimal number, a number with letters after it, or one beyond the range of a normal double), so LOSAT does not read it",
            PctWhy::MaxCountWraps => "names 4294967295 as its largest count, for which NCBI BLAST+ makes room for 0 percentages (the count plus one, in 32 bits) and reads past them (undefined behaviour in C++)",
            PctWhy::NoMemory => "names more percentages than the memory holds",
        }
    }
}

/// LOSAT's refusal of the `##pct:` line `line`, `why`.
fn refusal(line: &[u8], why: PctWhy) -> WinMaskError {
    WinMaskError::PctLine {
        line: String::from_utf8_lossy(line).into_owned(),
        why,
    }
}

/// The number of decimal digits at the start of `text`.
fn digits(text: &[u8]) -> usize {
    text.iter().take_while(|b| b.is_ascii_digit()).count()
}

/// One `>>` of a `double` from `text` at `*at`, which moves past the field read: `Ok(None)` when
/// the text ends at white space (nothing stored, the stream has failed); else the value stored
/// and whether the stream has failed. A stream that has failed, or that has reached the end of
/// the text, stores nothing more.
///
/// NCBI's Linux builds (libstdc++, the C++ standard's `num_get` stages 1 to 3 with `strtod`) take
/// into the field an optional sign, digits, a point, digits, and after a digit an `e` or `E`
/// with an optional sign and digits; they stop at any other character. A field without a digit
/// before its exponent, or with an exponent mark but no exponent digits, stores 0 and fails.
///
/// `Err` where libc++ (NCBI's macOS builds) reads otherwise: a field whose letters after the
/// sign start `inf` or `nan` (libc++ reads them as such), a number followed at once by a letter
/// that libc++ takes into its field (`a` to `f`, `i`, `n`, `p`, `x`, either case: `0x10`, `5x`,
/// `1e5e`; the field then does not convert whole and gives 0 there), and a number whose value is
/// not a normal double or an exact zero (the libraries flag overflow and underflow differently).
fn extract_double(text: &[u8], at: &mut usize) -> Result<Option<(f64, bool)>, ()> {
    while *at < text.len() && is_c_space(text[*at]) {
        *at += 1;
    }
    if *at == text.len() {
        return Ok(None);
    }
    let mut i = *at;
    let negative = text[i] == b'-';
    if matches!(text[i], b'+' | b'-') {
        i += 1;
    }
    let named = |word: &[u8]| text[i..].len() >= 3 && text[i..i + 3].eq_ignore_ascii_case(word);
    if named(b"inf") || named(b"nan") {
        return Err(());
    }
    let int = &text[i..i + digits(&text[i..])];
    i += int.len();
    let mut frac: &[u8] = &[];
    if text.get(i) == Some(&b'.') {
        i += 1;
        frac = &text[i..i + digits(&text[i..])];
        i += frac.len();
    }
    let mantissa = !int.is_empty() || !frac.is_empty();
    let mut exponent = None;
    if mantissa && matches!(text.get(i), Some(b'e' | b'E')) {
        i += 1;
        let sign = i;
        if matches!(text.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let count = digits(&text[i..]);
        exponent = Some(&text[sign..i + count]);
        i += count;
        if count == 0 {
            *at = i;
            return Ok(Some((0.0, true)));
        }
    }
    *at = i;
    if !mantissa {
        return Ok(Some((0.0, true)));
    }
    if text
        .get(i)
        .is_some_and(|c| b"abcdefinpxABCDEFINPX".contains(c))
    {
        return Err(());
    }
    let mut number = String::from(if negative { "-" } else { "" });
    for part in [int, &b"."[..], frac] {
        number.push_str(std::str::from_utf8(part).map_err(|_| ())?);
    }
    if let Some(exponent) = exponent {
        number.push('e');
        number.push_str(std::str::from_utf8(exponent).map_err(|_| ())?);
    }
    let value: f64 = number.parse().map_err(|_| ())?;
    let zero = int.iter().chain(frac).all(|&b| b == b'0');
    if value.is_normal() || (value == 0.0 && zero) {
        Ok(Some((value, false)))
    } else {
        Err(())
    }
}

/// The percentages of the `##pct:` lines of `metadata` (`count_map`), as the factory reads them
/// for a counts file of `stat_type` with `std::istringstream` and `>>`: the ascii formats a count
/// (`Uint4`) and a percentage from each line, the percentage 0 unless the count was read and the
/// stream goes on; the binary formats the largest count of the first line and, into a table of
/// that count plus one zeros, as many percentages as the stream gives before it ends or fails (a
/// failed read stores 0 and ends it). Words after those are not read.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:239-276
/// A line of the six letters alone has failed before (`check_pct_lines`). The counts are read as
/// [`extract_uint`] reads them, the percentages as [`extract_double`] does; a percentage that C++
/// libraries read differently is refused (`WinMaskError::PctLine`), and so are a binary line
/// whose largest count is `u32::MAX` (`check_pct_lines` refuses it first) and one whose table
/// does not fit the memory.
pub(super) fn count_map(
    metadata: &[Vec<u8>],
    stat_type: StatType,
) -> Result<Vec<f64>, WinMaskError> {
    let mut pct_lines = metadata.iter().filter(|line| line.starts_with(b"##pct:"));
    let mut map = Vec::new();
    let read = |line: &[u8], text: &[u8], at: &mut usize| {
        extract_double(text, at).map_err(|()| refusal(line, PctWhy::RuntimeWord))
    };
    match stat_type {
        StatType::Ascii | StatType::OAscii => {
            for line in pct_lines {
                let line = line.as_slice();
                let text = line.get(7..).unwrap_or_default();
                let mut at = 0;
                let going =
                    matches!(extract_uint(text, &mut at), Some((_, false))) && at < text.len();
                let pct = if going {
                    read(line, text, &mut at)?.map_or(0.0, |(value, _)| value)
                } else {
                    0.0
                };
                map.push(pct);
            }
        }
        StatType::Binary | StatType::OBinary => {
            if let Some(line) = pct_lines.next() {
                let line = line.as_slice();
                let text = line.get(7..).unwrap_or_default();
                let mut at = 0;
                let (max, mut going) = match extract_uint(text, &mut at) {
                    Some((max, failed)) => (max, !failed && at < text.len()),
                    None => (0, false),
                };
                if max == u32::MAX {
                    return Err(refusal(line, PctWhy::MaxCountWraps));
                }
                let size = max as usize + 1;
                map.try_reserve_exact(size)
                    .map_err(|_| refusal(line, PctWhy::NoMemory))?;
                map.resize(size, 0.0);
                for value in map.iter_mut() {
                    if !going {
                        break;
                    }
                    let Some((read_value, failed)) = read(line, text, &mut at)? else {
                        break;
                    };
                    *value = read_value;
                    going = !failed && at < text.len();
                }
            }
        }
        StatType::Unknown => {}
    }
    Ok(map)
}

/// The thresholds of `params` after the percentages `pct` of a counts file of `stat_type` with
/// the metadata lines `metadata`; `warn` gets the factory's warning when the file has no
/// percentages, as NCBI posts it (before the reader's warnings and errors).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:278-299
/// (The application's percentages are in 0..100 (`CArgAllow_Doubles`), so a given one is
/// `>= 0.0`. With no percentage given the map is not used, so its values are not checked.)
pub(super) fn apply_pct(
    params: &MaskerParams,
    metadata: &[Vec<u8>],
    stat_type: StatType,
    pct: &PctSettings,
    warn: &mut dyn FnMut(String),
) -> Result<MaskerParams, WinMaskError> {
    let mut params = *params;
    if !pct.any() {
        return Ok(params);
    }
    let map = count_map(metadata, stat_type)?;
    if map.is_empty() {
        let name = if pct.min_pct.is_some() {
            "min_pct"
        } else if pct.extend_pct.is_some() {
            "extend_pct"
        } else if pct.thres_pct.is_some() {
            "thres_pct"
        } else {
            "max_pct"
        };
        warn(format!(
            "counts file contains no percentage information; {name} settings will not be applied"
        ));
        return Ok(params);
    }
    if let Some(value) = pct.min_pct {
        params.min_score = count_by_pct(value, &map);
    }
    if let Some(value) = pct.extend_pct {
        params.textend = count_by_pct(value, &map);
    }
    if let Some(value) = pct.thres_pct {
        params.cutoff_score = count_by_pct(value, &map);
    }
    if let Some(value) = pct.max_pct {
        params.max_score = count_by_pct(value, &map);
    }
    Ok(params)
}

#[cfg(test)]
#[path = "pct_tests.rs"]
mod tests;

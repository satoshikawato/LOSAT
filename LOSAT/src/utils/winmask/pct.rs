//! The `-t_*_pct` settings of the windowmasker application: the percentages of the `##pct:`
//! metadata lines of a counts file, which `CSeqMaskerIstatFactory::create` reads before the
//! counts and which replace the masker's thresholds (`s_GetCountByPct`). blastn passes none
//! (`seq_masker.hpp:178-181`, -1 each), so its counts are read without this step.

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

/// The white space that `operator>>` skips in NCBI's `istringstream` (the "C" locale's
/// `isspace`).
fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
}

/// The words of `text` between white space.
fn words(text: &[u8]) -> impl Iterator<Item = &[u8]> {
    text.split(|&byte| is_space(byte))
        .filter(|word| !word.is_empty())
}

/// A `Uint4` that `operator>>` reads from `word`: decimal digits whose value fits. Any other
/// word is read by the C++ runtime's number parser (signs, overflow, a number that a letter
/// ends), which LOSAT does not reproduce: `None`.
fn unsigned_word(word: &[u8]) -> Option<u32> {
    if word.is_empty() || !word.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(word).ok()?.parse().ok()
}

/// A `double` that `operator>>` reads from `word`: a decimal number (`[+-]digits[.digits]` with
/// an optional exponent) whose value is a normal double or zero. A correctly rounded conversion
/// gives one value for such a word whatever the library (NCBI's `##pct:` values are written with
/// six significant digits). Any other word (NaN, infinity, a value out of the range of a normal
/// double, trailing letters) is read by the C++ runtime's rules, which LOSAT does not
/// reproduce: `None`.
fn double_word(word: &[u8]) -> Option<f64> {
    let mut i = 0;
    if matches!(word.first(), Some(b'+' | b'-')) {
        i += 1;
    }
    let int_digits = word[i..].iter().take_while(|b| b.is_ascii_digit()).count();
    i += int_digits;
    let mut frac_digits = 0;
    if word.get(i) == Some(&b'.') {
        i += 1;
        frac_digits = word[i..].iter().take_while(|b| b.is_ascii_digit()).count();
        i += frac_digits;
    }
    if int_digits + frac_digits == 0 {
        return None;
    }
    if matches!(word.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(word.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let exp_digits = word[i..].iter().take_while(|b| b.is_ascii_digit()).count();
        if exp_digits == 0 {
            return None;
        }
        i += exp_digits;
    }
    if i != word.len() {
        return None;
    }
    let value: f64 = std::str::from_utf8(word).ok()?.parse().ok()?;
    let nonzero_digits = word
        .iter()
        .take_while(|&&b| b != b'e' && b != b'E')
        .any(|b| (b'1'..=b'9').contains(b));
    let in_range = value.is_finite() && (value.is_normal() || (value == 0.0 && !nonzero_digits));
    in_range.then_some(value)
}

/// LOSAT's refusal of a `##pct:` line that NCBI reads by the C++ runtime's rules.
fn pct_line_refusal(line: &[u8]) -> WinMaskError {
    WinMaskError::PctLine {
        line: String::from_utf8_lossy(line).into_owned(),
    }
}

/// The percentages of the `##pct:` lines of `metadata` (`count_map`), as the factory reads them
/// for a counts file of `stat_type`.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:239-276
/// A line of the six letters alone has failed before (`check_pct_lines`). The words that the
/// stream reads must be the numbers NCBI's writer prints (`unsigned_word`, `double_word`), and a
/// binary line must hold its `max_map_count + 1` values: else LOSAT refuses the line
/// (`WinMaskError::PctLine`), as NCBI's reads of other words follow the C++ runtime's number
/// parser. Words after those the stream reads are not read.
pub(super) fn count_map(
    metadata: &[Vec<u8>],
    stat_type: StatType,
) -> Result<Vec<f64>, WinMaskError> {
    let mut pct_lines = metadata.iter().filter(|line| line.starts_with(b"##pct:"));
    let mut map = Vec::new();
    match stat_type {
        StatType::Ascii | StatType::OAscii => {
            for line in pct_lines {
                let mut read = words(&line[7..]);
                let count = read.next().and_then(unsigned_word);
                let pct = read.next().and_then(double_word);
                match (count, pct) {
                    (Some(_), Some(pct)) => map.push(pct),
                    _ => return Err(pct_line_refusal(line)),
                }
            }
        }
        StatType::Binary | StatType::OBinary => {
            if let Some(line) = pct_lines.next() {
                let mut read = words(&line[7..]);
                let max = read
                    .next()
                    .and_then(unsigned_word)
                    .ok_or_else(|| pct_line_refusal(line))?;
                for _ in 0..=u64::from(max) {
                    let value = read
                        .next()
                        .and_then(double_word)
                        .ok_or_else(|| pct_line_refusal(line))?;
                    map.push(value);
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
mod tests {
    use super::*;

    fn lines(text: &[&str]) -> Vec<Vec<u8>> {
        text.iter().map(|line| line.as_bytes().to_vec()).collect()
    }

    #[test]
    fn the_first_value_at_or_above_the_percentage_gives_the_count() {
        let map = [0.0, 50.0, 90.0, 99.5, 100.0];
        assert_eq!(count_by_pct(0.0, &map), 0);
        assert_eq!(count_by_pct(50.0, &map), 1);
        assert_eq!(count_by_pct(90.1, &map), 3);
        assert_eq!(count_by_pct(100.0, &map), 4);
        assert_eq!(count_by_pct(100.5, &[0.0, 99.0]), 1, "none: the last index");
    }

    #[test]
    fn text_formats_read_every_line_and_binary_ones_the_first() {
        let md = lines(&[
            "##note:x",
            "##pct: 0 0",
            "##pct: 1 89.7027",
            "##pct: 2 1e-05 extra",
        ]);
        assert_eq!(
            count_map(&md, StatType::Ascii).unwrap(),
            vec![0.0, 89.7027, 1e-05]
        );
        let md = lines(&["##pct: 2 0 50.5 100", "##pct: 1 7 8"]);
        assert_eq!(
            count_map(&md, StatType::OBinary).unwrap(),
            vec![0.0, 50.5, 100.0]
        );
        assert!(count_map(&lines(&["##note:x"]), StatType::Binary)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn words_the_runtime_would_parse_are_refused() {
        for line in [
            "##pct: 1 -nan",
            "##pct: -1 5",
            "##pct: 1 5x",
            "##pct: 1 1e999",
            "##pct: 1 1e-320",
            "##pct: 1",
            "##pct: 99999999999 1",
        ] {
            assert!(
                matches!(
                    count_map(&lines(&[line]), StatType::Ascii),
                    Err(WinMaskError::PctLine { .. })
                ),
                "{line}"
            );
        }
        assert!(count_map(&lines(&["##pct: 3 0 1 2"]), StatType::Binary).is_err());
    }

    #[test]
    fn a_file_without_percentages_names_the_first_setting_given() {
        let mut posts = Vec::new();
        let pct = PctSettings {
            extend_pct: Some(99.0),
            thres_pct: Some(99.5),
            ..PctSettings::default()
        };
        let params = MaskerParams::default();
        let got = apply_pct(&params, &[], StatType::Ascii, &pct, &mut |w| posts.push(w)).unwrap();
        assert_eq!(got, params);
        assert_eq!(
            posts,
            vec![
                "counts file contains no percentage information; extend_pct settings will not be applied"
                    .to_string()
            ]
        );
    }
}

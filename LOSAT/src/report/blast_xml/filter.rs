//! The filter string of `Parameters_filter` (`BlastFilteringOptionsToString`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:320-432
//! (and likewise `L;` / `S %d %1.1f %1.1f;` for SEG, `R;` / `R -d db;`, `W -t %d;` /
//! `W -d db;`, `m;` for mask at hash, `F` for an empty string.)

/// `kDustLevel`, `kDustWindow`, `kDustLinker` (core/blast_options.c:46-48).
const DUST_DEFAULTS: (i32, i32, i32) = (20, 64, 1);
/// `kSegWindow`, `kSegLocut`, `kSegHicut` (core/blast_seg.c:45-47).
const SEG_DEFAULTS: (i32, f64, f64) = (12, 2.2, 2.5);

/// What the window masker option writes (`windowMaskerOptions`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WindowMasker<'a> {
    /// `W -t <taxid>;`
    Taxid(i32),
    /// `W -d <database>;`
    Database(&'a str),
}

/// The filter string of a search: `dust` and `seg` are the options' (level, window, linker)
/// and (window, locut, hicut) when the filter is on, `repeats` the database of blastn's
/// `-filtering_db`, `window_masker` the window masker option, `mask_at_hash` whether the
/// filtering is applied at hashing (soft masking).
pub fn filter_string(
    dust: Option<(i32, i32, i32)>,
    seg: Option<(i32, f64, f64)>,
    repeats: Option<&str>,
    window_masker: Option<WindowMasker<'_>>,
    mask_at_hash: bool,
) -> Vec<u8> {
    let mut text = String::new();
    if let Some(dust) = dust {
        if dust == DUST_DEFAULTS {
            text.push_str("L;");
        } else {
            // `char buffer[24]`: snprintf cuts the text at 23 bytes.
            text.push_str(&truncated(format!("D {} {} {};", dust.0, dust.1, dust.2)));
        }
    }
    if let Some(seg) = seg {
        if seg == SEG_DEFAULTS {
            text.push_str("L;");
        } else {
            text.push_str(&truncated(format!(
                "S {} {} {};",
                seg.0,
                glibc_fixed1(seg.1),
                glibc_fixed1(seg.2)
            )));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:371-388 (the command
    // line always gives the database: `R -d <name>;`, never the plain `R;`.)
    if let Some(database) = repeats {
        text.push_str("R -d ");
        text.push_str(database);
        text.push(';');
    }
    match window_masker {
        Some(WindowMasker::Taxid(taxid)) if taxid != 0 => {
            text.push_str(&truncated(format!("W -t {taxid};")));
        }
        Some(WindowMasker::Taxid(_)) | None => {}
        Some(WindowMasker::Database(database)) => {
            text.push_str("W -d ");
            text.push_str(database);
            text.push(';');
        }
    }
    if mask_at_hash {
        text.push_str("m;");
    }
    if text.is_empty() {
        text.push('F');
    }
    text.into_bytes()
}

/// glibc's `%1.1f`: one decimal, the infinities as `inf` and `-inf`, NaN as `nan` or `-nan` by
/// its sign bit (Rust prints NaN as `NaN`). NCBI 2.17.0's blastp `-seg "12 2.2 -nan"` writes
/// `S 12 2.2 -nan;`, `-seg "12 +nan(1) 2.5"` `S 12 nan 2.5;` (sw8/nfix/probe/seg2).
fn glibc_fixed1(value: f64) -> String {
    match value {
        nan if nan.is_nan() && nan.is_sign_negative() => "-nan".into(),
        nan if nan.is_nan() => "nan".into(),
        value => format!("{value:.1}"),
    }
}

/// `snprintf(buffer, 24, ...)`: at most 23 bytes.
fn truncated(mut text: String) -> String {
    text.truncate(23);
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(bytes: Vec<u8>) -> String {
        String::from_utf8(bytes).expect("ASCII")
    }

    #[test]
    fn strings_observed_from_ncbi() {
        // blastn megablast / blastn / dc-megablast: DUST at its defaults, mask at hash.
        assert_eq!(
            s(filter_string(Some((20, 64, 1)), None, None, None, true)),
            "L;m;"
        );
        // blastn-short and `-dust no`: only mask at hash.
        assert_eq!(s(filter_string(None, None, None, None, true)), "m;");
        // `-dust no -soft_masking false` and blastp without SEG: empty.
        assert_eq!(s(filter_string(None, None, None, None, false)), "F");
        // blastp `-seg yes -soft_masking true`; blastx / tblastn / tblastx default SEG.
        assert_eq!(
            s(filter_string(None, Some((12, 2.2, 2.5)), None, None, true)),
            "L;m;"
        );
        assert_eq!(
            s(filter_string(None, Some((12, 2.2, 2.5)), None, None, false)),
            "L;"
        );
        // `-seg '10 1.8 2.1'`
        assert_eq!(
            s(filter_string(None, Some((10, 1.8, 2.1)), None, None, false)),
            "S 10 1.8 2.1;"
        );
    }

    #[test]
    fn custom_dust_and_seg_and_masker() {
        assert_eq!(
            s(filter_string(Some((25, 64, 1)), None, None, None, true)),
            "D 25 64 1;m;"
        );
        assert_eq!(
            s(filter_string(
                Some((20, 64, 1)),
                Some((12, 2.2, 2.5)),
                None,
                None,
                false
            )),
            "L;L;"
        );
        assert_eq!(
            s(filter_string(
                None,
                None,
                None,
                Some(WindowMasker::Taxid(9606)),
                true
            )),
            "W -t 9606;m;"
        );
        assert_eq!(
            s(filter_string(
                None,
                None,
                None,
                Some(WindowMasker::Database("db.wm")),
                false
            )),
            "W -d db.wm;"
        );
        assert_eq!(
            s(filter_string(
                None,
                None,
                None,
                Some(WindowMasker::Taxid(0)),
                false
            )),
            "F"
        );
    }

    #[test]
    fn seg_cutoffs_are_printed_as_glibc_rounds_them() {
        // %1.1f of an exact tie goes to the even digit: 2.25 -> 2.2, 0.35 (below the tie in
        // binary) -> 0.3, 0.45 (above) -> 0.5.
        assert_eq!(
            s(filter_string(
                None,
                Some((10, 2.25, 0.35)),
                None,
                None,
                false
            )),
            "S 10 2.2 0.3;"
        );
        assert_eq!(
            s(filter_string(
                None,
                Some((10, 0.45, 0.125)),
                None,
                None,
                false
            )),
            "S 10 0.5 0.1;"
        );
    }

    // NCBI 2.17.0 blastp `-outfmt 5` with `-seg` cuts that are not finite (sw8/nfix/probe/seg2).
    #[test]
    fn infinite_and_nan_seg_cutoffs_are_printed_as_glibc_prints_them() {
        let seg = |locut: f64, hicut: f64| {
            s(filter_string(
                None,
                Some((12, locut, hicut)),
                None,
                None,
                false,
            ))
        };
        assert_eq!(seg(f64::INFINITY, 2.5), "S 12 inf 2.5;");
        assert_eq!(seg(2.2, -f64::NAN), "S 12 2.2 -nan;");
        assert_eq!(
            seg(f64::from_bits(0x7ff8_0000_0000_0001), 2.5),
            "S 12 nan 2.5;"
        );
        assert_eq!(seg(f64::NEG_INFINITY, 2.5), "S 12 -inf 2.5;");
        assert_eq!(seg(-f64::NAN, -f64::NAN), "S 12 -nan -nan;");
    }

    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blastfilter_unit_test.cpp:1683-1697 FilterOptionsToStringLargeData
    #[test]
    fn a_long_repeat_database_is_written_whole() {
        // The database of blastn's `-filtering_db` is appended to the text whole, not through a
        // buffer of fixed size, after a custom DUST.
        let database = "X".repeat(4096);
        let text = s(filter_string(
            Some((20, 128, 1)),
            None,
            Some(&database),
            None,
            false,
        ));
        assert!(text.starts_with("D 20 128 1;R -d XXXXXXXXXXXXXXXXXXXX"));
        assert_eq!(text, format!("D 20 128 1;R -d {database};"));
    }

    #[test]
    fn a_long_custom_string_is_cut_like_the_24_byte_buffer() {
        let text = s(filter_string(
            Some((2147483647, 2147483647, 2147483647)),
            None,
            None,
            None,
            false,
        ));
        assert_eq!(text, "D 2147483647 2147483647");
        assert_eq!(text.len(), 23);
    }
}

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
            // `buffer` holds 24 bytes: snprintf cuts the text at 23 bytes.
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

/// `snprintf` into `buffer` of 24 bytes: at most 23 bytes.
fn truncated(mut text: String) -> String {
    text.truncate(23);
    text
}

#[cfg(test)]
#[path = "filter_tests.rs"]
mod tests;

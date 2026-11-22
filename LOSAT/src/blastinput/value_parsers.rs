//! Shared CLI value grammar; engine configuration remains typed.
use crate::utils::seg::SegParams;
use std::path::PathBuf;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-420
#[derive(Debug, Clone, PartialEq)]
pub enum DustSpec {
    No,
    Yes,
    Parameters {
        level: u32,
        window: usize,
        linker: usize,
    },
}

impl DustSpec {
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:46-48,62-64:
    pub fn params(&self) -> Option<(u32, usize, usize)> {
        match *self {
            Self::No => None,
            Self::Yes => Some((20, 64, 1)),
            Self::Parameters {
                level,
                window,
                linker,
            } => Some((level, window, linker)),
        }
    }
}

/// A `-dust` value as NCBI reads it: `no`, `yes`, or three signed decimal numbers separated
/// by single spaces. The error is NCBI's message; the masker replaces values out of its
/// ranges with its defaults (`CSymDustMasker`, `utils/dust.rs`), after NCBI's conversion
/// to `Uint4` (`as u32`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:375-384
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-426
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:195-198
/// `NStr::Split` keeps the empty fields between adjacent spaces; `NStr::StringToInt` (flags 0)
/// is `ncbi_str::string_to_int`.
pub fn parse_dust_filtering(value: &str) -> Result<DustSpec, String> {
    match value {
        "no" => return Ok(DustSpec::No),
        "yes" => return Ok(DustSpec::Yes),
        _ => {}
    }
    let tokens: Vec<_> = value.split(' ').collect();
    if tokens.len() != 3 {
        return Err("Invalid number of arguments to filtering option".into());
    }
    // `NStr::StringToInt(tokens[i])`, flags 0.
    let number = |token: &str| {
        crate::ncbi_str::string_to_int(token.as_bytes(), 0)
            .map_err(|_| "Invalid input for filtering parameters".to_string())
    };
    let (level, window, linker) = (number(tokens[0])?, number(tokens[1])?, number(tokens[2])?);
    Ok(DustSpec::Parameters {
        level: level as u32,
        window: window as u32 as usize,
        linker: linker as u32 as usize,
    })
}

/// An integer argument as NCBI reads it (`CArg_Integer`): a decimal number with an
/// optional sign, or else, for a value that starts with `0x` or `0X`, a hexadecimal number
/// (no digit after the prefix reads as 0), within the range of `int`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:383-390
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:788-791
/// `CArg_Integer` derives from `CArg_Int8`: `s_StringToInt8` is `ncbi_int8` (its messages),
/// then the range of `int`.
pub fn ncbi_integer(value: &str) -> Result<i32, String> {
    let int8 = ncbi_int8(value)?;
    i32::try_from(int8).map_err(|_| "Integer value is out of range".to_string())
}

/// A 64-bit integer argument as NCBI reads it (`CArg_Int8`, such as `-dbsize`): what
/// `ncbi_integer` reads, without the range of `int`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:351-360
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129
/// `s_StringToInt8` is `ncbi_app::args::string_to_int8` (the one reading of NCBI's integer
/// arguments). LOSAT's messages: "Integer value is out of range" for a hexadecimal number beyond
/// `Int8` (every byte after `0x` a hexadecimal digit), "expected an integer" for any other value
/// that does not read (a byte after `0x` that is not a hexadecimal digit included: NCBI's base
/// 16 conversion throws on an overflow before it reaches such a byte, ncbistr.cpp:852).
pub fn ncbi_int8(value: &str) -> Result<i64, String> {
    crate::ncbi_app::args::string_to_int8(value.as_bytes()).map_err(|_| {
        let bytes = value.as_bytes();
        let hexadecimal = bytes
            .get(..2)
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"0x"))
            && bytes[2..].iter().all(u8::is_ascii_hexdigit);
        if hexadecimal {
            "Integer value is out of range".to_string()
        } else {
            "expected an integer".to_string()
        }
    })
}

/// A real argument as NCBI reads it (`CArg_Double`): `NStr::StringToDouble(value,
/// fDecimalPosixOrLocal)` (`ncbi_str::string_to_double_posix_or_local`): the first byte is a
/// digit, a point, a comma or a sign, glibc's `strtod` reads the number (decimal, hexadecimal,
/// the infinities, NaN with its sign and `n-char-sequence`) or NCBI's own `StringToDoublePosix`
/// where it reads further (a mantissa that ends in `e` or `E`, in x87 extended precision, and an
/// `inf` after a point, `.inf`), and the number ends the value. LOSAT's message for NCBI's conversion error: "expected a number".
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
pub fn ncbi_double(value: &str) -> Result<f64, String> {
    crate::ncbi_str::string_to_double_posix_or_local(value.as_bytes())
        .map_err(|_| "expected a number".to_string())
}

/// An integer argument with a constraint: after the conversion, NCBI's constraint reads the
/// value again with `NStr::StringToDouble` (default flags, `ncbi_str::string_to_double`), which
/// does not read a `0x` prefix without digits (the only form that `ncbi_integer` reads and
/// `strtod` does not end at), so the value is illegal.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_constrained_integer(value: &str) -> Result<i32, String> {
    let n = ncbi_integer(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    Ok(n)
}

/// Whether the constraint's `NStr::StringToDouble(value)` (default flags) reads `value`
/// (`ncbi_constrained_integer`).
fn reads_as_double(value: &str) -> bool {
    crate::ncbi_str::string_to_double(value.as_bytes()).is_ok()
}

/// An integer argument of `CArg_Integer` with NCBI's lower bound (`at_least`,
/// `CArgAllowValuesGreaterThanOrEqual`).
pub fn ncbi_integer_at_least(value: &str, at_least: i32) -> Result<i32, String> {
    let n = ncbi_constrained_integer(value)?;
    if n < at_least {
        return Err(format!("expected an integer >= {at_least}"));
    }
    Ok(n)
}

/// A non-negative `Int8` argument (`CArg_Int8` with `CArgAllowValuesGreaterThanOrEqual(0)`):
/// `-searchsp`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
/// The constraint reads the value as `ncbi_constrained_integer`'s does.
pub fn nonnegative_ncbi_int8(value: &str) -> Result<i64, String> {
    let n = ncbi_int8(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    if n < 0 {
        return Err("expected an integer >= 0".into());
    }
    Ok(n)
}

/// A non-negative integer argument (`CArg_Integer` with `CArgAllowValuesGreaterThanOrEqual(0)`),
/// such as `-culling_limit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3302
pub fn nonnegative_ncbi_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 0)
}

/// A `-num_descriptions` or `-num_alignments` value: 0 or more, an `Int4` in NCBI's
/// integer forms (the sign and hexadecimal forms of `nonnegative_ncbi_integer`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685
pub fn report_count(value: &str) -> Result<usize, String> {
    nonnegative_ncbi_integer(value).map(|n| n as usize)
}

/// A discontiguous megablast template length: 16, 18 or 21.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:719-727
/// The constraint converts the value again with `NStr::StringToInt` in base 10, so a
/// hexadecimal value that the integer argument reads fails it.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-222,239
/// `NStr::StringToInt(value)` (flags 0) is `ncbi_str::string_to_int`.
pub fn blastn_template_length(value: &str) -> Result<u8, String> {
    ncbi_constrained_integer(value)?;
    match crate::ncbi_str::string_to_int(value.as_bytes(), 0) {
        Ok(n @ (16 | 18 | 21)) => Ok(n as u8),
        _ => Err("expected one of 16, 18, 21".into()),
    }
}
/// A BLASTN word size: 4 or more (the upper bound, 100, is an option check,
/// `blastn/scoring.rs`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:168-170
pub fn blastn_word_size(value: &str) -> Result<usize, String> {
    ncbi_integer_at_least(value, 4).map(|n| n as usize)
}
/// A BLASTN count of 1 or more (`-num_threads`, `-max_target_seqs`, `-max_hsps`), as
/// NCBI's arguments (c++/src/algo/blast/blastinput/blast_args.cpp:203-207,2731-2732,3162-3163).
pub fn blastn_count(value: &str) -> Result<usize, String> {
    ncbi_positive_integer(value).map(|n| n as usize)
}
/// `blastn_count` as an `int`, the type of BLASTX's fields.
pub fn ncbi_positive_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 1)
}
/// A BLASTN reward: 0 or more, as NCBI's argument (blast_args.cpp:658-659). A reward of 0
/// runs with a penalty of 0 (matrix-only scoring, `-task rmblastn`); with another penalty LOSAT
/// rejects it before the search (`blastn/scoring.rs`).
pub fn blastn_reward(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 0)
}
/// A BLASTN penalty: 0 or less, as NCBI's argument (blast_args.cpp:651-652), whose
/// constraint reads the value as the `>=` constraint does (blast_input_aux.hpp:135-138).
pub fn blastn_penalty(value: &str) -> Result<i32, String> {
    let n = ncbi_constrained_integer(value)?;
    if n > 0 {
        return Err("expected an integer <= 0".into());
    }
    Ok(n)
}
/// A BLASTN gap cost: any `int` (blast_args.cpp:175-182).
pub fn blastn_gap_cost(value: &str) -> Result<i32, String> {
    ncbi_integer(value)
}
/// A BLASTN e-value: any value that NCBI's argument reads (`ncbi_double`). NCBI's option
/// check rejects 0 or less (`blastn/scoring.rs`); an infinite or NaN e-value searches as in
/// NCBI (`+inf`, `-nan`, `1e999`; `blastn/ncbi_cutoffs.rs`, E2g).
pub fn blastn_evalue(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}
/// A BLASTN percent identity, in 0..=100 (NaN is outside).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:188-192
pub fn blastn_percentage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A percentage, in 0..=100 (NaN is outside): `-perc_identity`, and
/// `-qcov_hsp_perc` of every program (the `*_query_coverage` parsers).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201
/// The constraint reads the value with the argument's own conversion (`ncbi_double`), so it
/// tests the argument's number (`1e` is 1, `0x10` is 16).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
fn percentage(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(0.0..=100.0).contains(&n) {
        return Err("expected a percentage in 0..=100".into());
    }
    Ok(n)
}
/// A BLASTN `-qcov_hsp_perc` value (`percentage`).
pub fn blastn_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A BLASTP `-qcov_hsp_perc` value (`percentage`).
pub fn blastp_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A PSIBLAST `-qcov_hsp_perc` value (`percentage`).
pub fn psiblast_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// An RPSBLAST `-qcov_hsp_perc` value (`percentage`).
pub fn rpsblast_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// An RPSTBLASTN `-qcov_hsp_perc` value (`percentage`).
pub fn rpstblastn_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A BLASTX `-qcov_hsp_perc` value (`percentage`).
pub fn blastx_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A TBLASTN `-qcov_hsp_perc` value (`percentage`).
pub fn tblastn_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}
/// A TBLASTX `-qcov_hsp_perc` value (`percentage`).
pub fn tblastx_query_coverage(value: &str) -> Result<f64, String> {
    percentage(value)
}

/// A `-best_hit_overhang` or `-best_hit_score_edge` value: above 0 and below 0.5 (NaN is
/// outside), of every program that has the options.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/hspfilter_besthit.h:80-91
/// The constraint reads the value with the argument's own conversion (`ncbi_double`); NCBI
/// BLAST+ 2.17.0 refuses 0 and 0.5 with "Illegal value, expected (>0 and <0.5)" (line NP O6a's
/// probe).
pub fn best_hit_value(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(n > 0.0 && n < 0.5) {
        return Err("expected a number above 0 and below 0.5".into());
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:995-1006
pub fn genetic_code(value: &str) -> Result<u8, String> {
    let n = genetic_code_number(value).ok_or("invalid genetic code")?;
    if !crate::ncbi_app::args::is_genetic_code(i32::from(n)) {
        return Err("unsupported genetic code".into());
    }
    Ok(n)
}

/// The number of a genetic code argument as `CArgAllowGeneticCodeInteger` converts it
/// (`NStr::StringToInt(value)`, flags 0: `ncbi_str::string_to_int`), when it is 0 to 255 and
/// written without a minus sign; `None` for the other values, for which LOSAT's message says
/// "invalid" (`genetic_code`, TBLASTN's `tblastn_genetic_code`) where it says "unsupported" for
/// a number of that range outside NCBI's codes. NCBI refuses both alike.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:998-1007
pub(crate) fn genetic_code_number(value: &str) -> Option<u8> {
    match crate::ncbi_str::string_to_int(value.as_bytes(), 0) {
        Ok(n) if !value.starts_with('-') => u8::try_from(n).ok(),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-51
/// A BLASTN output file: a name shorter than 256 bytes, as NCBI's constraint; `-` is
/// standard output, and a file that cannot be created is reported when it is opened
/// (`blastn/blast_engine/run.rs`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
pub fn blastn_output_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().try_map(|value| {
        if crate::ncbi_app::args::file_name_length(value.as_encoded_bytes()) >= 256 {
            return Err("Illegal value, expected file name length < 256".to_string());
        }
        Ok(PathBuf::from(value))
    })
}

/// A BLASTN input file: any value, as NCBI's `eInputFile` argument; `-` is standard input
/// and a file that cannot be opened is reported when it is opened (`blastn/blast_engine/run.rs`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
pub fn blastn_input_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().map(PathBuf::from)
}

/// An input file of any program (`blastn_input_path`).
pub fn ncbi_input_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    blastn_input_path()
}

/// An output file of any program (`blastn_output_path`).
pub fn ncbi_output_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    blastn_output_path()
}

/// A real argument of BLASTP (`ncbi_double`).
pub fn blastp_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of PSIBLAST (`ncbi_double`).
pub fn psiblast_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of RPSBLAST (`ncbi_double`).
pub fn rpsblast_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of RPSTBLASTN (`ncbi_double`).
pub fn rpstblastn_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of TBLASTN (`ncbi_double`).
pub fn tblastn_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of BLASTX (`ncbi_double`).
pub fn blastx_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of TBLASTX (`ncbi_double`).
pub fn tblastx_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A real argument of BLASTN (`ncbi_double`).
pub fn blastn_real(value: &str) -> Result<f64, String> {
    ncbi_double(value)
}

/// A `-threshold` value: NCBI's constraint reads it again with `NStr::StringToDouble` and its
/// default flags (`ncbi_str::string_to_double`: glibc's `strtod` alone, so a value that only
/// `StringToDoublePosix` reads, such as `1e`, makes it throw), and requires at least 0 (a NaN
/// fails it); a value that the constraint reads is the same number as the argument's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_threshold(value: &str) -> Result<f64, String> {
    let threshold = ncbi_double(value)?;
    // `value >= m_MinValue` (a NaN fails it).
    let allowed =
        crate::ncbi_str::string_to_double(value.as_bytes()).is_ok_and(|checked| checked >= 0.0);
    if !allowed {
        return Err("Illegal value, expected greater or equal to 0".into());
    }
    Ok(threshold)
}

/// A BLASTP `-threshold` value (`ncbi_threshold`).
pub fn blastp_threshold(value: &str) -> Result<f64, String> {
    ncbi_threshold(value)
}

/// A BLASTX `-threshold` value (`ncbi_threshold`).
pub fn blastx_threshold(value: &str) -> Result<f64, String> {
    ncbi_threshold(value)
}

/// A PSIBLAST `-threshold` value (`ncbi_threshold`).
pub fn psiblast_threshold(value: &str) -> Result<f64, String> {
    ncbi_threshold(value)
}

/// A TBLASTN `-threshold` value (`ncbi_threshold`).
pub fn tblastn_threshold_value(value: &str) -> Result<f64, String> {
    ncbi_threshold(value)
}

/// A TBLASTX `-threshold` value (`ncbi_threshold`).
pub fn tblastx_threshold(value: &str) -> Result<f64, String> {
    ncbi_threshold(value)
}

/// A Boolean argument (`CArg_Boolean`, read with `NStr::StringToBool`:
/// `ncbi_str::string_to_bool`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:489-497
pub fn ncbi_boolean(value: &str) -> Result<bool, String> {
    crate::ncbi_str::string_to_bool(value)
        .ok_or_else(|| format!("Argument cannot be converted: `{value}'"))
}

/// A `-mt_mode` value: an integer from 0 to 2 (`CArgAllowValuesBetween(0, 2, true)`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3177
pub fn ncbi_mt_mode(value: &str) -> Result<u8, String> {
    let mode = ncbi_integer_at_least(value, 0)?;
    u8::try_from(mode)
        .ok()
        .filter(|mode| *mode <= 2)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 2".to_string())
}

/// A protein `-word_size` (`CArgAllowValuesGreaterThanOrEqual(2)` for a protein query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170
pub fn protein_word_size(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 2)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:386-407
#[derive(Debug, Clone, PartialEq)]
pub enum SegSpec {
    No,
    Yes,
    WindowLocutHicut { window: i32, locut: f64, hicut: f64 },
}

impl SegSpec {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    /// The options as NCBI's filter string prints them (`BlastFilteringOptionsToString`): the
    /// defaults for `yes`, the given window, locut and hicut otherwise (the substitution of
    /// the defaults for a value not above 0 is `params`', after the options are set).
    pub fn filter_options(&self) -> Option<(i32, f64, f64)> {
        match *self {
            Self::No => None,
            Self::Yes => Some((12, 2.2, 2.5)),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => Some((window, locut, hicut)),
        }
    }

    /// The SEG parameters of the query filter: a window, locut or hicut that is not above
    /// 0 keeps NCBI's default (`SegParametersNewAa`), then SEG checks the parameters
    /// (`SegParams::new`, blast_seg.c:2247-2258).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1147-1154
    #[inline]
    pub fn params(&self) -> Option<SegParams> {
        match *self {
            Self::No => None,
            Self::Yes => Some(SegParams::default()),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => {
                let defaults = SegParams::default();
                Some(SegParams::new(
                    usize::try_from(window)
                        .ok()
                        .filter(|&window| window > 0)
                        .unwrap_or(defaults.window),
                    if locut > 0.0 { locut } else { defaults.locut },
                    if hicut > 0.0 { hicut } else { defaults.hicut },
                ))
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-407
    pub fn to_ncbi_cli_string(&self) -> String {
        match self {
            Self::No => "no".to_string(),
            Self::Yes => "yes".to_string(),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => format!("{window} {locut} {hicut}"),
        }
    }
}

/// A `-seg` value as NCBI reads it: `no`, `yes`, or a window (an `int`), a locut and a
/// hicut separated by single spaces (`x_TokenizeFilteringArgs`, see `parse_dust_filtering`).
/// The errors are NCBI's messages; LOSAT also rejects a locut or hicut that is not finite.
/// The values are used as NCBI's query filter uses them (`SegSpec::params`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-406
#[cfg(any(test, target_arch = "wasm32"))]
pub fn parse_seg_filtering(value: &str) -> Result<SegSpec, String> {
    match value {
        "no" => return Ok(SegSpec::No),
        "yes" => return Ok(SegSpec::Yes),
        _ => {}
    }
    let tokens: Vec<&str> = value.split(' ').collect();
    if tokens.len() != 3 {
        return Err("Invalid number of arguments to filtering option".into());
    }
    let invalid = || "Invalid input for filtering parameters".to_string();
    let window = tokens[0].parse::<i32>().map_err(|_| invalid())?;
    let locut = tokens[1].parse::<f64>().map_err(|_| invalid())?;
    let hicut = tokens[2].parse::<f64>().map_err(|_| invalid())?;
    if !locut.is_finite() || !hicut.is_finite() {
        return Err("a SEG locut or hicut that is not finite is not supported by LOSAT".into());
    }
    Ok(SegSpec::WindowLocutHicut {
        window,
        locut,
        hicut,
    })
}

/// The TBLASTX `-outfmt` value of web ABI v1, which is frozen (plan TD-1): it keeps the
/// only format that it accepted before TBLASTX outfmt 0 and 7 were ported (session S08).
#[cfg(target_arch = "wasm32")]
pub fn tblastx_v1_outfmt(value: &str) -> Result<String, String> {
    if value.trim() != "6" {
        return Err(
            "unsupported TBLASTX outfmt: only 6 without custom fields is implemented".into(),
        );
    }
    Ok(value.into())
}

#[cfg(test)]
#[path = "value_parsers_tests.rs"]
mod tests;

/// The `-task` of a program: `accepted` names are the choices (help and `describe` list
/// them), a name in `rejected` is a task LOSAT does not implement and fails at the parse
/// with the given message (clap's error for a value its parser refuses, as every other
/// option's), and any other name is clap's invalid-value error.
///
/// NCBI's `-task` is a `CArgAllowStringSet` (`c++/src/algo/blast/blastinput/blast_args.cpp:105`); the rejected
/// names are NCBI tasks outside LOSAT's implementation.
#[derive(Clone)]
pub struct TaskParser {
    program: &'static str,
    accepted: &'static [&'static str],
    rejected: &'static [(&'static str, &'static str)],
}

/// The `-task` of `program` (its name in LOSAT's messages, `BLASTN`), `TaskParser`.
pub const fn task_parser(
    program: &'static str,
    accepted: &'static [&'static str],
    rejected: &'static [(&'static str, &'static str)],
) -> TaskParser {
    TaskParser {
        program,
        accepted,
        rejected,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:211-255
/// The tasks of NCBI's options factory (`GetTasks(eAll)`).
const NCBI_TASKS: &[&str] = &[
    "blastn",
    "blastn-short",
    "megablast",
    "dc-megablast",
    "vecscreen",
    "rmblastn",
    "blastp",
    "blastp-short",
    "blastp-fast",
    "psiblast",
    "phiblastp",
    "rpsblast",
    "rpstblastn",
    "blastx",
    "blastx-fast",
    "deltablast",
    "tblastn",
    "tblastn-fast",
    "psitblastn",
    "tblastx",
    "kblastp",
    "mapper",
    "mapr2g",
    "mapr2r",
    "mapg2g",
];

impl TaskParser {
    /// The check of a task that a caller of `run_local` set itself (the command line's parser
    /// makes it at parsing), as NCBI's options factory makes it: the name in lower case (the
    /// engine compares it exactly), then a task of `accepted`; the message of a `rejected` one;
    /// LOSAT's refusal for another task of NCBI's factory (another program's, which NCBI's API
    /// would set up); NCBI's error for a name NCBI does not know.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:336-342
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:662-673
    pub fn check(&self, task: &mut str) -> anyhow::Result<()> {
        task.make_ascii_lowercase();
        let task = &*task;
        if self.accepted.contains(&task) {
            return Ok(());
        }
        if let Some(&(_, message)) = self.rejected.iter().find(|(name, _)| *name == task) {
            anyhow::bail!("{message}");
        }
        if NCBI_TASKS.contains(&task) {
            let (last, rest) = self.accepted.split_last().expect("a program has a task");
            let tasks = if rest.is_empty() {
                last.to_string()
            } else {
                format!("{} or {last}", rest.join(", "))
            };
            anyhow::bail!(
                "the task {task} is not supported by LOSAT's {} (use {tasks})",
                self.program
            );
        }
        anyhow::bail!("'{task}' is not a supported task")
    }
}

impl clap::builder::TypedValueParser for TaskParser {
    type Value = String;

    fn parse_ref(
        &self,
        command: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<String, clap::Error> {
        let inner = clap::builder::PossibleValuesParser::new(self.accepted);
        if let Some(&(_, message)) = self
            .rejected
            .iter()
            .find(|(name, _)| value == std::ffi::OsStr::new(name))
        {
            let refuse = move |_: &str| -> Result<String, String> { Err(message.to_string()) };
            return clap::builder::TypedValueParser::parse_ref(&refuse, command, arg, value);
        }
        clap::builder::TypedValueParser::parse_ref(&inner, command, arg, value)
    }

    fn possible_values(
        &self,
    ) -> Option<Box<dyn Iterator<Item = clap::builder::PossibleValue> + '_>> {
        Some(Box::new(
            self.accepted
                .iter()
                .map(|name| clap::builder::PossibleValue::new(*name)),
        ))
    }
}

//! Shared CLI value grammar; engine configuration remains typed.
use crate::utils::seg::SegParams;
use std::path::PathBuf;

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:410-420
// ```c++
// if (dust_opts == kDfltArgNoFiltering) { opt.SetDustFiltering(false); }
// else if (dust_opts == kDfltArgApplyFiltering) { opt.SetDustFiltering(true); }
// else { x_TokenizeFilteringArgs(dust_opts, tokens);
//   opt.SetDustFilteringLevel(NStr::StringToInt(tokens[0]));
//   opt.SetDustFilteringWindow(NStr::StringToInt(tokens[1]));
//   opt.SetDustFilteringLinker(NStr::StringToInt(tokens[2])); }
// ```
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

    // NCBI core/blast_options.c:46-48,62-64:
    // const int kDustLevel = 20; const int kDustWindow = 64; const int kDustLinker = 1;
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

pub fn parse_dust_filtering(value: &str) -> Result<DustSpec, String> {
    match value {
        "no" => return Ok(DustSpec::No),
        "yes" => return Ok(DustSpec::Yes),
        _ => {}
    }
    let tokens: Vec<_> = value.split_whitespace().collect();
    if tokens.len() != 3 {
        return Err("DUST requires no, yes, or LEVEL WINDOW LINKER".into());
    }
    let level = nonnegative_i32(tokens[0])? as u32;
    let window = positive_usize(tokens[1])?;
    let linker = nonnegative_i32(tokens[2])? as usize;
    Ok(DustSpec::Parameters {
        level,
        window,
        linker,
    })
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:168-170,191-192,206-207,3162-3163
// ```c++
// new CArgAllowValuesGreaterThanOrEqual(1)
// new CArgAllow_Doubles(0.0, 100.0)
// ```
// CLI v2 additionally rejects all non-finite floating point values at entry.
pub fn positive_usize(value: &str) -> Result<usize, String> {
    let n = value
        .parse::<usize>()
        .map_err(|_| "expected a positive integer")?;
    if n == 0 || n > i32::MAX as usize {
        return Err("expected an integer in 1..=2147483647".into());
    }
    Ok(n)
}
pub fn nonnegative_f64(value: &str) -> Result<f64, String> {
    let n = value
        .parse::<f64>()
        .map_err(|_| "expected a finite number")?;
    if !n.is_finite() || n < 0.0 {
        return Err("expected a finite number >= 0".into());
    }
    Ok(n)
}
pub fn nonnegative_i32(value: &str) -> Result<i32, String> {
    let n = value
        .parse::<i32>()
        .map_err(|_| "expected an integer >= 0")?;
    if n < 0 {
        return Err("expected an integer >= 0".into());
    }
    Ok(n)
}
pub fn positive_i32(value: &str) -> Result<i32, String> {
    positive_usize(value).map(|n| n as i32)
}
/// An integer argument as NCBI reads it (`CArg_Integer`): a decimal number with an
/// optional sign, or else, for a value that starts with `0x` or `0X`, a hexadecimal number
/// (no digit after the prefix reads as 0), within the range of `int`.
///
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:118-129
/// ```c
/// Int8 s_StringToInt8(const string& value)
/// {
///     try {
///         return NStr::StringToInt8(value);
///     } catch (const CStringException&) {
///         if (NStr::StartsWith(value, "0x", NStr::eNocase)) {
///             return NStr::StringToInt8(value, 0, 16);
///         } else {
///             throw;
///         }
///     }
/// }
/// ```
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:383-390
/// ```c
/// inline CArg_Integer::CArg_Integer(const string& name, const string& value)
///     : CArg_Int8(name, value)
/// {
///     if (m_Integer < kMin_Int  ||  m_Integer > kMax_Int) {
///         NCBI_THROW(CArgException, eConvert, s_ArgExptMsg(GetName(),
///             "Integer value is out of range", value));
///     }
/// }
/// ```
/// NCBI reference: c++/src/corelib/ncbistr.cpp:788-791
/// ```c
///     // Remove leading '0x' for hex numbers
///     if ( base == 16 ) {
///         if (ch == '0'  &&  (next == 'x' || next == 'X')) {
///             pos += 2;
/// ```
/// `NStr::StringToInt8` in base 10 reads what `i64::from_str` reads (an optional sign and
/// ASCII digits); in base 16 it reads the hexadecimal digits after the prefix.
pub fn ncbi_integer(value: &str) -> Result<i32, String> {
    let int8 = match value.parse::<i64>() {
        Ok(n) => n,
        Err(_)
            if value
                .get(..2)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("0x")) =>
        {
            let digits = &value[2..];
            if !digits.bytes().all(|byte| byte.is_ascii_hexdigit()) {
                return Err("expected an integer".into());
            }
            if digits.is_empty() {
                0
            } else {
                i64::from_str_radix(digits, 16).map_err(|_| "Integer value is out of range")?
            }
        }
        Err(_) => return Err("expected an integer".into()),
    };
    i32::try_from(int8).map_err(|_| "Integer value is out of range".to_string())
}

/// A real argument as NCBI reads it (`CArg_Double`): the first character is a digit, a
/// point or a sign, and the rest is what `strtod` reads. LOSAT reads the decimal forms
/// and the signed infinities and NaN (`f64::from_str`), and rejects NCBI's other forms
/// (hexadecimal, an exponent mark without digits) with a message that names the program.
///
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:464
/// ```c
///         m_Double = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
/// ```
/// NCBI reference: c++/src/corelib/ncbistr.cpp:1313-1318
/// ```c
///     // Because strtod() may just skip such symbols.
///     if (!(flags & NStr::fAllowLeadingSymbols)) {
///         char c = str[pos];
///         if ( !isdigit((unsigned char)c)  &&  !s_IsDecimalPoint(c,flags)  &&  c != '-'  &&  c != '+') {
///             S2N_CONVERT_ERROR_INVAL(double);
///         }
/// ```
pub fn ncbi_double(value: &str, program: &str) -> Result<f64, String> {
    if !value
        .chars()
        .next()
        .is_some_and(|first| first.is_ascii_digit() || matches!(first, '.' | '-' | '+'))
    {
        return Err("expected a number".into());
    }
    value.parse::<f64>().map_err(|_| {
        format!("expected a decimal number (other forms that NCBI BLAST+ reads are not supported by LOSAT's {program})")
    })
}

/// A BLASTN integer argument of `CArg_Integer` with NCBI's lower bound (`at_least`).
fn blastn_integer_at_least(value: &str, at_least: i32) -> Result<i32, String> {
    let n = ncbi_integer(value)?;
    if n < at_least {
        return Err(format!("expected an integer >= {at_least}"));
    }
    Ok(n)
}

/// A BLASTN word size: 4 or more (the upper bound, 100, is an option check,
/// `blastn/scoring.rs`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:168-170
/// ```c
///         arg_desc.SetConstraint(kArgWordSize, m_QueryIsProtein
///                                ? new CArgAllowValuesGreaterThanOrEqual(2)
///                                : new CArgAllowValuesGreaterThanOrEqual(4));
/// ```
pub fn blastn_word_size(value: &str) -> Result<usize, String> {
    blastn_integer_at_least(value, 4).map(|n| n as usize)
}
/// A BLASTN count of 1 or more (`-num_threads`, `-max_target_seqs`, `-max_hsps`), as
/// NCBI's arguments (blast_args.cpp:203-207,2731-2732,3162-3163).
pub fn blastn_count(value: &str) -> Result<usize, String> {
    blastn_integer_at_least(value, 1).map(|n| n as usize)
}
/// A BLASTN reward: 0 or more, as NCBI's argument (blast_args.cpp:658-659). A reward of 0
/// is an option that LOSAT rejects before the search (`blastn/scoring.rs`).
pub fn blastn_reward(value: &str) -> Result<i32, String> {
    blastn_integer_at_least(value, 0)
}
/// A BLASTN penalty: 0 or less, as NCBI's argument (blast_args.cpp:651-652).
pub fn blastn_penalty(value: &str) -> Result<i32, String> {
    let n = ncbi_integer(value)?;
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
/// check rejects 0 or less, and LOSAT rejects infinity and NaN before the search
/// (`blastn/scoring.rs`).
pub fn blastn_evalue(value: &str) -> Result<f64, String> {
    ncbi_double(value, "BLASTN")
}
/// A BLASTN percent identity, in 0..=100 (NaN is outside).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:188-192
/// ```c
///         arg_desc.AddOptionalKey(kArgPercentIdentity, "float_value",
///                                 "Percent identity",
///                                 CArgDescriptions::eDouble);
///         arg_desc.SetConstraint(kArgPercentIdentity,
///                                new CArgAllow_Doubles(0.0, 100.0));
/// ```
pub fn blastn_percentage(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value, "BLASTN")?;
    if !(0.0..=100.0).contains(&n) {
        return Err("expected a percentage in 0..=100".into());
    }
    Ok(n)
}
pub fn blastp_word_size(value: &str) -> Result<usize, String> {
    let n = positive_usize(value)?;
    if n < 2 {
        return Err("protein word_size must be >= 2".into());
    }
    Ok(n)
}
// NCBI aa lookup word size is resolved before lookup construction; LOSAT's
// TBLASTX port currently constructs only 3-residue words (run_impl.rs).
pub fn tblastx_word_size(value: &str) -> Result<usize, String> {
    let n = positive_usize(value)?;
    if n != 3 {
        return Err("unsupported TBLASTX word_size: only 3 is implemented".into());
    }
    Ok(n)
}
// NCBI core/blast_options.c:1301-1308: threshold must be greater than zero.
pub fn positive_f64(value: &str) -> Result<f64, String> {
    let n = nonnegative_f64(value)?;
    if n == 0.0 {
        return Err("expected a finite number > 0".into());
    }
    Ok(n)
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:995-1006
// ```c++
// static const set<int> genetic_codes(gcs, gcs+sizeof(gcs)/sizeof(*gcs));
// return (genetic_codes.find(val) != genetic_codes.end());
// ```
pub fn genetic_code(value: &str) -> Result<u8, String> {
    let n = value.parse::<u8>().map_err(|_| "invalid genetic code")?;
    if !matches!(n, 1..=6 | 9..=16 | 21..=31 | 33) {
        return Err("unsupported genetic code".into());
    }
    Ok(n)
}

// NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-50
// ```c++
// const string kArgQuery("query");
// const string kArgSubject("subject");
// ```
/// A BLASTN input file: any value, as NCBI's `eInputFile` argument; `-` is standard input
/// and a file that cannot be opened is reported when it is opened (`blastn/blast_engine/run.rs`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
/// ```c
///     arg_desc.AddDefaultKey(kArgQuery, "input_file",
///                      "Input file name",
///                      CArgDescriptions::eInputFile, kDfltArgQuery);
/// ```
pub fn blastn_input_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().map(PathBuf::from)
}

// LOSAT CLI v2 capability boundary: both inputs are required file paths.
pub fn file_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().try_map(|value| {
        if value.is_empty() || value == "-" {
            return Err("a file path is required; stdin is not implemented".into());
        }
        Ok::<_, String>(PathBuf::from(value))
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:386-407
// ```c
// if (m_QueryIsProtein && args[kArgSegFiltering]) {
//     const string& seg_opts = args[kArgSegFiltering].AsString();
//     if (seg_opts == kDfltArgNoFiltering) {
//         opt.SetSegFiltering(false);
//     } else if (seg_opts == kDfltArgApplyFiltering) {
//         opt.SetSegFiltering(true);
//     } else {
//         x_TokenizeFilteringArgs(seg_opts, tokens);
//         opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
//         opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
//         opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
//     }
// }
// ```
#[derive(Debug, Clone, PartialEq)]
pub enum SegSpec {
    No,
    Yes,
    WindowLocutHicut {
        window: usize,
        locut: f64,
        hicut: f64,
    },
}

impl SegSpec {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    #[inline]
    pub fn params(&self) -> Option<SegParams> {
        match self {
            Self::No => None,
            Self::Yes => Some(SegParams::default()),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => Some(SegParams::new(*window, *locut, *hicut)),
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:396-407
    // ```c
    // if (seg_opts == kDfltArgNoFiltering) {
    //     opt.SetSegFiltering(false);
    // } else if (seg_opts == kDfltArgApplyFiltering) {
    //     opt.SetSegFiltering(true);
    // } else {
    //     x_TokenizeFilteringArgs(seg_opts, tokens);
    // }
    // ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:367-381,396-407
// ```c
// NStr::Split(filtering_args, " ", output);
// if (output.size() != 3) {
//     NCBI_THROW(CInputException, eInvalidInput,
//                "Invalid number of arguments to filtering option");
// }
// ...
// if (seg_opts == kDfltArgNoFiltering) {
//     opt.SetSegFiltering(false);
// } else if (seg_opts == kDfltArgApplyFiltering) {
//     opt.SetSegFiltering(true);
// } else {
//     x_TokenizeFilteringArgs(seg_opts, tokens);
//     opt.SetSegFilteringWindow(...);
//     opt.SetSegFilteringLocut(...);
//     opt.SetSegFilteringHicut(...);
// }
// ```
pub fn parse_seg_filtering(value: &str) -> Result<SegSpec, String> {
    if value == "no" {
        return Ok(SegSpec::No);
    }
    if value == "yes" {
        return Ok(SegSpec::Yes);
    }

    let tokens: Vec<&str> = value.split_whitespace().collect();
    if tokens.len() != 3 {
        return Err("invalid number of arguments to filtering option".to_string());
    }

    let window = tokens[0]
        .parse::<usize>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;
    let locut = tokens[1]
        .parse::<f64>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;
    let hicut = tokens[2]
        .parse::<f64>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;

    if window == 0 || window > i32::MAX as usize || !locut.is_finite() || !hicut.is_finite() {
        return Err("SEG requires window > 0 and finite locut/hicut".into());
    }
    // NCBI blast_seg.c:2247-2258 normalizes finite negative/inverted cutoffs.
    // Keep that behavior in SegParams::new at the typed configuration boundary.
    Ok(SegSpec::WindowLocutHicut {
        window,
        locut,
        hicut,
    })
}

// NCBI blast_args.cpp:2657-2660: AddDefaultKey(kArgOutputFormat, ..., eString, ...).
// CLI v2 exposes only the formatter capabilities actually implemented by LOSAT.
pub fn tblastx_outfmt(value: &str) -> Result<String, String> {
    if value.trim() != "6" {
        return Err(
            "unsupported TBLASTX outfmt: only 6 without custom fields is implemented".into(),
        );
    }
    Ok(value.into())
}

// NCBI blast_args.cpp:475-488: window size is an Int4, with zero selecting one-hit mode.
pub fn nonnegative_usize(value: &str) -> Result<usize, String> {
    nonnegative_i32(value).map(|n| n as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integers_are_read_as_ncbi_reads_them() {
        for (value, expected) in [
            ("7", 7),
            ("+7", 7),
            ("-07", -7),
            ("0x1F", 31),
            ("0XaB", 171),
            ("0x", 0),
            ("2147483647", i32::MAX),
            ("-2147483648", i32::MIN),
        ] {
            assert_eq!(ncbi_integer(value), Ok(expected), "{value}");
        }
        for value in [
            "",
            " 7",
            "7 ",
            "+",
            "-0x5",
            "0x+5",
            "0x5g",
            "1_000",
            "2147483648",
            "0x80000000",
            "00x5",
            "7.0",
        ] {
            assert!(ncbi_integer(value).is_err(), "{value:?}");
        }
    }

    #[test]
    fn reals_follow_ncbi_first_character_rule() {
        for (value, expected) in [("1e-5", 1e-5), (".5", 0.5), ("+2", 2.0), ("-0", 0.0)] {
            assert_eq!(ncbi_double(value, "BLASTN"), Ok(expected), "{value}");
        }
        assert!(ncbi_double("+inf", "BLASTN").unwrap().is_infinite());
        assert!(ncbi_double("-nan", "BLASTN").unwrap().is_nan());
        for value in ["inf", "nan", " 1", "", "e5"] {
            assert_eq!(
                ncbi_double(value, "BLASTN"),
                Err("expected a number".into()),
                "{value:?}"
            );
        }
        for value in ["0x10", "0x1p-3", "1e"] {
            assert!(ncbi_double(value, "BLASTN")
                .unwrap_err()
                .contains("not supported by LOSAT's BLASTN"));
        }
    }
}

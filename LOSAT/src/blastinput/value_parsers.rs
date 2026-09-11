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
pub fn percentage(value: &str) -> Result<f64, String> {
    let n = nonnegative_f64(value)?;
    if n > 100.0 {
        return Err("expected a percentage in 0..=100".into());
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
pub fn negative_i32(value: &str) -> Result<i32, String> {
    let n = value
        .parse::<i32>()
        .map_err(|_| "expected a negative integer")?;
    if n >= 0 {
        return Err("expected a negative integer".into());
    }
    Ok(n)
}
pub fn blastn_word_size(value: &str) -> Result<usize, String> {
    let n = positive_usize(value)?;
    // NCBI core/blast_options.c:1322-1334; include/algo/blast/core/blast_hits.h:192:
    // options->word_size < 4 || options->word_size > DBSEQ_CHUNK_OVERLAP
    // #define DBSEQ_CHUNK_OVERLAP 100
    if !(4..=100).contains(&n) {
        return Err("BLASTN word_size must be in 4..=100".into());
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
pub fn blastn_outfmt(value: &str) -> Result<String, String> {
    if value.trim().is_empty() {
        return Err("outfmt cannot be empty".into());
    }
    crate::algorithm::blastn::hsp::parse_blastn_output_format(value)?;
    Ok(value.into())
}
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

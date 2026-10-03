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

/// A `-dust` value as NCBI reads it: `no`, `yes`, or three signed decimal numbers separated
/// by single spaces. The error is NCBI's message; the masker replaces values out of its
/// ranges with its defaults (`CSymDustMasker`, `utils/dust.rs`), after NCBI's conversion
/// to `Uint4` (`as u32`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:375-384
/// ```c
/// CFilteringArgs::x_TokenizeFilteringArgs(const string& filtering_args,
///                                         vector<string>& output) const
/// {
///     output.clear();
///     NStr::Split(filtering_args, " ", output);
///     if (output.size() != 3) {
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "Invalid number of arguments to filtering option");
///     }
/// }
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:410-426
/// ```c
///         if ( !m_QueryIsProtein && args[kArgDustFiltering]) {
///             const string& dust_opts = args[kArgDustFiltering].AsString();
///             if (dust_opts == kDfltArgNoFiltering) {
///                 opt.SetDustFiltering(false);
///             } else if (dust_opts == kDfltArgApplyFiltering) {
///                 opt.SetDustFiltering(true);
///             } else {
///                 x_TokenizeFilteringArgs(dust_opts, tokens);
///                 opt.SetDustFilteringLevel(NStr::StringToInt(tokens[0]));
///                 opt.SetDustFilteringWindow(NStr::StringToInt(tokens[1]));
///                 opt.SetDustFilteringLinker(NStr::StringToInt(tokens[2]));
///             }
///         }
///     } catch (const CStringException& e) {
///         if (e.GetErrCode() == CStringException::eConvert) {
///             NCBI_THROW(CInputException, eInvalidInput,
///                        "Invalid input for filtering parameters");
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_objmgr_tools.cpp:195-198
/// ```c
///                 Blast_FindDustFilterLoc(*m_QueryVector,
///                     static_cast<Uint4>(m_Options->GetDustFilteringLevel()),
///                     static_cast<Uint4>(m_Options->GetDustFilteringWindow()),
///                     static_cast<Uint4>(m_Options->GetDustFilteringLinker()));
/// ```
/// `NStr::Split` keeps the empty fields between adjacent spaces, and `NStr::StringToInt`
/// reads what `i32::from_str` reads.
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
    let number = |token: &str| {
        token
            .parse::<i32>()
            .map_err(|_| "Invalid input for filtering parameters".to_string())
    };
    let (level, window, linker) = (number(tokens[0])?, number(tokens[1])?, number(tokens[2])?);
    Ok(DustSpec::Parameters {
        level: level as u32,
        window: window as u32 as usize,
        linker: linker as u32 as usize,
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
/// point or a sign, and the rest is what `strtod` reads. LOSAT reads the decimal forms,
/// the signed infinities and NaN (`f64::from_str`) and glibc's `nan(n-char-sequence)`
/// (letters, digits and `_`), and rejects NCBI's other forms (hexadecimal, an exponent
/// mark without digits) with a message that names the program.
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
    // glibc's strtod reads `nan(...)` as NaN; NCBI then checks that it ended at the end of
    // the string (ncbistr.cpp:1332-1376).
    let unsigned = value.trim_start_matches(['+', '-']);
    if unsigned.len() + 1 == value.len()
        && unsigned.len() >= 5
        && unsigned[..4].eq_ignore_ascii_case("nan(")
        && unsigned.ends_with(')')
        && unsigned[4..unsigned.len() - 1]
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
    {
        return Ok(f64::NAN);
    }
    value.parse::<f64>().map_err(|_| {
        format!("expected a decimal number (other forms, which NCBI BLAST+ may read, are not supported by LOSAT's {program})")
    })
}

/// An integer argument with a constraint: NCBI's constraint reads the value again with
/// `NStr::StringToDouble`, which does not read a `0x` prefix without digits (the only form
/// that `ncbi_integer` reads and `strtod` does not end at), so the value is illegal.
///
/// NCBI reference: c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// ```c
///     /// Overloaded method from CArgAllow
///     virtual bool Verify(const string& value) const {
///         return NStr::StringToDouble(value) >= m_MinValue;
///     }
/// ```
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:1231-1243
/// ```c
///     if ( m_Constraint ) {
///         bool err = false;
///         try {
///             bool check = m_Constraint->Verify(value);
///     ...
///         } catch (...) {
///             err = true;
///         }
/// ```
fn ncbi_constrained_integer(value: &str) -> Result<i32, String> {
    if value.eq_ignore_ascii_case("0x") {
        return Err("Illegal value".into());
    }
    ncbi_integer(value)
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

/// A non-negative integer argument (`CArg_Integer` with `CArgAllowValuesGreaterThanOrEqual(0)`),
/// such as `-culling_limit`.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3296-3302
/// ```c
///     arg_desc.AddOptionalKey(kArgCullingLimit, "int_value",
///                      "If the query range of a hit is enveloped by that of at "
///                      "least this many higher-scoring hits, delete the hit",
///                      CArgDescriptions::eInteger);
///     arg_desc.SetConstraint(kArgCullingLimit,
///     // best hit algorithm arguments
///                new CArgAllowValuesGreaterThanOrEqual(kDfltArgCullingLimit));
/// ```
pub fn nonnegative_ncbi_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 0)
}

/// A BLASTN task: NCBI's blastn tasks (case-sensitive), of which LOSAT implements
/// megablast, blastn, dc-megablast and blastn-short and rejects rmblastn (matrix scoring
/// and masklevel) explicitly.
///
/// NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:211-222
/// ```c
/// CBlastOptionsFactory::GetTasks(ETaskSets choice /* = eAll */)
/// {
///     set<string> retval;
///     if (choice == eNuclNucl || choice == eAll) {
///         retval.insert("blastn");
///         retval.insert("blastn-short");
///         retval.insert("megablast");
///         retval.insert("dc-megablast");
///         retval.insert("vecscreen");
///         // -RMH-
///         retval.insert("rmblastn");
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blastn_args.cpp:57-59
/// ```c
///     set<string> tasks
///         (CBlastOptionsFactory::GetTasks(CBlastOptionsFactory::eNuclNucl));
///     tasks.erase("vecscreen"); // vecscreen has its own program
/// ```
pub fn blastn_task(value: &str) -> Result<String, String> {
    match value {
        "megablast" | "blastn" | "dc-megablast" | "blastn-short" => Ok(value.to_string()),
        "rmblastn" => Err(format!(
            "the task {value} is not supported by LOSAT's BLASTN (use megablast, blastn, dc-megablast or blastn-short)"
        )),
        _ => Err("expected one of blastn, blastn-short, dc-megablast, megablast, rmblastn".into()),
    }
}
/// A discontiguous megablast template type: `coding`, `optimal` or `coding_and_optimal`.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:688-693,708-714
/// ```c
/// /// Value to specify coding template type
/// const char* kTemplType_Coding = "coding";
/// /// Value to specify optimal template type
/// const char* kTemplType_Optimal = "optimal";
/// /// Value to specify coding+optimal template type
/// const char* kTemplType_CodingAndOptimal = "coding_and_optimal";
/// ...
///     arg_desc.AddOptionalKey(kArgDMBTemplateType, "type",
///                  "Discontiguous MegaBLAST template type",
///                  CArgDescriptions::eString);
///     arg_desc.SetConstraint(kArgDMBTemplateType, &(*new CArgAllow_Strings,
///                                                   kTemplType_Coding,
///                                                   kTemplType_Optimal,
///                                                   kTemplType_CodingAndOptimal));
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:747-752
/// ```c
///         if (type == kTemplType_Coding) {
///             temp_type = eMBWordCoding;
///         } else if (type == kTemplType_Optimal) {
///             temp_type = eMBWordOptimal;
///         } else if (type == kTemplType_CodingAndOptimal) {
///             temp_type = eMBWordTwoTemplates;
/// ```
pub fn blastn_template_type(
    value: &str,
) -> Result<crate::algorithm::blastn::disc_lookup::DiscWordType, String> {
    use crate::algorithm::blastn::disc_lookup::DiscWordType;
    match value {
        "coding" => Ok(DiscWordType::Coding),
        "optimal" => Ok(DiscWordType::Optimal),
        "coding_and_optimal" => Ok(DiscWordType::TwoTemplates),
        _ => Err("expected one of coding, coding_and_optimal, optimal".into()),
    }
}
/// A discontiguous megablast template length: 16, 18 or 21.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:719-727
/// ```c
///     arg_desc.AddOptionalKey(kArgDMBTemplateLength, "int_value",
///                  "Discontiguous MegaBLAST template length",
///                  CArgDescriptions::eInteger);
///     set<int> allowed_values;
///     allowed_values.insert(16);
///     allowed_values.insert(18);
///     allowed_values.insert(21);
///     arg_desc.SetConstraint(kArgDMBTemplateLength,
///                            new CArgAllowIntegerSet(allowed_values));
/// ```
/// The constraint converts the value again with `NStr::StringToInt` in base 10, so a
/// hexadecimal value that the integer argument reads fails it.
///
/// NCBI reference: c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-222,239
/// ```c
///     virtual bool Verify(const string& value) const {                        \
///         DataType value2check = String2DataTypeFn(value);                    \
///         ITERATE(set<DataType>, itr, m_AllowedValues) {                      \
///             if (*itr == value2check) {                                      \
///                 return true;                                                \
///             }                                                               \
///         }                                                                   \
///         return false;                                                       \
///     }                                                                       \
/// ...
/// DEFINE_CARGALLOW_SET_CLASS(CArgAllowIntegerSet, int, NStr::StringToInt);
/// ```
/// `NStr::StringToInt` in base 10 reads what `i32::from_str` reads (an optional sign and
/// ASCII digits, `ncbi_integer`).
pub fn blastn_template_length(value: &str) -> Result<u8, String> {
    ncbi_constrained_integer(value)?;
    match value.parse::<i32>() {
        Ok(n @ (16 | 18 | 21)) => Ok(n as u8),
        _ => Err("expected one of 16, 18, 21".into()),
    }
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
    ncbi_integer_at_least(value, 4).map(|n| n as usize)
}
/// A BLASTN count of 1 or more (`-num_threads`, `-max_target_seqs`, `-max_hsps`), as
/// NCBI's arguments (blast_args.cpp:203-207,2731-2732,3162-3163).
pub fn blastn_count(value: &str) -> Result<usize, String> {
    ncbi_integer_at_least(value, 1).map(|n| n as usize)
}
/// A BLASTN reward: 0 or more, as NCBI's argument (blast_args.cpp:658-659). A reward of 0
/// is an option that LOSAT rejects before the search (`blastn/scoring.rs`).
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
/// A BLASTN output file: a name shorter than 256 bytes, as NCBI's constraint; `-` is
/// standard output, and a file that cannot be created is reported when it is opened
/// (`blastn/blast_engine/run.rs`).
///
/// NCBI reference: c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
/// ```c
///     static constexpr Uint4 kDfltMaxLength = 256;
///
///     CArgAllowMaximumFileNameLength(Uint4 max = kDfltMaxLength) : m_MaxLength(max) {}
///
/// protected:
///     /// Overloaded method from CArgAllow
///     virtual bool Verify(const string& value) const {
///         CFile fname(value);
///         return fname.GetName().size() < m_MaxLength;
/// ```
pub fn blastn_output_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().try_map(|value| {
        if ncbi_file_name_length(value.as_encoded_bytes()) >= 256 {
            return Err("Illegal value, expected file name length < 256".to_string());
        }
        Ok(PathBuf::from(value))
    })
}

/// The length of NCBI's `CDirEntry::GetName` of a path: the text after the last separator
/// once the trailing directory separators are removed (a path of one separator has none).
///
/// NCBI reference: c++/src/corelib/ncbifile.cpp:298-312
/// ```c
/// void CDirEntry::Reset(const string& path)
/// {
///     m_Path = path;
///     size_t len = path.length();
///     // Root dir
///     if ((len == 1)  &&  IsPathSeparator(path[0])) {
///         return;
///     }
///     // Disk name
/// #  if defined(DISK_SEPARATOR)
///     if ( (len == 2 || len == 3) && (path[1] == DISK_SEPARATOR) ) {
///         return;
///     }
/// #  endif
///     m_Path = DeleteTrailingPathSeparator(path);
/// ```
/// NCBI reference: c++/src/corelib/ncbifile.cpp:465-472
/// ```c
/// string CDirEntry::DeleteTrailingPathSeparator(const string& path)
/// {
///     size_t pos = path.find_last_not_of(DIR_SEPARATORS);
///     if (pos + 1 < path.length()) {
///         return path.substr(0, pos + 1);
///     }
///     return path;
/// }
/// ```
/// NCBI reference: c++/src/corelib/ncbifile.cpp:358-363
/// ```c
/// void CDirEntry::SplitPath(const string& path, string* dir,
///                           string* base, string* ext)
/// {
///     // Get file name
///     size_t pos = path.find_last_of(ALL_SEPARATORS);
///     string filename = (pos == NPOS) ? path : path.substr(pos+1);
/// ```
/// On Unix, `DIR_SEPARATORS` and `ALL_SEPARATORS` are `/` (ncbifile.cpp:110-115); on
/// Windows they are `/\` and `:/\`, with the disk name `C:`.
fn ncbi_file_name_length(path: &[u8]) -> usize {
    let (dir_separators, all_separators): (&[u8], &[u8]) = if cfg!(windows) {
        (b"/\\", b":/\\")
    } else {
        (b"/", b"/")
    };
    let root = path.len() == 1 && dir_separators.contains(&path[0]);
    let disk = cfg!(windows) && matches!(path.len(), 2 | 3) && path[1] == b':';
    let path = if root || disk {
        path
    } else {
        let end = path
            .iter()
            .rposition(|byte| !dir_separators.contains(byte))
            .map_or(0, |pos| pos + 1);
        &path[..end]
    };
    path.iter()
        .rposition(|byte| all_separators.contains(byte))
        .map_or(path.len(), |pos| path.len() - pos - 1)
}

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
    WindowLocutHicut { window: i32, locut: f64, hicut: f64 },
}

impl SegSpec {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    /// The SEG parameters of the query filter: a window, locut or hicut that is not above
    /// 0 keeps NCBI's default (`SegParametersNewAa`), then SEG checks the parameters
    /// (`SegParams::new`, blast_seg.c:2247-2258).
    ///
    /// NCBI reference: c++/src/algo/blast/core/blast_filter.c:1147-1154
    /// ```c
    ///         sparamsp = SegParametersNewAa();
    ///         sparamsp->overlaps = TRUE;
    ///         if (seg_options->window > 0)
    ///             sparamsp->window = seg_options->window;
    ///         if (seg_options->locut > 0.0)
    ///             sparamsp->locut = seg_options->locut;
    ///         if (seg_options->hicut > 0.0)
    ///             sparamsp->hicut = seg_options->hicut;
    /// ```
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

/// A `-seg` value as NCBI reads it: `no`, `yes`, or a window (an `int`), a locut and a
/// hicut separated by single spaces (`x_TokenizeFilteringArgs`, see `parse_dust_filtering`).
/// The errors are NCBI's messages; LOSAT also rejects a locut or hicut that is not finite.
/// The values are used as NCBI's query filter uses them (`SegSpec::params`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:396-406
/// ```c
///         if (m_QueryIsProtein && args[kArgSegFiltering]) {
///             const string& seg_opts = args[kArgSegFiltering].AsString();
///             if (seg_opts == kDfltArgNoFiltering) {
///                 opt.SetSegFiltering(false);
///             } else if (seg_opts == kDfltArgApplyFiltering) {
///                 opt.SetSegFiltering(true);
///             } else {
///                 x_TokenizeFilteringArgs(seg_opts, tokens);
///                 opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
///                 opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
///                 opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
///             }
/// ```
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

// NCBI blast_args.cpp:2657-2660: AddDefaultKey(kArgOutputFormat, ..., eString, ...).
// CLI v2 exposes only the formatter capabilities actually implemented by LOSAT.
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
// ```c
//     if (args[kArgOutputFormat]) {
//         string fmt_choice =
//             NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
// ```
/// A TBLASTX `-outfmt` value: the pairwise report (0) or the tabular formats (6, 7),
/// without a custom field list. The other formats are not implemented and fail.
pub fn tblastx_outfmt(value: &str) -> Result<String, String> {
    if !matches!(value.trim(), "0" | "6" | "7") {
        return Err(
            "unsupported TBLASTX outfmt: only 0, 6 and 7 without custom fields are implemented"
                .into(),
        );
    }
    Ok(value.into())
}

/// The TBLASTX `-outfmt` value of web ABI v1, which is frozen (plan TD-1): it keeps the
/// only format that it accepted before TBLASTX outfmt 0 and 7 were ported (session S08).
pub fn tblastx_v1_outfmt(value: &str) -> Result<String, String> {
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
        // NCBI 2.17.0 searches with each of these (-evalue, E2g).
        for value in ["+nan(1)", "-NaN()", "+nan(x_9)", "+Infinity", "1e999"] {
            let read = ncbi_double(value, "BLASTN").unwrap();
            assert!(read.is_nan() || read.is_infinite(), "{value}");
        }
        for value in [
            "+nan(", "+nan(1", "+nan(-)", "+nan(1)x", "+-nan(1)", "nan(1)",
        ] {
            assert!(ncbi_double(value, "BLASTN").is_err(), "{value}");
        }
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

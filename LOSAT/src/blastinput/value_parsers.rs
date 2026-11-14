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
/// ```c++
/// inline CArg_Int8::CArg_Int8(const string& name, const string& value)
///     : CArg_String(name, value)
/// {
///     try {
///         m_Integer = s_StringToInt8(value);
///     } catch (const CException& e) {
///         NCBI_RETHROW(e, CArgException, eConvert, s_ArgExptMsg(GetName(),
///             "Argument cannot be converted", value));
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129
/// ```c++
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
/// `NStr::StringToInt8(value)` and `NStr::StringToInt8(value, 0, 16)` are
/// `ncbi_str::string_to_int8` and `ncbi_str::string_to_int8_base`. LOSAT's messages: "expected
/// an integer" for a value that neither reads (a byte after `0x` that is not a hexadecimal
/// digit included, which is tested first: NCBI's base 16 conversion throws on an overflow
/// before it reaches such a byte, ncbistr.cpp:852), "Integer value is out of range" for a
/// hexadecimal number beyond `Int8`.
pub fn ncbi_int8(value: &str) -> Result<i64, String> {
    use crate::ncbi_str::{string_to_int8, string_to_int8_base};
    let value = value.as_bytes();
    match string_to_int8(value, 0) {
        Ok(n) => Ok(n),
        Err(_)
            if value
                .get(..2)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"0x")) =>
        {
            if !value[2..].iter().all(u8::is_ascii_hexdigit) {
                return Err("expected an integer".into());
            }
            string_to_int8_base(value, 0, 16)
                .map_err(|_| "Integer value is out of range".to_string())
        }
        Err(_) => Err("expected an integer".into()),
    }
}

/// A real argument as NCBI reads it (`CArg_Double`): `NStr::StringToDouble(value,
/// fDecimalPosixOrLocal)` (`ncbi_str::string_to_double_posix_or_local`): the first byte is a
/// digit, a point, a comma or a sign, glibc's `strtod` reads the number (decimal, hexadecimal,
/// the infinities, NaN with its sign and `n-char-sequence`) or NCBI's own `StringToDoublePosix`
/// where it reads further (a mantissa that ends in `e` or `E`, in x87 extended precision), and
/// the number ends the value. LOSAT's message for NCBI's conversion error: "expected a number".
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
/// ```c++
/// inline CArg_Double::CArg_Double(const string& name, const string& value)
///     : CArg_String(name, value)
/// {
///     try {
///         m_Double = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
///     } catch (const CException& e) {
///         NCBI_RETHROW(e,CArgException,eConvert,
///             s_ArgExptMsg(GetName(),"Argument cannot be converted",value));
///     }
/// }
/// ```
pub fn ncbi_double(value: &str) -> Result<f64, String> {
    crate::ncbi_str::string_to_double_posix_or_local(value.as_bytes())
        .map_err(|_| "expected a number".to_string())
}

/// An integer argument with a constraint: after the conversion, NCBI's constraint reads the
/// value again with `NStr::StringToDouble` (default flags, `ncbi_str::string_to_double`), which
/// does not read a `0x` prefix without digits (the only form that `ncbi_integer` reads and
/// `strtod` does not end at), so the value is illegal.
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
/// ```c
///     arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
///                             "Effective length of the search space",
///                             CArgDescriptions::eInt8);
///     arg_desc.SetConstraint(kArgEffSearchSpace,
///                            new CArgAllowValuesGreaterThanOrEqual(0));
/// ```
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

/// A `-num_descriptions` or `-num_alignments` value: 0 or more, an `Int4` in NCBI's
/// integer forms (the sign and hexadecimal forms of `nonnegative_ncbi_integer`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685
/// ```c
///     arg_desc.AddOptionalKey(kArgNumDescriptions, "int_value",
///                  "Number of database sequences to show one-line "
///                  "descriptions for\n"
///    		     "Not applicable for outfmt > 4\n"
///    			 "Default = `"+ NStr::IntToString(m_DfltNumDescriptions)+ "'",
///                  CArgDescriptions::eInteger);
///     arg_desc.SetConstraint(kArgNumDescriptions,
///                            new CArgAllowValuesGreaterThanOrEqual(0));
///
///     // number of alignments per DB sequence
///     arg_desc.AddOptionalKey(kArgNumAlignments, "int_value",
///                  "Number of database sequences to show alignments for\n"
///                  "Default = `" + NStr::IntToString(m_DfltNumAlignments) + "'",
///                  CArgDescriptions::eInteger );
///     arg_desc.SetConstraint(kArgNumAlignments,
///                            new CArgAllowValuesGreaterThanOrEqual(0));
/// ```
pub fn report_count(value: &str) -> Result<usize, String> {
    nonnegative_ncbi_integer(value).map(|n| n as usize)
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
    ncbi_positive_integer(value).map(|n| n as usize)
}
/// `blastn_count` as an `int`, the type of BLASTX's fields.
pub fn ncbi_positive_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 1)
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
/// check rejects 0 or less (`blastn/scoring.rs`); an infinite or NaN e-value searches as in
/// NCBI (`+inf`, `-nan`, `1e999`; `blastn/ncbi_cutoffs.rs`, E2g).
pub fn blastn_evalue(value: &str) -> Result<f64, String> {
    ncbi_double(value)
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
    percentage(value)
}
/// A percentage, in 0..=100 (NaN is outside): `-perc_identity`, and
/// `-qcov_hsp_perc` of every program (the `*_query_coverage` parsers).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:196-201
/// ```c
///         arg_desc.AddOptionalKey(kArgQueryCovHspPerc, "float_value",
///                             "Percent query coverage per hsp",
///                             CArgDescriptions::eDouble);
///         arg_desc.SetConstraint(kArgQueryCovHspPerc,
///                            new CArgAllow_Doubles(0.0, 100.0));
/// ```
/// The constraint reads the value with the argument's own conversion (`ncbi_double`), so it
/// tests the argument's number (`1e` is 1, `0x10` is 16).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
/// ```c++
/// bool CArgAllow_Doubles::Verify(const string& value) const
/// {
///     double val = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
///     ITERATE( set< TInterval >, pi, m_MinMax) {
///         if (pi->first <= val && val<= pi->second) {
///             return true;
///         }
///     }
///     return false;
/// }
/// ```
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

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:995-1006
// ```c++
// static const set<int> genetic_codes(gcs, gcs+sizeof(gcs)/sizeof(*gcs));
// return (genetic_codes.find(val) != genetic_codes.end());
// ```
pub fn genetic_code(value: &str) -> Result<u8, String> {
    let n = genetic_code_number(value).ok_or("invalid genetic code")?;
    if !matches!(n, 1..=6 | 9..=16 | 21..=31 | 33) {
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
/// ```c++
/// class CArgAllowGeneticCodeInteger : public CArgAllow
/// {
/// protected:
///      /// Overloaded method from CArgAllow
///      virtual bool Verify(const string& value) const {
///          static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33};
///          static const set<int> genetic_codes(gcs, gcs+sizeof(gcs)/sizeof(*gcs));
///          const int val = NStr::StringToInt(value);
///          return (genetic_codes.find(val) != genetic_codes.end());
///      }
/// ```
pub(crate) fn genetic_code_number(value: &str) -> Option<u8> {
    match crate::ncbi_str::string_to_int(value.as_bytes(), 0) {
        Ok(n) if !value.starts_with('-') => u8::try_from(n).ok(),
        _ => None,
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-51
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:578-583
/// ```c
///     arg_desc.AddOptionalKey(kArgWordScoreThreshold, "float_value",
///                  "Minimum word score such that the word is added to the "
///                  "BLAST lookup table",
///                  CArgDescriptions::eDouble);
///     arg_desc.SetConstraint(kArgWordScoreThreshold,
///                            new CArgAllowValuesGreaterThanOrEqual(0));
/// ```
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// ```c++
///     /// Overloaded method from CArgAllow
///     virtual bool Verify(const string& value) const {
///         return NStr::StringToDouble(value) >= m_MinValue;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
/// ```c++
///     if ( m_Constraint ) {
///         bool err = false;
///         try {
///             bool check = m_Constraint->Verify(value);
///     ...
///         } catch (...) {
///             err = true;
///         }
/// ```
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
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:489-497
/// ```c
/// inline CArg_Boolean::CArg_Boolean(const string& name, const string& value)
///     : CArg_String(name, value)
/// {
///     try {
///         m_Boolean = NStr::StringToBool(value);
///     } catch (const CException& e) {
///         NCBI_RETHROW(e,CArgException,eConvert, s_ArgExptMsg(GetName(),
///             "Argument cannot be converted",value));
///     }
/// ```
pub fn ncbi_boolean(value: &str) -> Result<bool, String> {
    crate::ncbi_str::string_to_bool(value)
        .ok_or_else(|| format!("Argument cannot be converted: `{value}'"))
}

/// A `-mt_mode` value: an integer from 0 to 2 (`CArgAllowValuesBetween(0, 2, true)`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3177
/// ```c++
///     	arg_desc.SetConstraint(kArgMTMode,
///                                new CArgAllowValuesBetween(0, 2, true));
/// ```
pub fn ncbi_mt_mode(value: &str) -> Result<u8, String> {
    let mode = ncbi_integer_at_least(value, 0)?;
    u8::try_from(mode)
        .ok()
        .filter(|mode| *mode <= 2)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 2".to_string())
}

/// A protein `-word_size` (`CArgAllowValuesGreaterThanOrEqual(2)` for a protein query).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-170
/// ```c
///         arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
///                                 CArgDescriptions::eInteger);
///         arg_desc.SetConstraint(kArgWordSize, m_QueryIsProtein
///                                ? new CArgAllowValuesGreaterThanOrEqual(2)
///                                : new CArgAllowValuesGreaterThanOrEqual(4));
/// ```
pub fn protein_word_size(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 2)
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

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129 (`s_StringToInt8`)
    // An Int8 argument (`-dbsize`) reads what an integer argument reads, in 64 bits.
    #[test]
    fn int8_arguments_are_read_as_ncbi_reads_them() {
        for (value, expected) in [
            ("7", 7),
            ("+7", 7),
            ("-07", -7),
            ("0", 0),
            ("0x1F", 31),
            ("0XaB", 171),
            ("0x", 0),
            ("3000000000", 3_000_000_000),
            ("-3000000000", -3_000_000_000),
            ("9223372036854775807", i64::MAX),
            ("-9223372036854775808", i64::MIN),
            ("0x7FFFFFFFFFFFFFFF", i64::MAX),
        ] {
            assert_eq!(ncbi_int8(value), Ok(expected), "{value}");
        }
        for value in [
            "",
            " 7",
            "7 ",
            "+",
            "-0x5",
            "0x5g",
            "1_000",
            "9223372036854775808",
            "-9223372036854775809",
            "0x8000000000000000",
            "7.0",
            "1e3",
            "abc",
        ] {
            assert!(ncbi_int8(value).is_err(), "{value:?}");
        }
    }

    // NCBI 2.17.0's `CArg_Double`: the bits of `NStr::StringToDouble(value,
    // fDecimalPosixOrLocal)` from NCBI's corelib (swebok-improve-20261105/sw8/nfix/probe/x87,
    // `s2dp`), and `blastn -evalue <value>` (exit 0, or 1 with the usage for the errors).
    #[test]
    fn reals_are_read_as_ncbis_carg_double() {
        for (value, bits) in [
            ("1e-5", 0x3ee4_f8b5_88e3_68f1),
            (".5", 0x3fe0_0000_0000_0000),
            ("+2", 0x4000_0000_0000_0000),
            ("-0", 0x8000_0000_0000_0000),
            // glibc's hexadecimal numbers
            ("0x10", 0x4030_0000_0000_0000),
            ("0x1p-3", 0x3fc0_0000_0000_0000),
            ("-0x.8", 0xbfe0_0000_0000_0000),
            // `StringToDoublePosix` reads a final exponent mark, in x87 extended precision
            ("1e", 0x3ff0_0000_0000_0000),
            ("2.5E", 0x4004_0000_0000_0000),
            ("-3e", 0xc008_0000_0000_0000),
            ("12345678901234567890e", 0x43e5_6a95_319d_63e1),
            ("0.0000000000000001e", 0x3c9c_d2b2_97d8_89bc),
            // the infinities, and NaN with glibc's sign and payload
            ("1e400", 0x7ff0_0000_0000_0000),
            ("+Infinity", 0x7ff0_0000_0000_0000),
            ("-nan", 0xfff8_0000_0000_0000),
            ("-nan(1)", 0xfff8_0000_0000_0001),
            ("+nan(5)", 0x7ff8_0000_0000_0005),
            ("-NaN()", 0xfff8_0000_0000_0000),
            ("+nan(x_9)", 0x7ff8_0000_0000_0000),
        ] {
            let read = ncbi_double(value).unwrap_or_else(|e| panic!("{value}: {e}"));
            assert_eq!(read.to_bits(), bits, "{value}");
        }
        for value in [
            "+nan(",
            "+nan(1",
            "+nan(-)",
            "+nan(1)x",
            "+-nan(1)",
            "nan(1)",
            "nan(5)",
            ",5",
            ",",
            "1x",
            "1e+",
            "1ex",
            "--5",
            "2,5",
            "inf",
            "nan",
            " 1",
            "",
            "e5",
            "0x",
            "1.2.3",
            // A character of several bytes across the fourth byte (S11 audit F-1).
            "-xx€",
            "-in\u{fffd}",
            "+\u{fffd}\u{fffd}",
            "-1éé",
        ] {
            assert_eq!(
                ncbi_double(value),
                Err("expected a number".into()),
                "{value:?}"
            );
        }
    }

    // NCBI 2.17.0: `-threshold`'s constraint reads the value with the default flags (`blastp
    // -threshold 1e` and `12345678901234567890e`: "Illegal value, expected >=0", exit 1; `0x10`
    // searches), `-perc_identity`'s with `fDecimalPosixOrLocal` (`1e` and `0x10` search, `0x1p7`
    // is illegal).
    #[test]
    fn each_constraint_reads_the_value_with_its_own_conversion() {
        assert_eq!(blastp_threshold("0x10"), Ok(16.0));
        assert_eq!(blastp_threshold("11"), Ok(11.0));
        assert!(blastp_threshold("1e400").unwrap().is_infinite());
        for value in [
            "1e",
            "2.5E",
            "12345678901234567890e",
            "-nan",
            "+nan(5)",
            "-1",
            "-0x1",
        ] {
            assert_eq!(
                blastp_threshold(value),
                Err("Illegal value, expected greater or equal to 0".into()),
                "{value}"
            );
        }
        assert_eq!(blastp_threshold("1x"), Err("expected a number".into()));
        assert_eq!(blastn_percentage("1e"), Ok(1.0));
        assert_eq!(blastn_percentage("0x10"), Ok(16.0));
        assert_eq!(blastp_query_coverage("99.5E"), Ok(99.5));
        for value in ["0x1p7", "1e400", "-nan", "-0.5e"] {
            assert_eq!(
                blastn_percentage(value),
                Err("expected a percentage in 0..=100".into()),
                "{value}"
            );
        }
    }
}

/// The `-task` of a program: `accepted` names are the choices (help and `describe` list
/// them), a name in `rejected` is a task LOSAT does not implement and fails at the parse
/// with the given message (clap's error for a value its parser refuses, as every other
/// option's), and any other name is clap's invalid-value error.
///
/// NCBI's `-task` is a `CArgAllowStringSet` (`blast_args.cpp:105`); the rejected
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
// ```c++
// set<string>
// CBlastOptionsFactory::GetTasks(ETaskSets choice /* = eAll */)
// {
//     set<string> retval;
//     if (choice == eNuclNucl || choice == eAll) {
//         retval.insert("blastn");
//         retval.insert("blastn-short");
//         retval.insert("megablast");
//         retval.insert("dc-megablast");
//         retval.insert("vecscreen");
//         // -RMH-
//         retval.insert("rmblastn");
//     }
// ...
// ```
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
    /// ```c++
    /// CBlastOptionsHandle*
    /// CBlastOptionsFactory::CreateTask(string task, EAPILocality locality)
    /// {
    ///     CBlastOptionsHandle* retval = NULL;
    ///
    ///     string lc_task(NStr::ToLower(task));
    ///     ThrowIfInvalidTask(lc_task);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:662-673
    /// ```c++
    /// void ThrowIfInvalidTask(const string& task)
    /// {
    ///     set<string> valid_tasks;
    ///     if (valid_tasks.empty()) {
    ///         valid_tasks = CBlastOptionsFactory::GetTasks();
    ///     }
    ///
    ///     if (valid_tasks.find(task) == valid_tasks.end()) {
    ///         ostringstream os;
    ///         os << "'" << task << "' is not a supported task";
    ///         NCBI_THROW(CBlastException, eInvalidArgument, os.str());
    ///     }
    /// }
    /// ```
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

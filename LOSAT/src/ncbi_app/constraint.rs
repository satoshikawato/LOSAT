//! The constraints on the values of NCBI's arguments (`CArgAllow` and its classes): the NCBI C++
//! Toolkit's (`CArgAllow_Strings`, `CArgAllow_Doubles`, `CArgAllow_Integers`, c++/src/corelib/
//! ncbiargs.cpp) and BLAST's (`CArgAllowValuesGreaterThanOrEqual` and the others of
//! c++/include/algo/blast/blastinput/blast_input_aux.hpp, `CArgAllowGeneticCodeInteger` of
//! blast_args.cpp). Each is NCBI's `Verify`, on the text of the value as the command line gives it,
//! and its `GetUsage`, the text of the reader's `Illegal value, expected ...` (`args::process_value`)
//! and of the XML's `<constraint>` (the tables of `args` are made from that XML).

use super::{Span, Text};

/// A constraint on a key's value (`CArgAllow`). The values of a set are rows of the tables of
/// `args` (`STRING_VALUES`, `INTEGER_VALUES`), so that a constraint holds no pointer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Constraint {
    /// `CArgAllow_Strings`: one of the values, compared with case, listed in the order of its
    /// `std::set`.
    Strings(Span),
    /// `CArgAllowValuesGreaterThanOrEqual` with `min`.
    AtLeast(f64),
    /// `CArgAllowValuesLessThanOrEqual` with `max`.
    AtMost(f64),
    /// `CArgAllowValuesBetween` with `min`, `max` and `inclusive`.
    Between { min: f64, max: f64, inclusive: bool },
    /// `CArgAllow_Doubles` with `min` and `max`.
    Doubles { min: f64, max: f64 },
    /// `CArgAllow_Integers` with `min` and `max`.
    Integers { min: i64, max: i64 },
    /// `CArgAllowStringSet`: one of the values, in the order of its `set<string>`.
    StringSet(Span),
    /// `CArgAllowIntegerSet`: one of the values, in the order of its `set<int>`.
    IntegerSet(Span),
    /// `CArgAllowGeneticCodeInteger`: a genetic code of NCBI's tables, or code 32 where LOSAT adds
    /// it (`with_32`, PD-TLOSAN-LOCAL-GENCODE-32: tblastn's `-db_gencode`). The usage text is
    /// NCBI's either way.
    GeneticCode { with_32: bool },
    /// `CArgAllowMaximumFileNameLength` with `max`: a file name shorter than `max` bytes.
    MaxFileNameLength(u32),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:998-1003
/// The genetic codes that `CArgAllowGeneticCodeInteger` allows.
const GENETIC_CODES: [i32; 26] = [
    1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33,
];

/// Whether `code` is a genetic code that `CArgAllowGeneticCodeInteger` allows: the one list of
/// them, for this reader and for BLASTX's check of a code that a library caller set
/// (`blastx::query_setup`).
pub(crate) fn is_genetic_code(code: i32) -> bool {
    GENETIC_CODES.contains(&code)
}

// No NCBI counterpart: PD-TLOSAN-LOCAL-GENCODE-32 (gc.prt ID 32, which LOSAT's TBLASTN translates
// with).
/// The genetic codes of `Constraint::GeneticCode { with_32 }`, in increasing order.
pub(crate) fn genetic_codes(with_32: bool) -> Vec<u32> {
    let codes = GENETIC_CODES.iter().copied().chain(with_32.then_some(32));
    let mut codes: Vec<u32> = codes.filter_map(|code| u32::try_from(code).ok()).collect();
    codes.sort_unstable();
    codes
}

impl Constraint {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1242
    // The reader calls `Verify` inside a `try`: an exception counts as a refused value.
    /// `Verify` of the constraint, on the value's bytes as the command line gives them; an
    /// exception that NCBI's `Verify` throws (a value that its conversion does not read) is a
    /// value it does not allow.
    pub(crate) fn verify(self, bytes: &[u8]) -> bool {
        let blast_double = || crate::ncbi_str::string_to_double(bytes).ok();
        match self {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4605-4609
            // A lookup of the value in the set, with case.
            Constraint::Strings(values) | Constraint::StringSet(values) => {
                strings(values).any(|allowed| allowed.as_bytes() == bytes)
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:111-113
            // BLAST's bounds read the value with `NStr::StringToDouble` (default flags), then
            // compare it; a value it does not read throws.
            Constraint::AtLeast(min) => blast_double().is_some_and(|number| number >= min),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:136-138
            Constraint::AtMost(max) => blast_double().is_some_and(|number| number <= max),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:163-172
            Constraint::Between {
                min,
                max,
                inclusive,
            } => blast_double().is_some_and(|number| {
                if inclusive {
                    number >= min && number <= max
                } else {
                    number > min && number < max
                }
            }),
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
            // Read with `fDecimalPosixOrLocal`, then inside the closed interval.
            Constraint::Doubles { min, max } => {
                crate::ncbi_str::string_to_double_posix_or_local(bytes)
                    .is_ok_and(|number| min <= number && number <= max)
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4686-4695
            // Read with `s_StringToInt8`, then inside the closed interval.
            Constraint::Integers { min, max } => {
                string_to_int8(bytes).is_ok_and(|number| min <= number && number <= max)
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-222, 239
            // Read with `NStr::StringToInt`, then looked up in the set.
            Constraint::IntegerSet(values) => crate::ncbi_str::string_to_int(bytes, 0)
                .is_ok_and(|number| integers(values).contains(&number)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1002-1007
            // Read with `NStr::StringToInt`, then looked up among the genetic codes (and code 32
            // where LOSAT adds it).
            Constraint::GeneticCode { with_32 } => crate::ncbi_str::string_to_int(bytes, 0)
                .is_ok_and(|code| is_genetic_code(code) || (with_32 && code == 32)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:85-88
            // The length of the file name part of the path, below the maximum.
            Constraint::MaxFileNameLength(max) => {
                u64::try_from(file_name_length(bytes)).is_ok_and(|length| length < u64::from(max))
            }
        }
    }

    /// `GetUsage()` of the constraint.
    pub(crate) fn usage(self) -> String {
        match self {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4612-4636
            // Each value in a backquote and a quote, joined by `, ` (the sets are case-sensitive,
            // so no `{case insensitive}`).
            Constraint::Strings(values) => strings(values)
                .map(|value| format!("`{value}'"))
                .collect::<Vec<_>>()
                .join(", "),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:116-118
            Constraint::AtLeast(min) => format!(">={}", double_to_string(min)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:141-143
            Constraint::AtMost(max) => format!("<={}", double_to_string(max)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:175-185
            Constraint::Between {
                min,
                max,
                inclusive,
            } => {
                let (min, max) = (double_to_string(min), double_to_string(max));
                if inclusive {
                    format!("(>={min} and =<{max})")
                } else {
                    format!("(>{min} and <{max})")
                }
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4833-4860
            // One value, one open end, no end, or `min..max`, each number as `DoubleToString`
            // writes it. `kMin_Double` is `DBL_MIN`, the least positive normal number
            // (ncbi_limits.h:143).
            Constraint::Doubles { min, max } => {
                let (low, high) = (min == f64::MIN_POSITIVE, max == f64::MAX);
                match (min == max, low, high) {
                    (true, _, _) => double_to_string(min),
                    (false, true, false) => format!("less or equal to {}", double_to_string(max)),
                    (false, false, true) => {
                        format!("greater or equal to {}", double_to_string(min))
                    }
                    (false, true, true) => String::new(),
                    (false, false, false) => {
                        format!("{}..{}", double_to_string(min), double_to_string(max))
                    }
                }
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4778-4794
            // As for `Doubles`, with the ends of `int` as the open ends; else
            // `CArgAllow_Int8s::GetUsage` of one interval that is not at the ends of `Int8`
            // (ncbiargs.cpp:4698-4726): `min..max`.
            Constraint::Integers { min, max } => {
                let (low, high) = (min == i64::from(i32::MIN), max == i64::from(i32::MAX));
                match (min == max, low, high) {
                    (true, _, _) => min.to_string(),
                    (false, true, false) => format!("less or equal to {max}"),
                    (false, false, true) => format!("greater or equal to {min}"),
                    (false, true, true) => String::new(),
                    (false, false, false) => format!("{min}..{max}"),
                }
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:224-231
            // `Permissible values: ` and each value in quotes, each followed by a space.
            Constraint::StringSet(values) => permissible_values(strings(values)),
            Constraint::IntegerSet(values) => permissible_values(integers(values).iter()),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1010-1012
            Constraint::GeneticCode { .. } => "values between: 1-6, 9-16, 21-31, 33".to_string(),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:91-93
            Constraint::MaxFileNameLength(max) => format!("file name length < {max}"),
        }
    }
}

/// The values of a set of strings (rows of `STRING_VALUES`).
pub(super) fn strings(values: Span) -> impl Iterator<Item = &'static str> {
    values
        .of(&super::STRING_VALUES)
        .iter()
        .map(|&value| Text::get(value))
}

/// The values of a set of integers (rows of `INTEGER_VALUES`).
pub(super) fn integers(values: Span) -> &'static [i32] {
    values.of(&super::INTEGER_VALUES)
}

/// `Permissible values: 'a' 'b' ` of BLAST's set classes.
fn permissible_values<T: std::fmt::Display>(values: impl Iterator<Item = T>) -> String {
    let mut usage = String::from("Permissible values: ");
    for value in values {
        usage += &format!("'{value}' ");
    }
    usage
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2339
/// `NStr::DoubleToString` of `value`: with no precision and the default flags, C's `%g`
/// (`report::sam::push_g`).
fn double_to_string(value: f64) -> String {
    let mut out = Vec::new();
    crate::report::sam::push_g(&mut out, value);
    String::from_utf8(out).expect("the digits of %g")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:115-129
/// `s_StringToInt8` of `value`, NCBI's reading of an integer argument (`CArg_Int8`, `CArg_Integer`)
/// and of `CArgAllow_Integers`: decimal, then, for a value that starts with `0x` or `0X`, base 16,
/// whose exception is then the one thrown (`0x` alone reads 0). The one reading of it, for this
/// reader and for LOSAT Web's ABI v1 (`web_v1_argv::ncbi_int8`).
pub(crate) fn string_to_int8(bytes: &[u8]) -> Result<i64, crate::ncbi_exception::NcbiException> {
    use crate::ncbi_str::{string_to_int8, string_to_int8_base};
    match string_to_int8(bytes, 0) {
        Err(_)
            if bytes
                .get(..2)
                .is_some_and(|head| head.eq_ignore_ascii_case(b"0x")) =>
        {
            string_to_int8_base(bytes, 0, 16)
        }
        converted => converted,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:298-312
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-363
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:465-472
/// The length of the `GetName()` of a `CFile` of `path`: the text after the last separator once the
/// trailing directory separators are removed (a path of one separator, or a disk name, keeps them).
/// On Unix `DIR_SEPARATORS` and `ALL_SEPARATORS` are `/` (ncbifile.cpp:110-115); on Windows they
/// are `/\` and `:/\`, with the disk name `C:`. The one reading of it, for this reader and for
/// LOSAT Web's ABI v1 (`web_v1_argv`, `CArgAllowMaximumFileNameLength` of `-out`).
pub(crate) fn file_name_length(path: &[u8]) -> usize {
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

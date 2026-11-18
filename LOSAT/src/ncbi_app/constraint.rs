//! The constraints on the values of NCBI's arguments (`CArgAllow` and its classes): the NCBI C++
//! Toolkit's (`CArgAllow_Strings`, `CArgAllow_Doubles`, `CArgAllow_Integers`, c++/src/corelib/
//! ncbiargs.cpp) and BLAST's (`CArgAllowValuesGreaterThanOrEqual` and the others of
//! c++/include/algo/blast/blastinput/blast_input_aux.hpp, `CArgAllowGeneticCodeInteger` of
//! blast_args.cpp). Each is NCBI's `Verify`, on the text of the value as the command line gives it,
//! and its `GetUsage`, the text of the reader's `Illegal value, expected ...` (`args::process_value`)
//! and of the XML's `<constraint>` (the tables of `args` are made from that XML).

/// A constraint on a key's value (`CArgAllow`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Constraint {
    /// `CArgAllow_Strings`: one of the values, compared with case, listed in the order of its
    /// `std::set`.
    Strings(&'static [&'static str]),
    /// `CArgAllowValuesGreaterThanOrEqual(min)`.
    AtLeast(f64),
    /// `CArgAllowValuesLessThanOrEqual(max)`.
    AtMost(f64),
    /// `CArgAllowValuesBetween(min, max, inclusive)`.
    Between { min: f64, max: f64, inclusive: bool },
    /// `CArgAllow_Doubles(min, max)`.
    Doubles { min: f64, max: f64 },
    /// `CArgAllow_Integers(min, max)`.
    Integers { min: i64, max: i64 },
    /// `CArgAllowStringSet`: one of the values, in the order of its `set<string>`.
    StringSet(&'static [&'static str]),
    /// `CArgAllowIntegerSet`: one of the values, in the order of its `set<int>`.
    IntegerSet(&'static [i32]),
    /// `CArgAllowGeneticCodeInteger`: a genetic code of NCBI's tables.
    GeneticCode,
    /// `CArgAllowMaximumFileNameLength(max)`: a file name shorter than `max` bytes.
    MaxFileNameLength(u32),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:998-1003
// ```c++
// class CArgAllowGeneticCodeInteger : public CArgAllow
// {
// protected:
//      /// Overloaded method from CArgAllow
//      virtual bool Verify(const string& value) const {
//          static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33};
// ```
/// The genetic codes that `CArgAllowGeneticCodeInteger` allows.
const GENETIC_CODES: [i32; 26] = [
    1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33,
];

impl Constraint {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1242
    // ```c++
    //     if ( m_Constraint ) {
    //         bool err = false;
    //         try {
    //             bool check = m_Constraint->Verify(value);
    //             if (m_NegateConstraint == CArgDescriptions::eConstraintInvert) {
    //                 err = check;
    //             } else {
    //                 err = !check;
    //             }
    //
    //         } catch (...) {
    //             err = true;
    //         }
    // ```
    /// `Verify(value)` of the constraint, on the value's text; an exception that NCBI's `Verify`
    /// throws (a value that its conversion does not read) is a value it does not allow.
    pub(crate) fn verify(self, value: &str) -> bool {
        let bytes = value.as_bytes();
        let blast_double = || crate::ncbi_str::string_to_double(bytes).ok();
        match self {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4605-4609
            // ```c++
            // bool CArgAllow_Strings::Verify(const string& value) const
            // {
            //     TStrings::const_iterator it = m_Strings.find(value);
            //     return it != m_Strings.end();
            // }
            // ```
            Constraint::Strings(values) | Constraint::StringSet(values) => values.contains(&value),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:111-113
            // ```c++
            //     virtual bool Verify(const string& value) const {
            //         return NStr::StringToDouble(value) >= m_MinValue;
            //     }
            // ```
            Constraint::AtLeast(min) => blast_double().is_some_and(|number| number >= min),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:136-138
            // ```c++
            //     virtual bool Verify(const string& value) const {
            //         return NStr::StringToDouble(value) <= m_MaxValue;
            //     }
            // ```
            Constraint::AtMost(max) => blast_double().is_some_and(|number| number <= max),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:163-172
            // ```c++
            //     virtual bool Verify(const string& value) const {
            //         double val = NStr::StringToDouble(value);
            //         bool retval = false;
            //         if ( !m_Inclusive ) {
            //             retval = (val > m_MinValue && val < m_MaxValue);
            //         } else {
            //             retval = (val >= m_MinValue && val <= m_MaxValue);
            //         }
            //         return retval;
            //     }
            // ```
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
            // ```c++
            // bool CArgAllow_Doubles::Verify(const string& value) const
            // {
            //     double val = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
            //     ITERATE( set< TInterval >, pi, m_MinMax) {
            //         if (pi->first <= val && val<= pi->second) {
            //             return true;
            //         }
            //     }
            //     return false;
            // }
            // ```
            Constraint::Doubles { min, max } => {
                crate::ncbi_str::string_to_double_posix_or_local(bytes)
                    .is_ok_and(|number| min <= number && number <= max)
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4686-4695
            // ```c++
            // bool CArgAllow_Int8s::Verify(const string& value) const
            // {
            //     Int8 val = s_StringToInt8(value);
            //     ITERATE( set< TInterval >, pi, m_MinMax) {
            //         if (pi->first <= val && val<= pi->second) {
            //             return true;
            //         }
            //     }
            //     return false;
            // }
            // ```
            Constraint::Integers { min, max } => {
                string_to_int8(value).is_ok_and(|number| min <= number && number <= max)
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-222, 239
            // ```c++
            //     virtual bool Verify(const string& value) const {                        \
            //         DataType value2check = String2DataTypeFn(value);                    \
            //         ITERATE(set<DataType>, itr, m_AllowedValues) {                      \
            //             if (*itr == value2check) {                                      \
            //                 return true;                                                \
            //             }                                                               \
            //         }                                                                   \
            //         return false;                                                       \
            //     }                                                                       \
            // ...
            // DEFINE_CARGALLOW_SET_CLASS(CArgAllowIntegerSet, int, NStr::StringToInt);
            // ```
            Constraint::IntegerSet(values) => crate::ncbi_str::string_to_int(bytes, 0)
                .is_ok_and(|number| values.contains(&number)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1002-1007
            // ```c++
            //      virtual bool Verify(const string& value) const {
            //          static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33};
            //          static const set<int> genetic_codes(gcs, gcs+sizeof(gcs)/sizeof(*gcs));
            //          const int val = NStr::StringToInt(value);
            //          return (genetic_codes.find(val) != genetic_codes.end());
            //      }
            // ```
            Constraint::GeneticCode => crate::ncbi_str::string_to_int(bytes, 0)
                .is_ok_and(|number| GENETIC_CODES.contains(&number)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:85-88
            // ```c++
            //     virtual bool Verify(const string& value) const {
            //         CFile fname(value);
            //         return fname.GetName().size() < m_MaxLength;
            //     }
            // ```
            Constraint::MaxFileNameLength(max) => {
                u64::try_from(file_name_length(bytes)).is_ok_and(|length| length < u64::from(max))
            }
        }
    }

    /// `GetUsage()` of the constraint.
    pub(crate) fn usage(self) -> String {
        match self {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4612-4636
            // ```c++
            //     string str;
            //     TStrings::const_iterator it = m_Strings.begin();
            //     for (;;) {
            //         str += "`";
            //         str += *it;
            //
            //         ++it;
            //         if (it == m_Strings.end()) {
            //             str += "'";
            //             if ( m_Strings.key_comp()("a", "A") ) {
            //                 str += "  {case insensitive}";
            //             }
            //             break;
            //         }
            //         str += "', ";
            //     }
            //     return str;
            // ```
            Constraint::Strings(values) => values
                .iter()
                .map(|value| format!("`{value}'"))
                .collect::<Vec<_>>()
                .join(", "),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:116-118
            // ```c++
            //     virtual string GetUsage(void) const {
            //         return ">=" + NStr::DoubleToString(m_MinValue);
            //     }
            // ```
            Constraint::AtLeast(min) => format!(">={}", double_to_string(min)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:141-143
            // ```c++
            //     virtual string GetUsage(void) const {
            //         return "<=" + NStr::DoubleToString(m_MaxValue);
            //     }
            // ```
            Constraint::AtMost(max) => format!("<={}", double_to_string(max)),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:175-185
            // ```c++
            //     virtual string GetUsage(void) const {
            //         string retval;
            //         if ( !m_Inclusive ) {
            //             retval = "(>" + NStr::DoubleToString(m_MinValue) + " and <"
            //                 + NStr::DoubleToString(m_MaxValue) + ")";
            //         } else {
            //             retval = "(>=" + NStr::DoubleToString(m_MinValue) + " and =<"
            //                 + NStr::DoubleToString(m_MaxValue) + ")";
            //         }
            //         return retval;
            //     }
            // ```
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
            // ```c++
            // string CArgAllow_Doubles::GetUsage(void) const
            // {
            //     if (m_MinMax.size() == 1) {
            //         double x_min = m_MinMax.begin()->first;
            //         double x_max = m_MinMax.begin()->second;
            //         if (x_min == x_max) {
            //             return NStr::DoubleToString(x_min);
            //         } else if (x_min == kMin_Double && x_max != kMax_Double) {
            //             return string("less or equal to ") + NStr::DoubleToString(x_max);
            //         } else if (x_min != kMin_Double && x_max == kMax_Double) {
            //             return string("greater or equal to ") + NStr::DoubleToString(x_min);
            //         } else if (x_min == kMin_Double && x_max == kMax_Double) {
            //             return kEmptyStr;
            //         }
            //     }
            //     string usage;
            //     ITERATE( set< TInterval >, pi, m_MinMax) {
            //         ...
            //             usage += NStr::DoubleToString(pi->first) + ".." + NStr::DoubleToString(pi->second);
            // ```
            // `kMin_Double` is `DBL_MIN`, the least positive normal number (ncbi_limits.h:143).
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
            // ```c++
            // string CArgAllow_Integers::GetUsage(void) const
            // {
            //     if (m_MinMax.size() == 1) {
            //         Int8 x_min = m_MinMax.begin()->first;
            //         Int8 x_max = m_MinMax.begin()->second;
            //         if (x_min == x_max) {
            //             return NStr::Int8ToString(x_min);
            //         } else if (x_min == kMin_Int && x_max != kMax_Int) {
            //             return string("less or equal to ") + NStr::Int8ToString(x_max);
            //         } else if (x_min != kMin_Int && x_max == kMax_Int) {
            //             return string("greater or equal to ") + NStr::Int8ToString(x_min);
            //         } else if (x_min == kMin_Int && x_max == kMax_Int) {
            //             return kEmptyStr;
            //         }
            //     }
            //     return CArgAllow_Int8s::GetUsage();;
            // }
            // ```
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
            // ```c++
            //     virtual string GetUsage(void) const {                                   \
            //         CNcbiOstrstream os;                                                 \
            //         os << "Permissible values: ";                                       \
            //         ITERATE(set<DataType>, itr, m_AllowedValues) {                      \
            //             os << "'" << *itr << "' ";                                      \
            //         }                                                                   \
            //         return CNcbiOstrstreamToString(os);                                 \
            //     }                                                                       \
            // ```
            Constraint::StringSet(values) => permissible_values(values.iter()),
            Constraint::IntegerSet(values) => permissible_values(values.iter()),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1010-1012
            // ```c++
            //      virtual string GetUsage(void) const {
            //          return "values between: 1-6, 9-16, 21-31, 33";
            //      }
            // ```
            Constraint::GeneticCode => "values between: 1-6, 9-16, 21-31, 33".to_string(),
            // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:91-93
            // ```c++
            //     virtual string GetUsage(void) const {
            //         return "file name length < " + NStr::IntToString(m_MaxLength);
            //     }
            // ```
            Constraint::MaxFileNameLength(max) => format!("file name length < {max}"),
        }
    }
}

/// `Permissible values: 'a' 'b' ` of BLAST's set classes.
fn permissible_values<T: std::fmt::Display>(values: impl Iterator<Item = T>) -> String {
    let mut usage = String::from("Permissible values: ");
    for value in values {
        usage += &format!("'{value}' ");
    }
    usage
}

/// `NStr::DoubleToString(value)`: C's `%g` (`report::sam::push_g`).
fn double_to_string(value: f64) -> String {
    let mut out = Vec::new();
    crate::report::sam::push_g(&mut out, value);
    String::from_utf8(out).expect("the digits of %g")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:115-129
// ```c++
// // Allow autodetection among decimal and hex, but NOT octal, in case
// // anyone's been relying on leading zeros being meaningless.
// inline
// Int8 s_StringToInt8(const string& value)
// {
//     try {
//         return NStr::StringToInt8(value);
//     } catch (const CStringException&) {
//         if (NStr::StartsWith(value, "0x", NStr::eNocase)) {
//             return NStr::StringToInt8(value, 0, 16);
//         } else {
//             throw;
//         }
//     }
// }
// ```
/// `s_StringToInt8(value)`, NCBI's reading of an integer argument (`CArg_Int8`, `CArg_Integer`)
/// and of `CArgAllow_Integers`: decimal, then, for a value that starts with `0x` or `0X`, base 16,
/// whose exception is then the one thrown (`0x` alone reads 0).
pub(crate) fn string_to_int8(value: &str) -> Result<i64, crate::ncbi_exception::NcbiException> {
    use crate::ncbi_str::{string_to_int8, string_to_int8_base};
    let bytes = value.as_bytes();
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
// ```c++
// void CDirEntry::Reset(const string& path)
// {
//     m_Path = path;
//     size_t len = path.length();
//     // Root dir
//     if ((len == 1)  &&  IsPathSeparator(path[0])) {
//         return;
//     }
//     // Disk name
// #  if defined(DISK_SEPARATOR)
//     if ( (len == 2 || len == 3) && (path[1] == DISK_SEPARATOR) ) {
//         return;
//     }
// #  endif
//     m_Path = DeleteTrailingPathSeparator(path);
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-363
// ```c++
// void CDirEntry::SplitPath(const string& path, string* dir,
//                           string* base, string* ext)
// {
//     // Get file name
//     size_t pos = path.find_last_of(ALL_SEPARATORS);
//     string filename = (pos == NPOS) ? path : path.substr(pos+1);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:465-472
// ```c++
// string CDirEntry::DeleteTrailingPathSeparator(const string& path)
// {
//     size_t pos = path.find_last_not_of(DIR_SEPARATORS);
//     if (pos + 1 < path.length()) {
//         return path.substr(0, pos + 1);
//     }
//     return path;
// }
// ```
/// The length of `CFile(path).GetName()`: the text after the last separator once the trailing
/// directory separators are removed (a path of one separator, or a disk name, keeps them). On Unix
/// `DIR_SEPARATORS` and `ALL_SEPARATORS` are `/` (ncbifile.cpp:110-115); on Windows they are `/\`
/// and `:/\`, with the disk name `C:`.
fn file_name_length(path: &[u8]) -> usize {
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

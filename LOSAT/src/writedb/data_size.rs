//! makeblastdb's file-local `Uint8ToString_DataSize(v, minprec = 10)` (app/blastdb/makeblastdb.cpp),
//! which writes the "Maximum file size: ..." log line of `-max_file_sz`, and the oracle pairs of that
//! option: the value `NStr::StringToUInt8_DataSize` (`ncbi_str::string_to_uint8_data_size`) gives each
//! input, or its error, and the log line.

/// makeblastdb's `Uint8ToString_DataSize(v, minprec = 10)`.
///
/// The loop in the C++ has no `break`, so it always runs all 8 iterations and `i` ends as `kMods.size()` (8)
/// whatever `v` is; `rv.append(kMods, 8, 1)` then appends nothing (pos == size is allowed), so only "B" is appended
/// and the K/M/G/... letter never appears. And a value below `minprec * 1024` (10240) is divided by 1024
/// in the first iteration, and the quotient (below 10) stays below 10240 and is divided again in each of the
/// remaining iterations, so it ends as 0; values from 10240 up are never divided.
pub fn uint8_to_string_data_size(value: u64) -> String {
    const MINPREC: u32 = 10;
    let mut v = value;
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:343-350
    // ```c++
    //     static string kMods = "KMGTPEZY";
    //
    //     size_t i(0);
    //     for(i = 0; i < kMods.size(); i++) {
    //         if (v < Uint8(minprec)*1024) {
    //             v /= 1024;
    //         }
    //     }
    // ```
    const K_MODS: &str = "KMGTPEZY";

    let mut i: usize = 0;
    while i < K_MODS.len() {
        if v < (MINPREC as u64) * 1024 {
            v /= 1024;
        }
        i += 1;
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:352-359
    // ```c++
    //     string rv = NStr::UInt8ToString(v);
    //
    //     if (i) {
    //         rv.append(kMods, i, 1);
    //         rv.append("B");
    //     }
    //
    //     return rv;
    // ```
    let mut rv = v.to_string();

    if i != 0 {
        // rv.append(kMods, i, 1): pos == size appends nothing (pos > size would throw, but i <= size)
        rv.push_str(K_MODS.get(i..i + 1).unwrap_or(""));
        rv.push('B');
    }

    rv
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncbi_str::string_to_uint8_data_size;

    // Oracle pairs: each input was given to makeblastdb 2.17.0 as `-max_file_sz <input>`
    // (`-blastdb_version 4 -dbtype prot`), and the expected strings are what it printed.

    /// Inputs whose log line is "Maximum file size: <expected>".
    const ORACLE_LOG: &[(&str, &str)] = &[
        ("0", "0B"),
        ("1", "0B"),
        ("5", "0B"),
        ("7", "0B"),
        ("10", "0B"),
        ("9", "0B"),
        ("10239", "0B"),
        ("10240", "10240B"),
        ("10241", "10241B"),
        ("10000", "0B"),
        ("10239B", "0B"),
        ("10240B", "10240B"),
        ("1023", "0B"),
        ("1024", "0B"),
        ("2048", "0B"),
        ("4096", "0B"),
        ("12345", "12345B"),
        ("300B", "0B"),
        ("300", "0B"),
        ("1KB", "0B"),
        ("1000", "0B"),
        ("1K", "0B"),
        ("1 KB", "0B"),
        ("1kb", "0B"),
        ("1KiB", "0B"),
        ("1.5KB", "0B"),
        ("1MiB", "1048576B"),
        ("3GB", "3000000000B"),
        ("4GB", "4000000000B"),
        ("4294967295", "4294967295B"),
        ("+5", "0B"),
        ("5.0", "0B"),
        ("5.5", "0B"),
        ("0.5", "0B"),
        ("0.4", "0B"),
        ("1.4999", "0B"),
        ("2.5", "0B"),
        ("10.5", "0B"),
        ("1.5", "0B"),
        ("1.50", "0B"),
        ("1.4", "0B"),
        ("00010", "0B"),
        ("0000000000000000000000000000000000012", "0B"),
        ("1B", "0B"),
        ("1b", "0B"),
        ("1 B", "0B"),
        ("1 b", "0B"),
        ("9999K", "9999000B"),
        ("10M", "10000000B"),
        ("10m", "10000000B"),
        ("10Mb", "10000000B"),
        ("10MB", "10000000B"),
        ("10MiB", "10485760B"),
        ("10mib", "10485760B"),
        ("10mIb", "10485760B"),
        ("10MIB", "10485760B"),
        ("4096M", "4096000000B"),
        ("4095M", "4095000000B"),
        ("4095MiB", "4293918720B"),
        ("4194303K", "4194303000B"),
        ("4194304K", "4194304000B"),
        ("4194303KiB", "4294966272B"),
        ("0T", "0B"),
        ("0TiB", "0B"),
        ("0EiB", "0B"),
        ("0.5KiB", "0B"),
        ("0.5K", "0B"),
        ("0.0004K", "0B"),
        ("0.0005K", "0B"),
        ("0.0004999K", "0B"),
        ("0.0005KiB", "0B"),
        ("0.0004999KiB", "0B"),
        ("1.0005K", "0B"),
        ("1.2345K", "0B"),
        ("1.2345678K", "0B"),
        ("1.2345678KiB", "0B"),
        ("0.9999M", "999900B"),
        ("0.0009999M", "0B"),
        ("1.5MiB", "1572864B"),
        ("1.5M", "1500000B"),
        ("2.5GiB", "2684354560B"),
        ("2.5G", "2500000000B"),
        ("2.5Gb", "2500000000B"),
        ("1 KiB", "0B"),
        ("10239.5", "10240B"),
        ("10239.4", "0B"),
        ("10240.4", "10240B"),
        ("9.9999K", "0B"),
        ("10.24K", "10240B"),
        ("10.2399K", "10240B"),
        ("10.2394K", "0B"),
        ("10KiB", "10240B"),
        ("10.0KiB", "10240B"),
        ("9.99KiB", "0B"),
        ("9.999KiB", "0B"),
        ("9.9999KiB", "10240B"),
        ("10.2400KiB", "10486B"),
        ("9.7656KiB", "0B"),
        ("9.7657KiB", "0B"),
        ("10239.4999", "0B"),
        ("0.00001KiB", "0B"),
        ("1k", "0B"),
        ("1m", "1000000B"),
        ("1g", "1000000000B"),
        ("1Kib", "0B"),
        ("1kIB", "0B"),
        ("2147483648", "2147483648B"),
        ("2147483648B", "2147483648B"),
        ("4294967295B", "4294967295B"),
        ("4294967295.4", "4294967295B"),
        ("4294967294.5", "4294967295B"),
        ("4294967295 B", "4294967295B"),
    ];

    /// Inputs that parse to a value of 4 GiB or more: makeblastdb 2.17.0 stops with
    /// "BLAST options error: max_file_sz must be < 4 GiB" before the log line, so the parse did not fail.
    const ORACLE_AT_LEAST_4_GIB: &[&str] = &[
        "4294967296",
        "4294967297",
        "5GiB",
        "4GiB",
        "4096MiB",
        "4194304KiB",
        "4194304 KiB",
        "1T",
        "1Tb",
        "1TB",
        "1TiB",
        "1P",
        "1PB",
        "1PiB",
        "1E",
        "1EB",
        "1EiB",
        "1e",
        "18446744073709551615",
        "18446744073709551614",
        "18446744073709551615.4",
        "18446744073709551615.0",
        "1844674407370955161",
        "1844674407370955162",
        "18446744073709551K",
        "18446744073709551.615K",
        "18446744073709551.615KB",
        "15EiB",
        "16E",
        "18E",
        "00000000000000000000018446744073709551615",
        "1t",
        "1p",
        "4294967295.5",
    ];

    /// Inputs for which makeblastdb 2.17.0 ended with the `CStringException` of `StringToUInt8_DataSize`:
    /// (input, text after "StringToUInt8_DataSize() - " and before " (m_Pos = ", the m_Pos, the line of
    /// ncbistr.cpp that its report names: probe of 2026-11-08, the official Linux x64 release).
    const ORACLE_ERRORS: &[(&str, &str, usize, u32)] = &[
        ("-5", "Cannot convert string '-5' to Uint8", 0, 1578),
        ("abc", "Cannot convert string 'abc' to Uint8", 0, 1622),
        ("", "Cannot convert string '' to Uint8", 0, 1622),
        ("+", "Cannot convert string '+' to Uint8", 1, 1622),
        ("-", "Cannot convert string '-' to Uint8", 0, 1578),
        ("+-5", "Cannot convert string '+-5' to Uint8", 1, 1622),
        ("++5", "Cannot convert string '++5' to Uint8", 1, 1622),
        (".5", "Cannot convert string '.5' to Uint8", 0, 1622),
        ("5.", "Cannot convert string '5.' to Uint8", 1, 1674),
        (" 1", "Cannot convert string ' 1' to Uint8", 0, 1622),
        (" 1KB", "Cannot convert string ' 1KB' to Uint8", 0, 1622),
        ("1 ", "Cannot convert string '1 ' to Uint8", 1, 1674),
        ("1KB ", "Cannot convert string '1KB ' to Uint8", 3, 1674),
        ("1  KB", "Cannot convert string '1  KB' to Uint8", 1, 1674),
        ("1 K B", "Cannot convert string '1 K B' to Uint8", 3, 1674),
        ("1K B", "Cannot convert string '1K B' to Uint8", 2, 1674),
        ("1KBB", "Cannot convert string '1KBB' to Uint8", 3, 1674),
        ("1KiBB", "Cannot convert string '1KiBB' to Uint8", 4, 1674),
        ("1KiB ", "Cannot convert string '1KiB ' to Uint8", 4, 1674),
        ("1iB", "Cannot convert string '1iB' to Uint8", 1, 1674),
        ("1B ", "Cannot convert string '1B ' to Uint8", 2, 1674),
        ("1BK", "Cannot convert string '1BK' to Uint8", 2, 1674),
        ("1 BK", "Cannot convert string '1 BK' to Uint8", 3, 1674),
        ("5Gi", "Cannot convert string '5Gi' to Uint8", 2, 1674),
        ("1Z", "Cannot convert string '1Z' to Uint8", 1, 1674),
        ("1Y", "Cannot convert string '1Y' to Uint8", 1, 1674),
        (
            "18446744073709551616",
            "Cannot convert string '18446744073709551616' to Uint8",
            19,
            1750,
        ),
        (
            "18446744073709551615.5",
            "Cannot convert string '18446744073709551615.5' to Uint8",
            20,
            1757,
        ),
        (
            "18446744073709551616.5",
            "Cannot convert string '18446744073709551616.5' to Uint8",
            19,
            1750,
        ),
        (
            "99999999999999999999",
            "Cannot convert string '99999999999999999999' to Uint8",
            19,
            1750,
        ),
        (
            "184467440737095516150",
            "Cannot convert string '184467440737095516150' to Uint8",
            20,
            1750,
        ),
        (
            "18446744073709552K",
            "Cannot convert string '18446744073709552K' to Uint8",
            19,
            1750,
        ),
        (
            "18446744073709551.616K",
            "Cannot convert string '18446744073709551.616K' to Uint8",
            19,
            1750,
        ),
        (
            "18446744073709551.6155K",
            "Cannot convert string '18446744073709551.6155K' to Uint8",
            20,
            1757,
        ),
        (
            "18446744073709551615KiB",
            "Cannot convert string '18446744073709551615KiB' to Uint8",
            19,
            1750,
        ),
        (
            "17179869184GiB",
            "Cannot convert string '17179869184GiB' to Uint8",
            19,
            1750,
        ),
        ("16EiB", "Cannot convert string '16EiB' to Uint8", 19, 1750),
        (
            "18446744073709551616E",
            "Cannot convert string '18446744073709551616E' to Uint8",
            19,
            1750,
        ),
        ("1.2.3", "Cannot convert string '1.2.3' to Uint8", 3, 1674),
        ("1..2", "Cannot convert string '1..2' to Uint8", 1, 1674),
        ("1.K", "Cannot convert string '1.K' to Uint8", 1, 1674),
        ("1.KiB", "Cannot convert string '1.KiB' to Uint8", 1, 1674),
        ("1. K", "Cannot convert string '1. K' to Uint8", 1, 1674),
        ("1 .5", "Cannot convert string '1 .5' to Uint8", 1, 1674),
        (".K", "Cannot convert string '.K' to Uint8", 0, 1622),
        (".", "Cannot convert string '.' to Uint8", 0, 1622),
        ("0x10", "Cannot convert string '0x10' to Uint8", 1, 1674),
        ("1e3", "Cannot convert string '1e3' to Uint8", 2, 1674),
        ("1,000", "Cannot convert string '1,000' to Uint8", 1, 1674),
        ("1_000", "Cannot convert string '1_000' to Uint8", 1, 1674),
        ("1KiBs", "Cannot convert string '1KiBs' to Uint8", 4, 1674),
        (
            "1\u{9}KB",
            "Cannot convert string '1\\tKB' to Uint8",
            1,
            1674,
        ),
        ("\u{9}1", "Cannot convert string '\\t1' to Uint8", 0, 1622),
        ("1\u{a}", "Cannot convert string '1\\n' to Uint8", 1, 1674),
        ("1\u{a}2", "Cannot convert string '1\\n2' to Uint8", 1, 1674),
        ("1?", "Cannot convert string '1?' to Uint8", 1, 1674),
        ("1??", "Cannot convert string '1\\?\\?' to Uint8", 1, 1674),
        (
            "1???",
            "Cannot convert string '1\\?\\?\\?' to Uint8",
            1,
            1674,
        ),
        ("?1", "Cannot convert string '?1' to Uint8", 0, 1622),
        ("??", "Cannot convert string '\\?\\?' to Uint8", 0, 1622),
        ("a&b", "Cannot convert string 'a&b' to Uint8", 0, 1622),
        ("1&", "Cannot convert string '1&' to Uint8", 1, 1674),
        (
            "\u{e9}1",
            "Cannot convert string '\\303\\2511' to Uint8",
            0,
            1622,
        ),
        (
            "1\u{e9}",
            "Cannot convert string '1\\303\\251' to Uint8",
            1,
            1674,
        ),
        ("1\"2", "Cannot convert string '1\\\"2' to Uint8", 1, 1674),
        ("1'2", "Cannot convert string '1\\'2' to Uint8", 1, 1674),
        ("1\\2", "Cannot convert string '1\\\\2' to Uint8", 1, 1674),
        ("1\u{7}", "Cannot convert string '1\\a' to Uint8", 1, 1674),
        ("1\u{1}", "Cannot convert string '1\\1' to Uint8", 1, 1674),
        (
            "1\u{7f}",
            "Cannot convert string '1\\177' to Uint8",
            1,
            1674,
        ),
        (
            "1\u{1}2",
            "Cannot convert string '1\\0012' to Uint8",
            1,
            1674,
        ),
        ("1KiB1", "Cannot convert string '1KiB1' to Uint8", 4, 1674),
        ("1 KiB ", "Cannot convert string '1 KiB ' to Uint8", 5, 1674),
        (
            "4294967295 ",
            "Cannot convert string '4294967295 ' to Uint8",
            10,
            1674,
        ),
    ];

    fn log_line(input: &str) -> String {
        uint8_to_string_data_size(string_to_uint8_data_size(input).unwrap())
    }

    #[test]
    fn oracle_log_lines() {
        for (input, want) in ORACLE_LOG {
            assert_eq!(&log_line(input), want, "input {:?}", input);
        }
    }

    #[test]
    fn oracle_at_least_4_gib() {
        for input in ORACLE_AT_LEAST_4_GIB {
            let v = string_to_uint8_data_size(input).unwrap();
            assert!(v >= 1 << 32, "input {:?} gave {}", input, v);
            assert_eq!(uint8_to_string_data_size(v), format!("{}B", v));
        }
    }

    #[test]
    fn oracle_errors() {
        for (input, message, pos, line) in ORACLE_ERRORS {
            let e = string_to_uint8_data_size(input).unwrap_err();
            assert_eq!(&e.message, message, "input {:?}", input);
            assert_eq!(e.pos, *pos, "input {:?}", input);
            assert_eq!(e.site.line, *line, "input {:?}", input);
        }
    }

    #[test]
    fn oracle_pair_count() {
        assert_eq!(
            ORACLE_LOG.len() + ORACLE_AT_LEAST_4_GIB.len() + ORACLE_ERRORS.len(),
            220
        );
    }

    #[test]
    fn data_size_formatting_thresholds() {
        // Rule 1: the letter index is always kMods.size(), so only "B" is appended.
        // Rule 2: values below 10 * 1024 are divided down to 0, larger ones are kept.
        assert_eq!(uint8_to_string_data_size(0), "0B");
        assert_eq!(uint8_to_string_data_size(1), "0B");
        assert_eq!(uint8_to_string_data_size(1023), "0B");
        assert_eq!(uint8_to_string_data_size(1024), "0B");
        assert_eq!(uint8_to_string_data_size(10239), "0B");
        assert_eq!(uint8_to_string_data_size(10240), "10240B");
        assert_eq!(uint8_to_string_data_size(10241), "10241B");
        assert_eq!(uint8_to_string_data_size(1 << 20), "1048576B");
        assert_eq!(uint8_to_string_data_size(4294967295), "4294967295B");
        assert_eq!(uint8_to_string_data_size(u64::MAX), "18446744073709551615B");
    }
}

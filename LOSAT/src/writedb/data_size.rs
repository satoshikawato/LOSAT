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
    const K_MODS: &str = "KMGTPEZY";

    let mut i: usize = 0;
    while i < K_MODS.len() {
        if v < (MINPREC as u64) * 1024 {
            v /= 1024;
        }
        i += 1;
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:352-359
    let mut rv = v.to_string();

    if i != 0 {
        // rv.append(kMods, i, 1): pos == size appends nothing (pos > size would throw, but i <= size)
        rv.push_str(K_MODS.get(i..i + 1).unwrap_or(""));
        rv.push('B');
    }

    rv
}

#[cfg(test)]
#[path = "data_size_tests.rs"]
mod tests;

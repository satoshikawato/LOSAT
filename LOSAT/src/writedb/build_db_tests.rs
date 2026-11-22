use super::*;

#[cfg(not(windows))]
#[test]
fn absolute_paths_are_normalized_as_ncbi_does() {
    assert_eq!(normalize_path("/a//b/./c/../d/"), "/a/b/d");
    assert_eq!(normalize_path("/a/b/.."), "/a");
    assert_eq!(normalize_path("/../a"), "/a");
    assert_eq!(normalize_path("/a/./"), "/a");
    assert_eq!(normalize_path("/"), "/");
    assert_eq!(normalize_path("/a/.."), "/");
}

// `CDirEntry::NormalizePath` of relative paths (ncbifile.cpp:820-971), the
// name that `CreateAbsolutePath` gives without a current directory.
#[cfg(not(windows))]
#[test]
fn relative_paths_are_normalized_as_ncbi_does() {
    assert_eq!(normalize_path(""), "");
    assert_eq!(normalize_path("n"), "n");
    assert_eq!(normalize_path("./n"), "n");
    assert_eq!(normalize_path("a/./b//c/"), "a/b/c");
    assert_eq!(normalize_path("a/../b"), "b");
    assert_eq!(normalize_path("a/.."), ".");
    assert_eq!(normalize_path("../a/../../b"), "../../b");
    assert_eq!(normalize_path("."), ".");
    assert_eq!(normalize_path(".//a"), ".//a");
}

// `CDirEntry::CreateAbsolutePath` relative to the current directory: the
// rule of every target but Windows, the WASI commands' included, whose
// current directory is the one the host gives in `PWD` (`/` without it).
// makeblastdb 2.17.0 printed "New DB name:   <cwd>/t" for `-out $'t\t'`
// (and wrote `t.nin`), and "<parent>/sp/m" for `-out m` in a directory
// `sp `.
#[cfg(not(windows))]
#[test]
fn a_relative_name_is_joined_to_the_current_directory_as_ncbi_does() {
    let cwd = |dir: &'static str| Some(Path::new(dir));
    assert_eq!(absolute_path_in(cwd("/w"), "n"), "/w/n");
    assert_eq!(absolute_path_in(cwd("/w/"), "n"), "/w/n");
    assert_eq!(absolute_path_in(cwd("/"), "n"), "/n");
    assert_eq!(absolute_path_in(cwd("/w"), "./sub/../n"), "/w/n");
    assert_eq!(absolute_path_in(cwd("/w"), "../n"), "/n");
    assert_eq!(absolute_path_in(cwd("/w"), "sub/"), "/w/sub");
    assert_eq!(absolute_path_in(cwd("/w"), "/abs//n/"), "/abs/n");
    // `ConcatPath` cuts the spaces at the ends of both parts, also those
    // of the current directory.
    assert_eq!(absolute_path_in(cwd("/w"), "t\t"), "/w/t");
    assert_eq!(absolute_path_in(cwd("/w"), "\t/x"), "/w/x");
    assert_eq!(absolute_path_in(cwd("/w/sp "), "m"), "/w/sp/m");
    // No current directory: the name, normalized.
    assert_eq!(absolute_path_in(None, "n"), "n");
    assert_eq!(absolute_path_in(None, "./a/../n\t"), "n");
    assert_eq!(absolute_path_in(None, "\t/x"), "/x");
    // The process's own directory, as `create_absolute_path` reads it.
    let here = std::env::current_dir().unwrap();
    assert_eq!(
        create_absolute_path("n"),
        format!("{}/n", here.to_string_lossy().trim_end_matches('/'))
    );
}

// Values makeblastdb 2.17.0 printed in "added N sequences in T seconds."
#[test]
fn seconds_print_as_six_significant_digits() {
    assert_eq!(format_seconds(0.000254154), "0.000254154");
    assert_eq!(format_seconds(0.00107908), "0.00107908");
    assert_eq!(format_seconds(0.0001719), "0.0001719");
    assert_eq!(format_seconds(0.000169992), "0.000169992");
    assert_eq!(format_seconds(5.5), "5.5");
    assert_eq!(format_seconds(12.3456789), "12.3457");
    assert_eq!(format_seconds(0.00001234567), "1.23457e-05");
    assert_eq!(format_seconds(1234567.0), "1.23457e+06");
    assert_eq!(format_seconds(100.0), "100");
    assert_eq!(format_seconds(0.0), "0");
}

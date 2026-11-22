use super::*;

// ncbifile.cpp:820-971 (cited at `ncbi_environment::normalize_path`).
#[test]
fn normalize_path_as_ncbi() {
    for (path, normalized) in [
        ("", ""),
        ("/", "/"),
        ("///", "/"),
        ("losat", "losat"),
        ("./losat", "losat"),
        ("bin/./losat", "bin/losat"),
        ("bin//losat", "bin/losat"),
        ("bin/../losat", "losat"),
        ("../bin/losat", "../bin/losat"),
        ("/../bin/losat", "/bin/losat"),
        ("bin/losat/", "bin/losat"),
        (".", "."),
        ("./", "."),
    ] {
        assert_eq!(
            normalize_path(path.as_bytes()),
            normalized.as_bytes(),
            "{path}"
        );
    }
}

// rpsbproc_main.cpp:1809,1836 (cited at `print_dom_line`, `print_cluster_line`): `m_iFrom` plus
// `COORDSBASE` for a hit that ends at 2^31 - 1 overflows `int` (both 0.5.1 builds print
// -2147483648: PS16c review oracle `plusmax`); its value, 2^31, is what NCBI prints for the hit
// moved 3 bases down plus 3 (PS16-Q1 oracle `ntplusmax` and its nearby input `ntplusmax3`;
// PD-LOSAT-NCBI-DEFECTS exception 12). The XML input's value one less than `INT_MIN` prints
// `INT_MIN`.
#[test]
fn printed_positions_past_int() {
    assert_eq!(printed_position(i32::MAX), 2_147_483_648);
    assert_eq!(printed_position(i32::MAX - 3), 2_147_483_645);
    assert_eq!(printed_position(i64::from(i32::MIN) - 1), -2_147_483_648);
    assert_eq!(printed_position(0), 1);
}

// Streaming a `double` of -0.0 to an `ostream` is glibc's `%g`, "-0".
#[test]
fn format_g_keeps_the_sign_of_zero() {
    assert_eq!(format_g(-0.0), "-0");
    assert_eq!(format_g(0.0), "0");
    assert_eq!(format_g(3.06789e-55), "3.06789e-55");
    assert_eq!(format_g(179.393), "179.393");
    assert_eq!(OVERLAPLEADING, 1_410_065_408);
}

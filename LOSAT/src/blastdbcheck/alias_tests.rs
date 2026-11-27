use super::*;

/// The lines of `bytes` as the test reads them (`alias_file_test`'s loop).
fn read_lines(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut reader = LineReader::new(FastaStream::from_bytes(bytes));
    let mut lines = Vec::new();
    loop {
        lines.push(reader.next_line().to_vec());
        if reader.at_eof() {
            return lines;
        }
    }
}

/// NCBI's `CStreamLineReader` learns the end of line from the first line; the other byte inside a
/// later line cuts it, and the rest after that byte runs into the next line (review F1 of row
/// RG-D2b, rows mix1 and mix2: `unknown keyword encountered: NSEQMAXLEN`).
#[test]
fn the_rest_after_the_other_end_of_line_runs_into_the_next_line() {
    let got = read_lines(b"DBLIST p5\nTITLE a\rNSEQ\nMAXLEN 5\n");
    assert!(got.contains(&b"NSEQMAXLEN 5".to_vec()), "{got:?}");
    let got = read_lines(b"A 1\r\nB 2\r\n");
    let lines: Vec<&[u8]> = got
        .iter()
        .map(Vec::as_slice)
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(lines, [&b"A 1"[..], b"B 2"]);
}

/// Approved exception (`Q-D2b-OV=A`): `a_ov1` (`FIRST_OID 65535`, `LAST_OID 65535`) and `a_ov2`
/// (`OIDLIST`, `MAXOID 65535`, `LENGTH 4294967295`, `NSEQ 65535`) compute exact products in `i128`
/// without false errors, while `a_ov3` (`FIRST_OID 2147483647`, `LAST_OID 0`, `NSEQ 2147483647`)
/// still reports the genuine NSEQ mismatch.
#[test]
fn widened_alias_checks_avoid_ncbi_overflow_false_errors() {
    use crate::ncbi_app::files::standard_output_file;
    use crate::report::run_files::{OpenFiles, RunFiles};

    let run_check = |values: Values| -> (u32, String) {
        let mut buf = Vec::new();
        let faults = {
            let files = OpenFiles::new(&RunFiles::Disk);
            let mut log = Log::new(standard_output_file(&mut buf, &files), DETAILS);
            check_values(&mut log, b"  p: ", &values).ok().unwrap()
        };
        (faults, String::from_utf8(buf).unwrap())
    };
    // Row a_ov1: (65535 + 1) * (65535 + 1) != 0 in i128, nseq == -1 -> 0 faults.
    let (faults, out) = run_check(Values {
        length: -1,
        oidlist: -1,
        gilist: -1,
        nseq: -1,
        first_oid: 65535,
        last_oid: 65535,
        maxoid: -1,
    });
    assert_eq!(faults, 0, "{out}");
    assert!(out.is_empty(), "{out}");

    // Row a_ov2: (65535 + 1) * (4294967295 + 1) * (65535 + 1) == 2^64 != 0 in i128 -> 0 faults.
    let (faults, out) = run_check(Values {
        length: 4_294_967_295,
        oidlist: 0,
        gilist: -1,
        nseq: 65535,
        first_oid: -1,
        last_oid: -1,
        maxoid: 65535,
    });
    assert_eq!(faults, 0, "{out}");
    assert!(out.is_empty(), "{out}");

    // Row a_ov3: genuinely inconsistent with NSEQ -> 1 fault.
    let (faults, out) = run_check(Values {
        length: -1,
        oidlist: -1,
        gilist: -1,
        nseq: i32::MAX,
        first_oid: i32::MAX,
        last_oid: 0,
        maxoid: -1,
    });
    assert_eq!(faults, 1, "{out}");
    assert!(
        out.contains("(FIRST_OID, LAST_OID) is not consistent with NSEQ"),
        "{out}"
    );
}

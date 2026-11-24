//! The rows of `tests/fixtures/blastdbcheck` whose output depends on the time-seeded random sample
//! (`seed` not `-` in `cases.tsv`): run in the process with NCBI's seed (`run_seeded`), their lines
//! of the sample test against NCBI's. The other lines of these rows depend on the run's directory
//! and environment, which the fixture test (`tests/blastdbcheck_fixtures.rs`) controls for every
//! other row.

use std::ffi::OsString;
use std::path::Path;

/// The directory that stands for `<RUN>`: its `in` is the fixtures' base databases.
fn run_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdb_aliastool")
}

/// The lines of the sample test (`/ Sample: `) of a report.
fn sample_lines(report: &[u8]) -> Vec<Vec<u8>> {
    report
        .split(|&byte| byte == b'\n')
        .filter(|line| line.windows(11).any(|w| w == b" / Sample: "))
        .map(<[u8]>::to_vec)
        .collect()
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:609-656
#[test]
fn the_seeded_rows_draw_ncbis_sample() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdbcheck");
    let table = std::fs::read_to_string(fixtures.join("cases.tsv")).unwrap();
    let run = run_dir();
    let run_text = run.to_str().unwrap();
    let mut checked = 0;
    for row in table.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = row.split('\t').collect();
        let (id, seed, args) = (fields[0], fields[3], fields[4]);
        if seed == "-" {
            continue;
        }
        let words: Vec<OsString> = args
            .split_whitespace()
            .map(|word| OsString::from(word.replace("../../in", &format!("{run_text}/in"))))
            .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        super::run_seeded(&words, seed.parse().unwrap(), &mut out, &mut err);
        let got: Vec<Vec<u8>> = sample_lines(&out)
            .into_iter()
            .map(|line| {
                String::from_utf8(line)
                    .unwrap()
                    .replace(run_text, "<RUN>")
                    .into_bytes()
            })
            .collect();
        let want = sample_lines(&std::fs::read(fixtures.join(format!("ncbi/{id}.out"))).unwrap());
        assert!(!want.is_empty(), "{id}: no sample lines in NCBI's report");
        assert_eq!(
            got,
            want,
            "{id}: LOSAT's report {}\nerrors {}",
            String::from_utf8_lossy(&out),
            String::from_utf8_lossy(&err)
        );
        checked += 1;
    }
    assert_eq!(checked, 4, "the seeded rows");
}

/// A `-db` or `-dir` that is not UTF-8 is LOSAT's refusal, which names the key, and nothing is
/// written to the standard output.
#[cfg(unix)]
#[test]
fn a_name_that_is_not_utf8_is_refused() {
    use std::os::unix::ffi::OsStringExt;
    for key in ["-db", "-dir"] {
        let words = vec![OsString::from(key), OsString::from_vec(b"d\xe9".to_vec())];
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = super::run_with(&words, &mut std::io::empty(), &mut out, &mut err);
        let err = String::from_utf8_lossy(&err);
        assert_eq!(status, 1, "{key}: {err}");
        assert!(
            err.contains(&format!("the value of {key} is not UTF-8")),
            "{key}: {err}"
        );
        assert!(out.is_empty(), "{key}");
    }
}

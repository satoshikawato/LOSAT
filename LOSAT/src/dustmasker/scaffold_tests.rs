//! TEMPORARY: runs the `-infmt blastdb` fixture rows through `work::run` until the command line
//! reader is merged.
use super::work::{report, run, Failure, Options, Report, SequenceSource};

fn fixture_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/dustmasker")
}

#[test]
fn database_rows_match_ncbi() {
    let dir = fixture_dir();
    let cases = std::fs::read_to_string(dir.join("cases.tsv")).unwrap();
    let mut checked = 0;
    let mut failed = Vec::new();
    for line in cases.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (name, args) = line.split_once('\t').unwrap();
        if !args.contains("-infmt blastdb")
            || args.contains('<')
            || name == "db_missing"
            || name == "db_prot"
        {
            continue;
        }
        let words: Vec<&str> = args.split_whitespace().collect();
        let value = |key: &str| words.iter().position(|w| *w == key).map(|i| words[i + 1]);
        let input = dir
            .join(value("-in").unwrap())
            .to_string_lossy()
            .into_owned();
        let options = Options {
            input: &input,
            infmt: value("-infmt").unwrap_or("fasta"),
            outfmt: value("-outfmt").unwrap_or("interval"),
            window: value("-window").map_or(64, |v| v.parse().unwrap()),
            level: value("-level").map_or(20, |v| v.parse().unwrap()),
            linker: value("-linker").map_or(1, |v| v.parse().unwrap()),
            parse_seqids: words.contains(&"-parse_seqids"),
            hard_masking: words.contains(&"-hard_masking"),
            blastdb_registry: None,
        };
        let mut out = Vec::new();
        let mut no_fasta = || -> Result<Box<dyn SequenceSource>, Failure> { unreachable!() };
        let (err, rc) = match run(&options, &mut no_fasta, &mut out) {
            Ok(()) => (String::new(), 0),
            Err(failure) => match report(failure) {
                Report::Posts { lines, status } => (lines, status),
                Report::Refused(reason) => (format!("REFUSED {reason}\n"), -1),
            },
        };
        let expected_out = std::fs::read(dir.join(format!("ncbi/{name}.out"))).unwrap();
        let expected_err = std::fs::read_to_string(dir.join(format!("ncbi/{name}.err"))).unwrap();
        let expected_rc: i32 = std::fs::read_to_string(dir.join(format!("ncbi/{name}.rc")))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        checked += 1;
        if out != expected_out || err != expected_err || rc != expected_rc {
            failed.push(format!(
                "{name}: out {} err {} rc {rc}/{expected_rc}\n  got err {err:?}\n  got out {:?}",
                out == expected_out,
                err == expected_err,
                String::from_utf8_lossy(&out[..out.len().min(300)])
            ));
        }
    }
    assert!(
        failed.is_empty(),
        "{} of {checked} differ:\n{}",
        failed.len(),
        failed.join("\n")
    );
}

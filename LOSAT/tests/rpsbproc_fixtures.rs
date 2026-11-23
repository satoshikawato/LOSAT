//! `losat rpsbproc` against NCBI's rpsbproc 0.5.1 (the binary of the 2025-05-13 package; bioconda's
//! 0.5.1 on the BLAST+ 2.17.0 toolkit gave the same bytes): every case of
//! `tests/fixtures/rpsbproc/CASES.tsv` (README.md there: the archives are NCBI rpsblast and
//! rpstblastn 2.17.0 `-outfmt 11` and `-outfmt 5` of E. coli K-12 MG1655 proteins and genome pieces
//! against CDD v3.22, the data files are CDD v3.22 cut to the rows the archives reach). A `compare` case gives NCBI's standard output, standard
//! error and exit status (and the `-o` file `out.txt` and the `-logfile` file `log.txt`); a
//! `refuse` case is refused by LOSAT (exit 1, the
//! last line of standard error `Error: ... is not supported by LOSAT's rpsbproc...`, standard
//! output cut at a line end) with the bytes of `refused/` (LOSAT's own; README.md); a `widened`
//! case gives the bytes of `widened/`, LOSAT's own where NCBI's `int` coordinate arithmetic
//! overflows (PD-LOSAT-NCBI-DEFECTS exception 12; each beside a `compare` case of its nearby
//! input, README.md); a `tied` case gives the bytes of `tied/`, LOSAT's own order of hits that
//! NCBI's `std::sort` leaves in an order of its own (PD-LOSAT-NCBI-DEFECTS 1.67), which differ
//! from NCBI's frozen bytes only as `tied/CLASSES.tsv` counts. NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rpsbproc")
}

/// A directory of its own for one test (the tests of a binary run in parallel) holding a copy
/// of the fixture's `arch/` and `data/` and of the files of `root/` (a copy, not a link:
/// Windows; a gzip-compressed archive `arch/<name>.asn.gz` or `arch/<name>.xml.gz` unpacked to
/// `arch/<name>.asn` or `arch/<name>.xml`), removed at the end.
struct RunDir(PathBuf);

impl RunDir {
    fn new(stem: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "losat_rpsbproc_{stem}_{}_{nanos}",
            std::process::id()
        ));
        for (from, to) in [
            ("arch", "arch"),
            ("arch/bad", "arch/bad"),
            ("data", "data"),
            ("root", ""),
        ] {
            copy_files(&fixtures().join(from), &dir.join(to));
        }
        RunDir(dir)
    }
}

/// Copies the files of `from` into `to` (made if missing); a file named `*.asn.gz` or `*.xml.gz`
/// is unpacked to `*.asn` or `*.xml` (the archives of nucleotide queries and the XML inputs are
/// stored compressed, README.md).
fn copy_files(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        if !entry.file_type().unwrap().is_file() {
            continue;
        }
        let name = entry.file_name();
        let packed = name.to_str().and_then(|name| {
            name.strip_suffix(".gz")
                .filter(|plain| plain.ends_with(".asn") || plain.ends_with(".xml"))
        });
        match packed {
            Some(plain) => {
                let mut reader =
                    flate2::read::GzDecoder::new(fs::File::open(entry.path()).unwrap());
                let mut out = fs::File::create(to.join(plain)).unwrap();
                std::io::copy(&mut reader, &mut out).unwrap();
            }
            None => {
                fs::copy(entry.path(), to.join(&name)).unwrap();
            }
        }
    }
}

impl Drop for RunDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One run of `losat rpsbproc words...` in `dir` with `stdin` (a file of the run directory, or
/// nothing): stdout, stderr, exit status. No LOSAT or NCBI variable, no registry file
/// (`NCBI_CONFIG_PATH` names the run directory, which has none).
fn rpsbproc(dir: &Path, words: &[&str], stdin: Option<&Path>) -> (Vec<u8>, Vec<u8>, i32) {
    let stdin = match stdin {
        Some(path) => Stdio::from(fs::File::open(path).unwrap()),
        None => Stdio::null(),
    };
    let output = command(dir).args(words).stdin(stdin).output().unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// `losat rpsbproc` in `dir`, without LOSAT or NCBI variables and registry files.
fn command(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .arg("rpsbproc");
    command
}

/// The rows of CASES.tsv: case, use, stdin, setup, arguments.
fn cases() -> Vec<[String; 5]> {
    let table = fs::read_to_string(fixtures().join("CASES.tsv")).unwrap();
    table
        .lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 5, "{line}");
            [0, 1, 2, 3, 4].map(|i| fields[i].to_string())
        })
        .collect()
}

/// Makes the files of a case's `setup` in `dir`, its parts separated by `,`: `data-with:<a>+<b>`,
/// the directory `<b>/` holding the data files with those of `data_with/<a>/` and then of
/// `data_with/<b>/` in their place; `data-without:<file>`, `data_missing/` without that file;
/// `conf:<file>`, `conf/<file>` as `c.ini` (for `-conffile c.ini`); `ini:<file>`, `conf/<file>`
/// as `rpsbproc.ini`; `ncbirc:<file>`, `conf/<file>` as the system-wide registry file that NCBI's
/// application reads on the platform (`SYSTEM_WIDE_REGISTRY`). The files that a run makes and the
/// registry files of the case before are removed first.
fn set_up(dir: &Path, setup: &str) {
    for made in [
        "out.txt",
        "log.txt",
        "c.ini",
        "rpsbproc.ini",
        SYSTEM_WIDE_REGISTRY,
    ] {
        let _ = fs::remove_file(dir.join(made));
    }
    let _ = fs::remove_dir_all(dir.join("data_missing"));
    for part in setup.split(',') {
        set_up_part(dir, part);
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:284-303 (cited at
// `ncbi_environment::SYSTEM_WIDE_REGISTRY`): the registry file that NCBI's application reads in
// the current directory is `.ncbirc`, and `ncbi.ini` on Windows (`eName_RcOrIni`). The fixtures were frozen on Linux
// with `.ncbirc`; a Windows run gives the file the name its NCBI reads.
const SYSTEM_WIDE_REGISTRY: &str = if cfg!(windows) { "ncbi.ini" } else { ".ncbirc" };

/// One part of `set_up`.
fn set_up_part(dir: &Path, setup: &str) {
    for (kind, name) in [
        ("conf:", "c.ini"),
        ("ini:", "rpsbproc.ini"),
        ("ncbirc:", SYSTEM_WIDE_REGISTRY),
    ] {
        if let Some(file) = setup.strip_prefix(kind) {
            fs::copy(fixtures().join("conf").join(file), dir.join(name)).unwrap();
        }
    }
    if let Some(with) = setup.strip_prefix("data-with:") {
        let layers: Vec<&str> = with.split('+').collect();
        let target = dir.join(layers[layers.len() - 1]);
        copy_files(&dir.join("data"), &target);
        for layer in layers {
            copy_files(&fixtures().join("data_with").join(layer), &target);
        }
    }
    if let Some(without) = setup.strip_prefix("data-without:") {
        fs::create_dir_all(dir.join("data_missing")).unwrap();
        for entry in fs::read_dir(dir.join("data")).unwrap() {
            let entry = entry.unwrap();
            if entry.file_name().to_str() != Some(without) {
                fs::copy(
                    entry.path(),
                    dir.join("data_missing").join(entry.file_name()),
                )
                .unwrap();
            }
        }
    }
}

/// `words`, NCBI's frozen bytes `ncbi` of a run in LOSAT's words (`losat::wording`, AGENTS.md
/// rule 8), with the version of the SESSION lines as NCBI's: rpsbproc copies the version of its
/// input into them, and the input of these cases is NCBI's archive or XML, which names NCBI's
/// version (`blastp\t2.17.0+`).
fn keep_the_sessions(ncbi: &[u8], mut words: Vec<u8>) -> Vec<u8> {
    for line in ncbi.split(|&byte| byte == b'\n') {
        if !line.starts_with(b"SESSION\t") {
            continue;
        }
        // `SESSION<tab>ordinal<tab>program<tab>version<tab>`: the rest of the line (a long
        // database name, an error that cuts it) is mapped with the other bytes.
        let Some(end) = line
            .iter()
            .enumerate()
            .filter(|&(_, &byte)| byte == b'\t')
            .map(|(at, _)| at)
            .nth(3)
        else {
            continue;
        };
        let head = &line[..=end];
        let mapped = losat::wording::losat_output_bytes(head).into_owned();
        if mapped == head {
            continue;
        }
        let find = |haystack: &[u8], needle: &[u8]| {
            haystack
                .windows(needle.len())
                .position(|window| window == needle)
        };
        // A file whose text is cut in the first line (see `file_in_losat_words`) keeps the head.
        let Some(at) = find(&words, &mapped) else {
            assert!(
                find(&words, head).is_some(),
                "the SESSION head is in the words"
            );
            continue;
        };
        words.splice(at..at + mapped.len(), head.iter().copied());
    }
    words
}

/// Compares a run's stdout, stderr and exit status with `<dir>/<case>.{out,err,rc}` by bytes
/// (the text is only shown). The frozen bytes of NCBI's are compared in LOSAT's words (NCBI's
/// bytes through `losat::wording`, AGENTS.md rule 8; LOSAT's own files are left as they are; the
/// SESSION lines keep the version of NCBI's input).
fn assert_same(dir: &Path, case: &str, out: &[u8], err: &[u8], status: i32) {
    let expected = |suffix: &str| {
        let bytes = fs::read(dir.join(format!("{case}.{suffix}"))).unwrap();
        keep_the_sessions(
            &bytes,
            losat::wording::losat_output_bytes(&bytes).into_owned(),
        )
    };
    for (name, got, want) in [
        ("stdout", out, expected("out")),
        ("stderr", err, expected("err")),
    ] {
        assert!(
            *got == want,
            "{case}: {name}\nLOSAT:\n{}\nexpected:\n{}",
            String::from_utf8_lossy(got),
            String::from_utf8_lossy(&want)
        );
    }
    let rc: i32 = String::from_utf8(expected("rc"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(status, rc, "{case}: exit status");
}

/// NCBI's frozen file in LOSAT's words (NCBI's bytes through `losat::wording`, AGENTS.md rule 8).
/// A file that `-o` and `-logfile` both name (the rows `e1c_*_same`) is NCBI's output whose end
/// the log's text has overwritten from the middle of its first line: the log's text (NCBI's two
/// lines, `Error: (<class>::<code>) <message>` and the application's `Error: (106.16) ...` line
/// naming the same exception) is rebuilt from the part that is left, mapped, and written from the
/// same place, so the part that the log keeps starts at the same byte of its text.
fn file_in_losat_words(frozen: &[u8]) -> Vec<u8> {
    const FAILED: &[u8] = b"Error: (106.16) Application's execution failed (";
    let mapped = |bytes: &[u8]| losat::wording::losat_output_bytes(bytes).into_owned();
    let second = frozen
        .windows(FAILED.len())
        .rposition(|window| window == FAILED);
    // The output (tab-delimited) ends in a tab; the log's text follows, cut in its first line.
    let tail_start = frozen
        .iter()
        .rposition(|&byte| byte == b'\t')
        .map(|at| at + 1);
    let (Some(second), Some(tail_start)) = (second, tail_start) else {
        return mapped(frozen);
    };
    if second < tail_start || frozen[tail_start..].starts_with(b"Error: ") {
        return mapped(frozen);
    }
    let tail = &frozen[tail_start..];
    let exception = &frozen[second + FAILED.len() - 1..];
    let first = [&b"Error: "[..], exception].concat();
    let text = [first, frozen[second..].to_vec()].concat();
    assert!(text.ends_with(tail), "the log's text ends the file");
    let kept = text.len() - tail.len();
    let mut file = frozen[..tail_start].to_vec();
    file.extend_from_slice(&mapped(&text)[kept..]);
    file
}

/// Runs the rows of CASES.tsv whose `use` is `usage` and compares them with `expected` (`oracle/`,
/// `widened/` or `tied/`): the number of rows.
fn compare_rows(stem: &str, usage: &str, expected: &str) -> usize {
    let run_dir = RunDir::new(stem);
    let dir = &run_dir.0;
    let oracle = fixtures().join(expected);
    let mut compared = 0;
    for [case, row_usage, stdin, setup, arguments] in cases() {
        if row_usage != usage {
            continue;
        }
        set_up(dir, &setup);
        let stdin = (stdin != "-").then(|| dir.join(&stdin));
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, stdin.as_deref());
        assert_same(&oracle, &case, &out, &err, status);
        for (suffix, made) in [("outfile", "out.txt"), ("log", "log.txt")] {
            let expected = oracle.join(format!("{case}.{suffix}"));
            let got = fs::read(dir.join(made)).ok();
            assert_eq!(
                got,
                fs::read(&expected)
                    .ok()
                    .map(|bytes| keep_the_sessions(&bytes, file_in_losat_words(&bytes))),
                "{case}: the file {made} (expected {})",
                expected.display()
            );
        }
        compared += 1;
    }
    compared
}

#[test]
fn rpsbproc_gives_ncbis_bytes() {
    assert_eq!(compare_rows("compare", "compare", "oracle"), 375);
}

// The rows where LOSAT's order of hits is its own (PD-LOSAT-NCBI-DEFECTS 1.67, the Owner's decision
// of 2026-11-20): hits equal in E-value and bit score keep their input order, where NCBI's
// `std::sort` leaves an order of its own; a NaN E-value or bit score sorts after every number of
// its key (Owner-delegated, 2026-11-20).
#[test]
fn rpsbproc_orders_tied_hits_by_its_rule() {
    assert_eq!(compare_rows("tied", "tied", "tied"), 9);
}

/// Where LOSAT's standard output of a `tied` row differs from NCBI's: the blocks and lines that
/// `tied/CLASSES.tsv` lists.
#[derive(Debug, PartialEq, Eq)]
struct TiedDifference {
    /// Query blocks whose DOMAINS lines differ only among hits of one E-value and bit score.
    tie_blocks: usize,
    /// Query blocks that hold a hit with a NaN E-value or bit score.
    nan_blocks: usize,
    /// NCBI's lines of the differing blocks that a longest common subsequence of the two blocks
    /// leaves out: the lines a diff removes.
    ncbi_lines: usize,
    /// LOSAT's lines of the differing blocks that the common subsequence leaves out: the lines a
    /// diff adds.
    losat_lines: usize,
}

/// The length of a longest common subsequence of two blocks' lines.
fn common_lines(a: &[&[u8]], b: &[&[u8]]) -> usize {
    let mut previous = vec![0; b.len() + 1];
    for line in a {
        let mut current = vec![0; b.len() + 1];
        for (j, other) in b.iter().enumerate() {
            current[j + 1] = if line == other {
                previous[j] + 1
            } else {
                previous[j + 1].max(current[j])
            };
        }
        previous = current;
    }
    previous[b.len()]
}

/// A line outside the query blocks, or a query block (`QUERY` to `ENDQUERY`).
#[derive(PartialEq)]
enum Segment<'a> {
    Line(&'a [u8]),
    Query(Vec<&'a [u8]>),
}

fn segments(text: &[u8]) -> Vec<Segment<'_>> {
    let mut segments = Vec::new();
    let mut block: Option<Vec<&[u8]>> = None;
    for line in text.split(|&byte| byte == b'\n') {
        match block.as_mut() {
            None if line.starts_with(b"QUERY\t") => block = Some(vec![line]),
            None => segments.push(Segment::Line(line)),
            Some(lines) => {
                lines.push(line);
                if line.starts_with(b"ENDQUERY") {
                    segments.push(Segment::Query(block.take().unwrap()));
                }
            }
        }
    }
    assert!(block.is_none(), "a QUERY block without its ENDQUERY");
    segments
}

fn fields(line: &[u8]) -> Vec<&[u8]> {
    line.split(|&byte| byte == b'\t').collect()
}

/// The section (`DOMAINS`, `SITES`, `MOTIFS`) of each line of a query block; `None` for the
/// markers and the lines outside the sections.
fn sections<'a>(block: &[&'a [u8]]) -> Vec<Option<&'a [u8]>> {
    let mut current = None;
    block
        .iter()
        .map(|&line| match line {
            b"DOMAINS" | b"SITES" | b"MOTIFS" => {
                current = Some(line);
                None
            }
            b"ENDDOMAINS" | b"ENDSITES" | b"ENDMOTIFS" => {
                current = None;
                None
            }
            _ => current,
        })
        .collect()
}

/// Whether a DOMAINS line of the block has a NaN E-value or bit score (`-nan`, NCBI's text).
fn holds_nan(block: &[&[u8]]) -> bool {
    let nan = |field: &[u8]| field == b"-nan" || field == b"nan";
    block.iter().zip(sections(block)).any(|(line, section)| {
        let line = fields(line);
        section == Some(&b"DOMAINS"[..]) && line.len() >= 8 && (nan(line[6]) || nan(line[7]))
    })
}

/// The query, E-value and bit score of a DOMAINS line (fields 2, 7 and 8): the key of a tie.
fn tie_key<'a>(line: &[&'a [u8]]) -> (&'a [u8], &'a [u8], &'a [u8]) {
    (line[1], line[6], line[7])
}

/// The hit types (field 3) of DOMAINS lines, sorted.
fn hit_types<'a>(lines: &[&'a [u8]]) -> Vec<&'a [u8]> {
    let mut types: Vec<&[u8]> = lines.iter().map(|&line| fields(line)[2]).collect();
    types.sort_unstable();
    types
}

/// Compares NCBI's and LOSAT's standard output of a `tied` row (README.md, "Tied hits"): equal
/// outside the query blocks, the blocks in pairs with the same `QUERY` line, and each block that
/// differs either a block with a NaN key, or one of the same length whose runs of differing lines
/// are DOMAINS lines of one query, E-value and bit score (fields 2, 7 and 8) with the same hit types
/// (field 3) on both sides. The lines of the differing blocks are counted as a diff counts them
/// (`common_lines`).
fn tied_difference(ncbi: &[u8], losat: &[u8]) -> Result<TiedDifference, String> {
    let (ncbi, losat) = (segments(ncbi), segments(losat));
    if ncbi.len() != losat.len() {
        return Err("the number of lines and query blocks".to_string());
    }
    let mut found = TiedDifference {
        tie_blocks: 0,
        nan_blocks: 0,
        ncbi_lines: 0,
        losat_lines: 0,
    };
    for (a, b) in ncbi.iter().zip(&losat) {
        let (a, b) = match (a, b) {
            _ if a == b => continue,
            (Segment::Query(a), Segment::Query(b)) => (a, b),
            _ => return Err("a line outside the query blocks".to_string()),
        };
        let name = String::from_utf8_lossy(a[0]).into_owned();
        if a[0] != b[0] {
            return Err(format!("{name}: the QUERY line"));
        }
        let differs = |i: usize| a.get(i) != b.get(i);
        let longer = a.len().max(b.len());
        let common = common_lines(a, b);
        found.ncbi_lines += a.len() - common;
        found.losat_lines += b.len() - common;
        if holds_nan(a) || holds_nan(b) {
            found.nan_blocks += 1;
            continue;
        }
        if a.len() != b.len() {
            return Err(format!("{name}: the number of lines, without a NaN"));
        }
        let sections = sections(a);
        let mut start = 0;
        while start < longer {
            if !differs(start) {
                start += 1;
                continue;
            }
            let end = (start..longer).find(|&i| !differs(i)).unwrap_or(longer);
            let run: Vec<Vec<&[u8]>> = (start..end)
                .flat_map(|i| [fields(a[i]), fields(b[i])])
                .collect();
            if sections[start..end]
                .iter()
                .any(|&section| section != Some(&b"DOMAINS"[..]))
            {
                return Err(format!("{name}: a difference outside DOMAINS"));
            }
            if run.iter().any(|line| line.len() < 8) {
                return Err(format!("{name}: a DOMAINS line of fewer than 8 fields"));
            }
            if run.iter().any(|line| tie_key(line) != tie_key(&run[0])) {
                return Err(format!(
                    "{name}: lines of different queries, E-values or bit scores"
                ));
            }
            if hit_types(&a[start..end]) != hit_types(&b[start..end]) {
                return Err(format!("{name}: the hit types of tied lines"));
            }
            start = end;
        }
        found.tie_blocks += 1;
    }
    Ok(found)
}

// Each `tied` row differs from NCBI's frozen bytes only where PD-LOSAT-NCBI-DEFECTS 1.67 says, in
// the blocks and lines that `tied/CLASSES.tsv` counts (no row is left out, and the table names each
// `tied` row once): the evidence of the approved difference kept as a check, so that `tied/` cannot
// drift from NCBI's bytes in any other way.
#[test]
fn tied_rows_differ_from_ncbi_only_in_the_order_of_tied_hits() {
    let read = |dir: &str, case: &str, suffix: &str| {
        fs::read(fixtures().join(dir).join(format!("{case}.{suffix}"))).unwrap()
    };
    // Both sides in LOSAT's words, as `assert_same` reads a case's files: NCBI's frozen bytes and
    // the `tied` files, which LOSAT wrote before it stated its own version (AGENTS.md rule 8).
    let in_losat_words = |bytes: Vec<u8>| {
        keep_the_sessions(
            &bytes,
            losat::wording::losat_output_bytes(&bytes).into_owned(),
        )
    };
    let table = fs::read_to_string(fixtures().join("tied").join("CLASSES.tsv")).unwrap();
    let mut listed = Vec::new();
    for line in table.lines().skip(1).filter(|line| !line.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 5, "{line}");
        let count = |i: usize| fields[i].parse::<usize>().unwrap();
        let (case, expected) = (
            fields[0],
            TiedDifference {
                tie_blocks: count(1),
                nan_blocks: count(2),
                ncbi_lines: count(3),
                losat_lines: count(4),
            },
        );
        for suffix in ["err", "rc"] {
            assert_eq!(
                in_losat_words(read("oracle", case, suffix)),
                in_losat_words(read("tied", case, suffix)),
                "{case}: {suffix}"
            );
        }
        let ncbi = in_losat_words(read("oracle", case, "out"));
        let losat = in_losat_words(read("tied", case, "out"));
        assert_eq!(tied_difference(&ncbi, &losat), Ok(expected), "{case}");
        listed.push(case.to_string());
    }
    let mut tied: Vec<String> = cases()
        .into_iter()
        .filter(|row| row[1] == "tied")
        .map(|row| row[0].clone())
        .collect();
    tied.sort();
    listed.sort();
    assert_eq!(listed, tied);
    assert_eq!(tied.len(), 9);
}

// The check of `tied_difference` itself: it takes a change of order, a moved hit type and another
// tied hit printed in a hit's place among hits of one E-value and bit score, and any change in a
// query with a NaN key; it refuses a move across E-values or bit scores, a changed hit type alone,
// and a change outside DOMAINS or outside the query blocks.
#[test]
fn the_tied_check_takes_only_the_order_of_tied_hits() {
    let report = |domains: &[&str], sites: &str| {
        let mut text = String::from("#DATA\nQUERY\tQuery_1\tPeptide\t9\tq\nDOMAINS\n");
        for line in domains {
            text.push_str(&line.replace(' ', "\t"));
            text.push('\n');
        }
        text.push_str("ENDDOMAINS\nSITES\n");
        text.push_str(&sites.replace(' ', "\t"));
        text.push_str("\nENDSITES\nENDQUERY\nENDDATA\n");
        text.into_bytes()
    };
    let a = "1 Query_1 Specific 1 1 9 1e-05 50 A a - -";
    let b = "1 Query_1 Non-specific 2 2 9 1e-05 50 B b - -";
    let b_specific = "1 Query_1 Specific 2 2 9 1e-05 50 B b - -";
    let a_non = "1 Query_1 Non-specific 1 1 9 1e-05 50 A a - -";
    let c = "1 Query_1 Non-specific 3 3 9 0.001 40 C c - -";
    let n = "1 Query_1 Specific 4 4 9 -nan 45 D d - -";
    let site = "1 Query_1 Specific site A1 1 1 1";
    let ncbi = report(&[a, b, c], site);
    let check = |domains: &[&str], sites: &str| tied_difference(&ncbi, &report(domains, sites));
    let counts = |tie_blocks, nan_blocks, ncbi_lines, losat_lines| {
        Ok(TiedDifference {
            tie_blocks,
            nan_blocks,
            ncbi_lines,
            losat_lines,
        })
    };
    assert_eq!(check(&[a, b, c], site), counts(0, 0, 0, 0));
    assert_eq!(check(&[b, a, c], site), counts(1, 0, 1, 1));
    assert_eq!(check(&[b_specific, a_non, c], site), counts(1, 0, 2, 2));
    assert_eq!(check(&[b_specific, b, c], site), counts(1, 0, 1, 1));
    assert!(check(&[a, c, b], site).is_err());
    assert!(check(&[a_non, b, c], site).is_err());
    assert!(check(&[a, b, c], "1 Query_1 Specific site A2 1 1 1").is_err());
    assert!(check(&[a, b], site).is_err());
    let with_nan = report(&[a, n, b, c], site);
    let nan_check = |domains: &[&str]| tied_difference(&with_nan, &report(domains, ""));
    assert_eq!(nan_check(&[a, b, c, n]), counts(0, 1, 2, 2));
    // A line changed after a moved one counts as well.
    let c_moved = "1 Query_1 Non-specific 3 3 8 0.001 40 C c - -";
    assert_eq!(nan_check(&[a, b, c_moved, n]), counts(0, 1, 3, 3));
    let mut header = ncbi.clone();
    header[1] = b'd';
    assert!(tied_difference(&ncbi, &header).is_err());
}

// The rows where NCBI's `int` coordinate arithmetic overflows (undefined; NCBI's 0.5.1 builds wrap
// it): LOSAT computes it as its value (PD-LOSAT-NCBI-DEFECTS exception 12), and its bytes equal
// NCBI's for the row's nearby input (a `compare` row) moved back by the shift that README.md names.
#[test]
fn rpsbproc_widens_ncbis_int_overflows() {
    assert_eq!(compare_rows("widened", "widened", "widened"), 9);
}

// Each `widened` row's bytes are NCBI's bytes of its nearby `compare` row moved back by the shift
// that `widened/PAIRS.tsv` names (the input file's name, and the printed numbers it lists, each of
// which must occur): the evidence of PD-LOSAT-NCBI-DEFECTS exception 12 kept as a check, so that
// `widened/` cannot drift from NCBI's nearby bytes.
#[test]
fn widened_rows_are_their_nearby_rows_moved_back() {
    let rows = cases();
    let row = |case: &str| {
        rows.iter()
            .find(|row| row[0] == case)
            .unwrap_or_else(|| panic!("{case}: no row"))
            .clone()
    };
    let input = |case: &str| {
        let arguments = row(case)[4].clone();
        let words: Vec<&str> = arguments.split(' ').collect();
        let at = words.iter().position(|&word| word == "-i").unwrap();
        words[at + 1].to_string()
    };
    let read = |dir: &str, case: &str, suffix: &str| {
        fs::read_to_string(fixtures().join(dir).join(format!("{case}.{suffix}"))).unwrap()
    };
    let table = fs::read_to_string(fixtures().join("widened").join("PAIRS.tsv")).unwrap();
    let mut pairs = 0;
    for line in table.lines().skip(1).filter(|line| !line.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (case, nearby) = (fields[0], fields[1]);
        assert_eq!(row(case)[1], "widened", "{case}");
        assert_eq!(row(nearby)[1], "compare", "{nearby}");
        let mut out = read("oracle", nearby, "out").replace(&input(nearby), &input(case));
        let mut err = read("oracle", nearby, "err").replace(&input(nearby), &input(case));
        for shift in &fields[2..] {
            let (old, new) = shift.split_once("=>").unwrap();
            let (old, new) = (old.replace("\\t", "\t"), new.replace("\\t", "\t"));
            assert!(
                out.contains(&old) || err.contains(&old),
                "{case}: {old:?} is not in {nearby}"
            );
            out = out.replace(&old, &new);
            err = err.replace(&old, &new);
        }
        assert_eq!(out, read("widened", case, "out"), "{case}: stdout");
        assert_eq!(err, read("widened", case, "err"), "{case}: stderr");
        assert_eq!(
            read("oracle", nearby, "rc"),
            read("widened", case, "rc"),
            "{case}: rc"
        );
        pairs += 1;
    }
    let widened = rows.iter().filter(|row| row[1] == "widened").count();
    assert_eq!((pairs, widened), (9, 9));
}

// The rows of `use` `linux`: `-logfile /dev/stdout` writes the diagnostics into the pipe of the
// standard output, where the report's part of a line that NCBI's `cout` keeps in glibc's buffer
// (of `st_blksize`, 4096 for a pipe on Linux) comes after them (`rpsbproc::NcbiStream::Cout`).
// NCBI's builds are Linux builds; the other platforms' buffers are not NCBI's.
#[cfg(target_os = "linux")]
#[test]
fn rpsbproc_gives_ncbis_bytes_through_a_pipe_on_linux() {
    assert_eq!(compare_rows("linux", "linux", "oracle"), 1);
}

// The refusals end LOSAT with exit 1 and the refusal as the last line of standard error, and
// leave no part of a line on standard output (the report is cut before the refused site; the
// files of `refused/` hold LOSAT's bytes, NCBI's output up to that line, README.md).
#[test]
fn rpsbproc_refuses_what_it_does_not_give() {
    let run_dir = RunDir::new("refuse");
    let dir = &run_dir.0;
    let expected = fixtures().join("refused");
    let mut refused = 0;
    for [case, usage, stdin, setup, arguments] in cases() {
        if usage != "refuse" {
            continue;
        }
        set_up(dir, &setup);
        let stdin = (stdin != "-").then(|| dir.join(&stdin));
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, stdin.as_deref());
        assert_eq!(status, 1, "{case}");
        assert!(out.is_empty() || out.ends_with(b"\n"), "{case}: stdout");
        let last = err
            .strip_suffix(b"\n")
            .and_then(|err| err.split(|&byte| byte == b'\n').next_back())
            .unwrap_or_default();
        assert!(
            last.starts_with(b"Error: ")
                && String::from_utf8_lossy(last).contains(" is not supported by LOSAT's rpsbproc"),
            "{case}: {}",
            String::from_utf8_lossy(&err)
        );
        assert_same(&expected, &case, &out, &err, status);
        refused += 1;
    }
    assert_eq!(refused, 14);
}

// The command line outside NCBI's report (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`): `-version`
// prints rpsbproc's version (NCBI: `rpsbproc: 0.5.1`), so does `-version-full` in NCBI's form
// with LOSAT's values (line RG item E3), `-help` NCBI's help in LOSAT's words, the argument
// errors are NCBI's reader's (status 1: `-dryrun=1` gives a value to a flag, `--i` is the unknown
// `-i`; line RG item K4i-2b), and a file name that is not UTF-8 is refused before anything is
// written.
#[test]
fn rpsbproc_command_line() {
    let run_dir = RunDir::new("cli");
    let dir = &run_dir.0;
    let (out, _, status) = rpsbproc(dir, &["-version"], None);
    assert_eq!(status, 0);
    let out = String::from_utf8(out).unwrap();
    assert!(
        out.starts_with(concat!(
            "rpsbproc: LOSAT ",
            env!("CARGO_PKG_VERSION"),
            "\n Package: LOSAT "
        )),
        "{out}"
    );
    let (out, _, status) = rpsbproc(dir, &["-help"], None);
    assert_eq!(status, 0);
    assert!(String::from_utf8(out)
        .unwrap()
        .contains("\n -i, --infile <String>\n"));
    let (out, err, status) = rpsbproc(
        dir,
        &["-i", "arch/one.asn", "-d", "data", "-version-full"],
        None,
    );
    assert_eq!((status, err), (0, Vec::new()));
    let out = String::from_utf8(out).unwrap();
    assert!(
        out.starts_with(concat!(
            "rpsbproc: LOSAT ",
            env!("CARGO_PKG_VERSION"),
            "\n Package: LOSAT "
        )) && out.ends_with(" Package-Config:  \n GI-64bit:  TRUE\n"),
        "{out}"
    );
    let (out, err, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-dryrun=1"], None);
    assert!(out.is_empty());
    let err = String::from_utf8(err).unwrap();
    assert_eq!(status, 1, "{err}");
    assert!(
        err.ends_with("\nError:  Invalid argument: -dryrun=1\n"),
        "{err}"
    );
    let (_, err, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-m", "xx"], None);
    let err = String::from_utf8(err).unwrap();
    assert_eq!(status, 1, "{err}");
    assert!(
        err.contains("Argument \"m\". Illegal value, expected"),
        "{err}"
    );
    let (_, err, status) = rpsbproc(dir, &["--i", "arch/one.asn"], None);
    let err = String::from_utf8(err).unwrap();
    assert_eq!(status, 1, "{err}");
    assert!(
        err.ends_with("\nError:  Unknown argument: \"-i\"\n"),
        "{err}"
    );
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let name = std::ffi::OsStr::from_bytes(b"x\xff.asn");
        for option in ["-i", "-o", "-d"] {
            let output = command(dir)
                .args(["-t", "doms", option])
                .arg(name)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(1), "{option}: {err}");
            assert!(output.stdout.is_empty());
            assert!(
                err.starts_with(&format!("Error: the -{} file name ", &option[1..]))
                    && err.contains("is not UTF-8")
                    && err.ends_with("not supported by LOSAT's rpsbproc\n"),
                "{option}: {err}"
            );
        }
    }
}

// The environment variables that change what NCBI's object stream reads or what its Seq-ids give
// are refused with the values that change it (`rpsbproc::reader_variable`; the gi parameters only
// for an archive with a gi query); the other values run.
#[test]
fn rpsbproc_refuses_reader_variables() {
    let run_dir = RunDir::new("env");
    let dir = &run_dir.0;
    let one = fs::read_to_string(dir.join("arch/one.asn")).unwrap();
    let with_gi = one.replacen(
        "local str \"Query_1\"",
        "local str \"Query_1\",\n gi 16128001",
        1,
    );
    assert_ne!(with_gi, one);
    fs::write(dir.join("gi.asn"), with_gi).unwrap();
    let run = |name: &str, value: &str, input: &str| {
        let output = command(dir)
            .env(name, value)
            .args(["-i", input, "-d", "data", "-t", "doms"])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        (
            output.status.code().unwrap(),
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    for (name, value, input) in [
        ("OBJECTS_SEQ_DESCR_ALLOW_EMPTY", "1", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_READ", "no", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_GET", "DefValue", "arch/one.asn"),
        ("SEQ_ID_AVOID_GI", "1", "gi.asn"),
        ("SEQ_ID_PREFER_ACCESSION_OVER_GI", "true", "gi.asn"),
    ] {
        let (status, err) = run(name, value, input);
        assert_eq!(status, 1, "{name}: {err}");
        assert!(
            err.contains(name) && err.contains("not supported by LOSAT's rpsbproc"),
            "{name}: {err}"
        );
    }
    for (name, value, input) in [
        ("SERIAL_VERIFY_DATA_READ", "yes", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_GET", "maybe", "arch/one.asn"),
        ("SEQ_ID_AVOID_GI", "1", "arch/one.asn"),
    ] {
        let (status, err) = run(name, value, input);
        assert_eq!(status, 0, "{name}={value}: {err}");
    }
    // A dry run reads no archive (NCBI rpsbproc 0.5.1 passes with these variables: line RG E1c
    // review 1, `r_dry_serial`).
    let output = command(dir)
        .env("SERIAL_VERIFY_DATA_READ", "no")
        .args(["-i", "arch/one.asn", "-d", "data", "-dryrun"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{err}");
    assert!(
        err.ends_with("Dryrun passed, all resources are present and loaded successfully.\n"),
        "{err}"
    );
}

// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 9: a failed write of the report stops LOSAT with
// its message and exit 4 (NCBI's rpsbproc drops the writes and exits 0), with no progress line
// after the failure.
#[cfg(target_os = "linux")]
#[test]
fn rpsbproc_stops_at_a_failed_write() {
    let run_dir = RunDir::new("full");
    let dir = &run_dir.0;
    let (out, err, status) = rpsbproc(
        dir,
        &[
            "-i",
            "arch/three.asn",
            "-d",
            "data",
            "-t",
            "doms",
            "-o",
            "/dev/full",
        ],
        None,
    );
    let err = String::from_utf8(err).unwrap();
    assert!(out.is_empty());
    assert_eq!(status, 4, "{err}");
    assert!(
        err.ends_with("For more information please refer to the README file\nrpsbproc failed to write output: No space left on device (os error 28)\n"),
        "{err}"
    );
    // A standard error whose reader has closed (NCBI: SIGPIPE): exit 4, no panic, also for a
    // refusal and for an error of the command line.
    for words in [
        &["-i", "arch/three.asn", "-d", "data", "-t", "doms"][..],
        &["-i", "arch/three.asn", "-d", "data", "-t", "feats"],
        &["-i", "arch/three.asn", "-x"],
    ] {
        let (reader, writer) = std::io::pipe().unwrap();
        drop(reader);
        let output = command(dir)
            .args(words)
            .stdin(Stdio::null())
            .stderr(Stdio::from(writer))
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(4), "{words:?}");
        assert!(output.stdout.is_empty(), "{words:?}");
    }
}

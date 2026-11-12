//! `losat dustmasker` against NCBI BLAST+ 2.17.0 dustmasker 1.0.0: each row of
//! `tests/dustmasker/cases.tsv` (its arguments as shell words, run in `tests/dustmasker`; `< file`
//! is the standard input) gives NCBI's standard output, standard error and exit status, frozen in
//! `tests/dustmasker/ncbi/<name>.{out,err,rc}` (probe of line RG, item K1, 2026-11-12 on ws). A
//! name the table repeats is compared on its last row, whose output the probe kept. The rows of
//! `NOT_YET` wait for dustmasker's masking; a row there that gives NCBI's bytes fails the test
//! until its name is removed from the list. NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The rows whose outcome needs dustmasker's work (the reader, the masking and the writers): not
/// compared yet.
const NOT_YET: &[&str] = &[
    "basic_int",
    "nruns_int",
    "ambig_int",
    "short_int",
    "deflines_int",
    "protein_int",
    "protfirst_int",
    "noheader_int",
    "gaps_int",
    "badchar_int",
    "crlf_int",
    "dup_int",
    "long_int",
    "width80_int",
    "basic_int_ps",
    "deflines_int_ps",
    "short_int_ps",
    "dup_int_ps",
    "nruns_int_ps",
    "basic_interval",
    "basic_interval_ps",
    "deflines_interval",
    "deflines_interval_ps",
    "short_interval",
    "short_interval_ps",
    "basic_fasta",
    "basic_fasta_ps",
    "deflines_fasta",
    "deflines_fasta_ps",
    "short_fasta",
    "short_fasta_ps",
    "basic_acclist",
    "basic_acclist_ps",
    "deflines_acclist",
    "deflines_acclist_ps",
    "short_acclist",
    "short_acclist_ps",
    "basic_maskinfo_asn1_text",
    "basic_maskinfo_asn1_text_ps",
    "deflines_maskinfo_asn1_text",
    "deflines_maskinfo_asn1_text_ps",
    "short_maskinfo_asn1_text",
    "short_maskinfo_asn1_text_ps",
    "basic_maskinfo_asn1_bin",
    "basic_maskinfo_asn1_bin_ps",
    "deflines_maskinfo_asn1_bin",
    "deflines_maskinfo_asn1_bin_ps",
    "short_maskinfo_asn1_bin",
    "short_maskinfo_asn1_bin_ps",
    "basic_maskinfo_xml",
    "basic_maskinfo_xml_ps",
    "deflines_maskinfo_xml",
    "deflines_maskinfo_xml_ps",
    "short_maskinfo_xml",
    "short_maskinfo_xml_ps",
    "basic_seqloc_asn1_text",
    "basic_seqloc_asn1_text_ps",
    "deflines_seqloc_asn1_text",
    "deflines_seqloc_asn1_text_ps",
    "short_seqloc_asn1_text",
    "short_seqloc_asn1_text_ps",
    "basic_seqloc_asn1_bin",
    "basic_seqloc_asn1_bin_ps",
    "deflines_seqloc_asn1_bin",
    "deflines_seqloc_asn1_bin_ps",
    "short_seqloc_asn1_bin",
    "short_seqloc_asn1_bin_ps",
    "basic_seqloc_xml",
    "basic_seqloc_xml_ps",
    "deflines_seqloc_xml",
    "deflines_seqloc_xml_ps",
    "short_seqloc_xml",
    "short_seqloc_xml_ps",
    "basic_fasta_hard",
    "ambig_fasta",
    "nruns_fasta",
    "nruns_fasta_hard",
    "gaps_fasta",
    "width80_fasta",
    "empty_mask_text",
    "empty_mask_bin",
    "empty_mask_xml",
    "nomask_mask_text",
    "win_0",
    "win_1",
    "win_2",
    "win_3",
    "win_4",
    "win_7",
    "win_8",
    "win_9",
    "win_32",
    "win_63",
    "win_64",
    "win_65",
    "win_100",
    "win_1000",
    "win_-1",
    "lev_0",
    "lev_1",
    "lev_2",
    "lev_4",
    "lev_5",
    "lev_10",
    "lev_19",
    "lev_21",
    "lev_64",
    "lev_65",
    "lev_100",
    "lev_1000",
    "lev_-1",
    "lnk_0",
    "lnk_1",
    "lnk_2",
    "lnk_10",
    "lnk_32",
    "lnk_33",
    "lnk_100",
    "lnk_-1",
    "win100_long",
    "win200_long",
    "lev10_long",
    "lnk5_long",
    "lnk5_nruns",
    "win20_nruns",
    "stdin_int",
    "stdin_dash",
    "infmt_bad",
    "nr_win_0",
    "nr_win_1",
    "nr_win_7",
    "nr_win_8",
    "nr_win_20",
    "nr_win_63",
    "nr_win_64",
    "nr_win_65",
    "nr_win_100",
    "nr_win_1000",
    "nr_win_-1",
    "nr_win_-64",
    "nr_win_2147483647",
    "nr_lnk_0",
    "nr_lnk_1",
    "nr_lnk_2",
    "nr_lnk_5",
    "nr_lnk_32",
    "nr_lnk_33",
    "nr_lnk_100",
    "nr_lnk_-1",
    "nr_lnk_-2",
    "nr_lnk_2147483647",
    "nr_win0_lnk0",
    "nr_win10_lnk3",
    "big_win_-2147483648",
    "big_win_+5",
    "big_win_0x10",
    "big_win_010",
    "opts_raw_text",
    "opts_raw_xml",
    "protein_fasta",
    "badchar_fasta",
    "lower_fasta_hard",
    "nruns_maskinfo_text",
    "protein_seqloc_text_ps",
    "deflines_seqloc_xml_ps2",
    "in_dir",
    "out_dash",
    "infmt_upper",
    "infmt_empty",
    "mdefl_interval",
    "mdefl_interval_ps",
    "mdefl_acclist",
    "mdefl_acclist_ps",
    "mdefl_maskinfo_asn1_text",
    "mdefl_maskinfo_asn1_text_ps",
    "mdefl_seqloc_xml",
    "mdefl_seqloc_xml_ps",
    "mdefl_maskinfo_asn1_bin",
    "mdefl_maskinfo_asn1_bin_ps",
];

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/dustmasker")
}

/// One row of the table: its name, the words after `dustmasker`, and the file of `< file`.
struct Case {
    name: String,
    words: Vec<String>,
    stdin: Option<String>,
}

/// The words of `line` as the probe's shell split them: blanks between words, double quotes
/// around a word with blanks or none, and `<` before the file of the standard input (no other
/// shell syntax is in the table).
fn shell_words(line: &str) -> (Vec<String>, Option<String>) {
    let mut words = Vec::new();
    let mut word: Option<String> = None;
    let mut quoted = false;
    for c in line.chars() {
        match c {
            '"' => {
                quoted = !quoted;
                word.get_or_insert_with(String::new);
            }
            ' ' | '\t' if !quoted => {
                if let Some(done) = word.take() {
                    words.push(done);
                }
            }
            _ => word.get_or_insert_with(String::new).push(c),
        }
    }
    assert!(!quoted, "unbalanced quotes: {line}");
    words.extend(word);
    let mut stdin = None;
    if let Some(at) = words.iter().position(|word| word == "<") {
        stdin = Some(words[at + 1].clone());
        words.drain(at..at + 2);
    }
    (words, stdin)
}

/// The rows of the table, the last of each name.
fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut cases: Vec<Case> = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (name, arguments) = line.split_once('\t').unwrap();
        let (words, stdin) = shell_words(arguments);
        cases.retain(|case| case.name != name);
        cases.push(Case {
            name: name.to_string(),
            words,
            stdin,
        });
    }
    cases
}

/// One run of `losat dustmasker` for `case` in `tests/dustmasker`: stdout, stderr, exit status.
/// No LOSAT, NCBI or BLAST variable, no registry file (`NCBI_CONFIG_PATH` names the fixture
/// directory, which has none).
fn run(case: &Case) -> (Vec<u8>, Vec<u8>, i32) {
    let dir = fixtures();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    let stdin = match &case.stdin {
        Some(file) => Stdio::from(fs::File::open(dir.join(file)).unwrap()),
        None => Stdio::null(),
    };
    let output = command
        .env("NCBI_CONFIG_PATH", &dir)
        .current_dir(&dir)
        .arg("dustmasker")
        .args(&case.words)
        .stdin(stdin)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// NCBI's stdout, stderr and exit status for the row `name`.
fn expected(name: &str) -> (Vec<u8>, Vec<u8>, i32) {
    let file = |extension: &str| fixtures().join("ncbi").join(format!("{name}.{extension}"));
    let status = fs::read_to_string(file("rc")).unwrap();
    (
        fs::read(file("out")).unwrap(),
        fs::read(file("err")).unwrap(),
        status.trim().parse().unwrap(),
    )
}

#[test]
fn dustmasker_gives_ncbis_bytes() {
    let cases = cases();
    for name in NOT_YET {
        assert!(
            cases.iter().any(|case| case.name == *name),
            "{name} of NOT_YET is no row of cases.tsv"
        );
    }
    let mut differ = Vec::new();
    let mut now_same = Vec::new();
    for case in &cases {
        let same = run(case) == expected(&case.name);
        match (NOT_YET.contains(&case.name.as_str()), same) {
            (false, false) => differ.push(case.name.as_str()),
            (true, true) => now_same.push(case.name.as_str()),
            _ => {}
        }
    }
    assert!(
        differ.is_empty() && now_same.is_empty(),
        "rows that differ from NCBI: {differ:?}; rows of NOT_YET that give NCBI's bytes (remove \
         them from the list): {now_same:?}"
    );
}

#[test]
fn the_shell_words_of_the_table_are_read_as_the_shell_reads_them() {
    let words = |line: &str| shell_words(line);
    assert_eq!(
        words("-in in/basic.fa -window \" 5\""),
        (
            vec![
                "-in".to_string(),
                "in/basic.fa".into(),
                "-window".into(),
                " 5".into()
            ],
            None
        )
    );
    assert_eq!(
        words("-outfmt \"\""),
        (vec!["-outfmt".to_string(), String::new()], None)
    );
    assert_eq!(
        words("-in - < in/basic.fa"),
        (
            vec!["-in".to_string(), "-".into()],
            Some("in/basic.fa".into())
        )
    );
}

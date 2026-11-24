//! `losat makembindex` against NCBI BLAST+ 2.17.0 makembindex: each row of
//! `tests/makembindex/cases.tsv` runs in its own copy of `tests/makembindex/in` (with `gen/`, the
//! generated inputs, beside it when the row names them), and gives NCBI's standard output,
//! standard error and exit status (`ncbi/<name>.{out,err,rc}`) and the size and FNV-1a digest of
//! each file the run made (`ncbi/<name>.files`), frozen by line RG item D5a (2026-11-21) with the
//! local official release. NCBI is not run here. Every row runs through the library API too
//! (`losat::api::run_makembindex`, by `examples/library_command_line`).
//!
//! `gen/` holds what the freeze made from the generator below (`gen_big.py` of the probe): two
//! FASTA files, the databases that `makeblastdb` builds of them (here LOSAT's, there NCBI's: an
//! index reads only their sequences and masks) and the lower case masks of `big.fa`. `spill/`
//! holds one more FASTA file of the generator, apart from `gen/` (28 MB, copied only for the row
//! that reads it).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::OnceLock;

use crate::library_examples::example_binary;
use crate::masker_fixtures::{shell_words, take_variables};
use crate::ncbi_windows_paths::{alias_list_path_on_windows, ncbi_paths_on_windows};

/// `tests/makembindex`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/makembindex")
}

/// One row of the table.
struct Case {
    name: String,
    words: Vec<String>,
    stdin: Option<String>,
    /// The variables of the row's first `NAME=value` words (`masker_fixtures::take_variables`).
    variables: Vec<(String, String)>,
    /// Whether the row reads `gen/`.
    generated: bool,
    /// Whether the row reads `spill/`.
    spill: bool,
}

/// The rows of the table.
fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut cases = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (name, arguments) = line.split_once('\t').unwrap_or((line, ""));
        let (mut words, stdin) = shell_words(arguments);
        let variables = take_variables(&mut words);
        cases.push(Case {
            name: name.to_string(),
            words,
            stdin,
            variables,
            generated: arguments.contains("gen/"),
            spill: arguments.contains("spill/"),
        });
    }
    cases
}

/// The letters of the generator: `x = x * 6364136223846793005 + 1442695040888963407` (mod
/// 2^64) from 1; a letter is `ACGT[x >> 62]`, and with `lower` a run of lower case letters starts
/// or ends where bits 40 to 47 of `x` are 0.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn letters(&mut self, n: usize, lower: bool) -> Vec<u8> {
        let mut low = false;
        (0..n)
            .map(|_| {
                let x = self.next();
                let letter = b"ACGT"[(x >> 62) as usize];
                if lower && (x >> 40) & 0xFF == 0 {
                    low = !low;
                }
                if low {
                    letter.to_ascii_lowercase()
                } else {
                    letter
                }
            })
            .collect()
    }
}

/// A FASTA record with lines of 80 letters.
fn record(out: &mut Vec<u8>, name: &str, letters: &[u8]) {
    out.extend_from_slice(format!(">{name}\n").as_bytes());
    for line in letters.chunks(80) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
}

/// `gen/big.fa` (5,100,000 letters with lower case runs, then 20,000) and `gen/many.fa` (600
/// records of 200 to 3,199 letters with lower case runs), each from a new generator.
fn generated_fasta() -> (Vec<u8>, Vec<u8>) {
    let mut lcg = Lcg(1);
    let mut big = Vec::new();
    record(&mut big, "big1 generated", &lcg.letters(5_100_000, true));
    record(&mut big, "big2", &lcg.letters(20_000, false));
    let mut lcg = Lcg(1);
    let mut many = Vec::new();
    for i in 0..600 {
        let n = 200 + (lcg.next() % 3000) as usize;
        let letters = lcg.letters(n, true);
        record(&mut many, &format!("m{i} generated {i}"), &letters);
    }
    (big, many)
}

/// `spill/spill.fa`: 2,200,000 records of the one letter `A` (one local id of `-stride 1` whose
/// offsets reach into the local id's bits), then `X` of 1,900,000 letters, `Y` of 100,000 and `R`
/// of 1,000,000 from a new generator (the volume of `-volsize 12` rolls back over `R`).
fn spill_fasta() -> Vec<u8> {
    let mut spill = Vec::new();
    for i in 0..2_200_000 {
        spill.extend_from_slice(format!(">t{i}\nA\n").as_bytes());
    }
    let mut lcg = Lcg(1);
    record(&mut spill, "X", &lcg.letters(1_900_000, false));
    record(&mut spill, "Y", &lcg.letters(100_000, false));
    record(&mut spill, "R", &lcg.letters(1_000_000, false));
    spill
}

/// The directory that holds `spill/`, made once.
fn spill() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let root = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join(format!("makembindex-spill-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("spill")).unwrap();
        fs::write(root.join("spill/spill.fa"), spill_fasta()).unwrap();
        root
    })
}

/// FNV-1a 64 of `bytes` (also the digest of `gi_mask_read`'s `.files`).
pub(crate) fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// A command of `binary` with no LOSAT, NCBI or BLAST variable, run in `dir`.
fn command(binary: &Path, dir: &Path) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command.env("NCBI_CONFIG_PATH", dir).current_dir(dir);
    command
}

/// The directory that holds `gen/`, made once: the FASTA files, then LOSAT's databases of them
/// with the commands of the freeze.
fn generated() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let root = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join(format!("makembindex-gen-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("gen")).unwrap();
        let (big, many) = generated_fasta();
        fs::write(root.join("gen/big.fa"), big).unwrap();
        fs::write(root.join("gen/many.fa"), many).unwrap();
        let losat = Path::new(env!("CARGO_BIN_EXE_losat"));
        for words in [
            &[
                "makeblastdb",
                "-in",
                "gen/big.fa",
                "-dbtype",
                "nucl",
                "-out",
                "gen/bigdb/big",
            ][..],
            &[
                "makeblastdb",
                "-in",
                "gen/many.fa",
                "-dbtype",
                "nucl",
                "-out",
                "gen/manydb/many",
            ],
            &[
                "convert2blastmask",
                "-in",
                "gen/big.fa",
                "-masking_algorithm",
                "other",
                "-masking_options",
                "lower case",
                "-outfmt",
                "maskinfo_asn1_bin",
                "-out",
                "gen/big.asnb",
            ],
            &[
                "makeblastdb",
                "-in",
                "gen/big.fa",
                "-dbtype",
                "nucl",
                "-mask_data",
                "gen/big.asnb",
                "-mask_id",
                "lower",
                "-out",
                "gen/mbigdb/big",
            ],
        ] {
            let output = command(losat, &root)
                .args(words)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            assert!(output.status.success(), "{words:?}: {output:?}");
        }
        root
    })
}

/// Copies the files of `from` into `to`, directories too.
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// The files below `dir`, as `/`-separated paths relative to it.
fn files_below(dir: &Path) -> Vec<String> {
    let mut found = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.push(parts.join("/"));
            }
        }
    }
    found.sort();
    found
}

/// What a run gave: stdout, stderr, exit status, and the files it made (`path\tsize\tdigest`).
type Run = (Vec<u8>, Vec<u8>, i32, String);

/// The run of `binary makembindex` for `case` in a fresh copy of the inputs (`tag` names the
/// copy); the directory, removed after the run unless `keep`.
fn run(binary: &Path, case: &Case, tag: &str) -> (Run, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "makembindex-{tag}-{}-{}",
        case.name,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    copy_tree(&fixtures().join("in"), &dir);
    if case.generated {
        copy_tree(&generated().join("gen"), &dir.join("gen"));
    }
    if case.spill {
        copy_tree(&spill().join("spill"), &dir.join("spill"));
    }
    let before = files_below(&dir);
    let stdin = match &case.stdin {
        Some(file) => Stdio::from(fs::File::open(dir.join(file)).unwrap()),
        None => Stdio::null(),
    };
    let output = command(binary, &dir)
        .envs(case.variables.iter().map(|(name, value)| (name, value)))
        .arg("makembindex")
        .args(&case.words)
        .stdin(stdin)
        .output()
        .unwrap();
    let mut files = String::new();
    for path in files_below(&dir) {
        if before.contains(&path) {
            continue;
        }
        let bytes = fs::read(dir.join(&path)).unwrap();
        files.push_str(&format!(
            "{path}\t{}\t{:016x}\n",
            bytes.len(),
            fnv1a64(&bytes)
        ));
    }
    let run = (
        output.stdout,
        output.stderr,
        output.status.code().unwrap(),
        files,
    );
    (run, dir)
}

/// The rows whose error quotes the missing `OIDLIST` of `in/olmiss.nal` or `in/olmissl.nal`
/// (`expected`).
const OIDLIST_ROWS: [&str; 3] = [
    "oidlist_new_miss_n8",
    "oidlist_old_miss_n8",
    "oidlist_old_missl_n8",
];

/// NCBI's run of `case` in `dir`: the frozen bytes with `@CWD@` read as `dir`, the standard
/// error in LOSAT's words (`losat::wording`).
fn expected(case: &Case, dir: &Path) -> Run {
    let file = |extension: &str| {
        fixtures()
            .join("ncbi")
            .join(format!("{}.{extension}", case.name))
    };
    let here = |bytes: Vec<u8>| -> Vec<u8> {
        let text = String::from_utf8(bytes).unwrap();
        let text = if cfg!(windows) {
            ncbi_paths_on_windows(&text, &["@CWD@"])
        } else {
            text
        };
        text.replace("@CWD@", &dir.to_string_lossy()).into_bytes()
    };
    let err = losat::wording::losat_bytes(&fs::read(file("err")).unwrap()).into_owned();
    // The rows whose `-input` is an alias file of `OIDLIST /losat-oidlist-fixture-missing/none.msk`
    // (`in/olmiss.nal`, `in/olmissl.nal`): on Windows NCBI joins that value to the alias file's
    // directory, the run's (`alias_list_path_on_windows`; its `DBLIST` names the alias `mvdb/mv`,
    // so the directory stays the alias file's: `seqdb::alias::compute_masks`).
    let err = if cfg!(windows) && OIDLIST_ROWS.contains(&case.name.as_str()) {
        let text = String::from_utf8(err).unwrap();
        alias_list_path_on_windows(&text, "@CWD@", "/losat-oidlist-fixture-missing/none.msk")
            .into_bytes()
    } else {
        err
    };
    (
        here(fs::read(file("out")).unwrap()),
        here(err),
        fs::read_to_string(file("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap(),
        fs::read_to_string(file("files")).unwrap(),
    )
}

/// The rows whose run by `binary` differs from NCBI's (`library`: the example, which ends every
/// failure with exit status 1).
fn rows_that_differ(binary: &Path, library: bool) -> Vec<String> {
    let tag = if library { "lib" } else { "cli" };
    let mut differ = Vec::new();
    for case in cases() {
        let (got, dir) = run(binary, &case, tag);
        let mut want = expected(&case, &dir);
        if library {
            want.2 = i32::from(want.2 != 0);
        }
        if got != want {
            let text = |bytes: &[u8]| {
                String::from_utf8_lossy(bytes)
                    .chars()
                    .take(800)
                    .collect::<String>()
            };
            eprintln!(
                "--- {}: LOSAT status {}, NCBI {}\nLOSAT stdout: {:?}\nNCBI stdout:  {:?}\n\
                 LOSAT stderr: {:?}\nNCBI stderr:  {:?}\nLOSAT files:\n{}NCBI files:\n{}",
                case.name,
                got.2,
                want.2,
                text(&got.0),
                text(&want.0),
                text(&got.1),
                text(&want.1),
                got.3,
                want.3,
            );
            differ.push(case.name.clone());
        } else {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    differ
}

#[test]
fn makembindex_gives_ncbis_bytes() {
    let differ = rows_that_differ(Path::new(env!("CARGO_BIN_EXE_losat")), false);
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn makembindex_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ(&example_binary("library_command_line"), true);
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

/// The generator writes the inputs of the freeze (the files of `gen_big.py big` and `gen_big.py
/// many`, by their sizes and FNV-1a digests).
#[test]
fn the_generated_inputs_are_the_freeze_s() {
    let (big, many) = generated_fasta();
    assert_eq!(
        (big.len(), fnv1a64(&big)),
        (5_184_022, 0x0a11_6bfc_0b6f_9996)
    );
    assert_eq!(
        (many.len(), fnv1a64(&many)),
        (1_066_698, 0x36e7_8e4b_02be_1c36)
    );
    let spill = spill_fasta();
    assert_eq!(
        (spill.len(), fnv1a64(&spill)),
        (28_326_399, 0xfdd8_1e8e_25df_89ba)
    );
}

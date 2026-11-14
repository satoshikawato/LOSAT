//! The NCBI application start of a run reads the registry files that the host registered as the
//! command line reads the files of its working directory (line RG item E1w-3a): `<program>.ini`
//! and `.ncbirc` too, not only the file of `-conffile`. A test process of its own: a registered
//! `.ncbirc` is read by every run of the process.

use losat_web_adapter::{files, run, store};

const SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

fn register(name: &str, bytes: &[u8]) -> u32 {
    let file = files::begin(files::KIND_REGISTRY, name, bytes.len()).unwrap();
    files::chunk(file, bytes).unwrap();
    files::end(file).unwrap();
    file
}

fn input(name: &str) -> Vec<u8> {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/fasta/outfmt11"
    );
    std::fs::read(format!("{root}/{name}")).unwrap()
}

fn run_blastp(words: &[&str]) -> Result<(), String> {
    let program = run::Program::parse("blastp").unwrap();
    let (query, _) = store::register("blastp", 0, &input("arch_bp_q.faa")).unwrap();
    let (subject, _) = store::register("blastp", 1, &input("arch_bp_s.faa")).unwrap();
    let result = store::with_inputs(program, query, subject, false, |queries, subjects| {
        run::run(words, queries, subjects)
    });
    store::release(query).unwrap();
    store::release(subject).unwrap();
    result.unwrap()
}

// NCBI BLAST+ 2.17.0 (g_reg probes of env parity, DECISIONS 13; line RG `E/SURVEY_E1.md`): a
// malformed `<program>.ini` ends the program with NCBI's report (exit status 2), unless `-conffile`
// names another file; a `.ncbirc` entry that LOSAT refuses is refused after the arguments, so an
// argument error comes first.
#[test]
fn the_start_reads_the_programs_ini_and_ncbirc() {
    let words = ["blastp", "-query", "q", "-subject", "s", "-outfmt", "6"];
    run::validate(&words).unwrap();
    let ini = register("blastp.ini", b"[x\n");
    let report = SECTION_REPORT.repeat(2).trim_end().to_string();
    assert_eq!(run::validate(&words).unwrap_err(), report);
    assert_eq!(run_blastp(&words).unwrap_err(), report);
    let good = register("e1w3a_other.ini", b"[BLAST]\n");
    let mut other = words.to_vec();
    other.extend(["-conffile", "e1w3a_other.ini"]);
    run::validate(&other).unwrap();
    run_blastp(&other).unwrap();
    files::release(ini).unwrap();
    files::release(good).unwrap();
    run::validate(&words).unwrap();

    let ncbirc = register(".ncbirc", b"[BLAST]\nE1W3A_UNKNOWN = 1\n");
    let refusal = run::validate(&words).unwrap_err();
    assert!(
        refusal.starts_with("the registry file ./.ncbirc sets [BLAST] E1W3A_UNKNOWN"),
        "{refusal}"
    );
    assert_eq!(run_blastp(&words).unwrap_err(), refusal);
    let mut wrong = words.to_vec();
    wrong.extend(["-evalue", "x"]);
    let error = run::validate(&wrong).unwrap_err();
    assert!(!error.starts_with("the registry file"), "{error}");
    files::release(ncbirc).unwrap();
    run::validate(&words).unwrap();
}

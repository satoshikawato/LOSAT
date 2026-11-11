//! The shared vector of web ABI v1's `-searchsp` and `-dbsize` (`tests/fixtures/web_v1_search_space/`):
//! its expected files are the bytes `losat blastn`, `losat tblastx` and `losat blastp` write for each
//! case, so `losat_web_run_pair` (`src/web_api.rs`, whose tests run on wasm32-wasip1 in the web
//! workflow) is held to the command line through them. When the vector was made (2026-11-10), every
//! case was also byte-identical to NCBI BLAST+ 2.17.0 on standard output.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/web_v1_search_space")
}

#[test]
fn the_cli_writes_the_bytes_of_the_web_v1_search_space_vector() {
    let dir = fixture();
    let root = std::env::temp_dir().join(format!("losat-web-v1-space-{}", std::process::id()));
    let cases = std::fs::read_to_string(dir.join("cases.tsv")).expect("cases");
    let mut count = 0;
    for line in cases.lines().filter(|line| !line.starts_with('#')) {
        let [name, program, query, subject, outfmt, extra] =
            line.split('\t').collect::<Vec<_>>()[..]
        else {
            panic!("six columns: {line:?}");
        };
        let work = root.join(name);
        std::fs::create_dir_all(&work).expect("work directory");
        std::fs::copy(dir.join(query), work.join("query")).expect("query");
        std::fs::copy(dir.join(subject), work.join("subject")).expect("subject");
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (variable, _) in std::env::vars_os() {
            if variable.to_string_lossy().starts_with("LOSAT_")
                && variable != "LOSAT_BLASTX_WORKER_LOG"
            {
                command.env_remove(variable);
            }
        }
        // The words of `losat_web_run_pair`: the extra words, then `-outfmt` and its value.
        let output = command
            .current_dir(&work)
            .args([program, "-query", "query", "-subject", "subject"])
            .args(extra.split_whitespace())
            .args(["-outfmt", outfmt])
            .output()
            .expect("losat");
        assert!(output.status.success(), "{name}: {output:?}");
        let expected = std::fs::read(dir.join(format!("{name}.out"))).expect(name);
        assert!(output.stdout == expected, "{name}: report differs");
        count += 1;
    }
    assert_eq!(count, 12);
    let _ = std::fs::remove_dir_all(&root);
}

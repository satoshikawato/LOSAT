//! The shared vector of web ABI v1's TBLASTN (`tests/fixtures/web_v1_tblastn/cases.tsv`): its
//! expected files are the bytes `losat tblastn` writes for each case, so `losat_web_run_pair` and
//! `losat_web_run_pair_handles` (`src/web_api.rs`, whose tests run on wasm32-wasip1 in the web
//! workflow) and `tblastn::web::run_pair` (`src/algorithm/tblastn/web.rs`) are held to the
//! command line through them. When the vector was made (2026-11-04), every case but the two
//! parser errors (LOSAT's own messages) was also byte-identical to NCBI BLAST+ 2.17.0 tblastn, on
//! standard output and standard error.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/web_v1_tblastn")
}

#[test]
fn the_cli_writes_the_bytes_of_the_web_v1_tblastn_vector() {
    let dir = fixture();
    let root = std::env::temp_dir().join(format!("losat-web-v1-tblastn-{}", std::process::id()));
    let cases = std::fs::read_to_string(dir.join("cases.tsv")).expect("cases");
    let mut count = 0;
    for line in cases.lines().filter(|line| !line.starts_with('#')) {
        let [name, query, subject, label, outfmt, extra, kind] =
            line.split('\t').collect::<Vec<_>>()[..]
        else {
            panic!("seven columns: {line:?}");
        };
        // The query is named `query` and the subject by its label, as `losat_web_run_pair`
        // (`subject`) and `losat_web_run_pair_handles` (the stored label) name them.
        let work = root.join(name);
        std::fs::create_dir_all(&work).expect("work directory");
        std::fs::copy(dir.join(query), work.join("query")).expect("query");
        std::fs::copy(dir.join(subject), work.join(label)).expect("subject");
        // The words of `losat_web_run_pair`: the extra words, then `-outfmt` and its value.
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (variable, _) in std::env::vars_os() {
            if variable.to_string_lossy().starts_with("LOSAT_")
                && variable != "LOSAT_BLASTX_WORKER_LOG"
            {
                command.env_remove(variable);
            }
        }
        let output = command
            .current_dir(&work)
            .args(["tblastn", "-query", "query", "-subject", label])
            .args(extra.split_whitespace())
            .args(["-outfmt", outfmt])
            .output()
            .expect("losat tblastn");
        let expected = std::fs::read(dir.join(format!("{name}.{kind}"))).expect(name);
        if kind == "out" {
            assert!(output.status.success(), "{name}: {output:?}");
            assert!(output.stdout == expected, "{name}: report differs");
        } else {
            assert!(!output.status.success(), "{name}: {output:?}");
            assert!(output.stdout.is_empty(), "{name}: wrote a report");
            assert!(
                output.stderr == expected,
                "{name}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
        count += 1;
    }
    assert_eq!(count, 11);
    let _ = std::fs::remove_dir_all(&root);
}

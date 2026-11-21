use super::*;

#[test]
fn commas_group_three_digits() {
    assert_eq!(with_commas(0), "0");
    assert_eq!(with_commas(999), "999");
    assert_eq!(with_commas(1_000), "1,000");
    assert_eq!(with_commas(186_063_409), "186,063,409");
}

// help.txt is NCBI's frozen `-help` text and is not edited or split: this keeps the version
// that it names (which LOSAT's output maps to its own, `crate::wording`) from drifting apart
// from NCBI's version that LOSAT's sources name elsewhere.
#[test]
fn the_frozen_help_names_the_version_that_losat_writes() {
    assert!(HELP.contains(&format!(
        "   BLAST database client, version {}\n",
        crate::NCBI_BLAST_VERSION
    )));
}

#[test]
fn the_short_usage_is_the_head_of_the_help() {
    let usage = crate::ncbi_app::short_usage(HELP);
    // NCBI's `2.17.0+` is LOSAT's `LOSAT <version>` (`crate::wording`), four bytes longer
    // than the 780 of NCBI's.
    assert_eq!(usage.len(), 784);
    assert!(usage.starts_with("USAGE\n  blastdbcmd [-h] [-help] [-db dbname]"));
    assert!(usage.ends_with(&format!(
        "   BLAST database client, version LOSAT {}\n\nUse '-help' to print detailed \
         descriptions of command line arguments\n",
        env!("CARGO_PKG_VERSION")
    )));
    assert_eq!(HELP.len(), 9450);
}

/// An output whose every write fails (a full disk).
struct FullDisk;

impl Write for FullDisk {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("no space left"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Err(std::io::Error::other("no space left"))
    }
}

/// A pipe whose reader has closed.
struct ClosedPipe;

impl Write for ClosedPipe {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The exit status and the diagnostics of `words` run with `input` as the standard input
/// and `out` as the standard output.
fn run_on(words: &[&str], input: &str, out: &mut dyn Write) -> (i32, String) {
    let words: Vec<OsString> = words.iter().map(OsString::from).collect();
    let mut err = Vec::new();
    let status = run_with(&words, &mut input.as_bytes(), out, &mut err);
    (status, String::from_utf8(err).unwrap())
}

fn fixture_db(name: &str) -> String {
    format!("{}/tests/{name}", env!("CARGO_MANIFEST_DIR"))
}

// NCBI 2.17.0 with `-out /dev/full` (review 2 of line RG item K1, F3): a failed write leaves
// the retrieval going on, so an id that is not found after it is reported as without the
// failure (`Skipped NOSUCHID`, exit status 1). The batch is the standard input (`-`).
#[test]
fn a_failed_write_leaves_the_later_entries_of_a_batch_and_their_errors() {
    let db = fixture_db("blastdbcmd/db/p5");
    let words = ["-db", db.as_str(), "-entry_batch", "-"];
    let batch = "XP_000001.1\nnosuchid\nP99999.2\n";
    let mut written = Vec::new();
    let expected = run_on(&words, batch, &mut written);
    assert_eq!(
        expected,
        (1, "Error: [blastdbcmd] Skipped NOSUCHID\n".to_string())
    );
    assert!(written.starts_with(b">"));
    assert_eq!(run_on(&words, batch, &mut FullDisk), expected);
    // A closed pipe ends NCBI's process by SIGPIPE (exit status 141) outside Windows before
    // the `Skipped` line, as `cerr`'s tie flushes the output first: the retrieval ends with
    // exit status 0 and no message (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10). On
    // Windows it is a failed write like the others.
    let closed = if cfg!(windows) {
        expected
    } else {
        (0, String::new())
    };
    assert_eq!(run_on(&words, batch, &mut ClosedPipe), closed);
}

/// Standard output and standard error into one file (`2>&1`).
#[derive(Clone)]
struct Merged(std::rc::Rc<std::cell::RefCell<Vec<u8>>>);

impl Write for Merged {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// `cerr` is tied to `cout`: NCBI flushes the output before each diagnostic, so with `2>&1`
// the sequences of the entries before a missing id come before its `Skipped` line, and an
// unfinished line before an `-outfmt` error (review of line RG item K1b, R2-F1).
#[test]
fn diagnostics_come_after_the_output_written_before_them() {
    let db = fixture_db("blastdbcmd/db/p5");
    for (words, batch, error) in [
        (
            vec!["-db", db.as_str(), "-entry_batch", "-", "-outfmt", "%a"],
            "XP_000001.1\nnosuchid\n",
            "Error: [blastdbcmd] Skipped NOSUCHID\n",
        ),
        (
            vec![
                "-db",
                db.as_str(),
                "-entry",
                "XP_000001.1",
                "-outfmt",
                "%a %Z",
            ],
            "",
            "Error: [blastdbcmd] Unrecognized format specification: '%Z'\n",
        ),
    ] {
        let mut alone = Vec::new();
        let (status, err) = run_on(&words, batch, &mut alone);
        assert_eq!(err, error, "{words:?}");
        let merged = Merged(Default::default());
        let words: Vec<OsString> = words.iter().map(OsString::from).collect();
        let status_merged = run_with(
            &words,
            &mut batch.as_bytes(),
            &mut merged.clone(),
            &mut merged.clone(),
        );
        assert_eq!(status_merged, status);
        let mut expected = alone;
        expected.extend_from_slice(error.as_bytes());
        assert_eq!(*merged.0.borrow(), expected, "{words:?}");
    }
    // A bad mask algorithm id of a batch line, between the outputs of its neighbours
    // (NCBI 2.17.0 `-ipg_batch - -outfmt "%a %g" 2>&1`, review of K1b, R3-F1).
    let db = fixture_db("blastdbcmd/db/m/pmg");
    let words: Vec<OsString> = ["-db", db.as_str(), "-ipg_batch", "-", "-outfmt", "%a %g"]
        .iter()
        .map(OsString::from)
        .collect();
    let merged = Merged(Default::default());
    let status = run_with(
        &words,
        &mut &b"1001 21\n1002 99\n1003 21\n"[..],
        &mut merged.clone(),
        &mut merged.clone(),
    );
    assert_eq!(status, 1);
    assert_eq!(
        String::from_utf8(merged.0.borrow().clone()).unwrap(),
        "XP_1001.1 1001\nError: [blastdbcmd] Invalid filtering algorithm ID: 99\nError: \
             [blastdbcmd] Skipped IPG : 1002\nXP_1003.1 1003\n"
    );
}

// NCBI 2.17.0 with `-out /dev/full` (F2): the serial stream of `-metadata` reports the
// failed write (five lines, exit status 6), which LOSAT refuses.
#[test]
fn a_failed_write_of_the_metadata_is_refused() {
    let db = fixture_db("dustmasker/db/basic");
    let words = ["-db", db.as_str(), "-metadata"];
    let mut written = Vec::new();
    assert_eq!(run_on(&words, "", &mut written), (0, String::new()));
    assert!(written.ends_with(b"}\n"));
    assert_eq!(
        run_on(&words, "", &mut FullDisk),
        (
            1,
            "Error: a failed write of the -metadata output (NCBI's serial stream reports it) \
                 is not supported by LOSAT's BLASTDBCMD\n"
                .to_string()
        )
    );
    // A closed pipe ends the program silently outside Windows, where SIGPIPE ends NCBI's.
    if !cfg!(windows) {
        assert_eq!(run_on(&words, "", &mut ClosedPipe), (0, String::new()));
    }
}

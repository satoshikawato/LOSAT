//! The `-culling_limit` option of TBLASTN as the parser reads it (the integration test
//! `tests/tblastn_culling.rs` runs the binary with the option).

use crate::cli::{try_parse_from, Commands};

fn limit(words: &[&str]) -> Result<i32, String> {
    let mut argv = vec!["losat", "tblastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(words);
    match try_parse_from(argv) {
        Ok(cli) => match cli.command {
            Commands::Tblastn(args) => Ok(args.culling_limit_value()),
            _ => unreachable!("tblastn"),
        },
        Err(error) => Err(error.to_string()),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp (CHspFilteringArgs): `-culling_limit` is an integer >= 0, default 0,
// hexadecimal accepted as NCBI's argument parser does.
#[test]
fn the_option_is_read_as_ncbi_reads_it() {
    assert_eq!(limit(&[]), Ok(0));
    assert_eq!(limit(&["-culling_limit", "0"]), Ok(0));
    assert_eq!(limit(&["-culling_limit", "1"]), Ok(1));
    assert_eq!(limit(&["-culling_limit", "0x2"]), Ok(2));
    assert_eq!(limit(&["-culling_limit", "2147483647"]), Ok(2147483647));
    for bad in ["-1", "1.5", "abc", "2147483648"] {
        let error = limit(&["-culling_limit", bad]).unwrap_err();
        assert!(!error.contains("is not supported by"), "{bad}: {error}");
    }
}

// The best-hit options are read (series NP O6a-3); with a limit the parser refuses the command
// line, as NCBI does, and not as an option that LOSAT lacks.
#[test]
fn the_best_hit_options_exclude_the_limit() {
    assert_eq!(limit(&["-best_hit_overhang", "0.1"]), Ok(0));
    let error = limit(&["-culling_limit", "1", "-best_hit_score_edge", "0.1"]).unwrap_err();
    assert!(!error.contains("is not supported by"), "{error}");
}

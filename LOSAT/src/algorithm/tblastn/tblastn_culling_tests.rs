//! The `-culling_limit` option of TBLASTN as the parser reads it (the integration test
//! `tests/tblastn_culling.rs` runs the binary with the option).

use crate::cli::{try_parse_from, Cli, Commands};

fn limit(words: &[&str]) -> Result<i32, String> {
    let mut argv = vec!["losat", "tblastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(words);
    match try_parse_from::<Cli, _, _>(argv) {
        Ok(cli) => match cli.command {
            Commands::Tblastn(args) => Ok(args.culling_limit_value()),
            _ => unreachable!("tblastn"),
        },
        Err(error) => Err(crate::cli::render_message(&error)),
    }
}

// NCBI blast_args.cpp (CHspFilteringArgs): `-culling_limit` is an integer >= 0, default 0,
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

// The best-hit options stay refused (they are not ported), with or without a limit.
#[test]
fn the_best_hit_options_stay_refused() {
    for words in [
        &["-best_hit_overhang", "0.1"][..],
        &["-culling_limit", "1", "-best_hit_score_edge", "0.1"],
    ] {
        let error = limit(words).unwrap_err();
        assert!(
            error.contains("is not supported by LOSAT's TBLASTN"),
            "{words:?}: {error}"
        );
    }
}

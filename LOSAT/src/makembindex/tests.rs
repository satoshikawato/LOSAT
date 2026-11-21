//! Unit tests of makembindex's command line (the index bytes are compared with NCBI's in
//! `tests/makembindex_fixtures.rs`).

use super::*;

/// `CArgAllow_Integers::GetUsage` of makembindex's constraints, as its `-help` prints them.
#[test]
fn integer_constraints_print_ncbi_usage() {
    let usage = |name: &str| {
        DESCS
            .iter()
            .find(|desc| desc.name() == name)
            .and_then(|desc| desc.constraint_usage())
            .unwrap()
    };
    assert_eq!(usage("nmer"), "8..15");
    for name in ["ws_hint", "volsize", "stride"] {
        assert_eq!(usage(name), "greater or equal to 1");
        assert!(HELP.contains(&format!(" -{name} <Integer, greater or equal to 1>\n")));
    }
    assert!(HELP.contains(" -nmer <Integer, 8..15>\n"));
}

//! The executable called by a program's name (`entry_argv`). The shims that `losat
//! --install-shims` links under these names are run in `tests/entry_names.rs`.

use super::{entry_argv, ENTRY_NAMES};

#[test]
fn the_name_of_a_program_is_its_subcommand() {
    let words = |argv: &[&str]| -> Vec<String> {
        entry_argv(argv.iter().copied())
            .iter()
            .map(|word| word.to_str().unwrap().to_string())
            .collect()
    };
    let losat_name = format!("losat{}", std::env::consts::EXE_SUFFIX);
    for name in ENTRY_NAMES {
        let expected = [losat_name.as_str(), name, "-query", "q"];
        assert_eq!(words(&[name, "-query", "q"]), expected);
        assert_eq!(
            words(&[&format!("/usr/bin/{name}"), "-query", "q"]),
            expected
        );
        assert_eq!(words(&[&format!("{name}.exe"), "-query", "q"]), expected);
        assert_eq!(words(&[&format!("{name}.EXE"), "-query", "q"]), expected);
        assert_eq!(words(&[name]), [losat_name.as_str(), name]);
    }
    // Every other name is `losat`, whatever its words are.
    for argv0 in [
        "losat",
        "/opt/losat",
        "losat.wasm",
        "blastn.1",
        "my-blastn",
        "blastn2",
    ] {
        assert_eq!(words(&[argv0, "blastn", "-h"]), [argv0, "blastn", "-h"]);
    }
    assert!(entry_argv(Vec::<String>::new()).is_empty());
}

//! The *validation* of an argv (docs/web/abi_v2.md §4 `losat_web2_validate`) gives the data loaders
//! of its registry also when an option of the argv is refused. A test process of its own: the
//! registry files that it registers are read by every start of the process.

use losat_web_adapter::{files, run};

fn register(name: &str, bytes: &[u8]) -> u32 {
    let file = files::begin(files::KIND_REGISTRY, name, bytes.len()).unwrap();
    files::chunk(file, bytes).unwrap();
    files::end(file).unwrap();
    file
}

// line WQA's D13 round 3, L1, 2026-11-20: the data loaders that the registry of `-conffile` sets
// are known as soon as the start has read it, also when an option of the argv is refused after it
// (`run::validate_with`, which `losat_web2_validate` emits from), so that the host reads the inputs
// with the setting that the user chose while the argv is still invalid.
#[test]
fn the_setting_of_the_data_loaders_is_known_when_an_option_is_refused() {
    let none = register("d13r3_none.ini", b"[BLAST]\nDATA_LOADERS=none\n");
    let words = [
        "blastp",
        "-query",
        "q",
        "-subject",
        "s",
        "-conffile",
        "d13r3_none.ini",
        "-evalue",
        "x",
    ];
    let mut known = Vec::new();
    let error = run::validate_with(&words, &mut |validation| known.push(validation)).unwrap_err();
    assert_eq!(error, run::validate(&words).unwrap_err());
    assert_eq!(
        known,
        [run::Validation {
            data_loaders: false,
            parse_deflines: None
        }]
    );
    // A valid argv gives the same setting, once.
    known.clear();
    assert_eq!(
        run::validate_with(&words[..7], &mut |validation| known.push(validation)).unwrap(),
        run::Validation {
            data_loaders: false,
            parse_deflines: Some(false)
        }
    );
    assert_eq!(
        known,
        [run::Validation {
            data_loaders: false,
            parse_deflines: Some(false)
        }]
    );
    // `-parse_deflines` is reported alongside `data_loaders` when an option check refuses the argv.
    known.clear();
    let with_pd_bad_evalue = [
        "blastp",
        "-query",
        "q",
        "-subject",
        "s",
        "-conffile",
        "d13r3_none.ini",
        "-parse_deflines",
        "-evalue",
        "0",
    ];
    run::validate_with(&with_pd_bad_evalue, &mut |validation| {
        known.push(validation)
    })
    .unwrap_err();
    assert_eq!(
        known,
        [run::Validation {
            data_loaders: false,
            parse_deflines: Some(true)
        }]
    );
    // Finding 3 (line RG-XT item X7-5 round R2): when an argv with `-parse_deflines` does not
    // parse (mid-edit option), `parse_deflines` is `None` while `data_loaders` is still reported,
    // so the host keeps its previous `-parse_deflines` setting without rescanning.
    for tail in [
        &["-evalue"][..],
        &["-evalu", "1"],
        &["-evalue", "1e-"],
        &["-word_size", "1"],
    ] {
        known.clear();
        let unparseable: Vec<&str> = [
            "blastp",
            "-query",
            "q",
            "-subject",
            "s",
            "-conffile",
            "d13r3_none.ini",
            "-parse_deflines",
        ]
        .into_iter()
        .chain(tail.iter().copied())
        .collect();
        run::validate_with(&unparseable, &mut |validation| known.push(validation)).unwrap_err();
        assert_eq!(
            known,
            [run::Validation {
                data_loaders: false,
                parse_deflines: None
            }],
            "{tail:?}"
        );
    }
    // `-parse_deflines` alone gives `Some(true)`; `-help -parse_deflines` gives `Some(false)`.
    known.clear();
    let pd_ok = [
        "blastp",
        "-query",
        "q",
        "-subject",
        "s",
        "-conffile",
        "d13r3_none.ini",
        "-parse_deflines",
    ];
    assert_eq!(
        run::validate_with(&pd_ok, &mut |validation| known.push(validation)).unwrap(),
        run::Validation {
            data_loaders: false,
            parse_deflines: Some(true)
        }
    );
    assert_eq!(
        known,
        [run::Validation {
            data_loaders: false,
            parse_deflines: Some(true)
        }]
    );
    known.clear();
    let help_pd = [
        "blastp",
        "-query",
        "q",
        "-subject",
        "s",
        "-conffile",
        "d13r3_none.ini",
        "-help",
        "-parse_deflines",
    ];
    assert_eq!(
        run::validate_with(&help_pd, &mut |validation| known.push(validation)).unwrap(),
        run::Validation {
            data_loaders: false,
            parse_deflines: Some(false)
        }
    );
    assert_eq!(
        known,
        [run::Validation {
            data_loaders: false,
            parse_deflines: Some(false)
        }]
    );
    // A start that fails (a -conffile that names no registered file) gives none.
    known.clear();
    let missing = [
        "blastp",
        "-query",
        "q",
        "-subject",
        "s",
        "-conffile",
        "d13r3_missing.ini",
    ];
    run::validate_with(&missing, &mut |validation| known.push(validation)).unwrap_err();
    assert!(known.is_empty());
    files::release(none).unwrap();
}

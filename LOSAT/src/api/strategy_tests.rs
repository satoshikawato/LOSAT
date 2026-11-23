use super::*;
use crate::api::{run, run_collect, run_to_file, validate, Inputs, Program};
use std::cell::Cell;

thread_local! {
    /// The calls of [`apply`] (the engine's import of a strategy file) made on this thread.
    pub(super) static APPLIED: Cell<usize> = const { Cell::new(0) };
}

fn strategy_bytes() -> Vec<u8> {
    std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/strategy/import/ss/blastp.sub.ss"
    ))
    .expect("a strategy file")
}

/// A strategy is imported once, when its search reads it (`Search::read_strategy`): every call
/// of that search (`run`, `run_collect`, `run_to_file`) takes a copy of the arguments of that
/// import, and each call writes the import's warnings as the command line writes them, once per
/// call (R2-N3 of review 2 of line RG item F2).
#[test]
fn a_strategy_is_imported_once_for_all_the_calls_of_its_search() {
    APPLIED.with(|count| count.set(0));
    let search = validate(Program::Blastp, &[]).unwrap();
    let strategy = search.read_strategy(&strategy_bytes()).unwrap();
    assert_eq!(APPLIED.with(Cell::get), 1, "the reading imports once");
    let mut reports = Vec::new();
    let mut warned = Vec::new();
    for _ in 0..3 {
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        run(
            &search,
            Inputs::strategy(&strategy),
            &mut [("0", &mut report)],
            &mut diagnostics,
        )
        .unwrap();
        reports.push(report);
        warned.push(diagnostics);
    }
    assert!(!reports[0].is_empty());
    assert!(
        !warned[0].is_empty(),
        "the strategy's own warning is written by each call"
    );
    assert!(
        reports.iter().all(|report| *report == reports[0])
            && warned.iter().all(|diagnostics| *diagnostics == warned[0]),
        "each call has the report and the warnings of the first"
    );
    let collected = run_collect(&search, Inputs::strategy(&strategy), &["0"]).unwrap();
    assert!(!collected.reports[0].bytes.is_empty());
    // `run_to_file`'s `-out` is the only word of the call that its arguments take.
    let path = std::env::temp_dir().join(format!("losat-n3-{}.txt", std::process::id()));
    let mut diagnostics = Vec::new();
    run_to_file(
        &search,
        Inputs::strategy(&strategy),
        "0",
        &path,
        &mut diagnostics,
    )
    .unwrap();
    let written = std::fs::read(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(!written.is_empty());
    assert!(
        !diagnostics.is_empty(),
        "the warning is written to the call's diagnostics"
    );
    assert_eq!(
        APPLIED.with(Cell::get),
        1,
        "no call of the search imports the strategy again"
    );
}

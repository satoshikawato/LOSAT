use super::*;

// Line RG series D, D-SO-5: the loop writes each batch's reports where `check` is `Ok`, holds them
// and writes them one batch at a time after the loop where only the debugging output keeps them
// there (the bytes of the stream, its refusals included), and writes them in one call for every
// other reason (outputs, formats and one batch that are not written per batch).
#[test]
fn the_reports_are_held_one_batch_at_a_time_only_for_the_debugging_output() {
    assert_eq!(Writes::of(Ok(())), Writes::Each);
    assert_eq!(Writes::of(Err(Off::Debug)), Writes::Held);
    for off in [Off::OneBatch, Off::Outputs, Off::Format, Off::CustomTabular] {
        assert_eq!(Writes::of(Err(off)), Writes::Whole, "{off:?}");
    }
    assert_eq!(Writes::default(), Writes::Whole);
    // The reasons `[X_STATS] TBXRPTSTREAM off=` can name: the refusals of the batches no longer
    // keep the reports after the loop.
    let names = [
        Off::OneBatch,
        Off::Outputs,
        Off::Format,
        Off::CustomTabular,
        Off::Debug,
    ]
    .map(Off::name);
    assert_eq!(
        names,
        ["one_batch", "outputs", "format", "custom_tabular", "debug"]
    );
}

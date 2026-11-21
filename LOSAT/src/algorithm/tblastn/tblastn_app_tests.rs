use super::*;

// A read that failed during the search fails its batch with
// NCBI's engine error (exit status 3) when the reader's error is NCBI's own, and stops the run
// with LOSAT's refusal or any other error, as the check's error does.
#[test]
fn a_read_failed_in_the_search_fails_as_the_check_s_error() {
    let engine = batch_failure(
        XReadError(PreliminaryReadError("Database memory map file error").into()).into(),
    )
    .expect("the engine error");
    let native = engine
        .downcast::<crate::cli::NativeError>()
        .expect("a NativeError");
    assert_eq!(native.exit, 3);
    assert_eq!(
        native.message,
        b"LOSAT engine error: Database memory map file error\n"
    );
    let refusal = batch_failure(XReadError(anyhow::anyhow!("refused")).into()).unwrap_err();
    assert_eq!(refusal.to_string(), "refused");
    let other = batch_failure(anyhow::anyhow!("other")).unwrap_err();
    assert_eq!(other.to_string(), "other");
}

// NCBI 2.17.0 tblastn with `-searchsp` or `OLD_FSC` and `-subject` of three empty records
// (env parity probe f_probe of 2026-11-12): the batch fails with NCBI's engine error.
#[test]
fn a_short_average_subject_length_fails_the_batch_with_the_engine_error() {
    let native = batch_failure(SetupEngineError::SUBJECT_LENGTH_INVALID.into())
        .expect("the engine error")
        .downcast::<crate::cli::NativeError>()
        .expect("a NativeError");
    assert_eq!(native.exit, 3);
    assert_eq!(
        native.message,
        b"LOSAT engine error: The average subject length is too short\n"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
// A statistics count of 0 without the Gumbel block (NCBI's division by 0, D1a-5 of line RG)
// fails the batch after the reports of the batches before, with LOSAT's refusal itself (not
// an engine error).
#[test]
fn a_setup_count_of_zero_fails_the_batch_with_the_refusal() {
    let refusal = crate::algorithm::common::database::check_setup_count(100, 0, false, "TBLASTN")
        .expect_err("the refusal");
    let failed = batch_failure(SetupRefusal(refusal).into()).expect("the batch's failure");
    assert!(failed.downcast_ref::<crate::cli::NativeError>().is_none());
    assert!(failed
        .to_string()
        .starts_with(crate::algorithm::common::database::STATISTICS_COUNT_ZERO));
    assert!(crate::algorithm::common::database::check_setup_count(100, 0, true, "TBLASTN").is_ok());
    assert!(crate::algorithm::common::database::check_setup_count(0, 0, false, "TBLASTN").is_ok());
}

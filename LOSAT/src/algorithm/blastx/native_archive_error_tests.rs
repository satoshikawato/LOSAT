use super::*;
use crate::report::serial::encoder::SerialError;

/// The encoder's refusal stays LOSAT's own error: no `BLAST engine error:` text and no
/// error archive. NCBI's engine errors still become both.
#[test]
fn the_encoders_refusal_gets_no_error_archive() {
    let refusal = archive_error(SerialError::NonVisibleChar(0xc3).into());
    assert!(refusal.is::<SerialError>());
    assert!(archive::error_lines(&refusal, false).is_none());
    assert!(archive::error_lines(&refusal, true).is_none());
    let engine_error = archive_error(anyhow::anyhow!("something else"));
    assert!(engine_error.downcast_ref::<NativeError>().is_some());
    assert!(archive::error_lines(&engine_error, false).is_some());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:209-214 (the
// subjects' `GetAllSeqs` turns `eNoDefline` into a `CInputException`; its text ends with the
// NUL of `ends`, which the error archive does not keep); a query keeps the reader's error.
#[test]
fn a_subject_without_a_first_defline_is_the_input_error_of_get_all_seqs() {
    let subject = subject_input_error(NoDeflineError(4).into());
    let native = subject.downcast_ref::<NativeError>().unwrap();
    assert_eq!(
        (native.exit, native.message.as_slice()),
        (
            1,
            b"LOSAT query/options error: Query input doesn't start with a defline or comment, \
                  line 4\0\nPlease refer to the LOSAT user manual.\n"
                .as_slice()
        )
    );
    let texts = |error: &anyhow::Error| -> Vec<Vec<u8>> {
        archive::error_lines(error, true)
            .unwrap()
            .into_iter()
            .map(|(_, text)| text)
            .collect()
    };
    assert_eq!(
        texts(&subject),
        [
            b"LOSAT query/options error: Query input doesn't start with a defline or comment, line 4".to_vec(),
            b"Please refer to the LOSAT user manual.".to_vec(),
        ]
    );
    let query = input_error(NoDeflineError(4).into());
    let native = query.downcast_ref::<NativeError>().unwrap();
    assert_eq!(
        (native.exit, native.message.as_slice()),
        (
            1,
            b"LOSAT query error: CFastaReader: Input doesn't start with a defline or comment \
                  around line 4\n"
                .as_slice()
        )
    );
    // Every other reader error of the subjects is the reader's (`BLAST query error`).
    let other = subject_input_error(FastaParseError("CFastaReader: x".into()).into());
    assert_eq!(
        other.downcast_ref::<NativeError>().unwrap().message,
        b"LOSAT query error: CFastaReader: x\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-944
// A BLASTX report with `-sorthits` whose one subject (the table's last row) totals more than
// 99999 bits: the report's `TotalScoreRefusal`, which `report_error` passes on. It is not
// NCBI's "BLAST failed to write output" (exit status 6), so `exit_on_native_error` leaves it
// and `main` writes "Error: " and the refusal's text and exits with status 1, as BLASTN,
// BLASTP, TBLASTN and TBLASTX do.
#[test]
fn the_sorthits_refusal_is_not_a_write_failure() {
    use crate::report::pairwise::{
        write_blastx_pairwise_report, BlastpPairwiseQuery, BlastpPairwiseReport,
        BlastxPairwiseOptions, PairwiseConfig, PairwiseHit,
    };
    let hit = PairwiseHit::from(crate::common::Hit {
        identity: 100.0,
        length: 50_000,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 150_000,
        s_start: 1,
        s_end: 50_000,
        e_value: 0.0,
        bit_score: 100_500.0,
        num_ident: 50_000,
        query_frame: 1,
        query_length: 150_000,
        q_idx: 0,
        s_idx: 0,
        raw_score: 260_000,
        sort_query_offset: 0,
        sort_query_end: 50_000,
        sort_subject_offset: 0,
        sort_subject_end: 50_000,
        has_sort_offsets: true,
        gap_info: None,
        num_positives: 50_000,
    });
    let query = BlastpPairwiseQuery {
        query_name: b"q".to_vec(),
        query_length: 150_000,
        ungapped_karlin: Default::default(),
        effective_search_space: 1,
        valid: true,
        batch_skipped: false,
    };
    let report = BlastpPairwiseReport {
        version: String::new(),
        database_name: Vec::new(),
        database_num_sequences: 1,
        database_total_letters: 50_000,
        matrix_name: Some("BLOSUM62".into()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 12.0,
        window_size: 40,
        gapped_karlin: Default::default(),
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: 500,
        num_alignments: 250,
    };
    let write = |hits: Option<u8>| {
        let mut out = Vec::new();
        let config = PairwiseConfig {
            program: "blastx".into(),
            sort: crate::report::hit_sort::ReportSort { hits, hsps: None },
            ..PairwiseConfig::default()
        };
        write_blastx_pairwise_report(
            std::slice::from_ref(&hit),
            &mut out,
            None,
            &config,
            std::slice::from_ref(&query),
            &[crate::report::pairwise::DescriptionLabel::title_only(
                b"s0".to_vec(),
            )],
            &report,
            &BlastxPairwiseOptions {
                word_threshold: 12.0,
                gapped: true,
                sum_stats: false,
                num_descriptions: 500,
                num_alignments: 0,
            },
            None,
            &[],
            &[],
            None,
            None,
            0,
            None,
            None,
            None,
            &crate::report::db_subject::HitLoads::default(),
        )
        .map(|()| out)
    };
    let error = report_error(write(Some(1)).unwrap_err().into());
    assert!(error.downcast_ref::<NativeError>().is_none(), "{error:?}");
    assert!(crate::report::pairwise::is_total_score_refusal(&error));
    let text = "-sorthits: the total score of the description table's last subject is wider \
                    than its column (more than 99999 bits), where NCBI BLAST+ 2.17.0 writes \
                    spaces without end (CAlignFormatUtil::AddSpace of a negative size); run \
                    without -sorthits";
    assert_eq!(error.to_string(), text);
    // No cause after the text (`main`'s "Caused by:" lines).
    assert_eq!(error.chain().count(), 1);
    // Without `-sorthits` the table is written (the score column takes the total's width).
    assert!(write(None).is_ok());
    // An error of the output stream stays NCBI's write failure.
    let failure = report_error(io::Error::other("disk full").into());
    assert_eq!(
        failure.downcast_ref::<NativeError>().map(|e| e.exit),
        Some(6)
    );
}

// Line DBID pfoff-1: `CATCH_ALL`'s report of a `CSeqDBException` that an outfmt 13 or 14
// report raises (the lookup of a reported subject that throws) stays as it is, not an
// engine error.
#[test]
fn a_database_error_of_catch_all_is_no_engine_error() {
    let thrown = crate::ncbi_exception::NcbiException::new(
        crate::seqdb::lookup_throw_site(crate::seqdb::tax_files::NO_ACCESSION_INFO.as_bytes())
            .unwrap(),
        crate::seqdb::tax_files::NO_ACCESSION_INFO,
    );
    let error = report_error(NativeError::exception(&thrown).into());
    let error = error.downcast_ref::<NativeError>().unwrap();
    assert_eq!(error.exit, 2);
    assert_eq!(
        error.message,
        b"LOSAT database error: DB contains no accession info.\n"
    );
}

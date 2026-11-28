use super::*;

#[test]
fn search_args_are_the_repeat_options() {
    let filter = RepeatFilter::new(
        std::ffi::OsString::from("rpt"),
        crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    );
    let args = filter.search_args("plus", 4).expect("valid arguments");
    assert_eq!(args.db.as_deref(), Some(std::ffi::OsStr::new("rpt")));
    assert_eq!(args.task, "blastn");
    assert_eq!((args.reward, args.penalty), (Some(1), Some(-1)));
    assert_eq!((args.gap_open, args.gap_extend), (Some(2), Some(1)));
    assert_eq!(args.word_size, Some(11));
    assert_eq!(args.min_raw_gapped_score, Some(26));
    assert_eq!(args.xdrop_ungap, Some(40.0));
    assert_eq!(args.xdrop_gap_final, Some(90.0));
    assert_eq!(args.strand, "plus");
    assert_eq!(args.num_threads, 4);
    assert!(args.one_query_batch);
    assert!(args.hits_only);
    assert!(args.db_reader_one_thread);
    assert!(args.filtering_db.is_none());
}

#[test]
fn lower_case_next_to_the_interval() {
    // The interval is letters 6..10 ("ACGT"); lower case 1 and 5 letters before it, and 3
    // after it.
    let record = b"AaCGTaACGTAAcGT";
    let outside = MasksOutside::lower_case(record, 6, 4);
    assert_eq!(outside.before, 0b10001);
    assert_eq!(outside.nearest_before(), Some(1));
    assert_eq!(outside.after, Some(3));
    assert_eq!(
        MasksOutside::lower_case(b"ACGT", 0, 4),
        MasksOutside::default()
    );
    assert_eq!(MasksOutside::before_part(Some(2)).before, 0b10);
    assert_eq!(MasksOutside::before_part(Some(6)).before, 0);
}

#[test]
fn hits_join_masks_within_four_letters() {
    let masks = |pairs: &[(usize, usize)]| -> Vec<MaskedInterval> {
        pairs
            .iter()
            .map(|&(start, end)| MaskedInterval::new(start, end))
            .collect()
    };
    // A gap of three letters joins (4 + 5 > 8), one of four does not (12 + 5 > 17 is false).
    assert_eq!(
        join(&[(8, 12)], &masks(&[(0, 5), (17, 20)]), 30, None),
        masks(&[(0, 13), (17, 20)])
    );
    // The lower case 3 letters before the interval joins a hit that starts at its letter 1,
    // which masks letter 0, and the letters it joins before the interval are recorded; the
    // lower case 3 letters after its last letter (26) joins a hit that ends on letter 25,
    // which masks letter 26.
    let mut outside = MasksOutside {
        before: 0b100,
        after: Some(3),
    };
    assert_eq!(
        join(&[(1, 10), (20, 25)], &[], 27, Some(&mut outside)),
        masks(&[(0, 11), (20, 27)])
    );
    assert_eq!(outside.before, 0b111);
    // Without a hit near them the masks outside add nothing.
    let mut outside = MasksOutside {
        before: 0b10000,
        after: Some(5),
    };
    assert_eq!(
        join(&[(3, 10)], &[], 20, Some(&mut outside)),
        masks(&[(3, 11)])
    );
    assert_eq!(outside.before, 0b10000);
}

#[test]
fn batch_error_repeats_the_text_for_each_query() {
    let error = RepeatError::Posted(b"Database name is required.".to_vec());
    let native = batch_error(&error, 2).expect("NCBI's error");
    assert_eq!(native.exit, 3);
    // NCBI's bytes in LOSAT's words (AGENTS.md rule 8).
    assert_eq!(
        native.message,
        crate::wording::losat_bytes(
            b"BLAST engine error: Warning: Database name is required. Warning: Database name is required. \n"
        )
        .as_ref()
    );
    assert!(batch_error(&RepeatError::Search(anyhow::anyhow!("x")), 2).is_none());
}

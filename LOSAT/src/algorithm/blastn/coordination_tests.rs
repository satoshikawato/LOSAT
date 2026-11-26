use super::*;

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:867-949
fn blastn_args(words: &[&str]) -> BlastnArgs {
    let argv = [
        "losat", "blastn", "-query", "q", "-subject", "s", "-outfmt", "6",
    ];
    let cli: crate::cli::Cli =
        crate::cli::try_parse_from(argv.iter().chain(words)).expect("valid arguments");
    let crate::cli::Commands::Blastn(args) = cli.command else {
        panic!("blastn")
    };
    args
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:344-380
// NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:55-84
// DI-6b: a value above 0 replaces the cutoff of the e-value; 0 and below are unset
// (blast_parameters.c:876), in every task.
#[test]
fn min_raw_gapped_score_above_zero_is_kept() {
    for task in ["megablast", "dc-megablast", "blastn", "blastn-short"] {
        let config = configure_task(&blastn_args(&["-task", task]));
        assert_eq!(config.min_raw_gapped_score, None, "{task}");
        for (value, expected) in [
            ("15", Some(15)),
            ("2147483647", Some(i32::MAX)),
            ("0", None),
            ("-5", None),
        ] {
            let config = configure_task(&blastn_args(&[
                "-task",
                task,
                "-min_raw_gapped_score",
                value,
            ]));
            assert_eq!(config.min_raw_gapped_score, expected, "{task} {value}");
        }
    }
}

// DI-6: the given options replace the task's values; the defaults stay the handles'
// (api/blast_nucl_options.cpp:163-193, blast_args.cpp:276-286,487-489,518-522,681-685).
#[test]
fn xdrop_window_and_greedy_options_replace_the_task_values() {
    let config = configure_task(&blastn_args(&[]));
    assert_eq!(
        (
            config.x_drop_gapped,
            config.x_drop_final,
            config.x_drop_ungapped
        ),
        (25.0, 100.0, 0.0)
    );
    assert_eq!(
        (config.window_size, config.scan_range, config.use_dp),
        (0, 0, false)
    );
    let config = configure_task(&blastn_args(&[
        "-xdrop_gap",
        "12.5",
        "-xdrop_gap_final",
        "150",
        "-xdrop_ungap",
        "7",
        "-window_size",
        "40",
        "-off_diagonal_range",
        "2",
        "-no_greedy",
        "-gapopen",
        "5",
        "-gapextend",
        "2",
    ]));
    assert_eq!(
        (
            config.x_drop_gapped,
            config.x_drop_final,
            config.x_drop_ungapped
        ),
        (12.5, 150.0, 7.0)
    );
    assert_eq!(
        (config.window_size, config.scan_range, config.use_dp),
        (40, 2, true)
    );
    let config = configure_task(&blastn_args(&[
        "-task",
        "dc-megablast",
        "-window_size",
        "0",
    ]));
    assert_eq!((config.window_size, config.x_drop_ungapped), (0, 20.0));
    assert_eq!(
        determine_window_size(&blastn_args(&["-task", "dc-megablast"])),
        40
    );
    assert_eq!(determine_window_size(&blastn_args(&["-task", "blastn"])), 0);
}

#[test]
fn rmblastn_is_blastn_with_matrix_only_scoring() {
    // NCBI blast_options_handle.cpp:365-370: rmblastn is the blastn handle with reward 0 and
    // penalty 0; blast_setup.c:359-368: the engine then scores with BLAST_REWARD and
    // BLAST_PENALTY (1/-3) while the options keep 0/0.
    assert_eq!(engine_scoring_pair(0, 0), (true, 1, -3));
    assert_eq!(engine_scoring_pair(2, -3), (false, 2, -3));
    assert_eq!(engine_scoring_pair(0, -2), (false, 0, -2));
    let blastn = configure_task(&blastn_args(&["-task", "blastn"]));
    for words in [
        &["-task", "rmblastn"][..],
        &["-task", "blastn", "-reward", "0", "-penalty", "0"],
    ] {
        let args = blastn_args(words);
        let config = configure_task(&args);
        assert_eq!((config.reward, config.penalty), (0, 0), "{words:?}");
        assert!(config.matrix_only, "{words:?}");
        assert_eq!((config.score_reward, config.score_penalty), (1, -3));
        assert_eq!(
            (
                config.effective_word_size,
                config.gap_open,
                config.gap_extend
            ),
            (11, 5, 2)
        );
        assert!(config.use_dp);
        assert_eq!(
            (
                config.x_drop_ungapped,
                config.x_drop_gapped,
                config.x_drop_final
            ),
            (20.0, 30.0, 100.0)
        );
        assert_eq!(
            (config.min_diag_separation, config.window_size),
            (blastn.min_diag_separation, blastn.window_size)
        );
        assert_eq!(determine_evalue(&args), 10.0);
    }
    assert!(task_dust_by_default("rmblastn"));
    assert!(!task_uses_megablast_chunks("rmblastn"));
    // The other tasks keep the option pair as the engine's.
    assert!(!blastn.matrix_only);
    assert_eq!((blastn.score_reward, blastn.score_penalty), (2, -3));
}

#[test]
fn task_defaults_follow_ncbi_handles() {
    let config = configure_task(&blastn_args(&["-task", "dc-megablast"]));
    assert_eq!(
        (config.effective_word_size, config.reward, config.penalty),
        (11, 2, -3)
    );
    assert_eq!((config.gap_open, config.gap_extend), (5, 2));
    assert!(config.use_dp);
    assert_eq!((config.x_drop_gapped, config.x_drop_final), (30.0, 100.0));
    assert_eq!(config.x_drop_ungapped, 20.0);
    assert_eq!(config.min_diag_separation, 6);
    assert_eq!(config.window_size, 40);
    assert_eq!(
        (config.mb_template_type, config.mb_template_length),
        (DiscWordType::Coding, 18)
    );

    let args = blastn_args(&["-task", "blastn-short"]);
    let config = configure_task(&args);
    assert_eq!(
        (config.effective_word_size, config.reward, config.penalty),
        (7, 1, -3)
    );
    assert_eq!((config.gap_open, config.gap_extend), (5, 2));
    assert!(config.use_dp);
    assert_eq!(config.min_diag_separation, 50);
    assert_eq!((config.window_size, config.mb_template_length), (0, 0));
    assert_eq!(determine_evalue(&args), 1000.0);
    assert!(!task_dust_by_default("blastn-short"));
    assert!(!task_uses_megablast_chunks("blastn-short"));
    assert!(task_uses_megablast_chunks("dc-megablast"));

    // Given options replace the task's values; the template applies to any task.
    let args = blastn_args(&[
        "-task",
        "megablast",
        "-word_size",
        "12",
        "-evalue",
        "5",
        "-template_type",
        "coding_and_optimal",
        "-template_length",
        "21",
    ]);
    let config = configure_task(&args);
    assert_eq!(config.effective_word_size, 12);
    assert!(!config.use_dp);
    assert_eq!(config.window_size, 0);
    assert_eq!(
        (config.mb_template_type, config.mb_template_length),
        (DiscWordType::TwoTemplates, 21)
    );
    assert_eq!(determine_evalue(&args), 5.0);
    for task in ["megablast", "blastn", "dc-megablast"] {
        assert!(task_dust_by_default(task));
        assert_eq!(determine_evalue(&blastn_args(&["-task", task])), 10.0);
    }
}

#[test]
fn test_collect_lowercase_masks() {
    let masks = collect_lowercase_masks(b"AAaaBBbC");
    assert_eq!(
        masks,
        vec![MaskedInterval::new(2, 4), MaskedInterval::new(6, 7)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:121-127
#[test]
fn test_merge_mask_intervals() {
    let merged = merge_mask_intervals(vec![
        MaskedInterval::new(0, 2),
        MaskedInterval::new(2, 5),
        MaskedInterval::new(6, 7),
    ]);
    assert_eq!(
        merged,
        vec![MaskedInterval::new(0, 5), MaskedInterval::new(6, 7)]
    );
}

// No NCBI counterpart: the name of LOSAT's debug and trace lines (`subject_trace_id`).
#[test]
fn subject_trace_id_is_the_title_up_to_its_first_space() {
    assert_eq!(
        subject_trace_id(b"NR_181915.1 Pseudomonas abyssi"),
        "NR_181915.1"
    );
    assert_eq!(subject_trace_id(b"no_space"), "no_space");
    assert_eq!(subject_trace_id(b""), "unknown");
    assert_eq!(subject_trace_id(b" leading space"), "unknown");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:171-180,91-117
// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:76-84
// A query with ambiguity codes (letters other than C, G, T read as A, N at random) has
// other DUST intervals on its reverse complement; `-strand minus` dusts the reverse
// complement and maps the intervals back. The query is the first line of
// O2a `results/dust_asym/asym_s12_2.txt`; the intervals are the ones that file lists for the
// plus strand and for the mapped reverse complement (inclusive there, [start, end) here).
// The NCBI binary's Query rows agree with the mapped reverse complement
// (`results/ncbi_dustasym`).
#[test]
fn dust_of_the_minus_strand_is_that_of_the_reverse_complement() {
    const QUERY: &str = "AGACCAANCCACCMCCCCCCCKNCCCYCCCCCCCCCCCNCCCCCCAMCCCCKCCCCCCCGMCCCCCTACCKCCNMKCYCCCNCCCCNGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGWGAAGAGAAGAGAAGAGATAAACTTGGTGCGCAACCGGGKCGGGASGGGYCGGGACAGMAWGGG";
    let fasta = format!(">dq\n{QUERY}\n");
    let records = crate::blastinput::fasta_reader::read_all(
        &mut crate::blastinput::fasta_reader::FastaInputSource::from_bytes(
            fasta.as_bytes(),
            crate::blastinput::fasta_reader::ReaderConfig::query("BLASTN", false, false),
        ),
        &mut |_| Ok(()),
    )
    .unwrap();
    let masks = |strand: &str| {
        let args = blastn_args(&["-strand", strand]);
        crate::utils::threading::with_search_pool(1, "blastn", |pool| {
            Ok(apply_dust_masking(&args, &records[..], pool))
        })
        .unwrap()
        .remove(0)
        .into_iter()
        .map(|mask| (mask.start, mask.end))
        .collect::<Vec<_>>()
    };
    let plus = masks("plus");
    assert_eq!(plus, [(2, 154)]);
    assert_eq!(masks("both"), plus);
    assert_eq!(masks("minus"), [(14, 86), (87, 154)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181 (cited at
// `subject_metadata_from_records`)
// The total of a subject set is NCBI's `Int8` sum: two records of 3,000,000,000 letters each
// (each length below 2^32, as on wasm32) give 6,000,000,000, which a 32-bit `usize` sum would
// overflow.
#[test]
fn subject_set_total_is_ncbis_int8_sum() {
    let records = [
        FastaRecord::new("Subject_1", b"", b"A"),
        FastaRecord::new("Subject_2", b"", b"C"),
    ];
    let metadata = subject_metadata_from_records(&records, |_| 3_000_000_000);
    assert_eq!(metadata.db_len_total, 6_000_000_000);
    assert!(metadata.db_len_total > i64::from(u32::MAX));
    assert_eq!(metadata.db_num_seqs, 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:479-497
// An imported search strategy carries any `Int4` `WindowSize`: the window stays signed from the
// import to the word finder.
#[test]
fn window_size_below_zero_is_kept_signed() {
    for task in ["megablast", "dc-megablast", "blastn", "blastn-short"] {
        for window in [-1, -40, i32::MIN, i32::MAX] {
            let mut args = blastn_args(&["-task", task]);
            args.window_size = Some(window);
            assert_eq!(determine_window_size(&args), window, "{task} {window}");
            assert_eq!(configure_task(&args).window_size, window, "{task} {window}");
        }
    }
}

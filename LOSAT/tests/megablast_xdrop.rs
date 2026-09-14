//! Regression for the zero-initialized traditional megablast X-drop option.
//! The complete output is from the fresh NCBI 2.17.0 oracle recorded in
//! docs/evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/.

// NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54
// m_InitWordOpts.Reset((BlastInitialWordOptions*)calloc(1, sizeof(BlastInitialWordOptions)));
// NCBI reference: c++/src/algo/blast/api/blast_nucl_options.cpp:163-174
// SetInitialWordOptionsDefaults() { SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL); ... }
// SetMBInitialWordOptionsDefaults() { SetWindowSize(BLAST_WINDOW_SIZE_NUCL); }
// NCBI reference: c++/src/algo/blast/core/blast_parameters.c:380-383
// if (curr_cutoffs->x_dropoff_init == 0) curr_cutoffs->x_dropoff = new_cutoff;
// else curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
// NCBI reference: c++/src/algo/blast/core/blast_gapalign.c:4012-4031
// init_hsp->offsets.qs_offsets.q_off =
//     init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
// BLAST_GreedyGappedAlignment(...);
// init_hsp->offsets.qs_offsets.q_off = gap_align->greedy_query_seed_start;
// This fixture crosses an ungapped score valley at raw X-drop 13; using
// blastn's 20-bit option (raw 11) changes the traceback seed and edit script.
// Keep the full input: shortening it can change the automatic word cutoff.
#[test]
fn megablast_automatic_xdrop_preserves_oracle_scripts() {
    let input =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/NZ_CP006932.fasta");
    let directory = std::env::temp_dir().join(format!(
        "losat-megablast-xdrop-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("output.tsv");
    let cli: LOSAT::cli::Cli = LOSAT::cli::try_parse_from([
        "LOSAT",
        "blastn",
        "-query",
        input.to_str().unwrap(),
        "-subject",
        input.to_str().unwrap(),
        "-task",
        "megablast",
        "-num_threads",
        "1",
        "-outfmt",
        "6",
        "-out",
        output.to_str().unwrap(),
    ])
    .unwrap();
    let LOSAT::cli::Commands::Blastn(args) = cli.command else {
        panic!("expected blastn command");
    };
    LOSAT::algorithm::blastn::run(args).unwrap();
    let actual = std::fs::read(&output).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
    assert_eq!(
        actual,
        include_bytes!("fixtures/megablast_automatic_xdrop.outfmt6")
    );
}

//! Regression for the zero-initialized traditional megablast X-drop option.
//! The complete output is from the fresh NCBI 2.17.0 oracle recorded in
//! docs/evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/.

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:163-174
// `SetInitialWordOptionsDefaults` sets the ungapped X-drop to `BLAST_UNGAPPED_X_DROPOFF_NUCL`,
// while `SetMBInitialWordOptionsDefaults` sets `BLAST_WINDOW_SIZE_NUCL` and leaves X-drop zero.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:380-383
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4031
// This fixture crosses an ungapped score valley at raw X-drop 13; using
// blastn's 20-bit option (raw 11) changes the traceback seed and edit script.
// Keep the full input: shortening it can change the automatic word cutoff.
#[test]
fn megablast_automatic_xdrop_preserves_oracle_scripts() {
    let input =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/NZ_CP006932.fasta");
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-megablast-xdrop-{}-{}",
        crate::utils::test_dir::process_id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("output.tsv");
    let cli: crate::cli::Cli = crate::cli::try_parse_from([
        "losat",
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
    let crate::cli::Commands::Blastn(args) = cli.command else {
        panic!("expected blastn command");
    };
    super::run(
        args,
        crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    )
    .unwrap();
    let actual = std::fs::read(&output).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
    assert_eq!(
        actual,
        include_bytes!("../../../../tests/fixtures/megablast_automatic_xdrop.outfmt6")
    );
}

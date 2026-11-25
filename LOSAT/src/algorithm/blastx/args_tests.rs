use super::*;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:288-300
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
// Item SC9: every word size of NCBI's runs; one above 4 is the compressed table with its own
// threshold, one of 4 or less the lookup table with the task's (blastx-fast's 20 too), and a
// given -threshold stays. IDENTITY is refused with the compressed table, where NCBI crashes.
#[test]
fn every_word_size_takes_ncbi_s_table_and_threshold() {
    for (word, threshold, compressed) in [
        ("2", 12.0, false),
        ("3", 12.0, false),
        ("4", 12.0, false),
        ("5", 19.3, true),
        ("6", 21.0, true),
        ("7", 20.25, true),
    ] {
        let resolved = parse(&["-word_size", word]).resolve().unwrap();
        assert_eq!(resolved.word_size.to_string(), word);
        assert_eq!(resolved.threshold, threshold, "{word}");
        assert_eq!(resolved.compressed_lookup(), compressed, "{word}");
        let fast = parse(&["-task", "blastx-fast", "-word_size", word])
            .resolve()
            .unwrap();
        let fast_threshold = if compressed { threshold } else { 20.0 };
        assert_eq!(fast.threshold, fast_threshold, "blastx-fast {word}");
        assert_eq!(fast.compressed_lookup(), compressed, "blastx-fast {word}");
        let given = parse(&["-word_size", word, "-threshold", "17"])
            .resolve()
            .unwrap();
        assert_eq!(given.threshold, 17.0, "{word}");
        let identity = parse(&[
            "-word_size",
            word,
            "-matrix",
            "IDENTITY",
            "-ungapped",
            "-comp_based_stats",
            "0",
        ])
        .resolve();
        if compressed {
            assert_eq!(
                identity.err().unwrap().to_string(),
                format!("-matrix IDENTITY with the compressed lookup table (-word_size {word}), with which NCBI BLAST+'s blastx crashes, is not supported by LOSAT's BLASTX")
            );
        } else {
            assert!(identity.is_ok(), "{word}");
        }
    }
    // The task's default word sizes: 3 is the lookup table, blastx-fast's 5 the compressed one.
    assert!(!parse(&[]).resolve().unwrap().compressed_lookup());
    assert!(parse(&["-task", "blastx-fast"])
        .resolve()
        .unwrap()
        .compressed_lookup());
}
fn parse(options: &[&str]) -> BlastxArgs {
    let mut argv = vec!["losat", "blastx", "-query", "q", "-subject", "s"];
    argv.extend_from_slice(options);
    let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
    match cli.command {
        crate::cli::Commands::Blastx(a) => a,
        _ => unreachable!(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:276-286 and
// core/blast_options.c:608-618 (cited at `resolve_with_format`): the X-drops reach the
// options (7, 15 and 25 bits by default, blastx-fast too), an ungapped X-drop of 0 or below is
// NCBI's options error, and a gapped X-drop that NCBI's kernels cannot size is refused
// (`core::x_drop`).
#[test]
fn the_x_drops_reach_the_options_or_ncbis_errors() {
    for task in [&[][..], &["-task", "blastx-fast"]] {
        let options = parse(task).resolve().unwrap();
        assert_eq!(
            (
                options.x_dropoff,
                options.gap_x_dropoff,
                options.gap_x_dropoff_final
            ),
            (7.0, 15.0, 25.0)
        );
    }
    let options = parse(&[
        "-xdrop_ungap",
        "3",
        "-xdrop_gap",
        "-1",
        "-xdrop_gap_final",
        "1e300",
    ])
    .resolve()
    .unwrap();
    assert_eq!(
        (
            options.x_dropoff,
            options.gap_x_dropoff,
            options.gap_x_dropoff_final
        ),
        (3.0, -1.0, 1e300)
    );
    for value in ["0", "-1", "-0.0", "-inf"] {
        let error = parse(&["-xdrop_ungap", value]).resolve().unwrap_err();
        assert!(
            format!("{error:#}").contains("x_dropoff must be greater than zero"),
            "{value}: {error:#}"
        );
    }
    // 4.2e8 bits are 1.09e9 raw at BLOSUM62 11/1's Lambda, past the kernels' 2^30 - 103: the
    // preliminary X-drop at every mode, the final one only in mode 0 (the redo of modes 1 to 3
    // converts its own at the scaled costs); an ungapped search has no gapped kernels.
    for extra in [
        &["-xdrop_gap", "4.2e8"][..],
        &["-xdrop_gap", "4.2e8", "-comp_based_stats", "2"],
        &["-xdrop_gap_final", "4.2e8", "-comp_based_stats", "0"],
    ] {
        let error = parse(extra).resolve().unwrap_err();
        assert!(
            format!("{error:#}").starts_with(&format!("an {} of 420000000.0 bits", extra[0])),
            "{extra:?}: {error:#}"
        );
    }
    for cbs in ["1", "2", "3"] {
        assert!(
            parse(&["-xdrop_gap_final", "4.2e8", "-comp_based_stats", cbs])
                .resolve()
                .is_ok()
        );
    }
    assert!(
        parse(&["-ungapped", "-comp_based_stats", "0", "-xdrop_gap", "4.2e8"])
            .resolve()
            .is_ok()
    );
}
#[test]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1300-1306
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
fn cli_acceptance_precedes_core_validation() {
    assert_eq!(parse(&["-threshold", "0"]).threshold, Some(0.0));
    assert!(parse(&["-threshold", "0"]).resolve().is_err());
    assert_eq!(parse(&["-evalue", "0"]).evalue, 0.0);
    assert!(parse(&["-evalue", "0"]).resolve().is_err());
}
#[test]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
fn omission_and_explicit_target_count_remain_distinct() {
    let omitted = parse(&[]).resolve().unwrap();
    let explicit = parse(&["-max_target_seqs", "500"]).resolve().unwrap();
    assert_eq!(
        (
            omitted.hitlist_size,
            omitted.num_descriptions,
            omitted.num_alignments
        ),
        (500, 500, 250)
    );
    assert_eq!(
        (
            explicit.hitlist_size,
            explicit.num_descriptions,
            explicit.num_alignments
        ),
        (500, 500, 500)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973 (`blastinput::report_counts`).
#[test]
fn num_descriptions_and_num_alignments_set_the_counts() {
    let counts = |extra: &[&str]| {
        let resolved = parse(extra).resolve().unwrap();
        (
            resolved.hitlist_size,
            resolved.num_descriptions,
            resolved.num_alignments,
        )
    };
    assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
    assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
    assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
    assert_eq!(
        counts(&["-outfmt", "6", "-num_descriptions", "2"]),
        (500, 500, 500)
    );
    let error = parse(&["-num_descriptions", "0", "-num_alignments", "0"])
        .resolve()
        .unwrap_err();
    assert_eq!(error.to_string(), "No hits are being saved");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2801-2855 and c++/src/objtools/align_format/tabular.cpp:70-99: BLASTX reads -outfmt and its
// fields as NCBI (a specification counts only for the tabular formats; unknown words are
// dropped, `-name` removes a field). Every tabular report goes to the shared writer
// (`tabular_report`): a list other than the `std` fields with a tab, and outfmt 10 and 20,
// with its spec (`tabular`), the default outfmt 6 and 7 with the `std` spec (`tabular` None).
#[test]
fn outfmt_and_fields_are_read_as_ncbi() {
    use crate::report::tabular::TabularKind;
    // The fields of the report the tabular writer writes (`tabular_report`), none for the
    // formats that are not tabular.
    let fields = |outfmt: &str| {
        let resolved = parse(&["-outfmt", outfmt]).resolve().unwrap();
        let names: Vec<String> = resolved.tabular_report().map_or(Vec::new(), |spec| {
            spec.fields.iter().map(|f| f.name().to_string()).collect()
        });
        (resolved.outfmt, names)
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40 (`kDfltArgTabularOutputFmt`, the `std` fields).
    let std: Vec<String> =
        "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            .split(' ')
            .map(String::from)
            .collect();
    assert_eq!(fields("0 std"), (0, Vec::new()));
    assert_eq!(fields("0 foo"), (0, Vec::new()));
    assert_eq!(fields(" 6 delim="), (6, std.clone()));
    assert_eq!(fields("7"), (7, std.clone()));
    assert_eq!(fields("6 foo"), (6, std.clone()));
    // The default reports' kind and delimiter (outfmt 6: rows; 7: with its comment lines).
    let report = |outfmt: &str| {
        parse(&["-outfmt", outfmt])
            .resolve()
            .unwrap()
            .tabular_report()
            .map(|spec| (spec.kind, spec.delimiter.clone()))
    };
    assert_eq!(report("6"), Some((TabularKind::Tabular, b"\t".to_vec())));
    assert_eq!(
        report("7 std"),
        Some((TabularKind::Commented, b"\t".to_vec()))
    );
    assert_eq!(report("5"), None);
    // Words between single spaces; an empty word names no field.
    let mut listed = vec!["qseqid".to_string(), "sseqid".to_string()];
    listed.extend(std.iter().cloned());
    assert_eq!(fields("7 qseqid  sseqid std"), (7, listed));
    // The XInclude formats to standard output fail later (`native::run`).
    assert_eq!(fields("13").0, 13);
    let tabular = |outfmt: &str| {
        parse(&["-outfmt", outfmt])
            .resolve()
            .unwrap()
            .tabular
            .map(|spec| (spec.kind, spec.fields.len(), spec.delimiter.clone()))
    };
    for outfmt in ["6", "7 std", "6 foo", "6 -qseqid", "0 qseqid", "5"] {
        assert_eq!(tabular(outfmt), None, "{outfmt:?}");
    }
    assert_eq!(
        tabular("6 std -evalue staxids"),
        Some((TabularKind::Tabular, 12, b"\t".to_vec()))
    );
    assert_eq!(
        tabular("7 delim=| std"),
        Some((TabularKind::Commented, 12, b"|".to_vec()))
    );
    assert_eq!(tabular("10"), Some((TabularKind::Csv, 12, b",".to_vec())));
    assert_eq!(
        tabular("20 qseqid"),
        Some((TabularKind::CsvWithHeader, 1, b",".to_vec()))
    );
}

use super::*;
use crate::algorithm::tblastn::stage_d_pipeline::{
    run_local_for_report, LocalStageDProfile, LocalStageDScoring,
};
use crate::api::local_blast::{FormatOutput, OutputSink};
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;
use crate::utils::seg::SegParams;

// The reports that read no letters of the reported database sequences: a value read wrongly
// here reaches `SubjectLetters::Lengths::nt`, which panics.
#[test]
fn only_the_reports_that_show_letters_read_the_subjects() {
    let reads_by = |hits: bool, recounted: bool, outfmts: &[&str]| {
        formats_read_subject_letters(hits, recounted, outfmts.iter().copied())
    };
    let reads = |outfmt: &str| reads_by(false, false, &[outfmt]);
    for outfmt in [
        "0",
        "5",
        "13",
        "14",
        "15",
        "16",
        "6 sseq",
        "6 qseq",
        "6 btop",
        "6 positive",
        "6 ppos",
        "7 std sseq",
        "10 sseq",
        "20 btop",
        // Alone in the counts' block, `nident` reads the rows of some HSPs.
        "6 nident",
        "not a format",
    ] {
        assert!(reads(outfmt), "{outfmt}");
    }
    for outfmt in [
        "6",
        "7",
        "10",
        "20",
        "6 std",
        "6 qseqid sseqid slen evalue",
        "6 nident pident",
        "8",
        "9",
        "11",
        "12",
    ] {
        assert!(!reads(outfmt), "{outfmt}");
    }
    assert!(reads_by(true, false, &["6"]));
    assert!(!reads_by(false, false, &[]));
    assert!(reads_by(false, false, &["6", "0"]));
    // A PSSM query whose reports count its HSPs again from their rows (lower-case iupacaa,
    // the gap residue: M4a-1c LG probe `ts_low_iupac.f1`, `ts_std0.g4`): every tabular
    // report reads the rows; the Seq-aligns and the archive keep the search's counts.
    for outfmt in ["6", "7", "10", "6 std", "6 qseqid sseqid slen evalue"] {
        assert!(reads_by(false, true, &[outfmt]), "{outfmt}");
    }
    for outfmt in ["8", "9", "11", "12", "18"] {
        assert!(!reads_by(false, true, &[outfmt]), "{outfmt}");
    }
}

// Outfmt 5 and the XML2 and JSON documents (16 and 15) take the structured-XML path of
// `render`; no other format does.
#[test]
fn the_xml_family_is_5_and_13_to_16() {
    use crate::report::blast_xml2::Xml2Kind;
    for (outfmt, kind, structured) in [
        ("5", None, true),
        ("13", Some(Xml2Kind::Files13), true),
        ("14", Some(Xml2Kind::Files14), true),
        ("15", Some(Xml2Kind::Single15), true),
        ("16", Some(Xml2Kind::Single16), true),
        (" 16 std", Some(Xml2Kind::Single16), true),
        ("0", None, false),
        ("6", None, false),
        ("12", None, false),
    ] {
        assert_eq!(xml2_kind(outfmt), kind, "{outfmt}");
        assert_eq!(is_structured_xml(outfmt), structured, "{outfmt}");
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288-301:
// Compare complete saved formatter bytes after the Stage D local search.
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467:
// All three formats are printed from one result.
#[test]
fn initial_tabular_oracle_bytes() {
    let query_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/query.faa"
    ));
    let subject_path = Path::new("docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna");
    let full_subject_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna"
    ));
    // The records of NCBI's reader, as TBLASTN reads the query and the subjects.
    let read = |path: &Path, prefix: &str, protein: bool| {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        let config = match prefix {
            "Query_" => ReaderConfig::query("TBLASTN", protein, false),
            _ => ReaderConfig::subject("TBLASTN", protein, false),
        };
        let file = std::fs::File::open(path).unwrap();
        read_all(&mut FastaInputSource::from_file(file, config), &mut |_| {
            Ok(())
        })
        .unwrap()
    };
    let queries = read(query_path, "Query_", true);
    let subjects = read(full_subject_path, "Subject_", false);
    let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
    let subject_seqs: Vec<&[u8]> = subjects.iter().map(|record| record.seq()).collect();
    let seg = SegParams::default();
    for mode in [0, 2] {
        let (mut results, parameters, ungapped_karlin, query_validity) = run_local_for_report(
            &query_seqs,
            &subject_seqs,
            LocalStageDProfile {
                max_intron_length: 0,
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
            },
            if mode == 2 {
                BlastCompoAdjustMode::CompositionMatrixAdjust
            } else {
                BlastCompoAdjustMode::NoCompositionBasedStats
            },
            true,
            LocalStageDScoring::default(),
        )
        .unwrap();
        const FORMATS: [&str; 3] = ["0", "6", "7"];
        let mut outputs_bytes = vec![Vec::new(); FORMATS.len()];
        let mut diagnostics = Vec::new();
        let mut outputs = ReportOutputs {
            other_outfmts: &[],
            html_posts: None,
            formats: FORMATS
                .iter()
                .zip(outputs_bytes.iter_mut())
                .map(|(outfmt, bytes)| FormatOutput {
                    outfmt,
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(bytes),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        render(
            &mut outputs,
            &queries,
            &subjects,
            crate::algorithm::tblastn::x_subject_store::SubjectLetters::Records(&subject_seqs),
            &ShownSubjects::Fasta {
                path: subject_path,
                parsed: &[],
            },
            SubjectSetNumbers::default(),
            &mut results,
            &parameters,
            &ungapped_karlin,
            &query_validity,
            &vec![false; queries.len()],
            LocalStageDScoring::default(),
            "BLOSUM62",
            Default::default(),
            Default::default(),
            1,
            Some(&seg),
            Some((12, 2.2, 2.5)),
            false,
            &[],
            &[],
            &ReportRanges {
                queries: &Default::default(),
                subjects: &Default::default(),
                epilog: true,
            },
            super::super::blast_xml::XmlValues {
                expect: 10.0,
                soft_masking: false,
                gapped: true,
                cbs: 2,
                db_gencode: 1,
                identity_warned: 0,
                filtered_empty: false,
            },
            Default::default(),
            None,
            None,
            &Default::default(),
            None,
            false,
            &[],
        )
        .unwrap();
        drop(outputs);
        for (fmt, actual) in FORMATS.iter().zip(&outputs_bytes) {
            let expected = std::fs::read(format!(
                "{}/../docs/evidence/tlosan_stage_e/initial_20260925/multi_hsp.mode{mode}.outfmt{fmt}",
                env!("CARGO_MANIFEST_DIR"),
            ))
            .unwrap();
            // The frozen files are NCBI's bytes; LOSAT writes them in its own words.
            let expected = crate::wording::losat_output_bytes(&expected).into_owned();
            if *actual != expected {
                // No NCBI counterpart: test plumbing; the platform's temporary directory
                // replaces the Unix-only "/tmp" that made a mismatch fail on Windows with
                // "path not found" instead of reporting the byte difference.
                let path = crate::utils::test_dir::temp_dir()
                    .join(format!("tblastn_stagee_mode{mode}_fmt{fmt}.actual"))
                    .display()
                    .to_string();
                std::fs::write(&path, actual).unwrap();
                let offset = actual
                    .iter()
                    .zip(&expected)
                    .position(|(a, b)| a != b)
                    .unwrap_or(actual.len().min(expected.len()));
                panic!("mode={mode} outfmt={fmt} first byte difference at {offset}, actual saved to {path}");
            }
        }
    }
}

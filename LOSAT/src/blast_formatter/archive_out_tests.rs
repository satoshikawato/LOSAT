//! Unit tests for `blast_formatter -outfmt 11` (`archive_out`).

use std::ffi::OsString;

use super::*;
use crate::report::archive::read::options::{BuiltOptions, EProgram, SeqRange, Setter};

fn fixture_archive(name: &str) -> Vec<u8> {
    std::fs::read(format!(
        "{}/tests/blast_formatter/in/{name}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn run_bytes(input: &[u8], words: &[&str]) -> (i32, Vec<u8>, String) {
    let mut line: Vec<OsString> = vec!["-archive".into(), "-".into()];
    line.extend(words.iter().map(OsString::from));
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status = super::super::run_with(&line, &mut &input[..], &mut out, &mut err);
    (status, out, String::from_utf8(err).unwrap())
}

fn empty_built(program: EProgram, setters: Vec<Setter>) -> BuiltOptions {
    BuiltOptions {
        program,
        setters,
        entrez_query: None,
        first_db_seq: None,
        final_db_seq: None,
        gi_list: None,
        negative_gi_list: None,
        taxid_list: None,
        negative_taxid_list: None,
        db_filtering_algorithm_id: None,
        db_filtering_algorithm_key: None,
        subject_masking_type: None,
        query_masks: None,
        query_range: SeqRange::EMPTY,
    }
}

#[test]
fn outfmt_11_reemits_subject_archive_and_is_idempotent() {
    let input = fixture_archive("bn.asn");
    let (status, out1, err1) = run_bytes(&input, &["-outfmt", "11"]);
    assert_eq!((status, err1.as_str()), (0, ""));
    let text1 = String::from_utf8(out1.clone()).unwrap();
    assert!(text1.starts_with("Blast4-archive ::= {\n"), "{text1}");
    assert!(
        text1.contains(&format!("ident \"{}\"", crate::wording::LOSAT_VERSION)),
        "{text1}"
    );
    assert!(text1.contains("program \"blastn\""), "{text1}");
    assert!(text1.contains("service \"megablast\""), "{text1}");
    // Feeding the emitted archive back into `blast_formatter -outfmt 11` yields identical bytes.
    let (status2, out2, err2) = run_bytes(&out1, &["-outfmt", "11"]);
    assert_eq!((status2, err2.as_str()), (0, ""));
    assert_eq!(out1, out2);
}

#[test]
fn outfmt_11_updates_hitlist_size_without_pruning_alignments() {
    let input = fixture_archive("bn.asn");
    let (status, out, err) = run_bytes(&input, &["-outfmt", "11", "-max_target_seqs", "1"]);
    assert_eq!(status, 0);
    assert!(err.contains("Examining 5 or more matches"), "{err}");
    let text = String::from_utf8(out).unwrap();
    assert!(
        text.contains("name \"HitlistSize\",\n          value integer 1"),
        "{text}"
    );
    assert!(text.contains("alignments {"), "{text}");
    assert!(text.contains("local str \"Subject_1\""), "{text}");

    let (status0, out0, err0) = run_bytes(&input, &["-outfmt", "11", "-num_alignments", "0"]);
    assert_eq!(status0, 0);
    assert!(err0.contains("Examining 5 or more matches"), "{err0}");
    let text0 = String::from_utf8(out0).unwrap();
    assert!(
        text0.contains("name \"HitlistSize\",\n          value integer 0"),
        "{text0}"
    );
    // Alignments are NOT pruned even with `-num_alignments 0`.
    assert!(text0.contains("local str \"Subject_1\""), "{text0}");
}

#[test]
fn outfmt_11_omits_messages_and_skips_error_archive() {
    let err_archive = fixture_archive("err.asn");
    let (status, out, err) = run_bytes(&err_archive, &["-outfmt", "11"]);
    assert_eq!((status, out.as_slice()), (0, &b""[..]));
    assert!(err.contains("gap existence and extension"), "{err}");

    let bn = String::from_utf8(fixture_archive("bn.asn")).unwrap();
    let with_msg = bn.replacen(
        "\n}\n",
        ",\n  messages {\n    {\n      code 1,\n      message \"Warning: sample note\"\n    }\n  }\n}\n",
        1,
    );
    let (status, out, err) = run_bytes(with_msg.as_bytes(), &["-outfmt", "11"]);
    assert_eq!(status, 0);
    assert_eq!(err, "Warning: sample note\n");
    let text = String::from_utf8(out).unwrap();
    assert!(!text.contains("messages {"), "{text}");
}

#[test]
fn filter_string_expansion_matches_ncbi_rules() {
    // Nucleotide `L;m;` -> MaskAtHash TRUE, DustFiltering TRUE + defaults, RepeatFiltering FALSE.
    let opts = AlgorithmOptions::from_built(&empty_built(
        EProgram::Blastn,
        vec![Setter::FilterString("L;m;".to_string())],
    ));
    let names: Vec<(&str, &ParamValue)> =
        opts.params().iter().map(|p| (p.name(), &p.value)).collect();
    assert_eq!(
        names,
        vec![
            ("EvalueThreshold", &ParamValue::Evalue(10.0)),
            ("DustFiltering", &ParamValue::Boolean(true)),
            ("RepeatFiltering", &ParamValue::Boolean(false)),
            ("MaskAtHash", &ParamValue::Boolean(true)),
            ("DustFilteringLevel", &ParamValue::Integer(20)),
            ("DustFilteringWindow", &ParamValue::Integer(64)),
            ("DustFilteringLinker", &ParamValue::Integer(1)),
        ]
    );

    // Custom DUST + Repeat DB on blastn.
    let opts = AlgorithmOptions::from_built(&empty_built(
        EProgram::Blastn,
        vec![Setter::FilterString(
            "D 15 32 2; R -d custom.lib".to_string(),
        )],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFilteringLevel" && p.value == ParamValue::Integer(15)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "RepeatFiltering" && p.value == ParamValue::Boolean(true)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "RepeatFilteringDB"
            && p.value == ParamValue::Str(b"custom.lib".to_vec())));

    // Protein `S` -> SegFiltering TRUE + defaults, no Dust/Repeat.
    let opts = AlgorithmOptions::from_built(&empty_built(
        EProgram::Blastp,
        vec![Setter::FilterString("S".to_string())],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "SegFiltering" && p.value == ParamValue::Boolean(true)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "SegFilteringWindow" && p.value == ParamValue::Integer(12)));
    assert!(!opts.params().iter().any(|p| p.name() == "DustFiltering"));

    // Invalid filter string falls back to cleared options.
    let opts = AlgorithmOptions::from_built(&empty_built(
        EProgram::Blastn,
        vec![Setter::FilterString("D 15 32".to_string())],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFiltering" && p.value == ParamValue::Boolean(false)));
    assert!(!opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFilteringLevel"));
}

#[test]
fn vecscreen_keeps_defaults_mode_stuck_except_evalue() {
    let mut opts = AlgorithmOptions::from_built(&empty_built(
        EProgram::VecScreen,
        vec![
            Setter::WordSize(11),
            Setter::EvalueThreshold(500.0),
            Setter::HitlistSize(10),
        ],
    ));
    opts.set_hitlist_size(25);
    assert_eq!(opts.program_name, "blastn");
    assert_eq!(opts.service_name, "vecscreen");
    assert_eq!(
        opts.params(),
        &[Param::new(
            fields::EVALUE_THRESHOLD,
            ParamValue::Evalue(500.0)
        )]
    );
}

#[test]
fn factory_defaults_for_disc_psi_psitblastn_and_delta() {
    let disc = AlgorithmOptions::from_built(&empty_built(EProgram::DiscMegablast, vec![]));
    let disc_names: Vec<&str> = disc.params().iter().map(Param::name).collect();
    assert_eq!(
        disc_names,
        vec![
            "EvalueThreshold",
            "MBTemplateType",
            "MBTemplateLength",
            "WordSize",
            "WindowSize",
        ]
    );

    let psi_tbn = AlgorithmOptions::from_built(&empty_built(EProgram::PsiTblastn, vec![]));
    let psi_tbn_names: Vec<&str> = psi_tbn.params().iter().map(Param::name).collect();
    assert_eq!(
        psi_tbn_names,
        vec![
            "EvalueThreshold",
            "InclusionThreshold",
            "PseudoCountWeight",
            "SmithWatermanMode",
            "DbGeneticCode",
            "WordThreshold",
        ]
    );

    let delta = AlgorithmOptions::from_built(&empty_built(EProgram::DeltaBlast, vec![]));
    assert_eq!(delta.service_name, "delta_blast");
    assert!(delta
        .params()
        .iter()
        .any(|p| p.name() == "InclusionThreshold" && p.value == ParamValue::Real(0.05)));
}

#[test]
fn taxid_and_negative_taxid_program_options_match_ncbi() {
    let mut built = empty_built(EProgram::Blastn, vec![]);
    built.taxid_list = Some(vec![9606, 562, 9606]);
    assert_eq!(
        build_program_options(&built, true),
        vec![Param::new(
            fields::TAXID_LIST,
            ParamValue::IntegerList(vec![562, 9606])
        )]
    );
    assert!(build_program_options(&built, false).is_empty());

    let mut neg = empty_built(EProgram::Blastn, vec![]);
    neg.negative_taxid_list = Some(vec![9606]);
    assert_eq!(
        build_program_options(&neg, true),
        vec![Param::new(
            fields::NEGATIVE_TAXID_LIST,
            ParamValue::IntegerList(Vec::new())
        )]
    );
}

#[test]
fn seq_loc_query_is_refused_in_outfmt_11() {
    let archive = br#"Blast4-archive ::= {
  request {
    ident "test",
    body queue-search {
      program "blastn",
      service "megablast",
      queries seq-loc-list {
        whole local str "Query_1"
      },
      subject sequences {
        {
          id {
            local str "Subject_1"
          },
          inst {
            repr raw,
            mol dna,
            length 4,
            seq-data ncbi2na '00'H
          }
        }
      },
      algorithm-options {
        {
          name "EvalueThreshold",
          value cutoff e-value { 1, 10, 1 }
        }
      }
    }
  },
  results {
    alignments {
    }
  }
}
"#;
    let (status, out, err) = run_bytes(archive, &["-outfmt", "11"]);
    assert_eq!(
        (status, out.as_slice(), err),
        (
            1,
            &b""[..],
            format!("Error: {}\n", Refusal::SeqLocQuery.reason())
        )
    );
}

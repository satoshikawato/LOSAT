//! Tests of the ASN.1 text PSSM reader (`pssm_text.rs`) and of the writer it pairs with
//! (`report::serial::pssm`).

use super::*;
use crate::report::serial::pssm::{asn_text, PssmWithParameters};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// A gzip-compressed fixture, inflated.
fn gunzip(name: &str) -> Vec<u8> {
    use std::io::Read as _;
    let file = std::fs::File::open(fixture(name)).unwrap();
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(file)
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Asserts that `written` is `expected`, naming the first line where they differ.
fn assert_same_text(written: &[u8], expected: &[u8], what: &str) {
    if written == expected {
        return;
    }
    let mut written_lines = written.split(|&byte| byte == b'\n');
    for (number, line) in expected.split(|&byte| byte == b'\n').enumerate() {
        let other = written_lines.next().unwrap_or_default();
        assert_eq!(
            String::from_utf8_lossy(other),
            String::from_utf8_lossy(line),
            "{what}: line {}",
            number + 1
        );
    }
    panic!(
        "{what}: {} bytes written, {} expected",
        written.len(),
        expected.len()
    );
}

fn real_of(text: &str) -> f64 {
    let document = format!(
        "PssmWithParameters ::= {{ pssm {{ numRows 1, numColumns 1, finalData {{ scores {{ 1 }}, lambda {text}, kappa {{ 0, 10, 0 }}, h {{ 0, 10, 0 }} }} }} }}"
    );
    read_pssm_with_parameters(document.as_bytes())
        .unwrap()
        .pssm
        .final_data
        .unwrap()
        .lambda
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (quoted at `Reader::real`):
// `mantissa * pow(10.0, exp)` with glibc's `pow`, 0 as read, the clamps to DBL_MIN and DBL_MAX
// of the magnitude, the sign put back after them, and the word forms. The expected bits are the
// products computed with glibc 2.39's `pow` (Python's `math.pow`).
#[test]
fn reals_are_read_with_ncbis_arithmetic() {
    for (text, bits) in [
        // kappa of NCBI 2.17.0 `rnd.p`; the frequency ratio 1 of `pssm_freq_ratios.asn`.
        ("{ 459597639216676, 10, -16 }", 0x3fa7_8809_c5f3_6d45u64),
        ("{ 576481676635514, 10, -16 }", 0x3fad_840f_8594_1c3f),
        ("{ 267, 10, -3 }", 0x3fd1_1687_2b02_0c4a),
        ("{ 765945647632464, 10, -15 }", 0x3fe8_82a0_7263_0a05),
        // glibc's pow(10, 23) is the double above 10^23's nearest.
        ("{ 3, 10, 23 }", 0x44cf_c384_2bd1_f072),
    ] {
        assert_eq!(real_of(text).to_bits(), bits, "{text}");
    }
    assert_eq!(real_of("{ -1, 10, 0 }"), -1.0);
    assert_eq!(real_of("{ 0, 10, 0 }").to_bits(), 0);
    assert_eq!(real_of("{ -0, 10, 0 }").to_bits(), (-0.0f64).to_bits());
    assert_eq!(real_of("{ 0, 10, 400 }").to_bits(), 0);
    assert_eq!(real_of("{ 1, 10, -400 }"), f64::MIN_POSITIVE);
    assert_eq!(real_of("{ -5, 10, -400 }"), -f64::MIN_POSITIVE);
    assert_eq!(real_of("{ 1, 10, 400 }"), f64::MAX);
    assert_eq!(real_of("{ -1, 10, 400 }"), -f64::MAX);
    assert_eq!(real_of("{ 3, 2, 4 }"), 48.0);
    assert_eq!(real_of("{ +7, 10, 1 }"), 70.0);
    assert_eq!(real_of("PLUS-INFINITY"), f64::INFINITY);
    assert_eq!(real_of("minus-infinity"), f64::NEG_INFINITY);
    // `HUGE_VAL/HUGE_VAL` of NCBI's x86-64 builds: the sign bit set (`asn_text_reader::NCBI_NAN`;
    // PS16d review oracle O2 `a_nan` prints `-nan`), also for the word form `nan`.
    assert!(real_of("NOT-A-NUMBER").is_nan() && real_of("NOT-A-NUMBER").is_sign_negative());
    assert!(real_of("nan").is_nan() && real_of("nan").is_sign_negative());
    assert_eq!(real_of("0.5"), 0.5);
    assert_eq!(real_of("-2.5e3"), -2500.0);
    assert_eq!(real_of("00.25"), 0.25);
    // -3 + 2147483647 is in range: far above DBL_MAX, which the word form clamps
    // (`fDecimalPosixFinite`).
    assert_eq!(real_of("0.005e2147483647"), f64::MAX);
    assert_eq!(real_of("-0.005e2147483647"), -f64::MAX);
    // -3 - 2147483647 is not (NCBI's `int` overflows): the reader's error, no panic.
    for bad in [
        "{ 1, 3, 0 }",
        "{ -, 10, 0 }",
        "{ 1, 10 }",
        "{ 1234567890123456789, 10, 0 }",
        "0.005e-2147483647",
        "-0.005e-2147483647",
    ] {
        let document = format!(
            "PssmWithParameters ::= {{ pssm {{ numRows 1, numColumns 1, finalData {{ scores {{ 1 }}, lambda {bad}, kappa {{ 0, 10, 0 }}, h {{ 0, 10, 0 }} }} }} }}"
        );
        assert!(
            read_pssm_with_parameters(document.as_bytes()).is_err(),
            "{bad}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:845-846 (quoted at `Reader::real`):
// `pow((double)10.0, exp)`, which the table of `libm_pow10` stands for.
#[test]
#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn libm_pow10_is_glibcs_pow() {
    use crate::report::serial::asn_text_reader::libm_pow10;

    for exp in -400..=400 {
        let pow = 10f64.powf(f64::from(exp));
        assert_eq!(libm_pow10(exp).to_bits(), pow.to_bits(), "10^{exp}");
    }
    assert_eq!(libm_pow10(i32::MIN), 0.0);
    assert_eq!(libm_pow10(i32::MAX), f64::INFINITY);
}

/// A PSSM of two rows and one column, with the members psiblast writes.
fn small_pssm(title: &[u8]) -> PssmWithParameters {
    use crate::report::serial::pssm::*;
    PssmWithParameters {
        pssm: Pssm {
            is_protein: Some(true),
            identifier: None,
            num_rows: 2,
            num_columns: 1,
            row_labels: None,
            by_row: Some(false),
            query: Some(PssmQuery {
                ids: vec![crate::seqid::ParsedSeqId::Local(ObjectId::Str(
                    b"Query_7".to_vec(),
                ))],
                titles: Some(vec![title.to_vec()]),
                inst: SeqInst {
                    repr: ("raw", 2),
                    mol: ("aa", 3),
                    length: Some(1),
                    seq_data: Some(SeqData::Ncbistdaa(vec![12])),
                },
            }),
            intermediate_data: Some(PssmIntermediateData {
                weighted_res_freqs_per_pos: Some(vec![0.0, 0.25]),
                freq_ratios: vec![0.0, 1.5],
                information_content: Some(vec![-1.0]),
                gapless_column_weights: Some(vec![0.0]),
                sigma: Some(vec![2.0]),
                interval_sizes: Some(vec![3]),
                num_matching_seqs: Some(vec![2]),
                ..PssmIntermediateData::default()
            }),
            final_data: Some(PssmFinalData {
                scores: vec![-32768, 5],
                lambda: 0.267,
                kappa: 0.041,
                h: 0.14,
                scaling_factor: None,
                lambda_ungapped: Some(-1.0),
                kappa_ungapped: Some(-1.0),
                h_ungapped: Some(-1.0),
                word_score_threshold: None,
            }),
        },
        params: Some(PssmParameters {
            pseudocount: Some(0),
            rpsdbparams: Some(FormatRpsDbParameters {
                matrix_name: b"BLOSUM62".to_vec(),
                gap_open: Some(11),
                gap_extend: Some(1),
            }),
            bit_score_thresh: None,
            bit_score_reporting_thresh: None,
        }),
    }
}

const SMALL_PSSM: &str = r#"PssmWithParameters ::= {
  pssm {
    isProtein TRUE,
    numRows 2,
    numColumns 1,
    byRow FALSE,
    query seq {
      id {
        local str "Query_7"
      },
      descr {
        title "a ""quoted"" title"
      },
      inst {
        repr raw,
        mol aa,
        length 1,
        seq-data ncbistdaa '0C'H
      }
    },
    intermediateData {
      weightedResFreqsPerPos {
        { 0, 10, 0 },
        { 25, 10, -2 }
      },
      freqRatios {
        { 0, 10, 0 },
        { 15, 10, -1 }
      },
      informationContent {
        { -1, 10, 0 }
      },
      gaplessColumnWeights {
        { 0, 10, 0 }
      },
      sigma {
        { 2, 10, 0 }
      },
      intervalSizes {
        3
      },
      numMatchingSeqs {
        2
      }
    },
    finalData {
      scores {
        -32768,
        5
      },
      lambda { 267, 10, -3 },
      kappa { 41, 10, -3 },
      h { 14, 10, -2 },
      lambdaUngapped { -1, 10, 0 },
      kappaUngapped { -1, 10, 0 },
      hUngapped { -1, 10, 0 }
    }
  },
  params {
    pseudocount 0,
    rpsdbparams {
      matrixName "BLOSUM62",
      gapOpen 11,
      gapExtend 1
    }
  }
}
"#;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:632-778 (quoted at
// `algorithm::psiblast::pssm::pssm_with_parameters`) and c++/src/objects/scoremat/scoremat.asn
// (the member order): the members that psiblast sets, in the module's order, the defaults that
// it sets written, and a quote of the title doubled.
#[test]
fn a_pssm_is_written_in_the_modules_order_and_read_back() {
    let value = small_pssm(b"a \"quoted\" title");
    let text = asn_text(&value).unwrap();
    assert_eq!(String::from_utf8(text.clone()).unwrap(), SMALL_PSSM);
    assert_eq!(read_pssm_with_parameters(&text).unwrap(), value);
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-986, 1477-1538 and 207-274
// (quoted at `Reader::string`, `Reader::octets` and `Reader::skip_white_space`): a string joins
// its lines, an OCTET STRING its hex digits across lines, and comments and white space between
// the values are skipped.
#[test]
fn strings_octets_and_comments_are_read_as_ncbi_reads_them() {
    let wrapped = SMALL_PSSM
        .replace(
            "title \"a \"\"quoted\"\" title\"",
            "title \"a \"\"quo\r\nted\"\" title\" -- a comment -- ",
        )
        .replace("'0C'H", "'0\nC'H")
        .replace("numRows 2,", "numRows 2, -- to the end of the line\n\t");
    let read = read_pssm_with_parameters(wrapped.as_bytes()).unwrap();
    assert_eq!(read, small_pssm(b"a \"quoted\" title"));
    // An odd number of hex digits: the last one is the high half of a byte.
    let odd = SMALL_PSSM.replace("'0C'H", "'0C1'H");
    let read = read_pssm_with_parameters(odd.as_bytes()).unwrap();
    let query = read.pssm.query.unwrap();
    assert_eq!(
        query.inst.seq_data,
        Some(crate::report::serial::pssm::SeqData::Ncbistdaa(vec![
            0x0C, 0x10
        ]))
    );
}

// The reader's errors name the line where it stopped.
#[test]
fn errors_name_their_line() {
    let error = read_pssm_with_parameters(b"Seq-entry ::= { }").unwrap_err();
    assert_eq!(error.line, 1);
    let bad = SMALL_PSSM.replace("byRow FALSE", "byRow MAYBE");
    let error = read_pssm_with_parameters(bad.as_bytes()).unwrap_err();
    assert_eq!(error.line, 6, "{error}");
    let truncated = &SMALL_PSSM[..SMALL_PSSM.len() / 2];
    assert!(read_pssm_with_parameters(truncated.as_bytes()).is_err());
    let out_of_order = SMALL_PSSM.replace(
        "    numRows 2,\n    numColumns 1,",
        "    numColumns 1,\n    numRows 2,",
    );
    let error = read_pssm_with_parameters(out_of_order.as_bytes()).unwrap_err();
    assert_eq!(error.line, 5, "{error}");
    let set = SMALL_PSSM.replace("query seq {", "query set {");
    assert!(read_pssm_with_parameters(set.as_bytes())
        .unwrap_err()
        .message
        .contains("not supported by LOSAT's PSSM reader"));
}

// The members the reader skips are named (NCBI's object keeps them and writes them back: psiblast
// `-in_pssm` with `-out_pssm`, NCBI 2.17.0 probes of PS5b), and an empty `descr` is NCBI's error
// (`CSeq_descr::PostRead`, quoted at `Reader::seq_entry`; NCBI: "Unrecognized format for PSSM").
#[test]
fn skipped_members_are_named_and_an_empty_descr_is_an_error() {
    let skips = |text: &str| {
        read_pssm_with_parameters_and_skips(text.as_bytes())
            .unwrap()
            .1
    };
    assert!(skips(SMALL_PSSM).is_empty());
    let edited = SMALL_PSSM
        .replace(
            "title \"a \"\"quoted\"\" title\"",
            "title \"a \"\"quoted\"\" title\",\n        comment \"a note\"",
        )
        .replace("length 1,", "length 1,\n        topology linear,")
        .replace(
            "      gapExtend 1\n    }\n  }",
            "      gapExtend 1\n    },\n    constraints {\n      nblocks 0,\n      blocks {\n      },\n      loops {\n      }\n    }\n  }",
        );
    let (value, skipped) = read_pssm_with_parameters_and_skips(edited.as_bytes()).unwrap();
    assert_eq!(skipped, ["Seqdesc", "topology", "constraints"]);
    assert_eq!(value, small_pssm(b"a \"quoted\" title"));
    let annot = SMALL_PSSM.replace(
        "        seq-data ncbistdaa '0C'H\n      }\n",
        "        seq-data ncbistdaa '0C'H\n      },\n      annot {\n        {\n          data ftable {\n          }\n        }\n      }\n",
    );
    assert_eq!(skips(&annot), ["annot"]);
    let molinfo = SMALL_PSSM.replace(
        "title \"a \"\"quoted\"\" title\"",
        "molinfo {\n          biomol peptide\n        },\n        title \"a \"\"quoted\"\" title\"",
    );
    assert_eq!(skips(&molinfo), ["molinfo"]);
    let empty = SMALL_PSSM.replace("        title \"a \"\"quoted\"\" title\"\n", "");
    let error = read_pssm_with_parameters(empty.as_bytes()).unwrap_err();
    assert_eq!(error.message, "empty Seq-descr is not allowed");
}

// NCBI's unit-test data c++/src/algo/blast/unit_tests/api/data/pssm_freq_ratios.asn (read by
// psiblast_unit_test.cpp:105): a query with a gi and a swissprot id and descriptors besides its
// title, frequency ratios of 26 rows, the parameters without gap costs.
#[test]
fn ncbis_pssm_with_frequency_ratios_is_read() {
    use crate::report::serial::pssm::SeqData;
    use crate::seqid::model::{ParsedSeqId, TextseqId};
    let value = read_pssm_with_parameters(&gunzip("pssm/pssm_freq_ratios.asn.gz")).unwrap();
    let pssm = &value.pssm;
    assert_eq!((pssm.num_rows, pssm.num_columns), (26, 232));
    assert_eq!((pssm.is_protein, pssm.by_row), (Some(true), Some(false)));
    let query = pssm.query.as_ref().unwrap();
    assert_eq!(
        query.ids,
        vec![
            ParsedSeqId::Gi(129295),
            ParsedSeqId::Swissprot(TextseqId {
                name: Some(b"OVALX_CHICK".to_vec()),
                accession: Some(b"P01013".to_vec()),
                release: None,
                version: None,
            }),
        ]
    );
    assert_eq!(
        query.titles,
        Some(vec![b"Gene X protein (Ovalbumin-related)".to_vec()])
    );
    assert_eq!(query.inst.repr, ("raw", 2));
    assert_eq!(query.inst.mol, ("aa", 3));
    assert_eq!(query.inst.length, Some(232));
    let Some(SeqData::Ncbistdaa(letters)) = &query.inst.seq_data else {
        panic!("ncbistdaa expected");
    };
    assert_eq!(letters.len(), 232);
    assert_eq!(&letters[..4], &[0x0F, 0x09, 0x0A, 0x04]);
    assert_eq!(letters[231], 0x0E);
    let data = pssm.intermediate_data.as_ref().unwrap();
    assert_eq!(data.freq_ratios.len(), 232 * 26);
    assert_eq!(data.freq_ratios[0], 0.0);
    assert_eq!(data.freq_ratios[1].to_bits(), 0x3fad_840f_8594_1c3f);
    assert!(data.weighted_res_freqs_per_pos.is_none());
    assert!(pssm.final_data.is_none());
    let params = value.params.as_ref().unwrap();
    assert_eq!(params.pseudocount, Some(9));
    let rps = params.rpsdbparams.as_ref().unwrap();
    assert_eq!(rps.matrix_name, b"BLOSUM62");
    assert_eq!((rps.gap_open, rps.gap_extend), (None, None));
    assert!(pssm.identifier.is_none() && pssm.row_labels.is_none());
}

// NCBI's unit-test data c++/src/algo/blast/unit_tests/api/data/pssm_zero_fratios.asn (read by
// psiblast_unit_test.cpp:1183-1205): no query, no parameters, the defaults not written, all
// frequency ratios 0, and the final data. Every member is one the reader keeps, so writing the
// value gives NCBI's bytes back, up to the empty line that ends the file (after the object and
// its separator, which the reader does not read).
#[test]
fn ncbis_pssm_with_zero_frequency_ratios_reads_and_writes_back() {
    let text = gunzip("pssm/pssm_zero_fratios.asn.gz");
    let value = read_pssm_with_parameters(&text).unwrap();
    let pssm = &value.pssm;
    assert_eq!((pssm.num_rows, pssm.num_columns), (28, 232));
    assert_eq!((pssm.is_protein, pssm.by_row), (None, None));
    assert!(pssm.query.is_none() && value.params.is_none());
    let data = pssm.intermediate_data.as_ref().unwrap();
    assert_eq!(data.freq_ratios.len(), 232 * 28);
    assert!(data.freq_ratios.iter().all(|&ratio| ratio.to_bits() == 0));
    let final_data = pssm.final_data.as_ref().unwrap();
    assert_eq!(final_data.scores.len(), 232 * 28);
    assert_eq!(final_data.lambda, 0.0);
    assert_eq!(final_data.kappa.to_bits(), 0x3fa4_9d6e_3994_5758);
    assert_eq!(final_data.h, 0.0);
    assert!(text.ends_with(b"}\n\n"));
    assert_same_text(
        &asn_text(&value).unwrap(),
        &text[..text.len() - 1],
        "pssm_zero_fratios.asn",
    );
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-396 (quoted at
// `report::serial::pssm`): psiblast's `-out_pssm` files from NCBI 2.17.0 (`db.ps.asn.*` rows of
// LOSAT/tests/db_search_manifest/psiblast.tsv, frozen when 64 KB or less), without the diagnostics of
// `-out_ascii_pssm` (`rnd`) and with them (`tkrnd`: the second query's), read and written again.
#[test]
fn psiblasts_pssm_files_read_and_write_back_byte_for_byte() {
    for name in ["db.ps.asn.rnd.6.asn", "db.ps.asn.tkrnd.6.asn"] {
        let text = std::fs::read(fixture("db_search").join(name)).unwrap();
        let value = read_pssm_with_parameters(&text).unwrap();
        assert_same_text(&asn_text(&value).unwrap(), &text, name);
        let data = value.pssm.intermediate_data.as_ref().unwrap();
        assert_eq!(
            data.information_content.is_some(),
            name.contains("tkrnd"),
            "{name}"
        );
    }
}

/// Reads every file of the directories in `PSSM_TEXT_DIRS` (separated by `:`) whose name ends in
/// `.asn`, `.smp` or `.pssm` and checks that writing it again gives its bytes (the PSSMs that NCBI
/// 2.17.0 psiblast wrote in step 1 of PS3 and the survey's profiles). Run by hand:
/// `PSSM_TEXT_DIRS=... cargo test --lib -- --ignored pssm_files_of_a_directory`.
#[test]
#[ignore = "reads the files of local directories (PSSM_TEXT_DIRS)"]
fn pssm_files_of_a_directory_read_and_write_back() {
    let dirs = std::env::var("PSSM_TEXT_DIRS").expect("PSSM_TEXT_DIRS");
    let mut paths = Vec::new();
    for dir in dirs.split(':') {
        let mut stack = vec![std::path::PathBuf::from(dir)];
        while let Some(dir) = stack.pop() {
            for entry in std::fs::read_dir(&dir).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    stack.push(path);
                } else if path
                    .extension()
                    .is_some_and(|ext| ext == "asn" || ext == "smp" || ext == "pssm")
                    || path
                        .file_name()
                        .is_some_and(|name| name.to_string_lossy().contains(".asn."))
                {
                    paths.push(path);
                }
            }
        }
    }
    paths.sort();
    let mut failures = Vec::new();
    for path in &paths {
        let text = std::fs::read(path).unwrap();
        if !text.starts_with(b"PssmWithParameters") {
            continue;
        }
        match read_pssm_with_parameters(&text) {
            Ok(value) => {
                if asn_text(&value).unwrap() != text {
                    failures.push(format!("{}: written back differently", path.display()));
                }
            }
            Err(error) => failures.push(format!("{}: {error}", path.display())),
        }
        eprintln!("read {}", path.display());
    }
    assert!(failures.is_empty(), "{failures:#?}");
    assert!(!paths.is_empty());
}

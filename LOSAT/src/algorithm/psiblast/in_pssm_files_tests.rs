use super::*;
use crate::report::serial::asn_types::Value;
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmIntermediateData, PssmParameters, PssmQuery,
    SeqInst, Seqdesc,
};
use crate::seqid::ParsedSeqId;

/// A PSSM of `columns` positions of the letter M whose query has the Seq-id `id`.
fn pssm(id: ObjectId, columns: usize, rows: usize) -> PssmWithParameters {
    PssmWithParameters {
        pssm: Pssm {
            is_protein: Some(true),
            identifier: None,
            num_rows: rows as i32,
            num_columns: columns as i32,
            row_labels: None,
            by_row: Some(false),
            query: Some(SeqEntry::Seq(Box::new(PssmQuery {
                ids: vec![ParsedSeqId::Local(id)],
                descr: None,
                inst: SeqInst::new(
                    ("raw", 2),
                    ("aa", 3),
                    Some(columns as u32),
                    Some(SeqData::Ncbistdaa(vec![12; columns])),
                ),
                annot: None,
            }))),
            intermediate_data: Some(PssmIntermediateData {
                freq_ratios: (0..columns * rows).map(|i| i as f64).collect(),
                ..PssmIntermediateData::default()
            }),
            final_data: None,
        },
        params: Some(PssmParameters {
            pseudocount: Some(0),
            rpsdbparams: Some(FormatRpsDbParameters {
                matrix_name: b"BLOSUM62".to_vec(),
                gap_open: None,
                gap_extend: None,
            }),
            ..PssmParameters::default()
        }),
    }
}

fn request(asn: bool, ascii: bool) -> PssmRequest {
    PssmRequest {
        asn_path: asn.then(|| std::path::PathBuf::from("f.asn")),
        ascii_path: ascii.then(|| std::path::PathBuf::from("f.txt")),
        save_each: false,
        pseudocount: 0,
        versions: Default::default(),
    }
}

// blastfmtutil.cpp:252-259, Pssm.cpp:88-93, sequtil_tables.cpp:1674-1707 (cited at
// `ascii_letters_defined`): NCBI's ASCII PSSM prints a letter of NCBIstdaa only for the letters
// of iupacaa (a lower-case one as its capital), the capitals, `*` and `-` of ncbieaa and the
// codes below 28.
#[test]
fn the_ascii_pssm_letters_are_defined_for_ncbistdaa_residues_only() {
    let with = |data: SeqData| {
        let mut object = pssm(ObjectId::Id(1), 3, 28);
        let Some(SeqEntry::Seq(query)) = object.pssm.query.as_mut() else {
            unreachable!("the test's PSSM has a Bioseq");
        };
        query.inst.seq_data = Some(data);
        ascii_letters_defined(&object)
    };
    assert!(with(SeqData::Ncbistdaa(vec![0, 12, 27])));
    assert!(!with(SeqData::Ncbistdaa(vec![1, 28, 2])));
    assert!(with(SeqData::Iupacaa(b"MKV".to_vec())));
    assert!(with(SeqData::Iupacaa(b"MkV".to_vec())));
    // A lower-case letter prints as its capital (`CIupacaaToStdaa`: M4a-1c LG probe
    // `tk_low_iupac.files`, the ASCII PSSM of the lower-case query is the capital one's).
    let codes = |data: SeqData| {
        let mut object = pssm(ObjectId::Id(1), 3, 28);
        let Some(SeqEntry::Seq(query)) = object.pssm.query.as_mut() else {
            unreachable!("the test's PSSM has a Bioseq");
        };
        query.inst.seq_data = Some(data);
        query_codes(&object)
    };
    assert_eq!(codes(SeqData::Iupacaa(b"mkV".to_vec())), [12, 10, 19]);
    assert_eq!(codes(SeqData::Iupacaa(b"MKV".to_vec())), [12, 10, 19]);
    for bad in [b"M*V", b"M-V", b"M#V"] {
        assert!(!with(SeqData::Iupacaa(bad.to_vec())), "{bad:?}");
    }
    assert!(with(SeqData::Ncbieaa(b"M*-".to_vec())));
    assert!(!with(SeqData::Ncbieaa(b"Mk#".to_vec())));
    let mut file_pssm = pssm(ObjectId::Id(1), 3, 28);
    if let Some(SeqEntry::Seq(query)) = file_pssm.pssm.query.as_mut() {
        query.inst.seq_data = Some(SeqData::Iupacaa(b"M#V".to_vec()));
    }
    let file = InPssmFile {
        object: file_pssm,
        computed: false,
        valid: true,
    };
    let error = written_back_refusal(&file, &request(false, true)).unwrap_err();
    assert!(error.to_string().contains("memory it never set"), "{error}");
    assert!(written_back_refusal(&file, &request(true, false)).is_ok());
}

const GAPPED: KarlinBlk = KarlinBlk {
    lambda: 0.267,
    k: 0.0425,
    h: 0.14,
};

// psiblast_aux_priv.cpp:440-502 (cited at `with_computed_scores`): the ratios of a PSSM of
// fewer rows get zeros up to 28 rows in each column, the scores are the engine's by column, the
// gapped values the engine's, and the parameters get the gap costs.
#[test]
fn a_ratio_only_pssm_gets_the_constructors_scores() {
    let mut value = pssm(ObjectId::Id(7), 2, 26);
    let rows = vec![[1; BLASTAA_SIZE], [2; BLASTAA_SIZE]];
    with_computed_scores(&mut value, &rows, GAPPED, 11, 1);
    let matrix = &value.pssm;
    assert_eq!(matrix.num_rows, 28);
    let ratios = &matrix.intermediate_data.as_ref().unwrap().freq_ratios;
    assert_eq!(ratios.len(), 56);
    assert_eq!(&ratios[24..30], &[24.0, 25.0, 0.0, 0.0, 26.0, 27.0]);
    assert_eq!(&ratios[54..], &[0.0, 0.0]);
    let data = matrix.final_data.as_ref().unwrap();
    assert_eq!(data.scores.len(), 56);
    assert_eq!((data.scores[27], data.scores[28]), (1, 2));
    assert_eq!((data.lambda, data.kappa, data.h), (0.267, 0.0425, 0.14));
    assert_eq!(data.lambda_ungapped, None);
    let rps = value.params.unwrap().rpsdbparams.unwrap();
    assert_eq!((rps.gap_open, rps.gap_extend), (Some(11), Some(1)));
    // Non-default search gap costs (PS7b): `PsiBlastAddAncillaryPssmData` writes the search's.
    let mut custom = pssm(ObjectId::Id(8), 2, 28);
    with_computed_scores(&mut custom, &rows, GAPPED, 9, 2);
    let rps = custom.params.unwrap().rpsdbparams.unwrap();
    assert_eq!((rps.gap_open, rps.gap_extend), (Some(9), Some(2)));
}

// The refusals of `rounds_input`, and what it keeps: the query's Seq-id (an integer stays
// one), the object only when a file is written.
#[test]
fn rounds_input_keeps_the_id_and_refuses_what_ncbi_writes_otherwise() {
    // A file whose PSSM holds the members `unread` that the search does not read.
    let file = |unread: Vec<&'static str>, computed: bool, params: bool| {
        let mut object = pssm(ObjectId::Id(7), 2, 28);
        if !params {
            object.params = None;
        }
        let query = object
            .pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap();
        for name in unread {
            match name {
                "Seqdesc" => {
                    let comment = crate::report::serial::asn_types::member_index(
                        crate::report::serial::asn_types::table::named::SEQDESC,
                        "comment",
                    );
                    query.descr = Some(vec![Seqdesc::Other(Value::Choice(
                        comment as u16,
                        Box::new(Value::Str(b"x".to_vec())),
                    ))]);
                }
                "topology" => query.inst.topology = Some(2),
                "annot" => query.annot = Some(Vec::new()),
                other => unreachable!("{other}"),
            }
        }
        InPssmFile {
            object,
            computed,
            valid: true,
        }
    };
    let refusal = |input: InPssmRounds| input.refuse_written_back().map(|()| input);
    let input = refusal(rounds_input(
        file(vec![], false, true),
        &request(false, false),
    ))
    .unwrap();
    assert!(input.object.is_none());
    let input = refusal(rounds_input(
        file(vec![], false, true),
        &request(true, false),
    ))
    .unwrap();
    assert!(input.object.is_some());
    // The members that the search does not read are written back with the others.
    let input = refusal(rounds_input(
        file(vec!["Seqdesc", "topology", "annot"], false, true),
        &request(true, false),
    ))
    .unwrap();
    let asn = input.files(&request(true, false)).unwrap().asn;
    let asn = String::from_utf8(asn.unwrap()).unwrap();
    assert!(
        asn.contains("comment \"x\"")
            && asn.contains("topology circular")
            && asn.contains("annot {\n      }"),
        "{asn}"
    );
    let error = refusal(rounds_input(
        file(vec![], true, false),
        &request(true, false),
    ))
    .unwrap_err()
    .to_string();
    assert!(error.contains("without params.rpsdbparams"), "{error}");
    assert!(error.contains("with -out_pssm"), "{error}");
    // The one owner of the condition, which the archive of `-outfmt 11` reads too: computed
    // scores without the file's own `rpsdbparams`.
    assert!(file(vec![], true, false).matrix_name_unwritable());
    assert!(!file(vec![], false, false).matrix_name_unwritable());
    assert!(!file(vec![], true, true).matrix_name_unwritable());
    let archive = matrix_name_refusal("in the archive of -outfmt 11");
    assert!(
        archive.contains("rpsdbparams in the archive of -outfmt 11 (NCBI")
            && archive.ends_with("is not supported by LOSAT's PSIBLAST"),
        "{archive}"
    );
    refusal(rounds_input(
        file(vec![], false, false),
        &request(true, false),
    ))
    .unwrap();
    // A PSSM that NCBI's checks refuse writes nothing.
    let mut invalid = file(vec!["annot"], true, false);
    invalid.valid = false;
    assert!(refusal(rounds_input(invalid, &request(true, true)))
        .unwrap()
        .object
        .is_none());
    // The long ASCII form needs every column it reads.
    let mut long = file(vec![], false, true);
    let data = long.object.pssm.intermediate_data.as_mut().unwrap();
    data.information_content = Some(vec![0.5, 0.5]);
    data.num_matching_seqs = Some(vec![2, 2]);
    data.gapless_column_weights = Some(vec![1.0]);
    data.weighted_res_freqs_per_pos = Some(vec![0.0; 56]);
    let error = refusal(rounds_input(long, &request(false, true)))
        .unwrap_err()
        .to_string();
    assert!(error.contains("informationContent"), "{error}");
}

// The review of PS5b, N1 and S1: a longer `weightedResFreqsPerPos` is read as NCBI's
// `Convert2Matrix` reads it (its first cells); a printed weighted frequency whose rounded
// percentage is not an `int`, or is NaN, is refused (NCBI's conversion is undefined behaviour).
#[test]
fn the_long_ascii_form_reads_the_first_cells_and_refuses_weights_outside_int() {
    let long = |weighted: Vec<f64>| {
        let mut object = pssm(ObjectId::Id(1), 1, 28);
        let mut row = [-1; BLASTAA_SIZE];
        row[0] = BLAST_SCORE_MIN;
        with_computed_scores(&mut object, &[row], GAPPED, 11, 1);
        let data = object.pssm.intermediate_data.as_mut().unwrap();
        data.information_content = Some(vec![0.5]);
        data.num_matching_seqs = Some(vec![1]);
        data.gapless_column_weights = Some(vec![0.0]);
        data.weighted_res_freqs_per_pos = Some(weighted);
        check_long_ascii(&object)
    };
    let mut weighted = vec![0.25; 28];
    long(weighted.clone()).unwrap();
    weighted.extend([7.0, 8.0]);
    long(weighted.clone()).unwrap();
    long(vec![0.25; 27]).unwrap_err();
    // Row 0 (`-`) is not printed: its value is not read.
    let mut unprinted = vec![0.25; 28];
    unprinted[0] = f64::NAN;
    long(unprinted).unwrap();
    for value in [1e8, -1e8, f64::NAN, f64::INFINITY] {
        let mut weighted = vec![0.25; 28];
        weighted[1] = value;
        let error = long(weighted).unwrap_err().to_string();
        assert!(error.contains("not an int or is NaN"), "{value}: {error}");
    }
    let mut edge = vec![0.25; 28];
    edge[1] = 21_474_836.47;
    long(edge).unwrap();
}

// psiblast_impl.cpp:150-168 (cited at `written_back`): the six values of `finalData`.
#[test]
fn the_write_back_sets_the_searched_blocks() {
    let mut object = pssm(ObjectId::Str(b"myq".to_vec()), 1, 28);
    with_computed_scores(&mut object, &[[0; BLASTAA_SIZE]], GAPPED, 11, 1);
    let input = InPssmRounds {
        object: Some(object),
        refusal: None,
    };
    let ungapped = KarlinBlk {
        lambda: 0.3176,
        k: 0.134,
        h: 0.401,
    };
    let text = input
        .written_back(&request(true, false), ungapped, GAPPED)
        .unwrap()
        .0
        .unwrap();
    let read = crate::blastinput::pssm_text::read_pssm_with_parameters(&text).unwrap();
    let data = read.pssm.final_data.unwrap();
    assert_eq!(
        (data.lambda_ungapped, data.kappa_ungapped, data.h_ungapped),
        (Some(0.3176), Some(0.134), Some(0.401))
    );
    assert_eq!((data.lambda, data.kappa, data.h), (0.267, 0.0425, 0.14));
    assert_eq!(
        input
            .written_back(&request(false, true), ungapped, GAPPED)
            .unwrap()
            .0,
        None
    );
}

// blastfmtutil.cpp:186-196 (cited at `ascii_pssm`): the short form of a PSSM without its
// information content, with the PSI blocks of its file (0 for those it lacks).
#[test]
fn a_file_pssm_without_information_content_prints_the_short_form() {
    let mut object = pssm(ObjectId::Id(1), 1, 28);
    let mut row = [-1; BLASTAA_SIZE];
    row[1] = 4;
    row[0] = BLAST_SCORE_MIN;
    with_computed_scores(&mut object, &[row], GAPPED, 11, 1);
    let text = String::from_utf8(ascii_pssm(&object)).unwrap();
    let expected = "\nLast position-specific scoring matrix computed\n           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n    1 M    4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 \n\n                      K         Lambda\nPSI Ungapped         0.0000     0.0000\nPSI Gapped           0.0425     0.2670\n";
    assert_eq!(text, expected);
}

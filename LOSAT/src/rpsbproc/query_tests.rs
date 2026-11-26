use super::*;
use crate::rpsbproc::archive::{StdAlign, StdLoc, StdPair};
use crate::rpsbproc::data::DataFiles;
use crate::seqid::model::{Dbtag, TextseqId};

/// The data files of the fixtures (`tests/fixtures/rpsbproc/data`).
fn fixture_data() -> CddData {
    let file = |name: &str| {
        format!(
            "{}/tests/fixtures/rpsbproc/data/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
        .into_bytes()
    };
    CddData::load(&DataFiles {
        cddids: file("cddid.tbl"),
        spfeats: file("cddannot.dat"),
        genfeats: file("cddannot_generic.dat"),
        bscores: file("bitscore_specific.txt"),
        cdtrack: file("cdtrack.txt"),
        famlinks: file("family_superfamily_links"),
    })
    .unwrap()
}

// The refusals of a nucleotide query's `Std-seg`s (cdalignproc_biodata.cpp:195-232,
// Seq_align.cpp:339-358 and rpsbproc_main.cpp:1809,1836, cited at `process_cd_query` and
// `query_row_id`): a first segment with fewer than 2 ids is read past its end (PS16c oracle
// `noids1`: segmentation fault in both 0.5.1 builds); a first location of a variant whose
// id LOSAT does not look up. A kept hit whose range starts or ends at 2^31 - 1, whose
// DOMAINS line overflows NCBI's `int` (PS16c review oracle `plusmax`, PS16-Q1 oracle
// `ntplusmax` and its nearby input), is kept whatever the strand, with its range
// (PD-LOSAT-NCBI-DEFECTS exception 12). A hit whose PSSM-ID is 0 is dropped before it is
// classified or printed.
#[test]
fn refusals_of_translated_hits() {
    let data = fixture_data();
    let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
    let cdd = |tag: i32| {
        Some(ParsedSeqId::General(Dbtag {
            db: b"CDD".to_vec(),
            tag: ObjectId::Id(tag),
        }))
    };
    let entries = [SeqEntry::Seq(Query {
        ids: vec![local.clone()],
        mol: 1,
        length: Some(6),
        seq_data: Some(("iupacna", b"ACGTAC".to_vec())),
        ..Query::default()
    })];
    let run =
        |query_id: QueryRowId, hit_id: Option<ParsedSeqId>, from: u32, to: u32, strand: i32| {
            let align = Align {
                segs: Segs::Std(Box::new(StdAlign {
                    query_id,
                    hit_id,
                    segs: vec![StdPair {
                        query: StdLoc::Int {
                            from,
                            to,
                            strand: Some(strand),
                        },
                        hit: StdLoc::Int {
                            from: 0,
                            to: 2,
                            strand: None,
                        },
                    }],
                })),
                score: Some(10),
                e_value: Some(1e-5),
                bit_score: Some(40.5),
                num_ident: None,
                score_of_other_type: None,
            };
            match process_blast_results(&entries, &[align], &data) {
                Ok(queries) => Ok(queries[0]
                    .aligns
                    .iter()
                    .map(|align| (align.from, align.to))
                    .collect::<Vec<_>>()),
                Err(ProcessError::Unsupported(message)) => Err(message),
                Err(other) => panic!("{other:?}"),
            }
        };
    let row = || QueryRowId::Id(local.clone());
    let refused = |message: &str| Err(format!("{message} is not supported by LOSAT's rpsbproc"));
    assert_eq!(
        run(row(), None, 0, 2, 1),
        refused("a Seq-align whose first Std-seg has fewer than 2 ids (NCBI reads past them)")
    );
    assert_eq!(
        run(QueryRowId::Location("packed-int"), cdd(5), 0, 2, 1),
        refused("a Std-seg without ids whose first location is a packed-int")
    );
    let end = 0x7fff_ffff;
    for strand in [0, 1, 3, 255] {
        assert_eq!(
            run(row(), cdd(5), 0, end, strand),
            Ok(vec![(0, i64::from(end))])
        );
    }
    let end = i64::from(end);
    assert_eq!(
        run(row(), cdd(5), 0x7fff_ffff, 0x7fff_ffff, 1),
        Ok(vec![(end, end)])
    );
    assert_eq!(run(row(), cdd(0), 0, 0x7fff_ffff, 1), Ok(vec![]));
}

// `Get1LetterSeqData` (biodata_core.cpp:199-253, cited at `one_letter_residues`): the
// literals of `ext delta` one after another (NCBIstdaa converted), cut to `inst.length`;
// the `seq-data` first when both are set; a literal in an encoding LOSAT does not convert is
// refused as a `seq-data` in it is.
#[test]
fn residues_of_delta_literals() {
    let query =
        |seq_data: Option<(&'static str, Vec<u8>)>, delta: Vec<(&'static str, Vec<u8>)>| Query {
            ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))],
            mol: 3,
            length: Some(5),
            seq_data,
            delta_literals: delta,
            ..Query::default()
        };
    let residues = |query: Query| fill_sequence(&query).map(|(_, _, _, residues)| residues);
    let delta = vec![
        ("ncbieaa", b"MK".to_vec()),
        ("ncbistdaa", vec![1, 25, 27, 0]),
    ];
    assert_eq!(
        residues(query(None, delta.clone())).ok(),
        Some(b"MKA*J".to_vec())
    );
    assert_eq!(
        residues(query(Some(("iupacaa", b"WW".to_vec())), delta)).ok(),
        Some(b"WW".to_vec())
    );
    let Err(ProcessError::Unsupported(message)) = residues(query(
        None,
        vec![("ncbieaa", b"MK".to_vec()), ("ncbi8aa", vec![1])],
    )) else {
        panic!("a literal of NCBI8aa is refused");
    };
    assert_eq!(
        message,
        "a query whose Seq-data of a delta literal is ncbi8aa is not supported by LOSAT's rpsbproc"
    );
}

// `SeqDataTo1LtrSeq` of the nucleotide encodings (biodata_core.cpp:149-197 and the
// `CSeqportUtil` tables, cited at `one_letter`): NCBI2na 4 bases per byte from the high
// bits, NCBI4na 2 (code 0 is `N`), NCBI8na 1 (codes from 16 the byte 0xff), IUPACna as it
// is; cut to `inst.length`; the literals of a delta of the kind of the first one only
// (`litType`); NCBIpna throws `CSeqportUtil`'s `runtime_error`.
#[test]
fn residues_of_nucleotide_encodings() {
    let query = |length: u32,
                 seq_data: Option<(&'static str, Vec<u8>)>,
                 delta: Vec<(&'static str, Vec<u8>)>| Query {
        ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))],
        mol: 1,
        length: Some(length),
        seq_data,
        delta_literals: delta,
        ..Query::default()
    };
    let residues = |query: Query| fill_sequence(&query).map(|(_, _, _, residues)| residues);
    // 0x1B = 00 01 10 11, 0xE4 = 11 10 01 00.
    assert_eq!(
        residues(query(7, Some(("ncbi2na", vec![0x1b, 0xe4])), vec![])).ok(),
        Some(b"ACGTTGC".to_vec())
    );
    assert_eq!(
        residues(query(4, Some(("ncbi4na", vec![0x01, 0x2f])), vec![])).ok(),
        Some(b"NACN".to_vec())
    );
    assert_eq!(
        residues(query(3, Some(("ncbi8na", vec![4, 15, 16])), vec![])).ok(),
        Some(vec![b'G', b'N', 0xff])
    );
    assert_eq!(
        residues(query(9, Some(("iupacna", b"acgTUx".to_vec())), vec![])).ok(),
        Some(b"acgTUx".to_vec())
    );
    let delta = vec![
        ("iupacna", b"AC".to_vec()),
        ("ncbieaa", b"MK".to_vec()),
        ("ncbi2na", vec![0xff]),
    ];
    assert_eq!(
        residues(query(10, None, delta)).ok(),
        Some(b"ACTTTT".to_vec())
    );
    let Err(ProcessError::Std(what)) = residues(query(4, Some(("ncbipna", vec![0; 5])), vec![]))
    else {
        panic!("NCBIpna throws");
    };
    assert_eq!(what, "Requested conversion not implemented");
}

// serial/exception.cpp:86-150 (cited at `choice_selection`); NCBI's rpsbproc printed the first
// one for a profile database with `lcl|` ids (PS16a probe 2).
#[test]
fn choice_selection_messages() {
    let message = |error| match error {
        ProcessError::Exception(exception) => {
            String::from_utf8(exception.get_msg().to_vec()).unwrap()
        }
        ProcessError::Unsupported(message) | ProcessError::Std(message) => {
            panic!("{message}")
        }
    };
    assert_eq!(
        message(choice_selection(
            "NCBI-Seqloc",
            "Seq-id",
            "general",
            "local"
        )),
        "CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local"
    );
    assert_eq!(
        message(choice_selection("NCBI-Blast4", "Blast4-cutoff", "e-value", "raw-score")),
        "CBlast4_cutoff::GetE_value(): Invalid choice selection: NCBI-Blast4::Blast4-cutoff.raw-score"
    );
    assert_eq!(
        message(choice_selection(
            "NCBI-Seqloc",
            "Seq-loc",
            "packed-int",
            "int"
        )),
        "CSeq_loc::GetPacked_int(): Invalid choice selection: NCBI-Seqloc::Seq-loc.int"
    );
}

// serialobject.cpp:266-306 (cited at `unassigned_member`); NCBI's rpsbproc printed it for a
// query without `inst.length` (PS16a review oracle O3-9).
#[test]
fn unassigned_member_message() {
    let ProcessError::Exception(exception) =
        unassigned_member("NCBI-Sequence", "Seq-inst", "length")
    else {
        panic!("an exception");
    };
    assert_eq!(exception.site().kind, "CUnassignedMember::eGet");
    assert_eq!(
        exception.get_msg(),
        b"CSeq_inst::GetLength(): Attempt to get unassigned member NCBI-Sequence::Seq-inst.length"
    );
}

// Seq_id.cpp:2292-2326 and biodata_core.cpp:1030-1046 (cited at `query_accession`): the
// label of the best id, after its first `:`, then after its first `|`.
#[test]
fn query_labels() {
    let query = |ids: Vec<ParsedSeqId>| Query {
        ids,
        mol: 3,
        ..Query::default()
    };
    let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
    assert_eq!(query_accession(&query(vec![local.clone()])), b"Query_1");
    let general = ParsedSeqId::General(Dbtag {
        db: b"CDD".to_vec(),
        tag: ObjectId::Id(5),
    });
    assert_eq!(query_accession(&query(vec![general])), b"5");
    let refseq = ParsedSeqId::Other(TextseqId {
        accession: Some(b"NP_414783".to_vec()),
        version: Some(1),
        ..TextseqId::default()
    });
    // FastaAARank prefers the RefSeq id to the local one.
    assert_eq!(query_accession(&query(vec![local, refseq])), b"NP_414783.1");
    assert_eq!(query_accession(&query(vec![ParsedSeqId::Gi(123)])), b"123");
    assert_eq!(query_accession(&query(Vec::new())), b"");
}

// `CBioseq::IsNa` (Bioseq.cpp:345-348, Seq_inst.hpp:90-95) and `GetDescrTitleFromBioseq`
// (biodata_core.cpp:297-366, ncbiobj.hpp:1047-1054): `dna` (1), `rna` (2) and `na` (4) are NA
// while `not-set` (0), `aa` (3) and `other` (255) are not; an explicit `title` descriptor (even
// empty) wins over `annot`, else `AnnotTitle::Title` is used, else empty; a `NULL` `Seq-annot` or
// `Seq-feat` reached before a matching feature is refused.
#[test]
fn query_title_and_is_na_mol() {
    assert!(!is_na_mol(0));
    assert!(is_na_mol(1));
    assert!(is_na_mol(2));
    assert!(!is_na_mol(3));
    assert!(is_na_mol(4));
    assert!(!is_na_mol(255));

    let base = |title: Option<Vec<u8>>, annot_title: AnnotTitle| Query {
        ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))],
        title,
        annot_title,
        mol: 3,
        length: Some(1),
        seq_data: Some(("ncbieaa", b"M".to_vec())),
        ..Query::default()
    };
    let title_of = |q: Query| fill_sequence(&q).map(|(_, t, _, _)| t);

    assert_eq!(
        title_of(base(
            Some(b"explicit".to_vec()),
            AnnotTitle::Title(b"feat".to_vec())
        ))
        .ok(),
        Some(b"explicit".to_vec())
    );
    assert_eq!(
        title_of(base(Some(Vec::new()), AnnotTitle::Title(b"feat".to_vec()))).ok(),
        Some(Vec::new())
    );
    assert_eq!(
        title_of(base(
            Some(b"ok".to_vec()),
            AnnotTitle::NullElement("a NULL Seq-annot in Bioseq.annot")
        ))
        .ok(),
        Some(b"ok".to_vec())
    );
    assert_eq!(
        title_of(base(None, AnnotTitle::Title(b"from prot".to_vec()))).ok(),
        Some(b"from prot".to_vec())
    );
    assert_eq!(
        title_of(base(None, AnnotTitle::None)).ok(),
        Some(Vec::new())
    );
    let Err(ProcessError::Unsupported(msg)) = title_of(base(
        None,
        AnnotTitle::NullElement("a NULL Seq-annot in Bioseq.annot"),
    )) else {
        panic!("expected refusal for NULL Seq-annot");
    };
    assert!(msg.contains("NULL Seq-annot"), "{msg}");
}

//! blast_formatter's tabular driver: the formatter's matrix, the options' matrix name and the
//! ids that the rows show for a local id.

use super::*;

#[test]
fn the_formatter_matrix_is_ncbis_standard_matrix_of_the_name() {
    assert_eq!(formatter_matrix(None), None);
    assert_eq!(formatter_matrix(Some(b"")), None);
    assert_eq!(
        formatter_matrix(Some(b"BLOSUM62")),
        Some(ScoringMatrix::Blosum62)
    );
    assert_eq!(
        formatter_matrix(Some(b"blosum45")),
        Some(ScoringMatrix::Blosum45)
    );
    assert_eq!(formatter_matrix(Some(b"PAM30")), Some(ScoringMatrix::Pam30));
    assert_eq!(
        formatter_matrix(Some(b"pam250")),
        Some(ScoringMatrix::Pam250)
    );
    assert_eq!(formatter_matrix(Some(b"BLOSUM62x")), None);
    assert_eq!(formatter_matrix(Some(b"PAM40")), None);
    assert_eq!(
        formatter_matrix(Some(b"IDENTITYxyz")),
        Some(ScoringMatrix::Identity)
    );
    assert_eq!(formatter_matrix(Some(b"other")), None);
}

#[test]
fn a_local_id_shows_the_first_word_of_the_title() {
    let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
    let word = |ids: Vec<ParsedSeqId>| match &ids[0] {
        ParsedSeqId::Local(ObjectId::Str(text)) => String::from_utf8(text.clone()).unwrap(),
        other => panic!("{other:?}"),
    };
    assert_eq!(word(replaced_ids([&local], b"mq1 a query", false)), "mq1");
    assert_eq!(
        word(replaced_ids([&local], b"mq1 a query", true)),
        "Query_1"
    );
    assert_eq!(word(replaced_ids([&local], b"", false)), "Query_1");
    assert_eq!(word(replaced_ids([&local], b" leading", false)), "Query_1");
    let numbered = ParsedSeqId::Local(ObjectId::Id(7));
    assert_eq!(word(replaced_ids([&numbered], b"", false)), "7");
}

#[test]
fn only_the_covered_programs_are_written() {
    for program in [
        EProgram::Blastn,
        EProgram::Megablast,
        EProgram::Blastp,
        EProgram::PsiBlast,
        EProgram::DeltaBlast,
        EProgram::PhiBlastp,
        EProgram::RpsBlast,
    ] {
        assert!(covered(program), "{program:?}");
    }
    for program in [
        EProgram::Blastx,
        EProgram::Tblastn,
        EProgram::Tblastx,
        EProgram::RpsTblastn,
        EProgram::PsiTblastn,
        EProgram::PhiBlastn,
    ] {
        assert!(!covered(program), "{program:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1153
// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:986-1117
/// The prefetch's `CanResolveRange` makes the map of each subject it reaches and stops at the
/// first throw; `GetBioseqLength` then reads a made map (probes of agent E3, 2026-11-27: t_huge1,
/// t_huge12, t_stop1, t_far1, t_dlen; probe L1 of agent E3L: t_dlen2).
#[test]
fn the_prefetch_makes_the_maps_up_to_the_first_throw() {
    use crate::report::serial::asn_input::Names;
    use crate::report::serial::asn_types::{read_document, table::named};
    // 0: raw, 4294967295 declared; 1: raw of 8; 2: a delta declaring 4294967295 whose literals
    // hold 8; 3: a delta of 12 whose second literal declares 4294967295.
    let insts = [
        "repr raw, mol na, length 4294967295, seq-data iupacna \"ACGTACGT\"",
        "repr raw, mol na, length 8, seq-data iupacna \"ACGTACGT\"",
        "repr delta, mol na, length 4294967295, ext delta { literal { length 4, seq-data \
         iupacna \"ACGT\" }, literal { length 4, seq-data iupacna \"ACGT\" } }",
        "repr delta, mol na, length 12, ext delta { literal { length 8, seq-data iupacna \
         \"ACGTACGT\" }, literal { length 4294967295, seq-data iupacna \"ACGT\" } }",
    ];
    let seqs: Vec<String> = insts
        .iter()
        .enumerate()
        .map(|(n, inst)| format!("seq {{ id {{ local str \"s{n}\" }}, inst {{ {inst} }} }}"))
        .collect();
    let text = format!("Bioseq-set ::= {{ seq-set {{ {} }} }}", seqs.join(", "));
    let (value, posts) = read_document(text.as_bytes(), named::BIOSEQ_SET, Names::STRATEGY);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.unwrap();
    let set = V::new(named::BIOSEQ_SET, &value);
    let bioseqs: Vec<Bioseq<'_>> = (0..4)
        .map(|n| Bioseq::of(set.field("seq-set").item(n).unwrap().choice().1))
        .collect();
    let bioseq_of = |key: BioseqKey| match key {
        BioseqKey::Subject(at) => bioseqs[at],
        BioseqKey::Query(_) => unreachable!("subjects only"),
    };
    let made_by = |order: &[usize]| {
        let intervals: Vec<_> = order
            .iter()
            .map(|&at| (Some(BioseqKey::Subject(at)), 0, 3))
            .collect();
        let mut made = MadeMaps::new(0, 4);
        resolve_maps(&intervals, &bioseq_of, &mut made).unwrap();
        made
    };
    let length = |made: &MadeMaps, at: usize| {
        made.map_length(BioseqKey::Subject(at), bioseqs[at])
            .map(|length| length.map_err(|_| ()))
    };
    // t_huge1: the first subject throws (made), the second is not reached.
    let made = made_by(&[0, 1]);
    assert_eq!(made.subjects, [true, false, false, false]);
    assert_eq!(length(&made, 0), Some(Err(())));
    // t_huge12: the second subject of 4294967295 keeps its declared length.
    let made = made_by(&[0, 0, 2]);
    assert_eq!(made.subjects, [true, false, false, false]);
    assert_eq!(length(&made, 2), None);
    // t_dlen: a delta that resolves gives the length of its literals; the rows after go on.
    let made = made_by(&[2, 1, 0]);
    assert_eq!(made.subjects, [true, true, true, false]);
    assert_eq!(length(&made, 2), Some(Ok(8)));
    assert_eq!(
        length(&made, 1),
        None,
        "a declared length other than 4294967295"
    );
    // t_stop1, t_far1: the move past the range's last segment resolves the next one, of
    // 4294967295 letters, whatever the range: the prefetch stops there.
    let made = made_by(&[3, 0]);
    assert_eq!(made.subjects, [false, false, false, true]);
    assert_eq!(length(&made, 3), None, "12 declared");
    // A row's read makes the map of its Bioseq: a later row reads the map's length.
    let mut made = MadeMaps::new(0, 4);
    assert_eq!(length(&made, 2), None);
    made.make(BioseqKey::Subject(2));
    assert_eq!(length(&made, 2), Some(Ok(8)));
}

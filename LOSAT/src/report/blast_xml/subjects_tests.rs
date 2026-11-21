use super::*;

// Line DBID pfoff-1 (DECISIONS 62; review 1, M1): the `Hit` of a database sequence whose
// lookup throws where a batch's prefetch is off has no names: by the one answer of the reports'
// loads (`HitLoads`), the cache keeps the lookup's exception for the `Hit` and raises it there
// (`CaughtException`); a subject shown only after it gets no names (nothing of it is made or
// refused); the others are named as before.
#[test]
fn the_hit_of_a_subject_whose_lookup_throws_raises_the_lookup_exception() {
    use crate::report::db_subject::{FoundInstead, HitLoads, LookupFailure};
    let set = |accession: &str| crate::seqid::DeflineSet {
        deflines: vec![crate::seqid::BlastDefline {
            title: Some(b"t".to_vec()),
            ids: vec![crate::seqid::ParsedSeqId::Genbank(
                crate::seqid::model::TextseqId {
                    accession: Some(accession.as_bytes().to_vec()),
                    version: Some(1),
                    ..Default::default()
                },
            )],
            ..Default::default()
        }],
        utf8: false,
    };
    let failure = LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO);
    let exception = failure.exception().unwrap();
    let thrown = ParsedSubject::new(set("AAA00001"))
        .unwrap()
        .found_instead(Some(FoundInstead::Failed(failure)));
    let found = |accession| Some(Arc::new(ParsedSubject::new(set(accession)).unwrap()));
    let parsed = vec![Some(Arc::new(thrown)), found("AAA00002"), found("AAA00003")];
    // One report read one by one: rank 1, then rank 0 (whose lookup throws), then rank 2.
    let loads = HitLoads::of(
        Some(crate::report::db_subject::RunThrow::every_row(
            crate::report::db_subject::AfterThrow::EndsReport,
        )),
        [(true, vec![Some(1), Some(0), Some(2)])],
        |rank: usize| parsed[rank].as_deref(),
        |_, _| 0,
    );
    let titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"t"[..]); 3];
    let subjects = XmlSubjects::Database {
        ordinals: &[0, 1, 2],
        oids: crate::seqdb::OidMap::unrestricted(),
        titles: XmlTitles::Shared(&titles),
        taxids: &[],
        taxonomy: &[],
        length: &|_| 10,
        parsed: &parsed,
        parsed_taxonomy: &[],
        loads: &loads,
    };
    // The names themselves do not decide the throw.
    assert!(subjects.names(0, "BLASTP").is_ok());
    let mut cache = SubjectNameCache::new(&subjects);
    // Made before the report (a refusal would be the error), raised at the `Hit`.
    let held = cache.get(0, "BLASTP").unwrap();
    assert!(matches!(held, HitNames::Throws(kept) if *kept == exception));
    let raised = held.names().err().unwrap();
    assert!(raised.is::<crate::report::blast_xml::CaughtException>());
    assert!(cache.get(1, "BLASTP").unwrap().names().is_ok());
    let unread = cache.get(2, "BLASTP").unwrap();
    assert!(matches!(unread, HitNames::Unread));
    assert!(!unread
        .names()
        .err()
        .unwrap()
        .is::<crate::report::blast_xml::CaughtException>());
    // A refusal of the names of a loaded subject passes through; nothing is made for one that
    // no report loads.
    assert!(HitNames::of(&loads, 1, || Err(anyhow::anyhow!("refused"))).is_err());
    assert!(matches!(
        HitNames::of(&loads, 2, || Err(anyhow::anyhow!("refused"))),
        Ok(HitNames::Unread)
    ));
    // Where the run writes another format, every subject is loaded and named.
    let every = HitLoads::of(
        None,
        [(true, vec![Some(0), Some(2)])],
        |rank: usize| parsed[rank].as_deref(),
        |_, _| 0,
    );
    assert!(matches!(
        HitNames::of(&every, 0, || subjects.names(0, "BLASTP")),
        Ok(HitNames::Names(_))
    ));
}

// Line DBID N3: following lists of joined subjects, or of one subject, make one `Hit`; a list
// without HSPs shows nothing and does not part them.
#[test]
fn lists_of_joined_subjects_make_one_hit() {
    let lists = [
        Some(4),
        Some(4),
        None,
        Some(9),
        Some(2),
        None,
        Some(3),
        Some(5),
    ];
    let joins = |previous: usize, s_idx: usize| matches!((previous, s_idx), (4, 9) | (3, 5));
    assert_eq!(hit_ranges(&lists, |list| *list, joins), [0..4, 4..5, 6..8]);
    assert_eq!(
        hit_ranges(&lists, |list| *list, |_, _| false),
        [0..2, 3..4, 4..5, 6..7, 7..8]
    );
    assert!(hit_ranges(&[None::<usize>], |list| *list, joins).is_empty());
}

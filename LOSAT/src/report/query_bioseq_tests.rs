use super::*;
use crate::blastinput::fasta_reader::believed::BelievedIds;
use crate::report::archive::bioseq::SeqLocStrand;

fn believed(text: &str) -> BelievedIds {
    BelievedIds::new(
        crate::seqid::parse::parse_ids(text.as_bytes(), 0, &mut |_| {})
            .unwrap()
            .into_iter()
            .flatten()
            .collect(),
        SeqLocStrand::Unknown,
    )
}

fn keys(text: &str) -> Vec<u64> {
    let mut keys = Vec::new();
    for id in &believed(text).ids {
        crate::seqid::defline_set::scope_keys(id, &mut keys);
    }
    keys
}

/// The run's searched queries `records` (ids, title, length) as `QueryBioseqs::note` reads them.
fn records<'a>(
    records: &'a [(BelievedIds, &'static [u8], usize)],
) -> impl Fn(usize) -> Option<QueryRecordRef<'a>> {
    move |h| {
        records
            .get(h)
            .map(|(believed, title, length)| QueryRecordRef {
                believed,
                title,
                length: *length,
            })
    }
}

#[test]
fn notes_the_query_of_the_span_that_holds_the_id() {
    let queries = [
        (believed("lcl|qq2"), &b"query two"[..], 270),
        (believed("gb|AAQ00001.1|"), &b"query one"[..], 179),
    ];
    let mut seen = QueryBioseqs::default();
    assert!(seen.is_empty());
    assert_eq!(
        seen.numbered_from(0).of(0, 7).map(QueryBioseq::length),
        None
    );
    seen.note(
        0,
        7,
        &keys("gb|AAQ00001.1|"),
        0..usize::MAX,
        0..usize::MAX,
        records(&queries),
        true,
    )
    .unwrap();
    let found = seen.of(0, 7).unwrap();
    assert_eq!(found.length(), 179);
    assert_eq!(found.tabular_ids()[0], b"gb|AAQ00001.1|".to_vec());
    assert_eq!(found.tabular_ids()[2], b"AAQ00001.1".to_vec());
    assert_eq!(found.taxonomy().taxid, 0);
    // Another query of the span reports the same sequence: the same holder.
    seen.note(
        1,
        7,
        &keys("gb|AAQ00001.1|"),
        0..usize::MAX,
        0..usize::MAX,
        records(&queries),
        true,
    )
    .unwrap();
    assert!(std::ptr::eq(seen.of(0, 7).unwrap(), seen.of(1, 7).unwrap()));
    // By the call's queries: a call that starts at the run's query 1.
    let from = seen.numbered_from(1);
    assert!(from.any());
    assert_eq!(from.of(0, 7).map(QueryBioseq::length), Some(179));
    assert_eq!(from.row(0, Some(7), 3), SEEN_ROW);
    assert_eq!(from.row(0, Some(8), 3), 3);
    assert_eq!(from.row(0, None, 3), 3);
}

#[test]
fn a_pdb_query_holds_the_chains_of_its_handle() {
    let queries = [(believed("pdb|1AB1|A"), &b"pdb query"[..], 46)];
    let mut seen = QueryBioseqs::default();
    // A sequence stored as pdb 1AB1 with chain 321 (no chain-id): its handle normalizes chain
    // 321 (0x141) to 'A' (65), matching the query `pdb|1AB1|A`.
    let mut chain_321_keys = Vec::new();
    let chain_321 = crate::seqid::model::ParsedSeqId::Pdb(crate::seqid::model::PdbSeqId {
        mol: b"1AB1".to_vec(),
        chain: Some(321),
        rel: None,
        chain_id: None,
    });
    crate::seqid::defline_set::scope_keys(&chain_321, &mut chain_321_keys);
    seen.note(
        0,
        3,
        &chain_321_keys,
        0..usize::MAX,
        0..usize::MAX,
        records(&queries),
        true,
    )
    .unwrap();
    assert_eq!(
        seen.of(0, 3).unwrap().tabular_ids()[0],
        b"pdb|1AB1|A".to_vec()
    );
    assert_eq!(seen.of(0, 3).unwrap().tax_deflines().len(), 1);
}

#[test]
fn keeps_the_refusal_where_the_holder_is_not_determined() {
    let twice = [
        (believed("gb|AAQ00001.1|"), &b"one"[..], 10),
        (believed("gb|AAQ00001.1|"), &b"again"[..], 12),
    ];
    let key = keys("gb|AAQ00001.1|");
    let mut seen = QueryBioseqs::default();
    // Two queries of the span hold the id.
    assert_eq!(
        seen.note(
            0,
            1,
            &key,
            0..usize::MAX,
            0..usize::MAX,
            records(&twice),
            true
        ),
        Err(NotNoted::NoHolder)
    );
    // The caller does not have every query of the span.
    let once = [(believed("gb|AAQ00001.1|"), &b"one"[..], 10)];
    let mut seen = QueryBioseqs::default();
    assert_eq!(
        seen.note(0, 1, &key, 0..usize::MAX, 0..1, records(&once), true),
        Err(NotNoted::NoHolder)
    );
    // No query of the span holds it (the next batch's query is not in the span).
    let mut seen = QueryBioseqs::default();
    assert_eq!(
        seen.note(
            1,
            1,
            &keys("gb|AAQ00009.1|"),
            1..2,
            0..usize::MAX,
            records(&twice),
            true
        ),
        Err(NotNoted::NoHolder)
    );
    assert!(seen.is_empty());
}

#[test]
fn tables_that_read_letters_keep_the_refusal() {
    use crate::blastinput::app::Program;
    let spec = |text: &str| {
        let choice = crate::blastinput::app::parse_formatting_string(text).unwrap();
        crate::report::tabular::TabularSpec::of_app(
            crate::blastinput::app::ReportFormat::from_number(choice.number).unwrap(),
            &choice,
            crate::report::archive::options::RequestKind::Blastp,
        )
        .expect("a field list")
    };
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid stitle slen staxids")),
        Program::Blastp
    ));
    assert!(tabular_shows(Some(&spec("18")), Program::Blastp));
    assert!(!tabular_shows(
        Some(&spec("6 qseqid sseqid sseq")),
        Program::Blastp
    ));
    assert!(!tabular_shows(
        Some(&spec("6 qseqid pident")),
        Program::Blastp
    ));
    assert!(tabular_shows(
        Some(&spec("6 qseqid pident")),
        Program::Blastn
    ));
    assert!(!tabular_shows(
        Some(&spec("6 qseqid sframe")),
        Program::Tblastx
    ));
    // The default fields: `pident` and `mismatch` read letters outside BLASTN.
    assert!(!tabular_shows(None, Program::Blastp));
    assert!(tabular_shows(None, Program::Blastn));
}

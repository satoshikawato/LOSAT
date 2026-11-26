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
    assert!(from.of(0, 8).is_none());
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
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid qseq qlen")),
        Program::Blastp
    ));
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
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid qseq qlen")),
        Program::Blastn
    ));
    // BLASTN `qseq` + `pident`/`nident`/`mismatch` enters the sequence block and recomputes
    // identities from the found query Bioseq.
    assert!(!tabular_shows(
        Some(&spec("6 qseqid qseq pident")),
        Program::Blastn
    ));
    // TBLASTX allows `qseq` (length-checked per HSP) and `sframe`/`frames`.
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid qseq qlen")),
        Program::Tblastx
    ));
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid sframe frames sstrand qframe")),
        Program::Tblastx
    ));
    // TBLASTN allows `sframe`/`frames` (using query Bioseq length) and refuses `qseq`.
    assert!(tabular_shows(
        Some(&spec("6 qseqid sseqid sframe frames sstrand")),
        Program::Tblastn
    ));
    assert!(!tabular_shows(
        Some(&spec("6 qseqid qseq")),
        Program::Tblastn
    ));
    // BLASTX refuses every field that enters the alignment block (`CreateTranslatedDensegFromNADenseg`
    // throws `Widths already exist`) and allows score/id/coverage fields and standalone `sframe`/`frames`.
    for field in [
        "qstart", "qend", "sstart", "send", "length", "gapopen", "gaps", "qframe", "sstrand",
        "qseq", "sseq", "pident", "nident", "mismatch", "positive", "ppos", "btop",
    ] {
        assert!(
            !tabular_shows(Some(&spec(&format!("6 qseqid {field}"))), Program::Blastx),
            "BLASTX should refuse {field}"
        );
    }
    for field in [
        "sseqid", "stitle", "slen", "evalue", "bitscore", "score", "qcovs", "qcovhsp", "qcovus",
        "sframe", "frames",
    ] {
        assert!(
            tabular_shows(Some(&spec(&format!("6 qseqid {field}"))), Program::Blastx),
            "BLASTX should show {field}"
        );
    }
    // The default fields: `pident` and `mismatch` read letters outside BLASTN.
    assert!(!tabular_shows(None, Program::Blastp));
    assert!(tabular_shows(None, Program::Blastn));
    assert!(!tabular_shows(None, Program::Blastx));
}

#[test]
fn ncbi_frame_matches_get_frame_formula() {
    // Plus strand: `(start0 % 3) + 1` regardless of Bioseq length.
    assert_eq!(ncbi_frame(true, 301, 50), 1);
    assert_eq!(ncbi_frame(true, 302, 50), 2);
    assert_eq!(ncbi_frame(true, 303, 50), 3);
    // Minus strand: `-(((len - start0 - 1) % 3) + 1)` with C++ signed remainder.
    assert_eq!(ncbi_frame(false, 659, 660), -2);
    assert_eq!(ncbi_frame(false, 659, 661), -3);
    assert_eq!(ncbi_frame(false, 659, 662), -1);
    assert_eq!(ncbi_frame(false, 659, 50), -1);
    assert_eq!(ncbi_frame(false, 659, 51), 1);
    assert_eq!(ncbi_frame(false, 659, 52), 0);
    assert_eq!(ncbi_frame(false, 360, 220), 1);
    assert_eq!(ncbi_frame(false, 360, 221), 0);
    assert_eq!(ncbi_frame(false, 360, 222), -1);
    assert_eq!(ncbi_frame(false, 194, 220), -3);
}

#[test]
fn hides_ordinal_requires_every_reporting_query() {
    let queries = [(believed("pdb|1AB1|A"), &b"q"[..], 46)];
    let mut seen = QueryBioseqs::default();
    assert!(!seen.hides_ordinal([0], 3));
    seen.note(
        0,
        3,
        &keys("pdb|1AB1|A"),
        0..1,
        0..1,
        records(&queries),
        true,
    )
    .unwrap();
    assert!(seen.hides_ordinal([0], 3));
    assert!(!seen.hides_ordinal([0, 1], 3));
    assert!(!seen.hides_ordinal(std::iter::empty(), 3));
    assert!(seen.numbered_from(0).hides_ordinal([0], 3));
    assert!(!seen.numbered_from(1).hides_ordinal([0], 3));
}

#[test]
fn tblastx_qseq_throws_matches_seqvector_clipping() {
    // Inside the found query Bioseq (`end <= len`): never throws.
    assert!(!tblastx_qseq_throws(true, 0, 300, 300));
    assert!(!tblastx_qseq_throws(false, 0, 300, 300));
    // Plus-strand HSP starting at or past `len` (`s3` shape): empty slice -> no throw.
    assert!(!tblastx_qseq_throws(true, 300, 450, 300));
    assert!(!tblastx_qseq_throws(true, 319, 469, 300));
    // Minus-strand HSP past or straddling `len` (`s6` shape): wrapped start -> empty -> no throw.
    assert!(!tblastx_qseq_throws(false, 330, 480, 300));
    assert!(!tblastx_qseq_throws(false, 292, 478, 300));
    assert!(!tblastx_qseq_throws(false, 250, 400, 300));
    // Plus-strand HSP straddling `len`: throws iff `(len - start0) % 3 != 0`.
    assert!(!tblastx_qseq_throws(true, 249, 399, 300)); // 51 % 3 == 0
    assert!(tblastx_qseq_throws(true, 250, 400, 300)); // 50 % 3 == 2
    assert!(tblastx_qseq_throws(true, 251, 401, 300)); // 49 % 3 == 1
}

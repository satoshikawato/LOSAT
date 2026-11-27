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

/// A reported database sequence with the ids `text` and no taxid (`ParsedSubject::of_record`).
fn subject(text: &str) -> ParsedSubject {
    ParsedSubject::of_record(&believed(text), b"db sequence")
}

// Line RG item X7-1 (probes `x71/probe2`, 2026-11-27): what each report does with a sequence that
// NCBI's scope finds by a generated id, and `generated_number`.
#[test]
fn a_generated_id_is_refused_where_a_report_differs() {
    use crate::blastinput::app::Program;
    let spec = |text: &str, kind| {
        let choice = crate::blastinput::app::parse_formatting_string(text).unwrap();
        crate::report::tabular::TabularSpec::of_app(
            crate::blastinput::app::ReportFormat::from_number(choice.number).unwrap(),
            &choice,
            kind,
        )
        .expect("a field list")
    };
    fn writers(program: Program, writers: Vec<FoundWriter<'_>>) -> FoundWriters<'_> {
        FoundWriters {
            program,
            writers,
            search_reads: false,
        }
    }
    let found = subject("lcl|Query_1");
    let differs = FoundMatch::default;
    let effect = |program, writer| {
        writers(program, vec![writer]).generated_effect(&found, 1, false, &differs)
    };
    let same = || FoundMatch {
        length: true,
        aligned: true,
        unstraddled: true,
    };
    let aligned = || FoundMatch {
        length: false,
        aligned: true,
        unstraddled: true,
    };
    let sam = |target_letters, found: &dyn Fn() -> FoundMatch| {
        writers(Program::Blastn, vec![FoundWriter::Sam { target_letters }]).generated_effect(
            &subject("lcl|Query_1"),
            1,
            false,
            found,
        )
    };
    // SAM reads the found query's length (the target's hard clip, `@SQ LN:` with `SR`) and with
    // `SQ` without `SR` its letters over the HSP (`bn.tail.f17`, `bn.sub.f17SQ` of `x71/fix2/ev2`).
    assert_eq!(sam(false, &same), None);
    assert_eq!(sam(true, &same), None);
    assert_eq!(sam(false, &aligned), Some(SAM_EFFECT));
    assert_eq!(
        sam(false, &|| FoundMatch {
            length: true,
            aligned: false,
            unstraddled: true,
        }),
        None
    );
    assert_eq!(
        sam(true, &|| FoundMatch {
            length: true,
            aligned: false,
            unstraddled: true,
        }),
        Some(SAM_EFFECT)
    );
    // Outfmt 0 to 5 and 13 to 16 show the query's title and letters (`bp.qid.f0`, `.f5`).
    for program in [Program::Blastp, Program::Blastn, Program::Tblastx] {
        assert_eq!(
            effect(program, FoundWriter::Pairwise),
            Some(PAIRWISE_EFFECT)
        );
        assert_eq!(
            effect(program, FoundWriter::Anchored),
            Some(QUERY_ANCHORED_EFFECT)
        );
        assert_eq!(effect(program, FoundWriter::Structured), Some(XML_EFFECT));
        // The Seq-annots and the archive write the Seq-align id (`bp.qid.f8`).
        assert_eq!(effect(program, FoundWriter::Annots), None);
        // No taxid: the organism report is NCBI's (`tx.bp.qid.f18` without a taxonomy database).
        assert_eq!(effect(program, FoundWriter::TaxReport), None);
    }
    // The default fields: BLASTN's read nothing of the found query (`bn.qid.f6`); `pident` does
    // elsewhere (`bp.qid.f6`); BLASTX's alignment block throws (`bx.qid.f6`).
    assert_eq!(effect(Program::Blastn, FoundWriter::Table(None)), None);
    assert_eq!(
        effect(Program::Blastp, FoundWriter::Table(None)),
        Some(TABLE_EFFECT)
    );
    assert_eq!(
        effect(Program::Blastx, FoundWriter::Table(None)),
        Some(TABLE_EFFECT)
    );
    // A found query with the sequence's letters over the HSP: the default fields are NCBI's
    // (`bp.tail.f6` of `x71/fix2/ev2`), not BLASTX's.
    let table = |program, found: &dyn Fn() -> FoundMatch| {
        writers(program, vec![FoundWriter::Table(None)]).generated_effect(
            &subject("lcl|Query_1"),
            1,
            false,
            found,
        )
    };
    assert_eq!(table(Program::Blastp, &aligned), None);
    assert_eq!(table(Program::Tblastx, &aligned), None);
    assert_eq!(table(Program::Blastx, &same), Some(TABLE_EFFECT));
    let kind = crate::report::archive::options::RequestKind::Blastn;
    let ids = spec("6 qseqid sseqid sacc evalue", kind);
    let titled = spec("6 qseqid sseqid stitle", kind);
    assert_eq!(
        effect(Program::Blastn, FoundWriter::Table(Some(&ids))),
        None
    );
    assert_eq!(
        effect(Program::Blastn, FoundWriter::Table(Some(&titled))),
        Some(TABLE_EFFECT)
    );
    // The first report that differs names the effect; psiblast's rounds that make a PSSM differ
    // in every format (`it2.ps.qid`).
    let both = writers(
        Program::Blastp,
        vec![FoundWriter::Annots, FoundWriter::Structured],
    );
    assert_eq!(
        both.generated_effect(&found, 1, false, &differs),
        Some(XML_EFFECT)
    );
    let rounds = FoundWriters {
        program: Program::Blastp,
        writers: vec![FoundWriter::Annots],
        search_reads: true,
    };
    assert_eq!(
        rounds.generated_effect(&found, 1, false, &differs),
        Some(SEARCH_EFFECT)
    );
    // The PSSM reads the found query's letters over the HSPs only (`ps1.ptail.save`).
    assert_eq!(rounds.generated_effect(&found, 1, false, &aligned), None);
    assert!(!rounds.shows_shown());
    // KB-12's substitution of a `-parse_deflines` id: the tables of ids and the Seq-annots.
    assert!(writers(
        Program::Blastn,
        vec![FoundWriter::Table(Some(&ids)), FoundWriter::Annots]
    )
    .shows_shown());
    assert!(!writers(
        Program::Blastn,
        vec![FoundWriter::Sam {
            target_letters: false
        }]
    )
    .shows_shown());
}

#[test]
fn a_generated_number_is_a_local_query_id_in_any_case() {
    let number =
        |text: &str| crate::blastinput::query_batch::generated_number(&believed(text).ids[0]);
    assert_eq!(number("lcl|Query_1"), Some(1));
    assert_eq!(number("lcl|query_27"), Some(27));
    assert_eq!(number("lcl|QUERY_3"), Some(3));
    assert_eq!(number("lcl|Query_01"), None);
    assert_eq!(number("lcl|Query_"), None);
    assert_eq!(number("lcl|Query_1a"), None);
    assert_eq!(number("lcl|Subject_1"), None);
    assert_eq!(number("ref|NP_1.1|"), None);
}

#[test]
fn the_generated_ids_of_a_batch_find_a_local_id_of_its_queries() {
    use crate::blastinput::query_batch::{GeneratedQueries, QueryScope};
    let make = || QueryScope::of_generated_ids(["Query_1", "Query_2"]);
    let batch = GeneratedQueries::new(&make, None);
    assert!(batch.finds(&subject("lcl|Query_2")));
    assert!(batch.finds(&subject("lcl|query_1")));
    assert!(!batch.finds(&subject("lcl|Query_3")));
    assert!(!batch.finds(&subject("gi|5|ref|NP_1.1|")));
    // Outfmt 5 with records not read yet: their numbers (`b1.bp.q2x.f5`).
    let kept = GeneratedQueries::new(&make, Some(3..5));
    assert!(kept.finds(&subject("lcl|Query_4")));
    assert!(!kept.finds(&subject("lcl|Query_5")));
}

/// `compare_found`: a query record has the database sequence's length where it is as long, and
/// its letters over the HSPs where its letters there (upper-cased, as NCBI's reader keeps them)
/// are the sequence's; a protein's `X` never is (a `>?` gap line), and a range past the record's
/// end is not (`tests/non_utf8_names` R3_7 and R7: `q` against its copy `lcl|Query_1` of
/// `q1_sq1`; `x71/fix2/ev2` `bp.tail`, `tx.tail`).
#[test]
fn a_query_record_has_the_database_sequence_letters_over_the_hsps() {
    use crate::seqdb::{SeqDb, SeqType};
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let protein = SeqDb::open(dir.join("q1_sq1").to_str().unwrap(), SeqType::Protein).unwrap();
    let q = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAWAHVDDMPNALSA\
LSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    let whole = [0..q.len()];
    let both = FoundMatch {
        length: true,
        aligned: true,
        unstraddled: true,
    };
    assert_eq!(compare_found(q, q.len(), &protein, 0, &whole), both);
    assert_eq!(
        compare_found(&q.to_ascii_lowercase(), q.len(), &protein, 0, &whole),
        both
    );
    // A record cut at its start: neither the length nor the letters over the whole sequence; its
    // letters over a range that it shares at the same place would be, here they are shifted.
    assert_eq!(
        compare_found(&q[1..], q.len() - 1, &protein, 0, &whole),
        FoundMatch {
            length: false,
            aligned: false,
            unstraddled: false,
        }
    );
    // A record that ends early: its letters over the HSPs before its end are the sequence's.
    assert_eq!(
        compare_found(&q[..100], 100, &protein, 0, &[3..40, 50..100]),
        FoundMatch {
            length: false,
            aligned: true,
            unstraddled: true,
        }
    );
    let past = compare_found(&q[..100], 100, &protein, 0, &[50..101]);
    assert!(!past.aligned && !past.unstraddled);
    // A range wholly past the record's end reads gaps (`tx.tail.f6`) but straddles nothing.
    let beyond = compare_found(&q[..100], 100, &protein, 0, &[110..120]);
    assert!(!beyond.aligned && beyond.unstraddled);
    assert!(!compare_found(q, q.len(), &protein, 1, &whole).aligned);
    // A record whose letters the caller dropped (BLASTN's `x_query_drop`) has the sequence's
    // length, and no letters over any HSP.
    assert_eq!(
        compare_found(&[], q.len(), &protein, 0, &whole),
        FoundMatch {
            length: true,
            aligned: false,
            unstraddled: true,
        }
    );
    let mut gap = q.to_vec();
    gap[10] = b'X';
    assert!(!compare_found(&gap, gap.len(), &protein, 0, &[5..20]).aligned);
    assert!(compare_found(&gap, gap.len(), &protein, 0, &[20..60]).aligned);
    let nucleotide =
        SeqDb::open(dir.join("of1/bnqid").to_str().unwrap(), SeqType::Nucleotide).unwrap();
    let copy = std::fs::read(dir.join("of1/bnqid.fa")).unwrap();
    let first: Vec<u8> = copy
        .split(|&byte| byte == b'>')
        .nth(1)
        .unwrap()
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flatten()
        .copied()
        .collect();
    assert_eq!(
        compare_found(&first, first.len(), &nucleotide, 0, &[0..first.len()]),
        both
    );
    assert_eq!(
        compare_found(
            &first[..first.len() - 1],
            first.len() - 1,
            &nucleotide,
            0,
            &[0..10]
        ),
        FoundMatch {
            length: false,
            aligned: true,
            unstraddled: true,
        }
    );
}

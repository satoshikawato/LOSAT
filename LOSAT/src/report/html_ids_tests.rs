// Expectations are NCBI BLAST+ 2.17.0's bytes of the probes in
// remaining-gaps-20261112/T/t4/out/ (bp.dm.f0.h, bn.dm.f0.h, bp.di.f0.rc_localid.h,
// bp.s.f0.rc_localid.h, bp.dp.f0.h).

use super::*;
use crate::report::formatter_registry::url_encode;

fn ids(text: &str) -> Vec<ParsedSeqId> {
    crate::seqid::parse::parse_ids(text.as_bytes(), 0, &mut |_| {})
        .unwrap()
        .into_iter()
        .flatten()
        .collect()
}

fn text(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).unwrap()
}

fn accession(ids: &[ParsedSeqId]) -> Vec<u8> {
    let mut out = Vec::new();
    description_accession(find_best(ids, text_score).unwrap(), &mut out);
    out
}

fn url(ids: &[ParsedSeqId], rank: usize, nucleotide: bool, place: LinkPlace) -> String {
    let accession = match place {
        LinkPlace::Description => accession(ids),
        LinkPlace::Heading => {
            let mut out = Vec::new();
            heading_label(find_best(ids, text_score).unwrap(), &mut out);
            out
        }
    };
    let info = SeqUrlInfo {
        accession: &accession,
        rank,
        nucleotide,
        place,
        taxid: 0,
    };
    let mut err = Vec::new();
    let url = id_url(ids, &info, &links(None, b""), &mut err).unwrap();
    assert!(err.is_empty());
    text(url)
}

/// The links of a run with the registry file `file` (`None`: none) and the `-db` text `db`.
fn links(file: Option<&[u8]>, db: &[u8]) -> HtmlLinks {
    HtmlLinks {
        registry: InitConfig::of(
            "blastp",
            file.map(|bytes| (".ncbirc".into(), bytes.to_vec())),
            false,
        )
        .unwrap(),
        db_name: db.to_vec(),
        nucleotide: false,
    }
}

#[test]
fn entrez_links_of_text_ids_pdb_and_patents() {
    assert_eq!(
        url(&ids("pir||A12345"), 1, false, LinkPlace::Description),
        "<a title=\"Show report for A12345\" href=\"https://www.ncbi.nlm.nih.gov/protein/A12345?report=genbank&log$=prottop&blast_rank=1&RID=\" >"
    );
    assert_eq!(
        url(&ids("pir||A12345"), 1, false, LinkPlace::Heading),
        "<a title=\"Show report for A12345\" href=\"https://www.ncbi.nlm.nih.gov/protein/A12345?report=genbank&log$=protalign&blast_rank=1&RID=0\" >"
    );
    assert_eq!(
        url(&ids("pat|US|5000000|1"), 2, false, LinkPlace::Description),
        "<a title=\"Show report for US5000000_1\" href=\"https://www.ncbi.nlm.nih.gov/protein/US5000000_1?report=genbank&log$=prottop&blast_rank=2&RID=\" >"
    );
    assert_eq!(
        url(&ids("pdb|1ABC|A"), 8, false, LinkPlace::Description),
        "<a title=\"Show report for 1ABC_A\" href=\"https://www.ncbi.nlm.nih.gov/protein/1ABC_A?report=genbank&log$=prottop&blast_rank=8&RID=\" >"
    );
    assert_eq!(
        url(&ids("emb|BN000001.1|"), 2, true, LinkPlace::Description),
        "<a title=\"Show report for BN000001.1\" href=\"https://www.ncbi.nlm.nih.gov/nucleotide/BN000001.1?report=genbank&log$=nucltop&blast_rank=2&RID=\" >"
    );
    assert_eq!(
        url(&ids("emb|BN000001.1|"), 2, true, LinkPlace::Heading),
        "<a title=\"Show report for BN000001.1\" href=\"https://www.ncbi.nlm.nih.gov/nucleotide/BN000001.1?report=genbank&log$=nuclalign&blast_rank=2&RID=0\" >"
    );
    // A GI with a text id: the link is the text id's; the score label is the GI.
    let gi_gb = ids("gi|1000|gb|AAA00001.1|");
    assert_eq!(
        url(&gi_gb, 7, false, LinkPlace::Description),
        "<a title=\"Show report for AAA00001.1\" href=\"https://www.ncbi.nlm.nih.gov/protein/AAA00001.1?report=genbank&log$=prottop&blast_rank=7&RID=\" >"
    );
    let mut score = Vec::new();
    score_link(1000, &accession(&gi_gb), &mut score);
    assert_eq!(text(score), "<a href=#1000>");
}

#[test]
fn trace_links_and_ids_without_a_link() {
    let trace = ids("gnl|TI|12345");
    assert_eq!(
        url(&trace, 6, false, LinkPlace::Description),
        "<a title=\"Show report for 12345\" href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=12345&RID=\">"
    );
    assert_eq!(
        url(&trace, 6, false, LinkPlace::Heading),
        "<a title=\"Show report for 12345\" href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=12345&RID=0\">"
    );
    assert_eq!(text(accession(&trace)), "TI:12345");
    for none in [
        "gi|999",
        "lcl|loc_x",
        "gnl|testdb|tag_s8",
        "gnl|BL_ORD_ID|4",
        "lcl|Subject_1",
    ] {
        assert_eq!(
            url(&ids(none), 1, false, LinkPlace::Description),
            "",
            "{none}"
        );
    }
    let mut score = Vec::new();
    score_link(0, &accession(&ids("gnl|BL_ORD_ID|4")), &mut score);
    assert_eq!(text(score), "<a href=#BL_ORD_ID:4>");
}

#[test]
fn local_id_template_tags() {
    // `T/t4/rc/localid/.ncbirc` (probes `bp.di.f0.rc_localid.h`, `bp.s.f0.rc_localid.h`).
    let file: &[u8] = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    let local = ids("lcl|s3_local_identifier");
    let accession = accession(&local);
    let info = SeqUrlInfo {
        accession: &accession,
        rank: 3,
        nucleotide: false,
        place: LinkPlace::Description,
        taxid: 9606,
    };
    let mut err = Vec::new();
    let db = links(Some(file), b"../../../t2/db/ps_ids");
    assert_eq!(
        text(id_url(&local, &info, &db, &mut err).unwrap()),
        "http://local.example/s3_local_identifier?db=..%2F..%2F..%2Ft2%2Fdb%2Fps_ids&taxid=9606&p=https:&t=s3_local_identifier"
    );
    let subject = ids("lcl|Subject_1");
    let info = SeqUrlInfo { taxid: 0, ..info };
    assert_eq!(
        text(id_url(&subject, &info, &links(Some(file), b""), &mut err).unwrap()),
        "http://local.example/Subject_1?db=&taxid=0&p=https:&t=Subject_1"
    );
    // Without the template a local id has no link.
    assert_eq!(
        text(id_url(&subject, &info, &links(None, b""), &mut err).unwrap()),
        ""
    );
    assert!(err.is_empty());
}

// D5: the description rows' links and score links of the ids that the formatter makes for a
// subject without parsed Seq-ids (NCBI BLAST+ 2.17.0 `bp.dp.f0.h`, `bp.s.f0.h`,
// `bp.s.f0.rc_localid.h`: `gnl|BL_ORD_ID|<oid>`, `lcl|Subject_<index + 1>`) and of a record with
// parsed ids (`bp.dm.f0.h` row 1).
#[test]
fn description_links_of_the_row_ids() {
    let file: &[u8] = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    let row = |report: &HtmlReport, parsed: Option<&ParsedSubject>, s_idx: u32, rank: usize| {
        let (link, score) = report.description_links(parsed, s_idx, rank).unwrap();
        (text(link), text(score))
    };
    let oids = [0, 2, 10].map(crate::seqdb::Oid::of_unrestricted);
    let database = HtmlReport::new(&Arc::new(links(None, b"../../db/ps_plain")), oids);
    assert_eq!(
        row(&database, None, 2, 3),
        (String::new(), "<a href=#BL_ORD_ID:10>".into())
    );
    let records = HtmlReport::new(&Arc::new(links(None, b"")), []);
    assert_eq!(
        row(&records, None, 0, 1),
        (String::new(), "<a href=#Subject_1>".into())
    );
    let records = HtmlReport::new(&Arc::new(links(Some(file), b"")), []);
    assert_eq!(
        row(&records, None, 2, 2),
        (
            "http://local.example/Subject_3?db=&taxid=0&p=https:&t=Subject_3".into(),
            "<a href=#Subject_3>".into()
        )
    );
    let believed = crate::blastinput::fasta_reader::believed::BelievedIds::new(
        ids("pir||A12345"),
        crate::report::archive::bioseq::SeqLocStrand::Unknown,
    );
    let record = ParsedSubject::of_record(&believed, b"pir subject");
    assert_eq!(
        row(&records, Some(&record), 0, 1),
        (
            "<a title=\"Show report for A12345\" href=\"https://www.ncbi.nlm.nih.gov/protein/A12345?report=genbank&log$=prottop&blast_rank=1&RID=\" >".into(),
            "<a href=#A12345>".into()
        )
    );
}

#[test]
fn anchors_with_and_without_version() {
    let id = &ids("ref|XP_000901.1|")[0];
    let mut out = Vec::new();
    anchor(0, id, true, &mut out);
    anchor(0, id, false, &mut out);
    anchor(1000, id, true, &mut out);
    anchor(0, &ids("gnl|BL_ORD_ID|0")[0], true, &mut out);
    assert_eq!(
        text(out),
        "<a name=XP_000901.1></a><a name=XP_000901></a><a name=1000></a><a name=BL_ORD_ID:0></a>"
    );
}

#[test]
fn url_encoding_and_templates() {
    // The registry's helpers (`formatter_registry`), the one owner of both.
    assert_eq!(
        text(url_encode(b"a b/c.d-e_f~g%h!$'()*,")),
        "a+b%2Fc.d-e_f%7Eg%25h!$'()*,"
    );
    // A value is not searched again; every occurrence is replaced.
    assert_eq!(
        text(map_template(b"<@a@>x<@a@>", "a", b"<@a@>")),
        "<@a@>x<@a@>"
    );
    assert_eq!(text(map_template(b"<@a@>x<@a@><@a", "a", b"1")), "1x1<@a");
}

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
    let mut posts = Vec::new();
    let url = id_url(ids, &info, &links(None, b""), &mut posts).unwrap();
    // Each id URL calls `InitConfig` (`HtmlPost::InitConfig`).
    assert_eq!(posts, [HtmlPost::InitConfig]);
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
        ncbirc: Default::default(),
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
    let mut posts = Vec::new();
    let db = links(Some(file), b"../../../t2/db/ps_ids");
    assert_eq!(
        text(id_url(&local, &info, &db, &mut posts).unwrap()),
        "http://local.example/s3_local_identifier?db=..%2F..%2F..%2Ft2%2Fdb%2Fps_ids&taxid=9606&p=https:&t=s3_local_identifier"
    );
    let subject = ids("lcl|Subject_1");
    let info = SeqUrlInfo { taxid: 0, ..info };
    assert_eq!(
        text(id_url(&subject, &info, &links(Some(file), b""), &mut posts).unwrap()),
        "http://local.example/Subject_1?db=&taxid=0&p=https:&t=Subject_1"
    );
    // Without the template a local id has no link.
    assert_eq!(
        text(id_url(&subject, &info, &links(None, b""), &mut posts).unwrap()),
        ""
    );
    assert_eq!(posts, [HtmlPost::InitConfig; 3]);
}

// D5: the description rows' links and score links of the ids that the formatter makes for a
// subject without parsed Seq-ids (NCBI BLAST+ 2.17.0 `bp.dp.f0.h`, `bp.s.f0.h`,
// `bp.s.f0.rc_localid.h`: `gnl|BL_ORD_ID|<oid>`, `lcl|Subject_<index + 1>`) and of a record with
// parsed ids (`bp.dm.f0.h` row 1).
#[test]
fn description_links_of_the_row_ids() {
    let file: &[u8] = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    let row = |report: &HtmlReport, parsed: Option<&ParsedSubject>, s_idx: u32, rank: usize| {
        let (link, score) = report
            .description_links(parsed, s_idx, rank, false, &mut Vec::new())
            .unwrap();
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

// F3 and O1 of line RG T4a2 review a: the row of NCBI's `catch` (the record's title cannot be
// decoded) names the Seq-align id (`gi|777`, not the display id `XP_000009.1`), links it and has
// no score link (NCBI BLAST+ 2.17.0, `blastp -subject us4.faa -parse_deflines -html`, probe
// `T/t4/fix1/out/unk4h`).
#[test]
fn the_catch_row_links_the_seq_align_id_without_a_score_link() {
    let believed = crate::blastinput::fasta_reader::believed::BelievedIds::new(
        ids("gi|777|ref|XP_000009.1|"),
        crate::report::archive::bioseq::SeqLocStrand::Unknown,
    );
    let record = ParsedSubject::of_record(&believed, b"bad \x81 title one");
    let (start, title, caught) = record.description_line(true);
    assert!(caught);
    assert_eq!(text(start), "gi|777 ");
    assert_eq!(text(title.text), "Unknown");
    let records = HtmlReport::new(&Arc::new(links(None, b"")), []);
    let mut posts = Vec::new();
    let (link, score) = records
        .description_links(Some(&record), 0, 1, true, &mut posts)
        .unwrap();
    assert_eq!(
        text(link),
        "<a title=\"Show report for 777\" href=\"https://www.ncbi.nlm.nih.gov/protein/777?report=genbank&log$=prottop&blast_rank=1&RID=\" >"
    );
    assert!(score.is_empty());
    assert_eq!(posts, [HtmlPost::InitConfig]);
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

/// The `-html` heading of a database sequence with the definition lines `lines` (ids, title) at
/// rank `rank` of a report of `report`.
fn heading(report: &HtmlReport, lines: &[(&str, &str)], rank: usize) -> String {
    let deflines = lines
        .iter()
        .map(|(id, title)| crate::seqid::BlastDefline {
            title: Some(title.as_bytes().to_vec()),
            ids: ids(id),
            ..Default::default()
        })
        .collect();
    let subject = ParsedSubject::new(crate::seqid::DeflineSet {
        deflines,
        utf8: false,
    })
    .unwrap();
    let posts = std::cell::RefCell::new(Vec::new());
    let html = HtmlHeading {
        report,
        s_idx: 0,
        rank,
        posts: &posts,
    };
    text(subject.heading(true, false, Some(html)).unwrap().unwrap())
}

fn entrez(accession: &str, rank: usize) -> String {
    format!("<a title=\"Show report for {accession}\" href=\"https://www.ncbi.nlm.nih.gov/protein/{accession}?report=genbank&log$=protalign&blast_rank={rank}&RID=0\" >")
}

// The alignment headings (D6, D7) of `bp.dm.f0.h` (lines 55, 117-127, 172, 187, 217, 232),
// `bp.dp.f0.h` (line 56), `bp.s.f0.h` and `bp.s.f0.rc_localid.h` (first headings).
#[test]
fn heading_lines_with_links_anchors_and_hidden_titles() {
    let report = HtmlReport::new(
        &Arc::new(links(None, b"../../db/pmix_ids")),
        [crate::seqdb::Oid::of_unrestricted(0)],
    );
    assert_eq!(
        heading(&report, &[("pir||A12345", "pir subject")], 1),
        format!(
            ">{}pir||A12345</a><a name=A12345></a> pir subject\n",
            entrez("A12345", 1)
        )
    );
    assert_eq!(
        heading(&report, &[("gnl|TI|12345", "trace subject")], 6),
        ">\
<a title=\"Show report for 12345\" href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=12345&RID=0\">\
TI:12345</a><a name=TI:12345></a> trace subject\n"
    );
    // The GI names the anchor; the title is HTML-encoded.
    assert_eq!(
        heading(
            &report,
            &[("gi|1000|gb|AAA00001.1|", "gi and gb <&> \"q\"")],
            7
        ),
        format!(
            ">{}AAA00001.1</a><a name=1000></a> gi and gb &lt;&amp;&gt; &quot;q&quot;\n",
            entrez("AAA00001.1", 7)
        )
    );
    // No link: the id, then the anchor.
    assert_eq!(
        heading(&report, &[("gi|999", "gi only subject")], 9),
        ">gi|999<a name=999></a> gi only subject\n"
    );
    assert_eq!(
        heading(&report, &[("lcl|loc_x", "local & id")], 10),
        ">loc_x<a name=loc_x></a> local &amp; id\n"
    );
    // Nine definition lines: the fourth to the ninth in the box of the third's label.
    let nine: Vec<(String, String)> = (1..=9)
        .map(|n| {
            (
                format!("ref|XP_00090{n}.1|"),
                format!("nine titles {n} <&>"),
            )
        })
        .collect();
    let lines: Vec<(&str, &str)> = nine
        .iter()
        .map(|(id, t)| (id.as_str(), t.as_str()))
        .collect();
    let mut expected = String::new();
    for n in 1..=9 {
        expected.push(if n == 1 { '>' } else { ' ' });
        if n == 4 {
            expected.push_str("<a href=\"#\" title=\"Other sequence titles\"  onmouseover=\"showInfo(this)\" class=\"resArrowLinkW mdl hiding\" id=\"XP_000903.1\">6 more sequence titles</a>\n <div id=\"info_XP_000903.1\" class=\"helpbox mdlbox hidden\">");
        }
        let accession = format!("XP_00090{n}.1");
        expected.push_str(&format!(
            "{}{accession}</a><a name={accession}></a> nine titles {n} &lt;&amp;&gt;\n",
            entrez(&accession, 3)
        ));
    }
    expected.push_str("</div>");
    assert_eq!(heading(&report, &lines, 3), expected);
    // Eight lines show every line.
    assert!(!heading(&report, &lines[..8], 3).contains("more sequence titles"));
    // A sequence without parsed Seq-ids, and a `-subject` record without them (with and without
    // the registry's `[LOCAL_ID] TOOL_URL`).
    let start = |report: &HtmlReport, s_idx: u32| {
        let mut out = Vec::new();
        let posts = std::cell::RefCell::new(Vec::new());
        let html = HtmlHeading {
            report,
            s_idx,
            rank: 1,
            posts: &posts,
        };
        html.unparsed_start(&mut out).unwrap();
        // The heading's one id link calls `InitConfig` once.
        assert_eq!(html.take_posts(), [HtmlPost::InitConfig]);
        text(out)
    };
    assert_eq!(start(&report, 0), "><a name=BL_ORD_ID:0></a> ");
    let records = HtmlReport::new(&Arc::new(links(None, b"")), []);
    assert_eq!(start(&records, 2), "><a name=Subject_3></a> ");
    let file: &[u8] = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    let records = HtmlReport::new(&Arc::new(links(Some(file), b"")), []);
    assert_eq!(
        start(&records, 0),
        ">http://local.example/Subject_1?db=&taxid=0&p=https:&t=Subject_1</a><a name=Subject_1></a> "
    );
}

// `cfg.fmtbad.nd0.*.h`: a malformed registry file throws at a heading's first id link, which
// NCBI catches for each alignment: the Seq-id text without the version.
#[test]
fn a_registry_that_throws_skips_the_heading() {
    let bad = HtmlLinks {
        registry: InitConfig::of(
            "blastp",
            Some((".ncbirc".into(), b"this line has no equals sign\n".to_vec())),
            false,
        )
        .unwrap(),
        db_name: Vec::new(),
        nucleotide: false,
        ncbirc: Default::default(),
    };
    let report = HtmlReport::new(&Arc::new(bad), []);
    let posts = std::cell::RefCell::new(Vec::new());
    let html = HtmlHeading {
        report: &report,
        s_idx: 0,
        rank: 1,
        posts: &posts,
    };
    let subject = ParsedSubject::new(crate::seqid::DeflineSet {
        deflines: vec![crate::seqid::BlastDefline {
            title: Some(b"t".to_vec()),
            ids: ids("ref|XP_000101.1|"),
            ..Default::default()
        }],
        utf8: false,
    })
    .unwrap();
    assert_eq!(
        subject
            .heading(true, false, Some(html))
            .unwrap()
            .unwrap_err(),
        b"XP_000101"
    );
    let error = html.unparsed_start(&mut Vec::new()).unwrap_err();
    assert_eq!(
        HtmlHeading::skips(error, || html.unparsed_skipped_id()).unwrap(),
        b"Subject_1"
    );
    let other = io::Error::other("write failed");
    assert!(HtmlHeading::skips(other, Vec::new).is_err());
}

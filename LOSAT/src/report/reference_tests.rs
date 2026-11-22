use super::*;

/// The references of `kPublication` hold `keyword` and its PubMed URL holds `pubmed_id`
/// (`s_DoTestReference`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:61-71
fn do_test_reference(keyword: &str, pubmed_id: &str, publication: Publication) {
    assert!(publication.string().contains(keyword), "{publication:?}");
    assert!(
        publication.pubmed_url().contains(pubmed_id),
        "{publication:?}"
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:73-79 testGappedBlastReference
#[test]
fn the_gapped_blast_reference_has_its_title_and_pubmed_id() {
    do_test_reference(
        "Gapped BLAST and PSI-BLAST",
        "9254694",
        Publication::GappedBlast,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:81-87 testMegaBlastReference
#[test]
fn the_megablast_reference_has_its_title_and_pubmed_id() {
    do_test_reference(
        "A greedy algorithm for aligning DNA",
        "10890397",
        Publication::MegaBlast,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:89-95 testCompositionBasedStatisticsReference
#[test]
fn the_composition_based_statistics_reference_has_its_title_and_pubmed_id() {
    do_test_reference(
        "with composition-based statistics",
        "11452024",
        Publication::CompBasedStats,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:97-103 testPhiBlastReference
#[test]
fn the_phi_blast_reference_has_its_title_and_pubmed_id() {
    do_test_reference("using patterns as seeds", "9705509", Publication::PhiBlast);
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:105-113 testInvalidReference
#[test]
fn the_sentinel_publication_has_no_text_and_no_url() {
    assert!(Publication::MaxPublications.string().is_empty());
    assert!(Publication::MaxPublications.pubmed_url().is_empty());
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/version_reference_unit_test.cpp:115-141 printAllReferencesHtml
#[test]
fn every_publication_has_one_www_url() {
    let publications = [
        Publication::GappedBlast,
        Publication::PhiBlast,
        Publication::MegaBlast,
        Publication::CompBasedStats,
        Publication::CompAdjustedMatrices,
        Publication::IndexedMegablast,
        Publication::DeltaBlast,
    ];
    // The publications before the sentinel, in the enum's order (the index cast to `EPublication`).
    for (index, publication) in publications.iter().enumerate() {
        assert_eq!(*publication as usize, index);
    }
    assert_eq!(publications.len(), Publication::MaxPublications as usize);
    let mut text = String::new();
    for publication in publications {
        text.push_str(publication.string());
        text.push('\n');
        text.push_str("<a href=");
        text.push_str(publication.pubmed_url());
        text.push_str(">Get abstract</a><br>\n");
    }
    // `text.find` of `kToken` from `i`, one past each match: overlapping matches count.
    let occurrences = (0..text.len())
        .filter(|&at| text[at..].starts_with("www"))
        .count();
    assert_eq!(occurrences, publications.len());
}

#[test]
fn the_html_free_text_writes_the_entity_as_a() {
    assert_eq!(
        Publication::GappedBlast.html_free_string(),
        "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    );
    // DELTA-BLAST's text has no entity; the sentinel's is empty.
    assert_eq!(
        Publication::DeltaBlast.html_free_parts(),
        [Publication::DeltaBlast.string(), "", ""]
    );
    assert_eq!(Publication::MaxPublications.html_free_string(), "");
    for index in 0..TABLE_LEN {
        assert!(!REFERENCES[index].replacen("&auml;", "a", 1).contains('&'));
    }
}

fn prolog_references(html: bool, protocol: &[u8], program: &str, megablast: bool) -> String {
    let mut out = Vec::new();
    write_prolog_references(&mut out, html, protocol, program, megablast, false).unwrap();
    String::from_utf8(out).unwrap()
}

/// The plain references of each program's prolog are the text that its intro wrote before the
/// table (frozen NCBI outfmt 0 of `outfmt0_manifest.tsv`, e.g. `e2e.blastp.many`).
#[test]
fn the_plain_references_are_ncbis_lines_for_each_program() {
    const GAPPED: &str = "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\n\
Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.\n\
Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of\n\
protein database search programs\", Nucleic Acids Res. 25:3389-3402.\n\n";
    let blastp = format!(
        "{GAPPED}\n\
Reference for composition-based statistics: Alejandro A. Schaffer,\n\
L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri\n\
I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001),\n\
\"Improving the accuracy of PSI-BLAST protein database searches with\n\
composition-based statistics and other refinements\", Nucleic Acids\n\
Res. 29:2994-3005.\n\n\n\n"
    );
    assert_eq!(
        prolog_references(false, DEFAULT_PROTOCOL, "blastp", false),
        blastp
    );
    let psiblast = format!(
        "{GAPPED}\n\
Reference for compositional score matrix adjustment: Stephen F.\n\
Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala,\n\
Aleksandr Morgulis, Alejandro A. Schaffer, and Yi-Kuo Yu (2005)\n\
\"Protein database searches using compositionally adjusted\n\
substitution matrices\", FEBS J. 272:5101-5109.\n\n\n\
Reference for composition-based statistics starting in round 2:\n\
Alejandro A. Schaffer, L. Aravind, Thomas L. Madden, Sergei\n\
Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and\n\
Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST\n\
protein database searches with composition-based statistics and\n\
other refinements\", Nucleic Acids Res. 29:2994-3005.\n\n\n\n"
    );
    assert_eq!(
        prolog_references(false, DEFAULT_PROTOCOL, "psiblast", false),
        psiblast
    );
    let gapped_only = format!("{GAPPED}\n\n");
    for program in [
        "rpsblast",
        "BLASTN",
        "TBLASTX",
        "TBLASTN",
        "PSITBLASTN",
        "BLASTX",
        "RPSTBLASTN",
    ] {
        assert_eq!(
            prolog_references(false, DEFAULT_PROTOCOL, program, false),
            gapped_only,
            "{program}"
        );
    }
    assert_eq!(
        prolog_references(false, DEFAULT_PROTOCOL, "blastn", true),
        "Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb\n\
Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J\n\
Comput Biol 2000; 7(1-2):203-14.\n\n\n\n"
    );
}

/// The `-html` references are NCBI BLAST+ 2.17.0's (line RG T4 probes, `T/t4/out`: `bp.s.f0.h`
/// lines 9-29 for blastp, `bn.s.f0.h` for megablast, `bn.s.f0.blastn.h` for task blastn).
#[test]
fn the_html_references_are_ncbis_lines() {
    assert_eq!(
        prolog_references(true, DEFAULT_PROTOCOL, "blastp", false),
        "<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=9254694&dopt=Citation\">Reference</a>:</b>\n\
Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer,\n\
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n\
\"Gapped BLAST and PSI-BLAST: a new generation of protein database\n\
search programs\", Nucleic Acids Res. 25:3389-3402.\n\
\n\
\n\
<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=11452024&dopt=Citation\">Reference for\n\
composition-based statistics</a>:</b>\n\
Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei\n\
Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and\n\
Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST\n\
protein database searches with composition-based statistics and\n\
other refinements\", Nucleic Acids Res. 29:2994-3005.\n\
\n\
\n\
\n"
    );
    assert_eq!(
        prolog_references(true, DEFAULT_PROTOCOL, "blastn", true),
        "<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=10890397&dopt=Citation\">Reference</a>:</b>\n\
Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000),\n\
\"A greedy algorithm for aligning DNA sequences\", J Comput Biol 2000;\n\
7(1-2):203-14.\n\
\n\
\n\
\n"
    );
    assert_eq!(
        prolog_references(true, DEFAULT_PROTOCOL, "BLASTN", false),
        "<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=9254694&dopt=Citation\">Reference</a>:</b>\n\
Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer,\n\
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n\
\"Gapped BLAST and PSI-BLAST: a new generation of protein database\n\
search programs\", Nucleic Acids Res. 25:3389-3402.\n\
\n\
\n\
\n"
    );
}

/// `[BLASTFMTUTIL] PROTOCOL` moves the wrap of the URL (probes `bp.di.f0.rc_p_http.h`,
/// `rc_p_empty.h`: `http:`, and an empty value, which `HasEntry` counts as set here).
#[test]
fn the_protocol_starts_the_pubmed_url() {
    let first_lines = |protocol: &[u8]| {
        let mut out = Vec::new();
        write_reference(&mut out, true, protocol, Publication::GappedBlast, false).unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .take(3)
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        first_lines(b"http:"),
        "<b><a\n\
href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Re\n\
trieve&list_uids=9254694&dopt=Citation\">Reference</a>:</b>"
    );
    assert_eq!(
        first_lines(b""),
        "<b><a\n\
href=\"//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retriev\n\
e&list_uids=9254694&dopt=Citation\">Reference</a>:</b>"
    );
}

#[test]
fn the_version_line_is_the_program_in_upper_case() {
    let line = |html: bool, program: &str| {
        let mut out = Vec::new();
        write_version_line(&mut out, html, program, "2.17.0+").unwrap();
        String::from_utf8(out).unwrap()
    };
    assert_eq!(line(false, "blastp"), "BLASTP 2.17.0+\n");
    assert_eq!(line(false, "TBLASTX"), "TBLASTX 2.17.0+\n");
    // NCBI BLAST+ 2.17.0 (probe `bp.s.f0.h`, line 6).
    assert_eq!(line(true, "blastp"), "<b>BLASTP 2.17.0+</b>\n");
}

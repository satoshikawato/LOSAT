use super::*;
use crate::report::archive::options::RequestKind;
use crate::seqdb::{SeqDb, SeqType};
use crate::seqid::model::ObjectId;
use crate::seqid::model::{Dbtag, PdbSeqId};

/// The subjects of `tests/db_search/ids/prot_v4` (`make_ids.py`), by the FASTA form of the
/// first id of their first definition line.
fn subject(first: &str) -> ParsedSubject {
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/prot_v4");
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    for ordinal in 0..db.num_seqs() {
        let set = db.defline_set(ordinal).unwrap();
        let mut fasta = Vec::new();
        label::fasta(&set.deflines[0].ids[0], &mut fasta);
        if fasta == first.as_bytes() {
            return ParsedSubject::new(set).unwrap();
        }
    }
    panic!("no sequence {first}");
}

fn text(bytes: &[u8]) -> &str {
    std::str::from_utf8(bytes).unwrap()
}

// Approved exception 2 of PD-LOSAT-NCBI-DEFECTS: BLASTN, TBLASTN and TBLASTX stop an outfmt 0
// title made only of punctuation at the end of the string (NCBI's `x_CleanAndCompress` reads
// past it and crashes); a nucleotide subject with parsed Seq-ids gets the title of the same
// stored title without them, and so does the tabular title word. A protein subject is refused.
#[test]
fn a_nucleotide_title_of_punctuation_stops_at_its_end() {
    use crate::algorithm::common::database::check_shown_parsed_subject;
    let subject = |id: ParsedSeqId| {
        let defline = BlastDefline {
            title: Some(b", ,".to_vec()),
            ids: vec![id],
            ..BlastDefline::default()
        };
        ParsedSubject::new(DeflineSet {
            deflines: vec![defline],
            utf8: false,
        })
        .unwrap()
    };
    let stored = crate::report::defline::generate_defline(b", ,", false, true).unwrap();
    assert!(stored.reads_past_end);
    let accession = subject(ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
        name: None,
        accession: Some(b"AF000050".to_vec()),
        release: None,
        version: Some(1),
    }));
    let (start, title) = accession.description(false);
    assert_eq!(text(&start), "AF000050.1 ");
    assert_eq!(title.text, stored.text);
    check_shown_parsed_subject(
        &accession,
        false,
        false,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTN",
    )
    .expect("a nucleotide title of punctuation is shown");
    assert!(check_shown_parsed_subject(
        &accession,
        true,
        false,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP"
    )
    .is_err());
    let local = subject(ParsedSeqId::Local(ObjectId::Str(b"Subject_51".to_vec())));
    assert!(local.tabular_ids(true).is_err());
    assert!(local.tabular_saccver(true).is_err());
    assert!(local.check_tabular_ids(true).is_err());
    assert!(local.check_tabular_ids(false).is_ok());
    let [id, _, _] = local.tabular_ids(false).unwrap().unwrap();
    let word = stored.text.split(|&byte| byte == b' ').next().unwrap();
    assert_eq!(id, word);
}

// NCBI BLAST+ 2.17.0 `blastn -db tests/db_search/ids/nucl_v5 -outfmt 17` names the subjects
// so (`LOSAT/tests/fixtures/db_search/db.id.blastn.17.v5.out`, ordinal order): the Seq-align
// id `gi|654321` and `gi|2020` (a GI in the second definition line) give the accession of
// the Bioseq's definition line, `AB000009` has no version, and `lcl|Subject_31` no Textseq-id.
#[test]
fn sam_names_the_subject_by_an_accession_of_its_bioseq() {
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/nucl_v5");
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Nucleotide).unwrap();
    let names: Vec<String> = (0..db.num_seqs())
        .map(|ordinal| {
            let subject = ParsedSubject::new(db.defline_set(ordinal).unwrap()).unwrap();
            let mut name = Vec::new();
            crate::report::seqid::SeqId::Db(subject.sam_id()).sam_name(&mut name);
            String::from_utf8(name).unwrap()
        })
        .collect();
    assert_eq!(
        names,
        [
            "NM_000001.1",
            "AF000002.3",
            "AJ000003.1",
            "AB000004.1",
            "NC_000005.11",
            "mynuc6",
            "mydb:nuc7",
            "NM_000008.2",
            "AB000009",
            "mynuc_10",
            "BK000012.1",
            "AF000013.1",
            "US1234567_14",
            "AF000015.1",
            "Subject_30",
            "Subject_31",
            "BK000018.1",
            "NC_000019.1",
            "AF000020.1",
            "AF000021.1",
            "AF000022.1",
            "NR_000023.1",
            "ti:123",
            "BN000025.1",
        ]
    );
}

fn description(first: &str) -> String {
    let (mut line, title) = subject(first).description(true);
    line.extend_from_slice(&title.text);
    String::from_utf8(line).unwrap()
}

// The strings are NCBI BLAST+ 2.17.0's `blastp -db` (outfmt 0 and 6) for these records
// (`LOSAT/tests/fixtures/db_search/db.id.blastp.0.v4.out`, `db.id.blastp.6c.v4.out`).
#[test]
fn titles_that_ncbi_makes_from_the_seq_ids() {
    assert_eq!(
        description("tpg|DAA10001.1|"),
        "DAA10001.1 tpa titled protein"
    );
    assert_eq!(
        description("tpg|DAA10002.1|"),
        "DAA10002.1 unnamed protein product"
    );
    assert_eq!(
        text(
            &subject("tpg|DAA10002.1|")
                .heading(true, false, None)
                .unwrap()
                .unwrap()
        ),
        ">DAA10002.1 TPA: unnamed protein product\n"
    );
    assert_eq!(description("pdb|2XYZ|B"), "2XYZ_B Chain B");
    assert_eq!(description("pdb|2XYY| "), "2XYY unnamed protein product");
    assert_eq!(
        description("pat|US|7654321|3"),
        "US7654321_3 Sequence 3 from Patent US 7654321"
    );
    assert_eq!(description("ref|NC_100013.1|"), "NC_100013.1 ");
    assert_eq!(
        text(
            &subject("ref|NC_100013.1|")
                .heading(true, false, None)
                .unwrap()
                .unwrap()
        ),
        ">NC_100013.1 \n"
    );
    assert_eq!(
        description("ref|WP_100014.1|"),
        "WP_100014.1 unnamed protein product"
    );
    assert_eq!(description("lcl|Subject_27"), "subject seven title");
    assert_eq!(
        text(
            &subject("lcl|Subject_30")
                .heading(true, false, None)
                .unwrap()
                .unwrap()
        ),
        ">Subject_30 unnamed protein product\n"
    );
    assert_eq!(description("pdb|5ABC|D"), "5ABC_D Chain D");
    assert_eq!(
        text(
            &subject("pdb|5ABC|D")
                .heading(true, false, None)
                .unwrap()
                .unwrap()
        ),
        ">5ABC_D Chain D\n XP_100015.1 second\n"
    );
    assert_eq!(
        text(
            &subject("lcl|Subject_36")
                .heading(true, false, None)
                .unwrap()
                .unwrap()
        ),
        ">Subject_36 unnamed protein product\n DAA10016.1 tpa second\n"
    );
    // A title from the source: the taxonomy database can name its organism.
    assert!(!subject("ref|NP_100006.1|").title_from_taxonomy());
    assert!(!subject("pdb|2XYZ|B").title_from_taxonomy());
    assert!(subject("gi|1818").title_from_taxonomy());
}

#[test]
fn tabular_ids_replace_subject_ids_by_the_title_word() {
    let ids = |first: &str| {
        subject(first)
            .tabular_ids(true)
            .unwrap()
            .unwrap()
            .map(|id| String::from_utf8(id).unwrap())
    };
    assert_eq!(ids("lcl|Subject_36"), ["unnamed", "unnamed", "unnamed"]);
    assert_eq!(ids("lcl|Subject_27"), ["subject", "subject", "subject"]);
    assert_eq!(
        ids("lcl|Subject_28"),
        ["tpg|DAA10008.1|", "DAA10008", "DAA10008.1"]
    );
    assert_eq!(ids("lcl|Subject_29"), ["pdb|3ABC|C", "3ABC_C", "3ABC_C"]);
    assert!(subject("lcl|Subject_29").check_tabular_ids(true).is_ok());
    // The `saccver`, made once.
    let once = subject("lcl|Subject_28");
    let first = once.tabular_saccver(true).unwrap().unwrap();
    assert_eq!(&*first, b"DAA10008.1");
    assert!(std::sync::Arc::ptr_eq(
        &first,
        &once.tabular_saccver(true).unwrap().unwrap()
    ));
}

fn one_defline(ids: Vec<ParsedSeqId>) -> DeflineSet {
    DeflineSet {
        deflines: vec![BlastDefline {
            title: Some(b"t".to_vec()),
            ids,
            ..BlastDefline::default()
        }],
        utf8: false,
    }
}

fn gb(accession: &str, version: i32) -> ParsedSeqId {
    ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
        accession: Some(accession.as_bytes().to_vec()),
        version: Some(version),
        ..Default::default()
    })
}

fn lines(lines: Vec<Vec<ParsedSeqId>>) -> ParsedSubject {
    let deflines = lines
        .into_iter()
        .map(|ids| BlastDefline {
            title: Some(b"t".to_vec()),
            ids,
            ..BlastDefline::default()
        })
        .collect();
    ParsedSubject::new(DeflineSet {
        deflines,
        utf8: false,
    })
    .unwrap()
}

fn strings(subject: &ParsedSubject) -> Vec<String> {
    subject
        .use_this_seqid()
        .iter()
        .map(|string| String::from_utf8(string.clone()).unwrap())
        .collect()
}

// Review 3 of line RG T3a, D1 (NCBI BLAST+ 2.17.0 probes `T/t3/review3_t3a`, `T/t3/fix3_t3a`):
// the description line shows the display id, the `WorstRank` best id of the line that holds
// the Seq-align id (`gi|777|lcl|Query_1`: the Seq-align id `gi|777`, the display id
// `lcl|Query_1`; `ref|NP_000777.1|lcl|Query_1`: `ref|NP_000777.1|` for both), and looks the
// title up again by it: its keys and those of the ids it matches weakly are `display`
// (`ref|NP_000781|` finds a query `ref|NP_000781.1|`, `ref|NP_000782.1|` not one
// `ref|NP_000782|`), apart from the Seq-align id's `own` and `weak`.
#[test]
fn the_display_id_has_its_own_scope_keys() {
    let local = |id: &str| ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec()));
    let refseq = |accession: &str, version: Option<i32>| {
        ParsedSeqId::Other(crate::seqid::model::TextseqId {
            accession: Some(accession.as_bytes().to_vec()),
            version,
            ..Default::default()
        })
    };
    let fasta = |id: Option<&ParsedSeqId>| {
        let mut out = Vec::new();
        label::fasta(id.unwrap(), &mut out);
        String::from_utf8(out).unwrap()
    };
    let keys = |id: &ParsedSeqId| {
        let mut keys = Vec::new();
        crate::seqid::defline_set::scope_keys(id, &mut keys);
        keys
    };
    let query = local("Query_1");
    let gi = lines(vec![vec![ParsedSeqId::Gi(777), local("Query_1")]]);
    assert_eq!(text(&gi.alignment_fasta()), "gi|777");
    assert_eq!(fasta(gi.display_id()), "lcl|Query_1");
    let scope = gi.query_scope_keys();
    assert!(keys(&query).iter().all(|key| scope.display.contains(key)));
    assert!(!keys(&query).iter().any(|key| scope.own.contains(key)));
    let refq = lines(vec![vec![refseq("NP_000777", Some(1)), local("Query_1")]]);
    assert_eq!(text(&refq.alignment_fasta()), "ref|NP_000777.1|");
    assert_eq!(fasta(refq.display_id()), "ref|NP_000777.1|");
    assert!(!keys(&query)
        .iter()
        .any(|key| refq.query_scope_keys().display.contains(key)));
    let versioned_query = |accession: &str, version: Option<i32>| {
        crate::seqid::id_scope::version_keys(&refseq(accession, version))
            .map(|(own, _)| own)
            .unwrap()
    };
    let weak = lines(vec![vec![ParsedSeqId::Gi(781), refseq("NP_000781", None)]]);
    assert_eq!(fasta(weak.display_id()), "ref|NP_000781|");
    assert!(weak
        .query_scope_keys()
        .display
        .contains(&versioned_query("NP_000781", Some(1))));
    let versioned = lines(vec![vec![
        ParsedSeqId::Gi(782),
        refseq("NP_000782", Some(1)),
    ]]);
    assert!(!versioned
        .query_scope_keys()
        .display
        .contains(&versioned_query("NP_000782", None)));
}

// Line DBID N1b: the titles whose bytes the taxids of the record can change, whatever the
// record's taxids are (`title_reads_source`).
#[test]
fn a_title_reads_the_source_where_it_is_made_from_it_or_its_end_is_rewritten() {
    let titled = |titles: &[&str]| {
        let deflines = titles
            .iter()
            .enumerate()
            .map(|(line, title)| BlastDefline {
                title: Some(title.as_bytes().to_vec()),
                ids: vec![gb(&format!("AF00000{line}"), 1)],
                ..BlastDefline::default()
            })
            .collect();
        ParsedSubject::new(DeflineSet {
            deflines,
            utf8: false,
        })
        .unwrap()
    };
    // A title made from the source names its organism with the taxonomy database only.
    let empty = titled(&[""]);
    assert!(empty.title_reads_source(true, true));
    assert!(empty.title_reads_source(false, true));
    assert!(!empty.title_reads_source(true, false));
    // A protein title ending in ']' where the heading generates the title (another line
    // without a title), with or without the taxonomy database; not for a nucleotide one.
    let bracket = titled(&["alpha, partial [Homo sapiens]"]);
    assert!(!bracket.title_reads_source(true, true));
    let generated = titled(&["alpha, partial [Homo sapiens]", ""]);
    assert!(generated.title_reads_source(true, false));
    assert!(!generated.title_reads_source(false, false));
    assert!(!titled(&["alpha", ""]).title_reads_source(true, true));
}

// Probes G01, G05 of line DI-8 (`GetFilteredRedundantSeqids`, `x_FillAlnDispParams`).
#[test]
fn a_positive_list_names_the_text_ids_or_the_gis() {
    let local = || ParsedSeqId::Local(ObjectId::Str(b"loc".to_vec()));
    // Text ids of every line, whichever comes first; a local id is not named.
    let subject = lines(vec![vec![local()], vec![gb("AF000001", 1)]]);
    assert!(subject.use_this_seqid().is_empty());
    let subject = subject.under_positive_list(true);
    assert_eq!(strings(&subject), ["seqid:AF000001.1"]);
    // The heading leaves out the line the strings do not name, and `>` goes to the next.
    assert_eq!(
        subject.heading(true, false, None).unwrap().unwrap(),
        b">AF000001.1 t\n"
    );
    // An ungapped Seq-align of two HSPs reaches the report without the strings.
    assert_eq!(
        subject.heading(true, true, None).unwrap().unwrap(),
        b">loc t\n AF000001.1 t\n"
    );
    // A local id alone: no string, every line in the heading.
    let subject = lines(vec![vec![local()]]).under_positive_list(true);
    assert!(subject.use_this_seqid().is_empty());
    assert_eq!(
        subject.heading(true, false, None).unwrap().unwrap(),
        b">loc t\n"
    );
    // A Seq-align id that is a GI: only GIs are named, and a line without one is left out.
    let subject = lines(vec![
        vec![ParsedSeqId::Gi(888), gb("GG000001", 1)],
        vec![gb("GG000002", 1)],
        vec![ParsedSeqId::Gi(890), gb("GG000003", 1)],
    ])
    .under_positive_list(true);
    assert_eq!(strings(&subject), ["gi:888", "gi:890"]);
    assert_eq!(
        subject.heading(true, false, None).unwrap().unwrap(),
        b">GG000001.1 t\n GG000003.1 t\n"
    );
    // A name without an accession is not named; a PDB id is.
    let name_only = ParsedSeqId::Prf(crate::seqid::model::TextseqId {
        name: Some(b"NAME".to_vec()),
        ..Default::default()
    });
    let pdb = ParsedSeqId::Pdb(crate::seqid::model::PdbSeqId {
        mol: b"3DEF".to_vec(),
        chain: None,
        rel: None,
        chain_id: Some(b"B".to_vec()),
    });
    let subject = lines(vec![vec![name_only], vec![pdb]]).under_positive_list(true);
    assert_eq!(strings(&subject), ["seqid:3DEF_B"]);
    // Without a positive list nothing changes.
    let subject = lines(vec![vec![local()], vec![gb("AF000001", 1)]]).under_positive_list(false);
    assert!(subject.use_this_seqid().is_empty());
    assert_eq!(
        subject.heading(true, false, None).unwrap().unwrap(),
        b">loc t\n AF000001.1 t\n"
    );
}

#[test]
fn the_kind_of_a_seq_align_id_refuses_nothing_before_its_lookup() {
    // NCBI's lookup decides whether its formatters find these (`found_instead`; probe of line
    // DBID N2b), not the kind of id.
    for id in [
        ParsedSeqId::Local(ObjectId::Id(123)),
        ParsedSeqId::Gibbsq(123),
        ParsedSeqId::NamedAnnotTrack(crate::seqid::model::TextseqId {
            accession: Some(b"NA000001".to_vec()),
            version: Some(1),
            ..Default::default()
        }),
        ParsedSeqId::General(Dbtag {
            db: b"PIG".to_vec(),
            tag: ObjectId::Id(123),
        }),
        ParsedSeqId::General(Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            tag: ObjectId::Id(5),
        }),
    ] {
        let subject = ParsedSubject::new(one_defline(vec![id])).unwrap();
        assert!(!subject.unfetchable(), "{subject:?}");
        assert!(subject.fetched(false, true, HitLoad::Loaded).is_ok());
    }
    // A local id that alone decides the definition line of the Bioseq: no lookup finds it.
    let shown = lines(vec![
        vec![ParsedSeqId::Local(ObjectId::Id(777))],
        vec![ParsedSeqId::Gi(777)],
    ]);
    assert!(shown.unfetchable());
    let error = shown.fetched(false, true, HitLoad::Loaded).unwrap_err();
    assert!(
        error.contains("Seq-align id gi|777, by which NCBI BLAST+'s formatters do not find"),
        "{error}"
    );
    // A GI in the set makes the GI the Seq-align id.
    let with_gi = one_defline(vec![ParsedSeqId::Gibbsq(123), ParsedSeqId::Gi(7)]);
    assert!(ParsedSubject::new(with_gi)
        .unwrap()
        .fetched(false, true, HitLoad::Loaded)
        .is_ok());
    let replaced = DeflineSet {
        deflines: vec![BlastDefline {
            ids: vec![ParsedSeqId::Gi(7)],
            replaced_title_bytes: 1,
            ..BlastDefline::default()
        }],
        utf8: false,
    };
    assert!(ParsedSubject::new(replaced).is_err());
    assert!(ParsedSubject::new(DeflineSet::default()).is_err());
}

#[test]
fn a_seq_align_id_by_which_ncbi_finds_another_sequence_or_none() {
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
    assert!(subject
        .clone()
        .found_instead(None)
        .fetched(false, true, HitLoad::Loaded)
        .is_ok());
    let error = subject
        .clone()
        .found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4))))
        .fetched(false, true, HitLoad::Loaded)
        .unwrap_err();
    assert!(
        error.contains("lookup gives database sequence 4 first"),
        "{error}"
    );
    let error = subject
        .clone()
        .found_instead(Some(FoundInstead::Nothing))
        .fetched(false, true, HitLoad::Loaded)
        .unwrap_err();
    assert!(error.contains("do not find the sequence"), "{error}");
    let error = subject
        .found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
            "DB contains no accession info.",
        ))))
        .fetched(false, true, HitLoad::Loaded)
        .unwrap_err();
    assert!(
        error.contains("whose lookup in the database fails (DB contains no accession info.;"),
        "{error}"
    );
}

// Line RG item T5d fixes 3 and 4: the sequence that the lookup gives first, another one whose
// header holds a UTF8String, is read by NCBI's scope with its warning, which the load of the
// id posts (`common::database::load_read_texts`), also where the batch reads one by one: no
// refusal.
#[test]
fn a_utf8string_header_that_the_scope_reads_is_no_refusal() {
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
    assert_eq!(subject.scope_utf8_read(), None);
    let read = subject.scope_utf8(Some(DbSequence::of_unrestricted(3)));
    assert_eq!(read.scope_utf8_read(), Some(DbSequence::of_unrestricted(3)));
    for one_by_one in [false, true] {
        assert_eq!(read.fetched(true, one_by_one, HitLoad::Loaded), Ok(()));
        assert_eq!(read.fetched(false, one_by_one, HitLoad::Loaded), Ok(()));
        assert_eq!(read.pssm_read(one_by_one), Ok(()));
    }
}

// Line DBID pfoff-1 (DECISIONS 62): where a batch's prefetch is off, the XML2 and JSON `Hit`
// loads the subject by its Seq-align id and the lookup throws NCBI's `CSeqDBException` (blastp
// 2.17.0, outfmt 15 of `pP pU`, `R/mixdb2/review1/out/ncbi/bp.lim1.U2.15.err`); only the text
// seen there gives that exception (review 1, M2).
#[test]
fn only_the_lookup_text_ncbi_was_seen_to_throw_gives_its_exception() {
    let failed = LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO);
    let exception = failed.exception().unwrap();
    assert_eq!(exception.get_msg(), b"DB contains no accession info.");
    assert_eq!(
        String::from_utf8(exception.what()).unwrap(),
        // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
        crate::wording::losat_text(&format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp\", line 142: Error: (CSeqDBException::eInvalid) BLASTDB::ncbi::CBlastLMDBManager::CBlastEnv::GetDbi() - DB contains no accession info.\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        ))
    );
    // LOSAT's own error, NCBI's text of a place that LOSAT does not know (the text of
    // `CATCH_ALL` for `eTooManyOpenFiles`), and the texts of NCBI's reader whose places LOSAT
    // knows but that no lookup was seen to throw (`seqdb::throw_site`): no exception.
    for failure in [
        LookupFailure::losat("byte range 0..9 passes the end of the file"),
        LookupFailure::ncbi("Too many open files, please raise the open file limit"),
        LookupFailure::ncbi("Error: File (x.pos) not found."),
        LookupFailure::ncbi("Database name is required."),
        LookupFailure::ncbi(crate::seqdb::tax_files::NO_VOLUME_INFO),
    ] {
        assert!(failure.exception().is_none(), "{failure:?}");
    }
}

// Line DBID pfoff-1 (DECISIONS 62; review 1, H1 and M1): the one answer of what the XML2 and
// JSON reports load (`HitLoads`): in a report read one by one, the `Hit`s before the first
// whose lookup throws load their subjects, that one throws, and those after it load nothing
// (NCBI's report function ends the report there); a subject is loaded where some report
// loads it; a report whose batch prefetches loads every subject; a run that writes another
// format loads every subject (which refuses the throw). The refusals take that answer.
#[test]
fn the_hits_after_one_whose_lookup_throws_load_nothing() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 0 and 4 throw; 2 and 3 are copies of one record (one `Hit`); 1, 5 and 6 are found.
    let copy = plain("AAA00003");
    let parsed: Vec<Option<ParsedSubject>> = vec![
        Some(plain("AAA00001").found_instead(throws())),
        Some(plain("AAA00002")),
        Some(copy.clone()),
        Some(copy),
        Some(plain("AAA00005").found_instead(throws())),
        Some(plain("AAA00006")),
        None,
    ];
    assert!(joins_as_copy(
        parsed[2].as_ref().unwrap(),
        parsed[3].as_ref().unwrap()
    ));
    let subject = |key: usize| parsed[key].as_ref();
    let none = |_: usize, _: usize| 0;
    let loads = HitLoads::of(
        Some(RunThrow::every_row(AfterThrow::EndsReport)),
        [
            // Read one by one: 1 loads, 0 throws, the copies 2 and 3, 4 and 6 load nothing.
            (
                true,
                vec![Some(1), None, Some(0), Some(2), Some(3), Some(4), Some(6)],
            ),
            // 6 loads in a later report; 5 only in a batch that prefetches.
            (true, vec![Some(6)]),
            (false, vec![Some(5)]),
        ],
        subject,
        none,
    );
    assert_eq!(loads.load(1), HitLoad::Loaded);
    assert_eq!(loads.load(0), HitLoad::Throws);
    assert_eq!(
        loads
            .thrown(0)
            .map(crate::ncbi_exception::NcbiException::get_msg),
        Some(b"DB contains no accession info.".as_slice())
    );
    for unread in [2, 3, 4] {
        assert_eq!(loads.load(unread), HitLoad::Unread, "{unread}");
        assert!(loads.thrown(unread).is_none());
    }
    assert!(loads.loaded(5) && loads.loaded(6));
    // A key the reports never show reads as loaded; renamed keys keep their answers.
    assert!(loads.loaded(7));
    let renamed = loads.renamed(|key| (key != 1).then_some(key + 10));
    assert_eq!(renamed.load(10), HitLoad::Throws);
    assert_eq!(renamed.load(14), HitLoad::Unread);
    assert_eq!(renamed.load(0), HitLoad::Loaded);
    // A report whose batch prefetches, and a run that writes another format, load every
    // subject.
    for (after, one_by_one) in [
        (Some(RunThrow::every_row(AfterThrow::EndsReport)), false),
        (None, true),
    ] {
        let every = HitLoads::of(after, [(one_by_one, vec![Some(0), Some(1)])], subject, none);
        assert!(every.loaded(0) && every.loaded(1));
    }
    // The refusals take the answer: a `Hit` that throws is not refused for what the lookup
    // finds, a subject that no report loads not at all, every other one as before (a lookup
    // that fails on a file too).
    let thrown = subject(0).unwrap();
    assert!(thrown.fetched(true, true, HitLoad::Throws).is_ok());
    assert!(thrown.fetched(true, true, HitLoad::Unread).is_ok());
    let error = thrown.fetched(false, true, HitLoad::Loaded).unwrap_err();
    assert!(
        error.contains("whose lookup in the database fails (DB contains no accession info.;"),
        "{error}"
    );
    let nothing = plain("AAA00007").found_instead(Some(FoundInstead::Nothing));
    assert!(nothing.fetched(false, true, HitLoad::Loaded).is_err());
    assert!(nothing.fetched(false, true, HitLoad::Unread).is_ok());
    // A weak match is refused where the `Hit` throws too (NCBI's scope may find it first).
    let weak = thrown
        .clone()
        .weak_match(Some(WeakMatch::Held(DbSequence::of_unrestricted(5))));
    assert!(weak.fetched(false, true, HitLoad::Throws).is_err());
}

// Line DBID pfoff-2 (DECISIONS 62 (12), probe `R/pfoff/p2`): what the tables load after a
// subject whose lookup throws (`AfterThrow`, `FormatThrow`): a table of a field of the subject
// block writes each row and goes on, every throwing subject throws and every other loads; a
// table that reads the subject after that block (`slen`, the letters) and outfmt 13 and 14
// end the run, so nothing after the first throw loads, in that report or a later one; a
// table of other fields reads nothing; the run's answer is its formats' (`FormatThrow::and`).
#[test]
fn the_formats_say_what_loads_after_a_throw() {
    use crate::blastinput::app::Program;
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 1 and 3 throw; 0, 2 and 4 are found.
    let parsed: Vec<ParsedSubject> = vec![
        plain("AAA00001"),
        plain("AAA00002").found_instead(throws()),
        plain("AAA00003"),
        plain("AAA00004").found_instead(throws()),
        plain("AAA00005"),
    ];
    let subject = |key: usize| parsed.get(key);
    let reports = || {
        [
            (true, vec![Some(0), Some(1), Some(2)]),
            (true, vec![Some(3), Some(4)]),
        ]
    };
    let every = |after| Some(RunThrow::every_row(after));
    let none = |_: usize, _: usize| 0;
    let rows = HitLoads::of(every(AfterThrow::WritesRow), reports(), subject, none);
    assert!(rows.loaded(0) && rows.loaded(2) && rows.loaded(4));
    assert_eq!(rows.load(1), HitLoad::Throws);
    assert_eq!(rows.load(3), HitLoad::Throws);
    assert!(rows.any_throws());
    let run = HitLoads::of(every(AfterThrow::EndsRun), reports(), subject, none);
    assert!(run.loaded(0));
    assert_eq!(run.load(1), HitLoad::Throws);
    // The run ended at 1: the next report reads nothing (3 does not throw).
    for unread in [2, 3, 4] {
        assert_eq!(run.load(unread), HitLoad::Unread, "{unread}");
    }
    let report = HitLoads::of(every(AfterThrow::EndsReport), reports(), subject, none);
    assert_eq!(report.load(2), HitLoad::Unread);
    assert_eq!(report.load(3), HitLoad::Throws);
    assert_eq!(report.load(4), HitLoad::Unread);
    // The formats' answers.
    let of = |outfmt: &str| FormatThrow::of_outfmt(outfmt, Program::Blastp, RequestKind::Blastp);
    let writes = |after| FormatThrow::Writes(RunThrow::every_row(after));
    for rows in [
        "6",
        "7 qseqid sseqid evalue",
        "10 qseqid stitle slen",
        "20",
        "20 qseqid sacc",
        "18",
    ] {
        assert_eq!(of(rows), writes(AfterThrow::WritesRow), "{rows}");
    }
    assert_eq!(of("6 qseqid slen"), writes(AfterThrow::EndsRun));
    assert_eq!(of("6 qseqid pident"), writes(AfterThrow::EndsRun));
    // `nident` alone reads the subject only at a row with identities.
    assert_eq!(
        of("6 qseqid nident"),
        FormatThrow::Writes(RunThrow {
            after: AfterThrow::EndsRun,
            rows: super::super::tabular::SubjectFetch::WithIdentities,
        })
    );
    assert_eq!(of("6 qseqid evalue bitscore"), FormatThrow::Unread);
    assert_eq!(of("13"), writes(AfterThrow::EndsRun));
    assert_eq!(of("15"), writes(AfterThrow::EndsReport));
    assert_eq!(of("8"), FormatThrow::Unread);
    for refused in ["5", "17"] {
        assert_eq!(
            of(refused),
            FormatThrow::Refused(ThrowRefusal::Format),
            "{refused}"
        );
    }
    // Pairwise: the first shown alignment ends the run, a later one is skipped; query-anchored:
    // any shown alignment ends the run.
    assert_eq!(of("0"), writes(AfterThrow::SkipsAlignment));
    for anchored in ["1", "2", "3", "4"] {
        assert_eq!(
            of(anchored),
            writes(AfterThrow::EndsAtAlignment),
            "{anchored}"
        );
    }
    // BLASTN's `pident` and `nident` read nothing of the subject.
    for unread in ["6 qseqid pident", "6 qseqid nident"] {
        assert_eq!(
            FormatThrow::of_outfmt(unread, Program::Blastn, RequestKind::Blastn),
            FormatThrow::Unread,
            "{unread}"
        );
    }
    // A run's answer: a refusal refuses, a format that reads nothing leaves the other's.
    let answer = |outfmts: &[&str]| {
        FormatThrow::of_outfmts(
            outfmts.iter().copied(),
            Program::Blastp,
            RequestKind::Blastp,
        )
    };
    let after = |outfmts: &[&str]| answer(outfmts).run();
    // The reason of each refusal: a format of its own (outfmt 5, 17) before a mix.
    assert_eq!(
        answer(&["6", "0"]),
        FormatThrow::Refused(ThrowRefusal::Mixed)
    );
    assert_eq!(
        answer(&["0", "5", "6"]),
        FormatThrow::Refused(ThrowRefusal::Format)
    );
    assert_eq!(
        after(&["6", "8"]),
        Some(RunThrow::every_row(AfterThrow::WritesRow))
    );
    // A pairwise format with another whose throws end elsewhere: refused.
    assert_eq!(after(&["6", "0"]), None);
    assert_eq!(after(&["0", "3"]), None);
    assert_eq!(
        after(&["0", "8"]),
        Some(RunThrow::every_row(AfterThrow::SkipsAlignment))
    );
    assert_eq!(after(&["6 qseqid evalue"]), None);
    assert_eq!(
        after(&["6 qseqid slen", "15"]),
        Some(RunThrow::every_row(AfterThrow::EndsReport))
    );
    // Rows of `nident` alone and rows of every HSP load at every row.
    assert_eq!(
        after(&["6 qseqid nident", "6 qseqid slen"]),
        Some(RunThrow::every_row(AfterThrow::EndsRun))
    );
}

// Review 2 of line DBID pfoff-2, M1 and L1 (probe `R/pfoff/p2/review2`, `zP`): where `nident`
// alone reads the subject (`SubjectFetch::WithIdentities`), a throwing subject whose rows have
// no identities reads nothing (its rows print the 0 of `x_ResetFields` and the run goes on);
// the run ends at the first throwing subject with a row of identities, and nothing after it
// loads, in that report, a later one, or one that a later call makes apart
// (`RunLoads::loads`).
#[test]
fn a_throw_ends_the_run_only_at_a_row_that_reads_the_subject() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 1 and 3 throw; 0, 2 and 4 are found.
    let parsed: Vec<ParsedSubject> = vec![
        plain("AAA00001"),
        plain("AAA00002").found_instead(throws()),
        plain("AAA00003"),
        plain("AAA00004").found_instead(throws()),
        plain("AAA00005"),
    ];
    let subject = |key: usize| parsed.get(key);
    let nident = Some(RunThrow {
        after: AfterThrow::EndsRun,
        rows: super::super::tabular::SubjectFetch::WithIdentities,
    });
    // The identities of each report's entries: 1's rows have none, 3's have some.
    let identities = |report: usize, entry: usize| [[5, 0, 7], [2, 4, 0]][report][entry];
    let reports = || {
        [
            (true, vec![Some(0), Some(1), Some(2)]),
            (true, vec![Some(3), Some(4)]),
        ]
    };
    let loads = HitLoads::of(nident, reports(), subject, identities);
    assert_eq!(loads.load(1), HitLoad::Unread);
    assert!(loads.thrown(1).is_none());
    assert!(loads.loaded(0) && loads.loaded(2));
    assert_eq!(loads.load(3), HitLoad::Throws);
    assert_eq!(loads.load(4), HitLoad::Unread);
    // Every row reads the subject: 1 throws and ends the run.
    let every = HitLoads::of(
        Some(RunThrow::every_row(AfterThrow::EndsRun)),
        reports(),
        subject,
        identities,
    );
    assert_eq!(every.load(1), HitLoad::Throws);
    assert_eq!(every.load(3), HitLoad::Unread);
    // The reports made apart, one call each, carry the end of the run.
    let run = RunLoads::new(nident);
    let [first, second] = reports();
    let first = run.loads([first], subject, |_, entry| [5, 0, 7][entry]);
    assert!(!run.ended());
    assert_eq!(first.load(1), HitLoad::Unread);
    let second = run.loads([second], subject, |_, entry| [2, 4][entry]);
    assert!(run.ended());
    assert_eq!(second.load(3), HitLoad::Throws);
    let later = run.loads([(true, vec![Some(0), Some(2)])], subject, |_, _| 9);
    assert_eq!(later.load(0), HitLoad::Unread);
    assert_eq!(later.load(2), HitLoad::Unread);
    // A run without a throw's answer loads every subject and never ends.
    let none = RunLoads::new(None);
    assert!(none.loads(reports(), subject, identities).loaded(3));
    assert!(!none.ended());
}

// Review 3 of line DBID pfoff-3, N1 (probe `R/pfoff/p3/fix4`, `probe_n1.py`): a report whose load
// ends the run ends NCBI's application there, so no report after it runs and no batch after it is
// prefetched (`HitLoads::reaches`): BLASTP and TBLASTN, which model every batch's prefetch before
// any byte, stop there, and a later batch's prefetch (its taxonomy reads, its stop) is never made.
// The run reaches every report where no load ends it, none where an earlier call ended it, and up
// to the report that ends it otherwise, at its last report too.
#[test]
fn the_run_reaches_no_report_after_the_one_that_ends_it() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 1 throws; 0 and 2 are found.
    let parsed: Vec<ParsedSubject> = vec![
        plain("AAA00001"),
        plain("AAA00002").found_instead(throws()),
        plain("AAA00003"),
    ];
    let subject = |key: usize| parsed.get(key);
    let none = |_: usize, _: usize| 0;
    let every = |after| Some(RunThrow::every_row(after));
    // Report 1 of three loads 1, which ends the run: reports 0 and 1 run, 2 does not.
    let reports = || {
        [
            (true, vec![Some(0)]),
            (true, vec![Some(2), Some(1)]),
            (true, vec![Some(0), Some(2)]),
        ]
    };
    let ended = HitLoads::of(every(AfterThrow::EndsRun), reports(), subject, none);
    assert!(ended.reaches(0) && ended.reaches(1));
    assert!(!ended.reaches(2) && !ended.reaches(3));
    // A format that writes the row and goes on, one that ends its report, a refusal and a run
    // whose formats read nothing: every report runs.
    for loads in [
        HitLoads::of(every(AfterThrow::WritesRow), reports(), subject, none),
        HitLoads::of(every(AfterThrow::EndsReport), reports(), subject, none),
        HitLoads::of(
            FormatThrow::Refused(ThrowRefusal::Mixed),
            reports(),
            subject,
            none,
        ),
        HitLoads::of(FormatThrow::Unread, reports(), subject, none),
    ] {
        assert!(loads.reaches(2) && loads.reaches(usize::MAX), "{loads:?}");
    }
    // The last report ends the run: every report runs, none after it.
    let last = HitLoads::of(
        every(AfterThrow::EndsRun),
        [(true, vec![Some(0)]), (true, vec![Some(1)])],
        subject,
        none,
    );
    assert!(last.reaches(1) && !last.reaches(2));
    // The reports made apart, one call each: a call after the end reaches none of its reports.
    let run = RunLoads::new(every(AfterThrow::EndsRun));
    let first = run.loads([(true, vec![Some(0), Some(1)])], subject, none);
    assert!(first.reaches(0) && !first.reaches(1));
    let later = run.loads([(true, vec![Some(2)])], subject, none);
    assert!(!later.reaches(0));
    // The keys renamed keep the reports reached.
    assert!(!later.renamed(Some).reaches(0));
}

// Review 3 of line DBID pfoff-3, N2: a call of the library writes several formats with one
// result, though `api::run` makes one search for each hit list size (`search_groups`): its
// refusal of a subject whose lookup throws is the call's, over every format of the call
// (`FormatThrow::of_search` with `ReportOutputs::other_outfmts`), so that the split by
// `-num_descriptions` and `-num_alignments` does not change it; a search keeps its own answer
// where the call does not refuse.
#[test]
fn the_call_refuses_a_throw_over_all_its_searches() {
    use crate::blastinput::app::Program;
    let search = |own: &[&str], others: &[&str]| {
        FormatThrow::of_outfmts(own.iter().copied(), Program::Blastp, RequestKind::Blastp)
            .of_search(others, Program::Blastp, RequestKind::Blastp)
    };
    let writes = |after| FormatThrow::Writes(RunThrow::every_row(after));
    // One search: its own answer.
    assert_eq!(search(&["6 qseqid slen"], &[]), writes(AfterThrow::EndsRun));
    assert_eq!(
        search(&["0", "6"], &[]),
        FormatThrow::Refused(ThrowRefusal::Mixed)
    );
    // The table's search and outfmt 0's search of one call: both refuse, as one search of both.
    assert_eq!(
        search(&["6 qseqid slen"], &["0"]),
        FormatThrow::Refused(ThrowRefusal::Mixed)
    );
    assert_eq!(
        search(&["0"], &["6 qseqid slen"]),
        FormatThrow::Refused(ThrowRefusal::Mixed)
    );
    // Outfmt 5 in another search: the call's refusal is outfmt 5's.
    assert_eq!(
        search(&["6"], &["5"]),
        FormatThrow::Refused(ThrowRefusal::Format)
    );
    // Formats that agree: each search keeps its own answer (a table of fields that read nothing
    // of the subject too).
    assert_eq!(search(&["0"], &["8"]), writes(AfterThrow::SkipsAlignment));
    assert_eq!(
        search(&["6 qseqid slen"], &["15"]),
        writes(AfterThrow::EndsRun)
    );
    assert_eq!(search(&["6 qseqid evalue"], &["6"]), FormatThrow::Unread);
}

// Review 3 of line DBID pfoff-3, L-2: the reason of a refusal of a subject whose lookup throws
// says what NCBI does and why LOSAT refuses, one text for each situation
// (`ThrowRefusal::words`): outfmt 5 and 17, which the command line meets too, say nothing of the
// library; a call that mixes the formats says the library's reason.
#[test]
fn a_refused_throw_gives_the_reason_of_its_situation() {
    let thrown = ParsedSubject::new(one_defline(vec![gb("AAA00002", 1)]))
        .unwrap()
        .found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        ))));
    let reason = |load| thrown.fetched(false, true, load).unwrap_err();
    let format = reason(HitLoad::Refused(ThrowRefusal::Format));
    assert!(format.contains(LookupFailure::NCBI_THROWS), "{format}");
    assert!(format.contains(ThrowRefusal::Format.words()), "{format}");
    assert!(!format.contains("library"), "{format}");
    let mixed = reason(HitLoad::Refused(ThrowRefusal::Mixed));
    assert!(mixed.contains(ThrowRefusal::Mixed.words()), "{mixed}");
    assert!(!mixed.contains(ThrowRefusal::Format.words()), "{mixed}");
    // A throw that the formats write as NCBI does is not refused for what the lookup finds.
    assert!(thrown.fetched(false, true, HitLoad::Throws).is_ok());
    // The loads of a run whose formats refuse give each subject the reason.
    let loads = HitLoads::of(
        FormatThrow::Refused(ThrowRefusal::Format),
        [(true, vec![Some(0)])],
        |_| Some(&thrown),
        |_, _| 0,
    );
    assert_eq!(loads.load(0), HitLoad::Refused(ThrowRefusal::Format));
    assert!(loads.loaded(0));
}

// Review 1 of line DBID pfoff-2, H1 and M1 (probe `R/pfoff/p2/review1`, `*.H.*`, `*.H2.*`):
// the tables' taxonomy, from which their names are looked up, is made only of the subjects
// that a report loads (`common::tabular::tabular_taxonomy` with `HitLoads`): a throwing
// subject's taxid (9606 here), which NCBI never reads, opens no taxonomy database and is not
// looked up, so its record is never refused and the texts of opening wait for a loaded row.
#[test]
fn the_names_of_the_tables_take_only_the_loaded_subjects() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let with_taxid = |accession: &str, taxid: i32| {
        let mut set = one_defline(vec![gb(accession, 1)]);
        set.deflines[0].taxid = Some(taxid);
        ParsedSubject::new(set).unwrap()
    };
    // 0 throws (taxid 9606); 1 is found (taxid 99990).
    let parsed: Vec<Option<std::sync::Arc<ParsedSubject>>> = vec![
        Some(std::sync::Arc::new(
            with_taxid("NM_000001", 9606).found_instead(throws()),
        )),
        Some(std::sync::Arc::new(with_taxid("XM_000002", 99990))),
    ];
    let loads = HitLoads::of(
        Some(RunThrow::every_row(AfterThrow::WritesRow)),
        [(true, vec![Some(0), Some(1)])],
        |rank| parsed[rank].as_deref(),
        |_, _| 0,
    );
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/prot_v4");
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    let taxonomy = crate::algorithm::common::tabular::tabular_taxonomy(
        &db,
        &[0, 1],
        &parsed,
        [0, 1],
        &loads,
        "BLASTP",
    )
    .unwrap();
    assert_eq!(
        taxonomy[0],
        super::super::tabular::TabularTaxonomy::default()
    );
    assert_eq!(taxonomy[1].taxid, 99990);
    // Every subject loaded (no report reads one by one): both taxids.
    let every = crate::algorithm::common::tabular::tabular_taxonomy(
        &db,
        &[0, 1],
        &parsed,
        [0, 1],
        &HitLoads::default(),
        "BLASTP",
    )
    .unwrap();
    assert_eq!(every[0].taxid, 9606);
    // The names of the rows: no lookup asks a taxid of a subject that is not loaded, so the
    // damaged taxonomy database (`bad_magic`, whose open writes a text) is not opened for it.
    let spec = super::super::tabular::TabularSpec::of_app(
        crate::blastinput::app::ReportFormat::Tabular,
        &crate::blastinput::app::FormatChoice {
            number: 6,
            spec: "qseqid sseqid ssciname".to_string(),
            delimiter: String::new(),
        },
        RequestKind::Blastp,
    )
    .unwrap();
    let damaged =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/bad_magic");
    let names = |taxonomy: &super::super::tabular::TabularTaxonomy| {
        let search =
            crate::seqdb::SearchPath::new(&damaged, Some(damaged.as_os_str()), None).unwrap();
        super::super::tabular::TabularNames::look_up(
            &spec,
            &crate::seqdb::taxdb::LazyTaxDb::new(search),
            [taxonomy],
            "BLASTP",
        )
        .unwrap()
    };
    assert!(!names(&taxonomy[0]).looks_up());
    assert!(names(&every[0]).looks_up());
}

// Line DBID pssmread (probe `R/pssmread`, reviews 1 and 2): LOSAT's refusals of what
// psiblast's PSSM reads by the Seq-align id (what the scope gets is
// `common::database::scope_bioseq`'s).
#[test]
fn the_pssm_refuses_reads_that_losat_does_not_make_as_ncbi() {
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
    let read = |subject: &ParsedSubject| subject.pssm_read(false);
    assert_eq!(read(&subject), Ok(()));
    assert!(subject.bioseq_holds_alignment_id());
    assert_eq!(subject.other_holder(), None);
    let other = subject
        .clone()
        .found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4))));
    assert_eq!(read(&other), Ok(()));
    for found in [
        FoundInstead::Nothing,
        FoundInstead::ThrowsOnId("Error: ... is not a valid ordinal id".into()),
        FoundInstead::OtherTaxids(DbSequence::of_unrestricted(2)),
    ] {
        assert_eq!(read(&subject.clone().found_instead(Some(found))), Ok(()));
    }
    // Refused: another holder besides what the lookup gives, a lookup failing on a file, the
    // subject first with another holder.
    let held = other
        .clone()
        .also_held(Some(DbSequence::of_unrestricted(6)));
    assert_eq!(held.other_holder(), Some(DbSequence::of_unrestricted(6)));
    assert!(read(&held)
        .unwrap_err()
        .contains("which database sequence 6 also holds"));
    let none_held = subject
        .clone()
        .found_instead(Some(FoundInstead::Nothing))
        .also_held(Some(DbSequence::of_unrestricted(6)));
    assert!(read(&none_held).is_err());
    for found in [
        FoundInstead::AlsoHeld(DbSequence::of_unrestricted(6)),
        FoundInstead::Failed(LookupFailure::ncbi("DB contains no accession info.")),
    ] {
        assert!(read(&subject.clone().found_instead(Some(found))).is_err());
    }
    // A weak match refuses where the batch reads one by one, whatever the lookup gives.
    let weak = other
        .clone()
        .weak_match(Some(WeakMatch::Held(DbSequence::of_unrestricted(7))));
    assert_eq!(weak.pssm_read(false), Ok(()));
    assert!(weak.pssm_read(true).is_err());
    // A definition line that the formatters do not find: its Bioseq does not hold the id, and
    // the reasons are what the lookup finds (review 2, M2).
    let shown = lines(vec![
        vec![ParsedSeqId::Local(ObjectId::Id(777))],
        vec![ParsedSeqId::Gi(777)],
    ]);
    assert!(shown.unfetchable());
    assert!(!shown.bioseq_holds_alignment_id());
    assert_eq!(read(&shown), Ok(()));
    let also = shown
        .clone()
        .found_instead(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(2))));
    assert_eq!(also.other_holder(), Some(DbSequence::of_unrestricted(2)));
    let error = read(&also).unwrap_err();
    assert!(
        error.contains("which database sequence 2 also holds"),
        "{error}"
    );
    assert!(!error.contains("formatters do not find"), "{error}");
    let failed = shown.found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
        "no accession info",
    ))));
    assert!(read(&failed)
        .unwrap_err()
        .contains("whose lookup in the database fails"));
}

// Line DBID weakbatch (probe w0): NCBI's scope matches a Seq-align id weakly only where the
// formatters read the batch's subjects one by one, its prefetch off; the prefetch loads each
// subject by its own id.
#[test]
fn a_weak_match_is_refused_only_where_the_batch_reads_one_by_one() {
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
    for weak in [
        WeakMatch::Held(DbSequence::of_unrestricted(4)),
        WeakMatch::Unlisted,
        WeakMatch::Failed("no accession info".into()),
    ] {
        let subject = subject.clone().weak_match(Some(weak.clone()));
        assert!(
            subject.fetched(true, false, HitLoad::Loaded).is_ok(),
            "{weak:?}"
        );
        let error = subject.fetched(false, true, HitLoad::Loaded).unwrap_err();
        assert!(
            error.contains("in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off"),
            "{error}"
        );
        assert_eq!(subject.prefetched(), Prefetched::Itself, "{weak:?}");
        // What the lookup finds is refused whatever the batch: a copy with other taxids in
        // the reports of taxids, before the weak match.
        let copy = subject.found_instead(Some(FoundInstead::OtherTaxids(
            DbSequence::of_unrestricted(2),
        )));
        assert!(copy.fetched(false, false, HitLoad::Loaded).is_ok());
        assert!(copy
            .fetched(true, false, HitLoad::Loaded)
            .unwrap_err()
            .contains("other taxids"));
        assert!(copy
            .fetched(false, true, HitLoad::Loaded)
            .unwrap_err()
            .contains("weakly"));
    }
    let other = subject.found_instead(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(4))));
    assert!(other.fetched(false, false, HitLoad::Loaded).is_err());
}

// Line DBID P1: what NCBI's prefetch reads for a subject that no report shows (psiblast's
// traceback hit list beyond -max_target_seqs), from what its lookup finds.
#[test]
fn the_prefetch_of_an_unshown_subject_reads_what_the_lookup_finds() {
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
    let prefetched =
        |found: Option<FoundInstead>| subject.clone().found_instead(found).prefetched();
    for itself in [
        None,
        Some(FoundInstead::OtherTaxids(DbSequence::of_unrestricted(2))),
        Some(FoundInstead::CopyApart(
            DbSequence::of_unrestricted(2),
            CopyApart::Title,
        )),
    ] {
        assert_eq!(prefetched(itself.clone()), Prefetched::Itself, "{itself:?}");
    }
    assert_eq!(prefetched(Some(FoundInstead::Nothing)), Prefetched::Nothing);
    // A file's error stops NCBI's prefetch; the exception of the id is its warning.
    assert!(matches!(
        prefetched(Some(FoundInstead::Failed(LookupFailure::ncbi("no accession info")))),
        Prefetched::Failed(reason) if reason.contains("prefetch stops with a BLAST Database error")
    ));
    let throws = Some(FoundInstead::ThrowsOnId("Invalid choice selection".into()));
    assert!(matches!(
        prefetched(throws.clone()),
        Prefetched::Failed(reason) if reason.contains("writes Error pre-fetching sequence data")
    ));
    assert!(subject.clone().found_instead(throws).prefetch_fails());
    assert!(!subject
        .clone()
        .found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
            "no accession info"
        ))))
        .prefetch_fails());
    assert!(matches!(
        prefetched(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4)))),
        Prefetched::Another {
            other: Some(other),
            or_itself: false,
            ..
        } if other == DbSequence::of_unrestricted(4)
    ));
    assert!(matches!(
        prefetched(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(4)))),
        Prefetched::Another {
            other: Some(other),
            or_itself: true,
            ..
        } if other == DbSequence::of_unrestricted(4)
    ));
    assert!(matches!(
        prefetched(Some(FoundInstead::CopyApart(
            DbSequence::of_unrestricted(4),
            CopyApart::TaxonomyTexts
        ))),
        Prefetched::Another { other: None, .. }
    ));
    // Other holders that LOSAT cannot list (line DBID v4late): what the lookup gives, or
    // nothing where such a holder's Bioseq read before stands for it.
    assert!(matches!(
        prefetched(Some(FoundInstead::HoldersUnknown {
            instead: Some(DbSequence::of_unrestricted(4)),
            why: String::new(),
        })),
        Prefetched::Another {
            other: Some(other),
            or_itself: false,
            or_none: true,
            ..
        } if other == DbSequence::of_unrestricted(4)
    ));
    assert!(matches!(
        prefetched(Some(FoundInstead::HoldersUnknown {
            instead: None,
            why: String::new(),
        })),
        Prefetched::ItselfOrNothing { .. }
    ));
    // A subject whose definition line the formatters do not find: what its lookup gives, as
    // any other (review 2 of line DBID pssmread, H3).
    let shown = lines(vec![
        vec![ParsedSeqId::Local(ObjectId::Id(777))],
        vec![ParsedSeqId::Gi(777)],
    ]);
    assert_eq!(shown.prefetched(), Prefetched::Itself);
    assert!(matches!(
        shown.clone().found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(3)))).prefetched(),
        Prefetched::Another { other: Some(other), or_itself: false, or_none: false, reason }
            if other == DbSequence::of_unrestricted(3) && reason.contains("gives database sequence 3 first")
    ));
    assert_eq!(
        shown
            .found_instead(Some(FoundInstead::Nothing))
            .prefetched(),
        Prefetched::Nothing
    );
}

// Line DBID N1b L1: a copy with other taxids from which NCBI can make the subject's title is
// refused by the reports that show a generated title (outfmt 0, the tabular word of a
// replaced id) and by those of taxids, not by the others (tabular ids, `stitle`, XML).
#[test]
fn a_title_made_from_a_copy_s_taxids_is_refused_where_a_generated_title_is_shown() {
    use crate::algorithm::common::database::check_shown_parsed_subject;
    let title_apart = Some(FoundInstead::CopyApart(
        DbSequence::of_unrestricted(3),
        CopyApart::Title,
    ));
    let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
        .unwrap()
        .found_instead(title_apart.clone());
    assert!(subject.fetched(false, true, HitLoad::Loaded).is_ok());
    let error = subject.fetched(true, true, HitLoad::Loaded).unwrap_err();
    assert!(error.contains("database sequence 3 also holds"), "{error}");
    assert!(error.contains("taxids and their names for it"), "{error}");
    assert!(subject.check_tabular_ids(true).is_ok());
    assert!(subject.all_tabular_ids(true).is_ok());
    // Outfmt 0 where it shows a title that can depend on the source: an empty title with
    // the taxonomy database, a heading generated of a title ending in `]`.
    assert!(check_shown_parsed_subject(
        &subject,
        true,
        true,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP"
    )
    .is_ok());
    let titled = |titles: &[&str]| {
        let deflines = titles
            .iter()
            .enumerate()
            .map(|(line, title)| BlastDefline {
                title: Some(title.as_bytes().to_vec()),
                ids: vec![gb(&format!("AF00000{line}"), 1)],
                ..BlastDefline::default()
            })
            .collect();
        ParsedSubject::new(DeflineSet {
            deflines,
            utf8: false,
        })
        .unwrap()
        .found_instead(title_apart.clone())
    };
    let error = check_shown_parsed_subject(
        &titled(&[""]),
        false,
        true,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTN",
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("can make that sequence's title"), "{error}");
    assert!(check_shown_parsed_subject(
        &titled(&[""]),
        false,
        false,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTN"
    )
    .is_ok());
    let heading = titled(&["alpha, partial [Homo sapiens]", ""]);
    assert!(check_shown_parsed_subject(
        &heading,
        true,
        false,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP"
    )
    .is_err());
    let bracket = titled(&["alpha, partial [Homo sapiens]"]);
    assert!(check_shown_parsed_subject(
        &bracket,
        true,
        true,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP"
    )
    .is_ok());
    // A local id `Subject_N` that the title's first word replaces, of the target or of
    // another definition line: the word can depend on the source for a nucleotide and for a
    // protein title ending in `]`; refused where a printed id shows the word.
    use crate::report::tabular::PrintedIds;
    let (subject_ids, all_ids) = (
        PrintedIds {
            subject: true,
            all: false,
        },
        PrintedIds {
            subject: false,
            all: true,
        },
    );
    let local = || ParsedSeqId::Local(ObjectId::Str(b"Subject_1".to_vec()));
    let replaced = ParsedSubject::new(one_defline(vec![local()]))
        .unwrap()
        .found_instead(title_apart.clone());
    assert!(replaced.check_tabular_ids(false).is_ok());
    let error = replaced.check_printed_word(false, subject_ids).unwrap_err();
    assert!(
        error.contains("whose first word replaces a printed id"),
        "{error}"
    );
    assert!(replaced.check_printed_word(false, all_ids).is_err());
    assert!(replaced
        .check_printed_word(false, PrintedIds::default())
        .is_ok());
    assert!(replaced.check_printed_word(true, subject_ids).is_ok());
    let second = ParsedSubject::new(DeflineSet {
        deflines: vec![
            one_defline(vec![gb("AAA00001", 1)]).deflines.remove(0),
            one_defline(vec![local()]).deflines.remove(0),
        ],
        utf8: false,
    })
    .unwrap()
    .found_instead(title_apart.clone());
    assert!(second.check_printed_word(false, subject_ids).is_ok());
    assert!(second.check_printed_word(false, all_ids).is_err());
    assert!(second.check_printed_word(true, all_ids).is_ok());
    // A replaced id beside a better one is not printed (`GetSeqIdListString`).
    let beside = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1), local()]))
        .unwrap()
        .found_instead(title_apart.clone());
    assert!(beside.check_printed_word(false, subject_ids).is_ok());
    assert!(beside.check_printed_word(false, all_ids).is_ok());
    let bracket = ParsedSubject::new(DeflineSet {
        deflines: vec![BlastDefline {
            title: Some(b"alpha, partial [Homo sapiens]".to_vec()),
            ids: vec![local()],
            ..BlastDefline::default()
        }],
        utf8: false,
    })
    .unwrap()
    .found_instead(title_apart);
    assert!(bracket.check_printed_word(true, subject_ids).is_err());
    // The other copies are refused as before: the taxonomy texts by every report that reads
    // the Bioseq, the taxids by those that write them.
    let texts = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
        .unwrap()
        .found_instead(Some(FoundInstead::CopyApart(
            DbSequence::of_unrestricted(3),
            CopyApart::TaxonomyTexts,
        )));
    assert!(texts.fetched(false, true, HitLoad::Loaded).is_err());
    assert!(texts.generated_title_fetched().is_ok());
    let taxids = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
        .unwrap()
        .found_instead(Some(FoundInstead::OtherTaxids(
            DbSequence::of_unrestricted(3),
        )));
    assert!(taxids.fetched(false, true, HitLoad::Loaded).is_ok());
    assert!(taxids.generated_title_fetched().is_ok());
}

// Line DBID N1b L1 probe: of a nucleotide `lcl|Subject_1` without a title and with taxid 9606,
// NCBI's tabular id is the first word of the title it makes from the source where the
// taxonomy database names the taxid (`Homo`); a protein's is `unnamed` either way (probe of
// line RG T5, `T/t5/out/ncbi/s.*.ids`: `Homo`, `Uncultured`, `unnamed`). A source that
// LOSAT refuses (a lookup that crashes NCBI, a named RefSeq record) refuses the word.
#[test]
fn a_tabular_word_of_a_title_made_from_the_source_names_the_organism() {
    use crate::report::tabular::PrintedIds;
    let subject = |first: &str, taxid: Option<i32>, source: SourceTaxname| {
        ParsedSubject::new(DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![match first.strip_prefix("ref|") {
                    Some(accession) => ParsedSeqId::Other(crate::seqid::model::TextseqId {
                        accession: Some(accession.as_bytes().to_vec()),
                        version: Some(1),
                        ..Default::default()
                    }),
                    None => ParsedSeqId::Local(ObjectId::Str(first.as_bytes().to_vec())),
                }],
                taxid,
                ..BlastDefline::default()
            }],
            utf8: false,
        })
        .unwrap()
        .with_source(true, |taxid| {
            assert_eq!(taxid, 9606);
            source
        })
    };
    let named = || SourceTaxname::Named(std::sync::Arc::from(&b"uncultured bacterium"[..]));
    let printed = PrintedIds::of(None);
    let local = subject("Subject_1", Some(9606), named());
    assert!(local.check_printed_word(false, printed).is_ok());
    let [id, _, _] = local.tabular_ids(false).unwrap().unwrap();
    assert_eq!(text(&id), "Uncultured");
    let [id, _, _] = local.tabular_ids(true).unwrap().unwrap();
    assert_eq!(text(&id), "unnamed");
    assert_eq!(
        text(&local.bioseq_title(true, true).unwrap().text),
        "unnamed protein product [uncultured bacterium]"
    );
    assert_eq!(
        text(&local.bioseq_title(false, false).unwrap().text),
        "Uncultured bacterium"
    );
    let unnamed = subject("Subject_1", Some(9606), SourceTaxname::Unnamed);
    assert!(unnamed.tabular_ids(false).is_ok());
    assert_eq!(text(&unnamed.bioseq_title(false, false).unwrap().text), "");
    let failed = subject(
        "Subject_1",
        Some(9606),
        SourceTaxname::Refused("a failed lookup"),
    );
    assert_eq!(failed.tabular_ids(false).unwrap_err(), "a failed lookup");
    assert_eq!(failed.tabular_ids(true).unwrap_err(), "a failed lookup");
    // Without a taxid the Bioseq has no source: nothing is looked up.
    let none = ParsedSubject::new(DeflineSet {
        deflines: vec![BlastDefline {
            ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Subject_1".to_vec()))],
            ..BlastDefline::default()
        }],
        utf8: false,
    })
    .unwrap()
    .with_source(true, |_| panic!("no source to name"));
    assert_eq!(none.source(), &SourceTaxname::Unnamed);
    // A named RefSeq chromosome gets the name (`x_SetTitleFromNC` makes nothing of a database
    // Bioseq), a protein's not upper-cased and not `unnamed protein product`; an unnamed one
    // keeps its empty title.
    let refseq = subject("ref|NC_000001", Some(9606), named());
    assert!(matches!(refseq.source(), SourceTaxname::Named(_)));
    assert_eq!(
        text(&refseq.bioseq_title(false, false).unwrap().text),
        "Uncultured bacterium"
    );
    assert_eq!(
        text(&refseq.bioseq_title(true, false).unwrap().text),
        "uncultured bacterium"
    );
    let refseq = subject("ref|NC_000001", Some(9606), SourceTaxname::Unnamed);
    assert_eq!(refseq.source(), &SourceTaxname::Unnamed);
}

// NCBI writes `Sequence with id AB000010 no longer exists in database...alignment skipped` in
// place of each alignment of a sequence whose heading's generated title it cannot decode (probe
// of line RG T5, `pt.100009.*.f0`); a title it can decode, or a stored one, is shown.
#[test]
fn a_heading_that_does_not_decode_names_the_seq_align_id() {
    let subject = |title: Option<&[u8]>, name: &[u8]| {
        ParsedSubject::new(DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![gb("AB000010", 1)],
                title: title.map(<[u8]>::to_vec),
                taxid: Some(100009),
                ..BlastDefline::default()
            }],
            utf8: false,
        })
        .unwrap()
        .with_source(false, |_| SourceTaxname::Named(std::sync::Arc::from(name)))
    };
    let bad = subject(None, b"Bad \x81 name");
    assert_eq!(
        bad.heading_skipped_id(false, false).as_deref(),
        Some(&b"AB000010"[..])
    );
    // The heading's own making of the title gives the same id (one GenerateDefline).
    assert_eq!(
        bad.heading(false, false, None).unwrap().unwrap_err(),
        bad.heading_skipped_id(false, false).unwrap()
    );
    assert_eq!(text(&bad.description(false).1.text), "Unknown");
    assert_eq!(
        subject(None, b"Homo sapiens").heading_skipped_id(false, false),
        None
    );
    assert_eq!(
        text(
            &subject(None, b"Homo sapiens")
                .heading(false, false, None)
                .unwrap()
                .unwrap()
        ),
        ">AB000010.1 Homo sapiens\n"
    );
    assert_eq!(
        subject(Some(b"stored"), b"Bad \x81 name").heading_skipped_id(true, false),
        None
    );
    assert!(subject(Some(b"stored"), b"Bad \x81 name")
        .heading(true, false, None)
        .unwrap()
        .is_ok());
}

#[test]
fn pdb_ids_without_a_chain_match_any_chain_of_the_molecule() {
    let pdb = |chain_id: Option<&str>| {
        let id = ParsedSeqId::Pdb(PdbSeqId {
            mol: b"1ABC".to_vec(),
            chain: chain_id.map(|id| i32::from(id.as_bytes()[0])),
            rel: None,
            chain_id: chain_id.map(|id| id.as_bytes().to_vec()),
        });
        ParsedSubject::new(one_defline(vec![id])).unwrap()
    };
    assert!(joined_in_reports(&pdb(None), &pdb(Some("AA"))));
    assert!(!joined_in_reports(&pdb(Some("A")), &pdb(Some("B"))));
}

// Line DBID N3: the copies of one record in the databases of an alias (three of r1, two of
// r3, the b ones with taxid 9606 in pjoin_abt) are joined; r6 (r1's letters under another id)
// and a Seq-align id that matches another one's without being its copy (pfound_v5's
// gb|DDD00004| and gb|DDD00004.1|, pdb|2XYZ| and pdb|2XYZ|AA) are not.
#[test]
fn copies_of_one_record_are_joined_and_other_matching_ids_are_not() {
    let open = |name: &str| {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
        (0..db.num_seqs())
            .map(|ordinal| ParsedSubject::new(db.defline_set(ordinal).unwrap()).unwrap())
            .collect::<Vec<_>>()
    };
    // a: r1 r2 r3, b: r1 r3 r4 r6, c: r1.
    for name in ["pjoin_abc", "pjoin_abt"] {
        let s = open(name);
        assert_eq!(s.len(), 8);
        for (earlier, later) in [(0, 3), (3, 0), (3, 7), (7, 0), (2, 4)] {
            assert!(
                joins_as_copy(&s[earlier], &s[later]),
                "{name} {earlier} {later}"
            );
        }
        for (earlier, later) in [(7, 6), (6, 3), (0, 1), (2, 5)] {
            assert!(
                !joined_in_reports(&s[earlier], &s[later]),
                "{name} {earlier} {later}"
            );
            assert!(
                !joins_as_copy(&s[earlier], &s[later]),
                "{name} {earlier} {later}"
            );
        }
    }
    let s = open("pfound_v5");
    for prefix in ["gb|DDD00004", "pdb|2XYZ|"] {
        let pair: Vec<&ParsedSubject> = s
            .iter()
            .filter(|subject| subject.alignment_fasta().starts_with(prefix.as_bytes()))
            .collect();
        assert_eq!(pair.len(), 2, "{prefix}");
        assert!(joined_in_reports(pair[0], pair[1]), "{prefix}");
        assert!(!joins_as_copy(pair[0], pair[1]), "{prefix}");
    }
    let parsed: Vec<Option<std::sync::Arc<ParsedSubject>>> = open("pjoin_abc")
        .into_iter()
        .map(|subject| Some(std::sync::Arc::new(subject)))
        .chain([None])
        .collect();
    assert!(rank_joins_as_copy(&parsed, 3, 0));
    assert!(!rank_joins_as_copy(&parsed, 0, 8));
    assert!(!rank_joins_as_copy(&parsed, 8, 0));
}

#[test]
fn refseq_accessions_whose_empty_titles_come_from_the_source() {
    use RefseqKind::{ChromosomeOrNcrna, Mrna};
    for (accession, kind) in [
        ("NC_000001", Some(ChromosomeOrNcrna)),
        ("AC_000001", Some(ChromosomeOrNcrna)),
        ("NM_000001", Some(Mrna)),
        ("NM_000000001", Some(Mrna)),
        ("XM_000000001", Some(Mrna)),
        ("NR_000001", Some(ChromosomeOrNcrna)),
        ("NR_000000001", None),
        ("NC_0000001", None),
        ("NP_000001", None),
        ("WP_000000001", None),
        ("XR_000001", None),
    ] {
        assert_eq!(
            refseq_accession_kind(accession.as_bytes()),
            kind,
            "{accession}"
        );
    }
}
/// The subject of a `-subject` record read with `-parse_deflines` (`ParsedSubject::of_record`).
fn record_subject(defline: &str, protein: bool) -> ParsedSubject {
    record_subject_of(defline.as_bytes(), protein)
}

/// `record_subject` of a defline of any bytes.
fn record_subject_of(defline: &[u8], protein: bool) -> ParsedSubject {
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
    let mut fasta = defline.to_vec();
    fasta.extend_from_slice(if protein { b"\nMKV\n" } else { b"\nACGT\n" });
    let config = ReaderConfig::subject("BLASTN", protein, false).with_parse_deflines(true);
    let mut source = FastaInputSource::from_bytes(&fasta, config);
    let record = read_all(&mut source, &mut |_| Ok(())).unwrap().remove(0);
    ParsedSubject::of_record(record.believed().unwrap(), &record.title)
}

// NCBI BLAST+ 2.17.0 -parse_deflines -subject (probe `T/t7/review` of line RG, `C.lclS9`,
// `J.pS9tpa`): a parsed `lcl|Subject_9` keeps its id (`s_ReplaceLocalId` with `parse_local`),
// and `GetSeqIdList` puts the first word of the title generated from the rest of the defline
// in its place: `s7` of `TPA: s7 third ~`; where that title is in no encoding NCBI guesses the
// lookup throws (`sseqid` Unknown, `sall*` the Seq-align id `Subject_9`), a record's
// definition line is never looked up.
#[test]
fn parsed_subject_ids_of_a_title_ncbi_cannot_decode_fail_the_lookup() {
    let decodable = record_subject(">lcl|Subject_9 TPA: s7 third ~", false);
    assert_eq!(
        decodable.tabular_ids(false).unwrap().unwrap()[2],
        b"s7".to_vec()
    );
    let bad = record_subject_of(b">lcl|Subject_9 bad x\x81y title", false);
    assert_eq!(bad.tabular_ids(false).unwrap(), None);
    assert_eq!(bad.tabular_saccver(false).unwrap(), None);
    assert!(bad.check_tabular_ids(false).is_ok());
    assert_eq!(
        bad.lookup_failure(false).unwrap(),
        Some(crate::report::tabular::LookupFailure {
            subject_id: true,
            all_ids: false,
            bioseq: false,
            length: false,
        })
    );
    assert_eq!(
        bad.alignment_tabular_ids(),
        [b"Subject_9".to_vec(), b"Subject_9".to_vec()]
    );
    assert_eq!(bad.skipped_id(), b"Subject_9");
    let all = bad.all_tabular_ids(false).unwrap();
    assert!(!all.failed);
    assert_eq!(text(&all.ids[0].id), "Subject_9");
    // A record of a decodable title has no failure.
    assert_eq!(decodable.lookup_failure(false).unwrap(), None);
}

// NCBI BLAST+ 2.17.0 -parse_deflines -subject (probe `DI7/ws/pd_sub.py`, 2026-11-01): `sseqid`,
// `sacc`, `saccver`, `sallseqid`, the outfmt 0 heading and description of a record.
#[test]
fn records_read_with_parse_deflines_follow_the_ncbi_oracle() {
    for (defline, protein, ids, all, heading, local) in [
        (
            ">NM_000003.1 x y",
            false,
            ["ref|NM_000003.1|", "NM_000003", "NM_000003.1"],
            "ref|NM_000003.1|",
            ">NM_000003.1 x y\n",
            false,
        ),
        (
            ">plain1 a",
            false,
            ["plain1", "plain1", "plain1"],
            "plain1",
            ">plain1 a\n",
            false,
        ),
        (
            ">tpg|BK000001.1| y",
            false,
            ["tpg|BK000001.1|", "BK000001", "BK000001.1"],
            "tpg|BK000001.1|",
            ">BK000001.1 TPA: y\n",
            false,
        ),
        (">s1", false, ["s1", "s1", "s1"], "s1", ">s1 \n", false),
        (
            ">lcl|Subject_1 x y",
            false,
            ["x", "x", "x"],
            "Subject_1",
            "",
            true,
        ),
        (
            ">",
            false,
            ["Subject_1", "Subject_1", "Subject_1"],
            "Subject_1",
            "",
            true,
        ),
        (
            ">",
            true,
            ["unnamed", "unnamed", "unnamed"],
            "Subject_1",
            "",
            true,
        ),
    ] {
        let subject = record_subject(defline, protein);
        let shown = subject.tabular_ids(protein).unwrap().unwrap();
        assert_eq!(
            shown.map(|id| String::from_utf8(id).unwrap()),
            ids,
            "{defline}"
        );
        let all_ids = subject.all_tabular_ids(protein).unwrap().ids;
        assert_eq!(text(&all_ids[0].id), all, "{defline}");
        assert_eq!(subject.shown_as_record(), local, "{defline}");
        if !local {
            assert_eq!(
                text(&subject.heading(protein, false, None).unwrap().unwrap()),
                heading,
                "{defline}"
            );
            assert!(subject.first_title().is_empty() || defline.contains(' '));
        }
    }
}

// NCBI BLAST+ 2.17.0 `blastp -num_threads 1` on the version 5 database `psfx` of line RG
// T5c (`tests/db_search/t5/make_t5.py --suffix-only`; rows `db.t5.sfx.blastp.*`): the heading
// of an untitled line, outfmt 18 and the tabular word of `lcl|Subject_9nn` show the target's
// title with its end rewritten from the first source (`x_AdjustProteinTitleSuffix`); the
// description (`fLeavePrefixSuffix`) keeps the stored end. Only those subjects ask for the
// source's name.
#[test]
fn a_protein_title_ending_in_a_bracket_gets_its_end_from_the_source() {
    use crate::algorithm::common::database::check_shown_parsed_subject;
    let set = |lines: &[(ParsedSeqId, &str, Option<i32>)]| DeflineSet {
        deflines: lines
            .iter()
            .map(|(id, title, taxid)| BlastDefline {
                ids: vec![id.clone()],
                title: (!title.is_empty()).then(|| title.as_bytes().to_vec()),
                taxid: *taxid,
                ..BlastDefline::default()
            })
            .collect(),
        utf8: false,
    };
    let bacillus = || SourceTaxname::Named(std::sync::Arc::from(&b"Bacillus sp."[..]));
    let named = |lines: &[(ParsedSeqId, &str, Option<i32>)], protein| {
        ParsedSubject::new(set(lines))
            .unwrap()
            .with_source(protein, |taxid| {
                assert_eq!(taxid, 100002);
                bacillus()
            })
    };
    // The first source comes from the second line.
    let heading = named(
        &[
            (gb("AC400001", 1), "alpha, partial [Bacillus sp.]", None),
            (gb("AD400001", 1), "", Some(100002)),
        ],
        true,
    );
    assert_eq!(heading.source(), &bacillus());
    assert_eq!(
        text(&heading.heading(true, false, None).unwrap().unwrap()),
        ">AC400001.1 alpha, partial [Bacillus sp.]\n AD400001.1 alpha [Bacillus sp.]\n"
    );
    assert_eq!(
        text(&heading.description(true).1.text),
        "alpha, partial [Bacillus sp.]"
    );
    assert!(check_shown_parsed_subject(
        &heading,
        true,
        true,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP"
    )
    .is_ok());
    let deflines = heading.tax_deflines(true);
    assert!(matches!(
        &deflines[1].title,
        crate::report::tax_report::TaxTitle::Text(title) if title == b"alpha [Bacillus sp.]"
    ));
    // A nucleotide is not adjusted, and asks nothing.
    let nucleotide = ParsedSubject::new(set(&[
        (gb("AC400001", 1), "alpha, partial [Bacillus sp.]", None),
        (gb("AD400001", 1), "", Some(100002)),
    ]))
    .unwrap()
    .with_source(false, |_| panic!("a nucleotide title is not adjusted"));
    assert_eq!(
        text(&nucleotide.heading(false, false, None).unwrap().unwrap()),
        ">AC400001.1 alpha, partial [Bacillus sp.]\n AD400001.1 alpha, partial [Bacillus sp.]\n"
    );
    // The tabular word of `lcl|Subject_9nn`.
    let local = |tag: &str| ParsedSeqId::Local(ObjectId::Str(tag.as_bytes().to_vec()));
    let word = named(
        &[(
            local("Subject_901"),
            "alpha, partial [Bacillus sp.]",
            Some(100002),
        )],
        true,
    );
    let [id, accession, version] = word.tabular_ids(true).unwrap().unwrap();
    assert_eq!([text(&id), text(&accession), text(&version)], ["alpha"; 3]);
    // Without a name (no taxonomy database) a ` [ ]` end goes.
    let nameless = ParsedSubject::new(set(&[(local("Subject_903"), "a, [ ]", Some(77777))]))
        .unwrap()
        .with_source(true, |_| SourceTaxname::Unnamed);
    assert_eq!(text(&nameless.tabular_ids(true).unwrap().unwrap()[0]), "a");
    // A title whose generated form no report shows (no untitled line, no replaced id), a
    // title not ending in `]` and a Bioseq without a source ask nothing.
    for lines in [
        vec![(
            gb("AC400002", 1),
            "alpha, partial [Bacillus sp.]",
            Some(100002),
        )],
        vec![
            (gb("AC400003", 1), "alpha protein", Some(100002)),
            (gb("AD400003", 1), "", Some(100002)),
        ],
        vec![
            (gb("AC400031", 1), "gg, partial [Bacillus sp.]", None),
            (gb("AD400031", 1), "", None),
        ],
    ] {
        let subject = ParsedSubject::new(set(&lines))
            .unwrap()
            .with_source(true, |_| panic!("no name is asked for {lines:?}"));
        assert_eq!(subject.source(), &SourceTaxname::Unnamed);
    }
    // An adjusted title of punctuation is read past its end (NCBI BLAST+ 2.17.0 crashes:
    // probe of line RG T5c, `pedge_e3`): refused.
    let edge = ParsedSubject::new(set(&[
        (gb("AG400003", 1), ", , [ ]", Some(77777)),
        (gb("AH400003", 1), "", None),
    ]))
    .unwrap()
    .with_source(true, |_| SourceTaxname::Unnamed);
    let error = check_shown_parsed_subject(
        &edge,
        true,
        true,
        || crate::seqdb::Oid::of_unrestricted(0),
        "BLASTP",
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("x_CleanAndCompress reads past the end"),
        "{error}"
    );
}

/// A lookup whose file cannot be read (line WQA D18 review r2 L1, 2026-11-21): NCBI crashes only
/// where it reads bytes of a file it has mapped and the system no longer gives them (a later read,
/// `SeqDbError::Read`: SIGBUS). A first byte that cannot be read where the file is opened
/// (`SeqDbError::Map`) is NCBI's failed mapping, an error it reports, also where LOSAT has no
/// text for it (a code without NCBI's text on WASI, Windows): LOSAT's text, as a file that does
/// not open.
#[test]
fn only_a_later_read_of_a_lookup_file_is_ncbis_crash() {
    use crate::seqdb::SeqDbError;
    let path = std::path::PathBuf::from("/d/db.psi");
    let later = LookupFailure::of_error(&SeqDbError::Read {
        path: path.clone(),
        source: std::io::Error::from_raw_os_error(5),
    });
    assert_eq!(later.what_ncbi_does(None), LookupFailure::NCBI_CRASHES);
    assert_eq!(later.ncbi_message(), None);
    let without_text = |error: SeqDbError| {
        let failure = LookupFailure::of_error(&error);
        assert_eq!(
            failure.what_ncbi_does(None),
            LookupFailure::NCBI_READS,
            "{error}"
        );
        assert_eq!(failure.ncbi_message(), None, "{error}");
        assert_eq!(failure.to_string(), error.to_string());
    };
    without_text(SeqDbError::Map {
        path: path.clone(),
        len: 64,
        source: std::io::Error::other("no system code"),
    });
    without_text(SeqDbError::Open {
        path,
        source: std::io::Error::other("no system code"),
    });
}

// DBID pfoff-3 (DECISIONS 62 (13); review 2 of the item, L7): where a pairwise report ends the
// run is decided on the alignments that it shows, after the sort and the cut: NCBI's probe
// `R/pfoff/p3` `bp.U1.b1.dflt.0` ends (exit 2) where `bp.U1.b1.sh4.0`, whose `-sorthits 4` puts
// a loaded subject first, writes skipped lines (exit 0), and `bp.U2.b1.dflt.3` (outfmt 3) ends
// where `bp.U2.b1.dflt.0` does not; a local Seq-align id is not looked up there, and
// `-num_alignments 0` shows none. The model of the loads (`HitLoads::of` of the reads that
// `QueryGroup::read_order` makes before the rows) and the writer (`aligned_end` on its group
// and the subjects' own lookups, `lookup_throws`, not the model's loads; review 3 of the item,
// L-4) end at the same alignment: the subjects aligned before it load, it throws, those after it
// are never read (with `-num_descriptions 0`, no description line reads them), and no report
// after it runs (`HitLoads::reaches`, N1).
#[test]
fn the_model_and_the_writer_end_the_run_at_the_same_alignment() {
    use crate::report::hit_sort::ReportSort;
    use crate::report::x_pairwise_par::{aligned_end, writer_reads, QueryGroup};
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    let local = ParsedSubject::new(one_defline(vec![ParsedSeqId::Local(ObjectId::Str(
        b"ENST1".to_vec(),
    ))]))
    .unwrap();
    // Subject 0 throws and comes first by E-value; 1 loads and comes first by bit score; 2
    // throws. With `local_first`, 0's Seq-align id is a local one.
    let parsed = |local_first: bool| {
        vec![
            if local_first {
                local.clone().found_instead(throws())
            } else {
                plain("AAA00001").found_instead(throws())
            },
            plain("AAA00002"),
            plain("AAA00003").found_instead(throws()),
        ]
    };
    let rows: Vec<crate::report::pairwise::PairwiseHit> =
        [(0, 30.0, 1e-20), (1, 90.0, 1e-10), (2, 20.0, 1e-5)]
            .into_iter()
            .map(|(s_idx, bit_score, e_value)| {
                crate::report::pairwise::PairwiseHit::from(crate::common::Hit {
                    identity: 100.0,
                    length: 10,
                    mismatch: 0,
                    gapopen: 0,
                    q_start: 1,
                    q_end: 10,
                    s_start: 1,
                    s_end: 10,
                    e_value,
                    bit_score,
                    num_ident: 10,
                    query_frame: 1,
                    query_length: 0,
                    q_idx: 0,
                    s_idx,
                    raw_score: 25,
                    sort_query_offset: 0,
                    sort_query_end: 0,
                    sort_subject_offset: 0,
                    sort_subject_end: 0,
                    has_sort_offsets: false,
                    gap_info: None,
                    num_positives: 10,
                })
            })
            .collect();
    // Where the model and the writer end the run, as the position among the shown alignments.
    let end = |after: AfterThrow, hits: Option<u8>, alignments: usize, local_first: bool| {
        let parsed = parsed(local_first);
        // What the writer's lookups throw, from the subjects themselves.
        let lookups: Vec<Option<crate::ncbi_exception::NcbiException>> =
            parsed.iter().map(ParsedSubject::lookup_throws).collect();
        let reading = crate::report::pairwise::PairwiseReading {
            sort: ReportSort { hits, hsps: None },
            shown: Shown {
                descriptions: 0,
                alignments,
            },
            html: false,
        };
        let reads = QueryGroup::read_order(
            rows.iter().map(|row| (row.hit.s_idx, row)),
            |row| crate::report::hit_sort::HspScores::of(row),
            |_, _| false,
            reading,
            false,
            None,
        );
        let run = RunLoads::new(Some(RunThrow::every_row(after)));
        let loads = run.loads(
            [(true, ReportKeys::Pairwise(reads))],
            |key| parsed.get(key),
            |_, _| 0,
        );
        let writer = writer_reads(&rows, |_, _| false, reading, false, None);
        let aligned: Vec<u32> = (0..writer.aligned)
            .map(|hit| writer.hits[hit][0] as u32)
            .collect();
        let writer_end = aligned_end(after, &aligned, |s_idx| {
            lookups[s_idx as usize]
                .as_ref()
                .map(|exception| (exception, parsed[s_idx as usize].alignment_id_is_local()))
        })
        .map(|(at, _)| at);
        assert_eq!(run.ended(), writer_end.is_some());
        // The report after it runs where the run did not end there.
        assert_eq!(loads.reaches(1), writer_end.is_none());
        for (position, &s_idx) in aligned.iter().enumerate() {
            let expected = match writer_end {
                Some(at) if position > at => HitLoad::Unread,
                _ if parsed[s_idx as usize].lookup_throws().is_some() => HitLoad::Throws,
                _ => HitLoad::Loaded,
            };
            assert_eq!(loads.load(s_idx as usize), expected, "{s_idx}");
        }
        writer_end
    };
    let every = usize::MAX;
    // Outfmt 0 ends at its first alignment where it throws, not after `-sorthits` moves a
    // loaded one first, nor with no alignment shown.
    assert_eq!(end(AfterThrow::SkipsAlignment, None, every, false), Some(0));
    assert_eq!(end(AfterThrow::SkipsAlignment, Some(1), every, false), None);
    assert_eq!(end(AfterThrow::SkipsAlignment, None, 0, false), None);
    assert_eq!(end(AfterThrow::SkipsAlignment, None, every, true), None);
    // Outfmt 1 to 4 end at the first shown alignment that throws, a local one too, unless
    // the cut leaves none.
    assert_eq!(
        end(AfterThrow::EndsAtAlignment, Some(1), every, false),
        Some(1)
    );
    assert_eq!(end(AfterThrow::EndsAtAlignment, None, every, true), Some(0));
    assert_eq!(end(AfterThrow::EndsAtAlignment, Some(1), 1, false), None);
    for other in [
        AfterThrow::EndsRun,
        AfterThrow::EndsReport,
        AfterThrow::WritesRow,
    ] {
        assert!(!other.is_pairwise());
        assert_eq!(other.alignment_end([0], |_| true, |_| false), None);
    }
}

// DBID pfoff-3: the description lines load each shown subject in order (one that throws
// loads nothing), then the alignments load theirs up to the one that ends the run; after the
// run's end nothing loads, in a later report either; a batch that prefetches loads all; a
// subject that neither the description lines nor the alignments show is never read (review 1
// of pfoff-3, H1: psiblast's later round shows the subjects of its first HSPs).
#[test]
fn a_pairwise_report_loads_its_descriptions_then_its_alignments() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 1 and 3 throw.
    let parsed: Vec<ParsedSubject> = vec![
        plain("AAA00001"),
        plain("AAA00002").found_instead(throws()),
        plain("AAA00003"),
        plain("AAA00004").found_instead(throws()),
        plain("AAA00005"),
        plain("AAA00006"),
    ];
    let subject = |key: usize| parsed.get(key);
    let none = |_: usize, _: usize| 0;
    let run = |after| Some(RunThrow::every_row(after));
    // A report of one subject per hit whose description lines show its first `descriptions`
    // hits and whose alignments its first `alignments`; a batch that prefetches gives its
    // subjects alone.
    let report = |one_by_one: bool, keys: &[usize], descriptions: usize, alignments: usize| {
        let keys = if one_by_one {
            ReportKeys::Pairwise(PairwiseReads {
                hits: keys.iter().map(|&key| vec![key]).collect(),
                described: (0..keys.len().min(descriptions)).collect(),
                aligned: keys.len().min(alignments),
                html: false,
            })
        } else {
            ReportKeys::Entries(keys.iter().map(|&key| Some(key)).collect())
        };
        (one_by_one, keys)
    };
    let loads = |run, reports: Vec<(bool, ReportKeys)>| {
        let loads = HitLoads::of(run, reports, subject, none);
        (0..6).map(|key| loads.load(key)).collect::<Vec<_>>()
    };
    use HitLoad::{Loaded as L, Throws as T, Unread as U};
    let every = usize::MAX;
    let three = |descriptions, alignments, keys: [&[usize]; 3]| {
        keys.iter()
            .map(|keys| report(true, keys, descriptions, alignments))
            .collect::<Vec<_>>()
    };
    // Outfmt 0: the second report's first alignment throws; its description lines loaded 4
    // before; the third report loads nothing.
    let ends_second: [&[usize]; 3] = [&[0, 1, 2], &[3, 4], &[5]];
    assert_eq!(
        loads(
            run(AfterThrow::SkipsAlignment),
            three(every, every, ends_second)
        ),
        [L, T, L, T, L, U]
    );
    // A first alignment that loads: the throwing ones are skipped and the run goes on.
    let skips: [&[usize]; 3] = [&[0, 1, 2], &[4, 3], &[5]];
    assert_eq!(
        loads(run(AfterThrow::SkipsAlignment), three(every, every, skips)),
        [L, T, L, T, L, L]
    );
    // `-num_alignments 0`: no alignment section, nothing ends.
    assert_eq!(
        loads(
            run(AfterThrow::SkipsAlignment),
            three(every, 0, ends_second)
        ),
        [L, T, L, T, L, L]
    );
    // Outfmt 1 to 4: the first report's second alignment ends the run.
    assert_eq!(
        loads(run(AfterThrow::EndsAtAlignment), three(every, every, skips)),
        [L, T, L, U, U, U]
    );
    // `-num_alignments 1`: only the first alignment is shown, which loads; 3 throws in the
    // description lines alone.
    assert_eq!(
        loads(run(AfterThrow::EndsAtAlignment), three(every, 1, skips)),
        [L, T, L, T, L, L]
    );
    // `-num_descriptions 0`: only the alignments up to the end load, and 2 is never reached.
    assert_eq!(
        loads(run(AfterThrow::EndsAtAlignment), three(0, every, skips)),
        [L, T, U, U, U, U]
    );
    // The description lines that a report shows are its writer's (`PairwiseReads::described`,
    // in the order of the lines): psiblast's later round shows those of its first HSPs, here
    // 4 alone (the part of the subjects found before) then 2, so 0, 5, 1 and 3, which neither
    // the lines nor the alignments (`-num_alignments 0`) show, are never read (1 and 3, whose
    // lookups throw, throw nothing).
    let two_part = (
        true,
        ReportKeys::Pairwise(PairwiseReads {
            hits: vec![vec![0], vec![2], vec![4], vec![5], vec![1], vec![3]],
            described: vec![2, 1],
            aligned: 0,
            html: false,
        }),
    );
    assert_eq!(
        loads(run(AfterThrow::SkipsAlignment), vec![two_part]),
        [U, U, L, U, L, U]
    );
    // A batch that prefetches loads every subject it shows, where `-num_descriptions 0
    // -num_alignments 1` reads the first alignment alone of a batch that does not.
    assert_eq!(
        loads(
            run(AfterThrow::SkipsAlignment),
            vec![
                report(false, &[0, 2, 4], 0, 1),
                report(true, &[5, 1, 3], 0, 1)
            ]
        ),
        [L, U, L, U, L, L]
    );
    // After the run's end a batch that prefetches loads nothing either.
    assert_eq!(
        loads(
            run(AfterThrow::SkipsAlignment),
            vec![
                report(true, &[1, 0], every, every),
                report(false, &[2, 3, 4, 5], every, every)
            ]
        ),
        [L, T, U, U, U, U]
    );
    // A hit that shows copies of its record: a thrown first subject loads none of them (1 and
    // 2, in hits that are not shown, are not read).
    let copies = (
        true,
        ReportKeys::Pairwise(PairwiseReads {
            hits: vec![vec![4, 5], vec![3, 0], vec![1], vec![2]],
            described: vec![0, 1],
            aligned: 0,
            html: false,
        }),
    );
    assert_eq!(
        loads(run(AfterThrow::SkipsAlignment), vec![copies]),
        [U, U, U, T, L, L]
    );
    // A later call's reports read nothing after the end (`RunLoads`).
    let calls = RunLoads::new(run(AfterThrow::SkipsAlignment));
    let first = calls.loads([report(true, &[3, 4], every, every)], subject, none);
    assert_eq!(first.load(3), HitLoad::Throws);
    assert!(calls.ended());
    let later = calls.loads([report(true, &[5], every, every)], subject, none);
    assert_eq!(later.load(5), HitLoad::Unread);
}

// Line RG T4b review b F2 (NCBI BLAST+ 2.17.0, probes `nb.pfo.bpe.na0.f{0,1,3}.rep.h` and
// `nb2.pfo.bpe.na0.f0.rep.p`): with `-html` the description line of a subject whose load throws
// ends the run (its `catch` makes the line's link, which loads the subject again outside any
// `catch`): the lines before it load, no later line, alignment, report or batch's prefetch loads.
// Without `-html` the line is `Unknown` and the run goes on, with `-num_alignments 0` into the
// next batch, whose prefetch loads every subject; with `-html` too where no shown line throws.
#[test]
fn an_html_description_line_whose_load_throws_ends_the_run() {
    let throws = || {
        Some(FoundInstead::Failed(LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )))
    };
    let plain = |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
    // 1 and 3 throw.
    let parsed: Vec<ParsedSubject> = vec![
        plain("AAA00001"),
        plain("AAA00002").found_instead(throws()),
        plain("AAA00003"),
        plain("AAA00004").found_instead(throws()),
        plain("AAA00005"),
        plain("AAA00006"),
    ];
    let subject = |key: usize| parsed.get(key);
    let none = |_: usize, _: usize| 0;
    let every = usize::MAX;
    // A batch whose prefetch is off reports 0 to 3 (its description lines show the first
    // `descriptions`, its alignments the first `alignments`), then a batch that prefetches 4, 5.
    let run = |after, descriptions: usize, alignments: usize, html| {
        let keys = [0, 1, 2, 3];
        let first = ReportKeys::Pairwise(PairwiseReads {
            hits: keys.iter().map(|&key| vec![key]).collect(),
            described: (0..keys.len().min(descriptions)).collect(),
            aligned: keys.len().min(alignments),
            html,
        });
        let second = ReportKeys::Entries(vec![Some(4), Some(5)]);
        let loads = HitLoads::of(
            Some(RunThrow::every_row(after)),
            [(true, first), (false, second)],
            subject,
            none,
        );
        let reaches_second = loads.reaches(1);
        (
            (0..6).map(|key| loads.load(key)).collect::<Vec<_>>(),
            reaches_second,
        )
    };
    use HitLoad::{Loaded as L, Throws as T, Unread as U};
    for after in [AfterThrow::SkipsAlignment, AfterThrow::EndsAtAlignment] {
        // `-num_alignments 0`: the plain lines of 1 and 3 are `Unknown`, the next batch is
        // prefetched; with `-html` the line of 1 ends the run.
        assert_eq!(
            run(after, every, 0, false),
            (vec![L, T, L, T, L, L], true),
            "{after:?}"
        );
        assert_eq!(
            run(after, every, 0, true),
            (vec![L, T, U, U, U, U], false),
            "{after:?}"
        );
        // With alignments too: the line ends it before them.
        assert_eq!(
            run(after, every, every, true),
            (vec![L, T, U, U, U, U], false),
            "{after:?}"
        );
        // No shown line throws (`-num_descriptions 1`, `-num_descriptions 0 -num_alignments 1`):
        // the run goes on as without `-html`.
        assert_eq!(
            run(after, 1, 0, true),
            (vec![L, U, U, U, L, L], true),
            "{after:?}"
        );
        assert_eq!(
            run(after, 0, 1, true),
            (vec![L, U, U, U, L, L], true),
            "{after:?}"
        );
    }
}

// DBID pfoff-3, NCBI 2.17.0 (probe `R/pfoff/p3`, `bp.U2.b1.dflt.0`, `bp.U2.b1.gis.0`): the
// description line of a subject whose load throws is the Seq-align id and `Unknown`
// (`x_GetDeflineInfo`'s catch), the id's FASTA form with `-show_gis` (no space), and its
// alignments' line names the accession without its version.
#[test]
fn a_thrown_description_is_the_alignment_id_and_unknown() {
    let refseq = ParsedSubject::new(one_defline(vec![ParsedSeqId::Other(
        crate::seqid::model::TextseqId {
            accession: Some(b"NP_100001".to_vec()),
            version: Some(1),
            ..Default::default()
        },
    )]))
    .unwrap();
    let line = |subject: &ParsedSubject| {
        let (mut start, title) = subject.thrown_description();
        start.extend_from_slice(&title.text);
        String::from_utf8(start).unwrap()
    };
    assert_eq!(line(&refseq), "NP_100001.1 Unknown");
    assert_eq!(
        line(&refseq.clone().showing_gis(true)),
        "ref|NP_100001.1|Unknown"
    );
    assert_eq!(refseq.skipped_id(), b"NP_100001");
    assert!(!refseq.alignment_id_is_local());
    // A GI Seq-align id is a reading of the model, not NCBI's bytes: no database makes a GI's
    // lookup throw (NCBI finds a GI by its ISAM index, not by the accession table whose
    // absence throws; review 1 of pfoff-3, L3), so no probe or fixture shows the line. What
    // can be tested is that it is made as the line of a loaded subject of that id
    // (`description_start`, which the probes check), then `Unknown`, after two spaces where
    // `-show_gis` wrote `gi|N|`.
    let gi = ParsedSubject::new(one_defline(vec![ParsedSeqId::Gi(1234)])).unwrap();
    for subject in [gi.clone(), gi.showing_gis(true)] {
        let (start, title) = subject.thrown_description();
        assert_eq!(start, subject.description_start());
        let gap: &[u8] = if subject.show_gis { b"  " } else { b"" };
        assert_eq!(title.text, [gap, b"Unknown"].concat());
    }
}

fn pdb_id(mol: &str, chain: Option<i32>, chain_id: Option<&str>) -> ParsedSeqId {
    ParsedSeqId::Pdb(PdbSeqId {
        mol: mol.as_bytes().to_vec(),
        chain,
        rel: None,
        chain_id: chain_id.map(|c| c.as_bytes().to_vec()),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:386-431
// The data loader looks up the Seq-id of the Seq-align id's handle and keeps the sequence only
// where one of its ids compares equal to it: a record whose only id is `pdb` with chain 321
// (`char` 65 for the handle) is not served, and NCBI BLAST+ 2.17.0 writes the `catch` of each
// read (line RG-KB item 12, `i12/ncbi`: `1AB1_A Unknown`, `Sequence with id 1AB1_A no longer
// exists in database...alignment skipped`).
#[test]
fn a_pdb_chain_that_the_handle_changes_finds_no_bioseq() {
    let unserved = lines(vec![vec![pdb_id("1AB1", Some(321), None)]]);
    assert!(unserved.scope_finds_none());
    assert_eq!(
        unserved.loader_id().into_owned(),
        pdb_id("1AB1", Some(65), None)
    );
    let (start, title, caught) = unserved.description_line(true);
    assert_eq!(
        (text(&start), text(&title.text), caught),
        ("1AB1_A ", "Unknown", true)
    );
    assert_eq!(
        unserved.heading(true, false, None).unwrap(),
        Err(b"1AB1_A".to_vec())
    );
    assert_eq!(
        unserved.heading_skipped_id(true, false).as_deref(),
        Some(&b"1AB1_A"[..])
    );
    assert!(matches!(unserved.prefetched(), Prefetched::Nothing));
    // The record is not refused for what the lookup finds: the formats write the `catch`.
    assert_eq!(unserved.fetched(true, false, HitLoad::Loaded), Ok(()));
    // Chains that the handle keeps, and a second id that the handle's Seq-id matches, are served.
    for served in [
        lines(vec![vec![pdb_id("1AB1", Some(65), None)]]),
        lines(vec![vec![pdb_id("1AB1", Some(-5), None)]]),
        lines(vec![vec![pdb_id("1AB1", None, None)]]),
        lines(vec![vec![pdb_id("1AB1", Some(66), Some("A"))]]),
        lines(vec![vec![
            gb("AB600001", 1),
            pdb_id("1AB1", Some(321), None),
        ]]),
    ] {
        assert!(!served.scope_finds_none(), "{:?}", served.alignment_id());
        assert!(!served.description_line(true).2);
    }
    // Where the handle's Seq-id is held by another definition line than the Seq-align id's, the
    // Bioseq NCBI builds is made of that line: its title is refused only where the description
    // line or the organism report shows it; the tables, the heading and the XML documents show
    // the sequence's own lines (review 1 of line RG-KB item 12, R2: probe `e15_twoline`).
    use crate::algorithm::common::database::description_line_refusal;
    let other_line = ParsedSubject::new(DeflineSet {
        deflines: vec![
            BlastDefline {
                title: Some(b"t".to_vec()),
                ids: vec![pdb_id("1AB1", Some(321), None)],
                ..BlastDefline::default()
            },
            BlastDefline {
                title: Some(b"u".to_vec()),
                ids: vec![pdb_id("1AB1", Some(65), Some("A"))],
                ..BlastDefline::default()
            },
        ],
        utf8: false,
    })
    .unwrap();
    assert!(!other_line.scope_finds_none());
    assert!(other_line
        .handle_title_refusal()
        .unwrap()
        .contains("a definition line other than the Seq-align id's"));
    assert!(description_line_refusal(&other_line, HitLoad::Loaded)
        .unwrap()
        .contains("description line and organism report"));
    assert_eq!(other_line.fetched(true, false, HitLoad::Loaded), Ok(()));
    assert!(matches!(other_line.prefetched(), Prefetched::Itself));
    // One line that holds both chains: NCBI's Bioseq reads the line's first pdb id, which is not
    // the handle's, so its load throws: refused where a format reads the Bioseq or the prefetch
    // reads it; the Seq-annots and the archive read nothing of it (outfmt 8 is NCBI's bytes).
    let one_line = ParsedSubject::new(DeflineSet {
        deflines: vec![BlastDefline {
            title: Some(b"t".to_vec()),
            ids: vec![
                pdb_id("1AB1", Some(321), None),
                pdb_id("1AB1", Some(65), None),
            ],
            ..BlastDefline::default()
        }],
        utf8: false,
    })
    .unwrap();
    assert!(one_line
        .fetched(false, false, HitLoad::Loaded)
        .unwrap_err()
        .contains("its load throws"));
    assert!(matches!(
        one_line.prefetched(),
        Prefetched::Failed(reason) if reason.contains("oid headers do not contain target gi/seq_id")
    ));
    assert!(one_line.pssm_read(false).is_err());
    assert_eq!(description_line_refusal(&one_line, HitLoad::Loaded), None);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:886-892
// The description line of a sequence whose Bioseq the scope does not get prints, with
// `-show_gis`, a GI that NCBI never sets; a description id longer than 65 bytes pads the line
// without end. Both are refused only for a sequence whose line a report prints (the first
// `-num_descriptions`; review 1 of line RG-KB item 12, R1).
#[test]
fn a_description_line_is_refused_only_where_a_report_prints_it() {
    use crate::algorithm::common::database::{description_line_refusal, refuse_described_lines};
    let unserved = lines(vec![vec![pdb_id("1AB1", Some(321), None)]]);
    assert_eq!(description_line_refusal(&unserved, HitLoad::Loaded), None);
    let unserved = unserved.showing_gis(true);
    assert!(description_line_refusal(&unserved, HitLoad::Loaded)
        .unwrap()
        .contains("prints a GI that NCBI BLAST+ 2.17.0 never sets"));
    let long = lines(vec![vec![ParsedSeqId::Local(ObjectId::Str(
        "x".repeat(70).into_bytes(),
    ))]]);
    assert!(description_line_refusal(&long, HitLoad::Loaded)
        .unwrap()
        .contains("longer than 65 bytes"));
    let served = lines(vec![vec![pdb_id("1AB4", Some(65), None)]]);
    let parsed = vec![
        Some(std::sync::Arc::new(served)),
        Some(std::sync::Arc::new(unserved)),
    ];
    let refuse = |described: Vec<u32>| {
        refuse_described_lines(
            1,
            |_| vec![0, 1],
            |_| Ok(described.clone()),
            &parsed,
            |_| HitLoad::Loaded,
            crate::seqdb::Oid::of_unrestricted,
            "BLASTP",
        )
    };
    // `-num_descriptions 1`: rank 1 prints no line.
    assert!(refuse(vec![0]).is_ok());
    assert!(refuse(vec![]).is_ok());
    let error = refuse(vec![0, 1]).unwrap_err().to_string();
    assert!(error.contains("never sets"), "{error}");
}

#[test]
fn the_scope_keys_of_a_pdb_id_are_its_handle() {
    // `chain` 321 and 65, with or without the `chain-id` "A", are one handle; "B" is another.
    let keys = |id: &ParsedSeqId| {
        let mut keys = Vec::new();
        crate::seqid::defline_set::scope_keys(id, &mut keys);
        keys
    };
    let a = keys(&pdb_id("1AB1", Some(65), Some("A")));
    assert_eq!(keys(&pdb_id("1AB1", Some(321), None)), a);
    assert_eq!(keys(&pdb_id("1AB1", Some(65), None)), a);
    assert_eq!(keys(&pdb_id("1ab1", Some(65), Some("A"))), a);
    assert_ne!(keys(&pdb_id("1AB1", Some(66), Some("B"))), a);
    assert_ne!(keys(&pdb_id("1AB1", Some(200), None)), a);
}

use super::*;

fn list(ids: &[&str]) -> TextMatchIds {
    let mut set = TextMatchIds::default();
    for id in ids {
        set.insert(id.as_bytes()).unwrap();
    }
    set
}

#[test]
fn split_gives_the_word_starts_and_the_end() {
    assert_eq!(split(b"gb|AB000001.1|"), vec![0, 3, 14]);
    assert_eq!(split(b">rec"), vec![1, 5]);
    assert_eq!(split(b"a||b"), vec![0, 2, 3, 5]);
    assert_eq!(split(b""), vec![1]);
    assert_eq!(split(b"|"), vec![1]);
    assert_eq!(split(b">"), vec![2]);
}

#[test]
fn ids_of_no_words_are_refused() {
    let mut set = TextMatchIds::default();
    for id in ["", ">", "|", ">|"] {
        assert_eq!(
            set.insert(id.as_bytes()),
            Err(EmptyId {
                id: id.as_bytes().to_vec()
            })
        );
    }
}

#[test]
fn runs_of_words_of_the_first_word_match() {
    let set = list(&["AB000001.1", "P12345|ABC_HUMAN", "NC_000099.1|", "rec2"]);
    assert!(set.find(b"gb|AB000001.1| GenBank style"));
    assert!(set.find(b"sp|P12345|ABC_HUMAN protein-like id"));
    assert!(set.find(b"ref|NC_000099.1|"));
    assert!(set.find(b"lcl|rec2 local"));
    assert!(!set.find(b"rec22"));
    assert!(!set.find(b"x AB000001.1"));
    let two = list(&["lcl|rec2"]);
    assert!(two.find(b"lcl|rec2 local"));
    assert!(!two.find(b"rec2"));
}

fn seq_id(text: &str) -> ParsedSeqId {
    seq_id_from_text(text.as_bytes(), parse_flags::ANY_RAW, &mut |_| {}).unwrap()
}

fn seq_id_list(ids: &[&str]) -> SeqIdSet {
    let mut set = SeqIdSet::default();
    for id in ids {
        set.add(&seq_id(id));
    }
    set
}

fn record(ids: Vec<ParsedSeqId>) -> Record {
    Record {
        ids,
        title: Some(b"title".to_vec()),
        taxname: None,
        has_source: false,
        get_title: None,
        iupac: b"ACGT".to_vec(),
        texts: Vec::new(),
    }
}

// NCBI BLAST+ 2.17.0 windowmasker -text_match F on a database made with -parse_seqids (probe
// s1 of line NP M5, 2026-11-21: tmf_wm_ps_ids_case, _unver, _forms, _lcl).
#[test]
fn seq_ids_match_as_ncbi_s_handles() {
    let set = seq_id_list(&["ab000001.1", "nc_000099.1", "lcl|REC2", "lcl|rec4."]);
    assert!(set.find(&[seq_id("gb|AB000001.1|")]));
    assert!(set.find(&[seq_id("ref|NC_000099.1|")]));
    assert!(set.find(&[seq_id("lcl|rec2")]));
    assert!(set.find(&[seq_id("lcl|rec4.")]));
    assert!(!set.find(&[seq_id("lcl|rec1")]));
    let bare = seq_id_list(&["BK000001.1", "P12345"]);
    assert!(bare.find(&[seq_id("tpg|BK000001.1|")]));
    let unversioned = seq_id_list(&["AB000001", "NC_000099", "AB000001.2"]);
    assert!(!unversioned.find(&[seq_id("gb|AB000001.1|")]));
    assert!(!unversioned.find(&[seq_id("ref|NC_000099.1|")]));
}

// tmf_wm_ids_ord of the same probe: `gnl|BL_ORD_ID|0` names OID 0 by a string tag "0"
// (`CObject_id::SetStrOrId` takes 1-9 first), which is not the integer tag of the database.
#[test]
fn an_ordinal_id_of_zero_is_a_string_tag() {
    let set = seq_id_list(&["gnl|BL_ORD_ID|0", "gnl|BL_ORD_ID|3"]);
    let ord = |oid| {
        ParsedSeqId::General(crate::seqid::model::Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            tag: crate::seqid::model::ObjectId::Id(oid),
        })
    };
    assert!(!set.find(&[ord(0)]));
    assert!(set.find(&[ord(3)]));
    assert!(set.find(&[ord(1), ord(3)]));
}

// The giim tree finds the first handle whose Giimport-id is not `Equals` (seq_id_tree.cpp
// 2926-2941): an equal giim id of a sequence gets a new handle, another one of the same
// number the list's.
#[test]
fn a_giim_id_finds_the_handle_of_another_giimport_id() {
    let set = seq_id_list(&["gim|5", "gim|5"]);
    assert!(!set.is_empty());
    assert!(!set.find(&[seq_id("gim|5")]));
    assert!(!set.find(&[seq_id("gim|6")]));
    let other = ParsedSeqId::Giim(crate::seqid::model::GiimportId {
        id: 5,
        db: Some(b"x".to_vec()),
        release: None,
    });
    assert!(set.find(&[other]));
}

fn printed(lists: &IdLists, ids: &[&str]) -> Option<String> {
    let ids: Vec<ParsedSeqId> = ids.iter().map(|id| seq_id(id)).collect();
    lists.handle_seq_id(&ids).map(|id| {
        let mut out = Vec::new();
        crate::seqid::label::fasta(&id, &mut out);
        String::from_utf8(out).unwrap()
    })
}

// The list's handles live for the run, so the info of a database id's handle is the one a
// list id made: a packed GenBank/EMBL/DDBJ accession takes its alternative (another number of
// the same prefix, digit count and version too), a name keeps the list's case, a local string
// gets the sequence's case back. The access handle is the first gi, else Textseq-id.
#[test]
fn the_writers_print_the_seq_id_that_a_list_id_made() {
    let lists = |ids: &[&str]| IdLists {
        ids: Some(IdSet::SeqId(seq_id_list(ids))),
        exclude_ids: None,
    };
    let ddbj = lists(&["AB000001.1", "sp|P12345|abc_human", "lcl|REC1"]);
    assert_eq!(
        printed(&ddbj, &["gb|AB000001.1|"]).as_deref(),
        Some("dbj|AB000001.1|")
    );
    assert_eq!(
        printed(&ddbj, &["gb|AB000003.1|"]).as_deref(),
        Some("dbj|AB000003.1|")
    );
    assert_eq!(
        printed(&ddbj, &["lcl|x", "emb|AB000004.1|"]).as_deref(),
        Some("dbj|AB000004.1|")
    );
    assert_eq!(printed(&ddbj, &["gb|AB000001.2|"]), None);
    assert_eq!(printed(&ddbj, &["gb|AB0000001.1|"]), None);
    assert_eq!(printed(&ddbj, &["gi|5", "gb|AB000001.1|"]), None);
    assert_eq!(
        printed(&ddbj, &["sp|P12345|ABC_HUMAN"]).as_deref(),
        Some("sp|P12345|abc_human")
    );
    assert_eq!(printed(&ddbj, &["lcl|rec1"]).as_deref(), Some("lcl|rec1"));
    let excluded = IdLists {
        ids: None,
        exclude_ids: Some(IdSet::SeqId(seq_id_list(&["emb|AB000002.1|"]))),
    };
    assert_eq!(
        printed(&excluded, &["gb|AB000001.1|"]).as_deref(),
        Some("emb|AB000001.1|")
    );
    let text = IdLists {
        ids: Some(IdSet::Text(TextMatchIds::default())),
        exclude_ids: None,
    };
    assert_eq!(printed(&text, &["gb|AB000001.1|"]), None);
}

// A list of which no id was read is empty: `consider` keeps every sequence.
#[test]
fn an_empty_seq_id_list_keeps_every_sequence() {
    let lists = IdLists {
        ids: Some(IdSet::SeqId(SeqIdSet::default())),
        exclude_ids: None,
    };
    assert!(lists.consider(&record(vec![seq_id("lcl|x")])));
    let lists = IdLists {
        ids: Some(IdSet::SeqId(seq_id_list(&["lcl|y"]))),
        exclude_ids: None,
    };
    assert!(!lists.consider(&record(vec![seq_id("lcl|x")])));
    assert!(lists.consider(&record(vec![seq_id("lcl|Y")])));
    let lists = IdLists {
        ids: None,
        exclude_ids: Some(IdSet::SeqId(seq_id_list(&["lcl|y"]))),
    };
    assert!(!lists.consider(&record(vec![seq_id("lcl|y")])));
}

// A list matched as text reads the title the source made (`Record::get_title`), and only a
// list that is not empty asks for it.
#[test]
fn a_text_list_matches_the_title_the_source_made() {
    let lists = IdLists {
        ids: Some(IdSet::Text(list(&["Homo"]))),
        exclude_ids: None,
    };
    assert!(lists.match_titles());
    let mut titled = record(vec![seq_id("lcl|x")]);
    titled.get_title = Some(b"Homo sapiens".to_vec());
    assert!(lists.consider(&titled));
    titled.get_title = Some(b"No definition line found".to_vec());
    assert!(!lists.consider(&titled));
    let empty = IdLists {
        ids: None,
        exclude_ids: Some(IdSet::Text(TextMatchIds::default())),
    };
    assert!(!empty.match_titles());
    let seq_ids = IdLists {
        ids: Some(IdSet::SeqId(seq_id_list(&["lcl|y"]))),
        exclude_ids: None,
    };
    assert!(!seq_ids.match_titles());
}

// tmf_wm_ids_bad of probe s1: NCBI's post (its build tree aside), here in LOSAT's words; the
// id keeps its bytes.
#[test]
fn an_id_that_cseq_id_cannot_read_is_posted_with_ncbi_s_report() {
    let parse = |id: &[u8]| seq_id_from_text(id, parse_flags::ANY_RAW, &mut |_| {});
    let error = parse(b"gb|").unwrap_err();
    assert_eq!(
        String::from_utf8(error.ncbi_what_bytes()).unwrap(),
        format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/objects/seq/../seqloc/Seq_id.cpp\", line \
                 2501: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::Set() - \
                 Malformatted ID gb|\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    );
    assert_eq!(
        cannot_read(b"gb|", &error),
        b"CWinMaskConfig::FillIdList(): can't understand id: gb|: Malformatted ID gb|\n: ignoring"
    );
    let error = parse(b"caf\xe9\r").unwrap_err();
    assert_eq!(
        cannot_read(b"caf\xe9\r", &error),
        b"CWinMaskConfig::FillIdList(): can't understand id: caf\xe9\r: Malformatted ID caf\xe9\n: ignoring"
    );
    let error = parse(b"").unwrap_err();
    assert_eq!(
        cannot_read(b"", &error),
        b"CWinMaskConfig::FillIdList(): can't understand id: : Empty bare accession supplied\n: ignoring"
    );
}

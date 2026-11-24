use super::*;
use crate::seqid::model::{ObjectId, TextseqId};

fn parsed(defline: &[u8], ids: Vec<ParsedSeqId>, title: Option<&[u8]>) -> ParsedIds {
    // The DefLine field encodes the defline as an ASN.1 `PrintableString`; for these lines (no
    // backslash, no control byte except ^A and tab) that is the line with \1 and \t escapes.
    let mut field = Vec::new();
    for &c in defline {
        match c {
            0x01 => field.extend_from_slice(b"\\1"),
            b'\t' => field.extend_from_slice(b"\\t"),
            _ => field.push(c),
        }
    }
    ParsedIds {
        ids,
        hist_old_id: None,
        defline_field: field,
        title_descr: title.map(<[u8]>::to_vec),
    }
}

fn refseq(acc: &str, version: i32) -> ParsedSeqId {
    ParsedSeqId::Other(TextseqId {
        accession: Some(acc.as_bytes().to_vec()),
        version: Some(version),
        ..TextseqId::default()
    })
}

fn titles_and_ids(set: &DeflineSet) -> Vec<(Vec<u8>, Vec<ParsedSeqId>)> {
    set.deflines
        .iter()
        .map(|d| (d.title.clone().unwrap(), d.ids.clone()))
        .collect()
}

// survey B fixture p02_multi (makeblastdb 2.17.0, blastdbcmd %i|%t).
#[test]
fn the_second_parse_splits_at_ctrl_a_and_keeps_the_spaces_of_titles() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    let p = parsed(
        b">ref|NP_100003.1| first title\x01ref|NP_100004.1| second title\x01gi|1005| third title",
        vec![refseq("NP_100003", 1)],
        Some(b"first title"),
    );
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(
        titles_and_ids(&set),
        [
            (b"first title".to_vec(), vec![refseq("NP_100003", 1)]),
            (b"second title".to_vec(), vec![refseq("NP_100004", 1)]),
            (b"third title".to_vec(), vec![ParsedSeqId::Gi(1005)]),
        ]
    );
    let p = parsed(
        b">ref|NP_100011.1|\ttab right after id",
        vec![refseq("NP_100011", 1)],
        None,
    );
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(titles_and_ids(&set)[0].0, b"  tab right after id");
    // A space right after '>': the id is empty, NCBI's parse throws.
    let p = parsed(
        b">  ref|NP_100008.1|  leading spaces",
        vec![refseq("NP_100008", 1)],
        Some(b"leading spaces"),
    );
    assert!(second_parse(&p, true, &mut no_warn).is_none());
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: true,
        length: 1,
        seq_data: None,
        parsed: Some(p),
        delta_fix: None,
    };
    let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
    assert_eq!(
        titles_and_ids(&set),
        [(b"leading spaces".to_vec(), vec![refseq("NP_100008", 1)])]
    );
}

// survey B fixture p04_molmismatch and the probe of a nucleotide range.
#[test]
fn ids_of_the_other_molecule_become_local() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    let p = parsed(b">NC_000913.3 raw nucleotide", vec![], None);
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(
        set.deflines[0].ids,
        [ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))]
    );
    let mut ids = vec![refseq("NC_000913", 3), ParsedSeqId::Gi(5)];
    localize_other_molecule(&mut ids, true, &mut no_warn);
    assert_eq!(
        ids[0],
        ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))
    );
    assert_eq!(ids[1], ParsedSeqId::Gi(5));
    // A nucleotide id with a range: the Bioseq has lcl|3 and the replaced ref id in
    // its history, so the local second parse is not used.
    let mut p = parsed(
        b">NC_000913.3:1-5 t",
        vec![ParsedSeqId::Local(ObjectId::Id(3))],
        Some(b"t"),
    );
    p.hist_old_id = Some(refseq("NC_000913", 3));
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: false,
        length: 1,
        seq_data: None,
        parsed: Some(p),
        delta_fix: None,
    };
    let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
    assert_eq!(
        titles_and_ids(&set),
        [(b"t".to_vec(), vec![ParsedSeqId::Local(ObjectId::Id(3))])]
    );
}

#[test]
fn large_gis_skip_the_sequence_and_log_the_last_gi() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    assert_eq!(
        large_gi_to_skip(&[ParsedSeqId::Gi(4_294_967_295)], &mut no_warn),
        None
    );
    assert_eq!(
        large_gi_to_skip(
            &[ParsedSeqId::Gi(4_294_967_296), ParsedSeqId::Gi(7)],
            &mut no_warn
        ),
        Some(7)
    );
}

// makeprofiledb's query with a descriptor stated as NULL (`Seqdesc::Empty`, line RG item D3c-3b's
// reader): NCBI's `x_GetBioseqBinaryHeader` dereferences each descriptor until it finds the binary
// lines, so without them the definition lines stop at the empty one wherever it stands
// (makeprofiledb 2.17.0, probe D3c-3c `probeB` B01-B03), and so does the taxid loop before the
// first `org` (`first_org_taxid`).
#[test]
fn a_descriptor_stated_as_null_makes_no_definition_lines() {
    for descr in [
        vec![Seqdesc::Empty],
        vec![Seqdesc::Title(b"t".to_vec()), Seqdesc::Empty],
        vec![Seqdesc::Empty, Seqdesc::Title(b"t".to_vec())],
    ] {
        assert_eq!(
            bioseq_deflines(&[], &descr),
            Err(BioseqDeflinesError::NullDescriptor),
            "{descr:?}"
        );
        assert_eq!(
            first_org_taxid(&descr),
            Err(BioseqDeflinesError::NullDescriptor),
            "{descr:?}"
        );
    }
    assert_eq!(first_org_taxid(&[Seqdesc::Title(b"t".to_vec())]), Ok(None));
}

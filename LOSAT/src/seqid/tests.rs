//! Tests of the Seq-id model, decoder and printings against NCBI BLAST+ 2.17.0: the
//! databases, the definition lines blastdbcmd reads back and the id fields of blastp and
//! blastn reports in `LOSAT/tests/seqid` (`make_fixtures.py` there).

use std::collections::HashSet;
use std::path::PathBuf;

use super::ber::decode_full_set;
use super::defline_set::{
    alignment_line_id, compare, description_id, seq_id_list, tabular_id, Sic, TabularId,
};
use super::label::{bare, content, fasta, first_gi, format_label};
use super::model::{
    Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId,
};
use super::rank::{best_rank_score, find_best, text_score};
use super::{BlastDefline, DeflineSet};
use crate::report::defline::generate_defline;
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::Encoder;
use crate::report::serial::json::JsonEncoder;
use crate::seqdb::{SeqDb, SeqDbError, SeqType};

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("seqid")
}

fn unescape(field: &str) -> String {
    let mut out = String::new();
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            other => panic!("bad escape {other:?} in {field:?}"),
        }
    }
    out
}

fn read_tsv(name: &str) -> Vec<Vec<String>> {
    let path = fixture_dir().join("expected").join(name);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|line| line.split('\t').map(unescape).collect())
        .collect()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).expect("ASCII")
}

fn with<F: FnOnce(&mut Vec<u8>)>(f: F) -> String {
    let mut out = Vec::new();
    f(&mut out);
    text(&out)
}

fn open(name: &str, seqtype: SeqType) -> SeqDb {
    let base = fixture_dir().join("db").join(name);
    SeqDb::open(base.to_str().unwrap(), seqtype).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn sets(db: &SeqDb) -> Vec<DeflineSet> {
    (0..db.num_seqs())
        .map(|oid| db.defline_set(oid).unwrap_or_else(|e| panic!("{oid}: {e}")))
        .collect()
}

const DATABASES: [(&str, SeqType); 4] = [
    ("prot_ids", SeqType::Protein),
    ("nucl_ids", SeqType::Nucleotide),
    ("prot_extra", SeqType::Protein),
    ("variants", SeqType::Protein),
];

/// The fields of blastdbcmd `-outfmt "%o %i %a %g %T %e %t"` for one definition line.
fn blastdbcmd_fields(oid: u32, defline: &BlastDefline) -> Vec<String> {
    let best = defline.display_id().expect("an id");
    let gi = first_gi(&defline.ids);
    vec![
        oid.to_string(),
        with(|out| fasta(best, out)),
        with(|out| bare(best, out)),
        if gi == 0 {
            "N/A".to_string()
        } else {
            gi.to_string()
        },
        defline.taxid.unwrap_or(0).to_string(),
        defline
            .memberships
            .iter()
            .flatten()
            .map(|&m| i64::from(m))
            .sum::<i64>()
            .to_string(),
        defline.title.as_deref().map_or("N/A".to_string(), text),
    ]
}

fn check_deflines(name: &str, sets: &[(u32, DeflineSet)], expected: &[Vec<String>]) {
    let got: Vec<Vec<String>> = sets
        .iter()
        .flat_map(|(oid, set)| set.deflines.iter().map(|d| blastdbcmd_fields(*oid, d)))
        .collect();
    assert_eq!(got.len(), expected.len(), "{name}: definition lines");
    for (got, expected) in got.iter().zip(expected) {
        assert_eq!(got, expected, "{name}");
    }
}

#[test]
fn every_definition_line_reads_as_ncbis_blastdbcmd_reads_it() {
    for (name, seqtype) in DATABASES {
        let db = open(name, seqtype);
        let error = db
            .check_ordinal_ids()
            .expect_err("a database with parsed Seq-ids");
        assert!(
            matches!(error, SeqDbError::Unsupported { .. }),
            "{name}: {error}"
        );
        let sets: Vec<(u32, DeflineSet)> = sets(&db)
            .into_iter()
            .zip(0..)
            .map(|(s, o)| (o, s))
            .collect();
        check_deflines(name, &sets, &read_tsv(&format!("{name}.deflines.tsv")));
        for (oid, set) in &sets {
            for defline in &set.deflines {
                assert_eq!(
                    defline.replaced_title_bytes + defline.replaced_id_bytes,
                    0,
                    "{name} {oid}"
                );
            }
        }
    }
}

/// The records of `expected/18S_sample.ber`: ordinal id and the stored header record.
fn sample_18s() -> Vec<(u32, Vec<u8>)> {
    let bytes = std::fs::read(fixture_dir().join("expected").join("18S_sample.ber")).unwrap();
    let mut records = Vec::new();
    let mut rest = &bytes[..];
    while !rest.is_empty() {
        let oid = u32::from_be_bytes(rest[0..4].try_into().unwrap());
        let len = u32::from_be_bytes(rest[4..8].try_into().unwrap()) as usize;
        records.push((oid, rest[8..8 + len].to_vec()));
        rest = &rest[8 + len..];
    }
    records
}

#[test]
fn the_18s_sample_reads_as_ncbis_blastdbcmd_reads_it() {
    let sets: Vec<(u32, DeflineSet)> = sample_18s()
        .into_iter()
        .map(|(oid, record)| (oid, decode_full_set(&record).unwrap()))
        .collect();
    assert_eq!(sets.len(), 26);
    check_deflines("18S sample", &sets, &read_tsv("18S_sample.deflines.tsv"));
    let first = &sets[0].1.deflines[0];
    assert_eq!(first.memberships.as_deref(), Some(&[4][..]));
    assert_eq!(first.taxid, Some(4956));
    assert_eq!(
        sets.iter()
            .find(|(oid, _)| *oid == 1479)
            .unwrap()
            .1
            .deflines
            .len(),
        23
    );
}

/// The whole database 18S_fungal_sequences that NCBI distributes (3,905 sequences, 4,035
/// definition lines) against blastdbcmd's dump: run with `LOSAT_SEQID_18S_DB` (the
/// database's path without extension) and `LOSAT_SEQID_18S_DEFLINES` (the dump of
/// `blastdbcmd -entry all -outfmt '%o%%%i%%%a%%%g%%%T%%%e%%%t'`).
#[test]
#[ignore = "needs NCBI's 18S_fungal_sequences database (LOSAT_SEQID_18S_DB)"]
fn the_whole_distributed_18s_database_reads_as_ncbis_blastdbcmd_reads_it() {
    let base = std::env::var("LOSAT_SEQID_18S_DB").expect("LOSAT_SEQID_18S_DB");
    let dump = std::env::var("LOSAT_SEQID_18S_DEFLINES").expect("LOSAT_SEQID_18S_DEFLINES");
    let expected: Vec<Vec<String>> = std::fs::read_to_string(dump)
        .unwrap()
        .lines()
        .map(|line| line.splitn(7, '%').map(str::to_string).collect())
        .collect();
    let db = SeqDb::open(&base, SeqType::Nucleotide).unwrap();
    assert!(db.check_ordinal_ids().is_err());
    let sets: Vec<(u32, DeflineSet)> = sets(&db)
        .into_iter()
        .zip(0..)
        .map(|(s, o)| (o, s))
        .collect();
    assert_eq!(sets.len(), 3905);
    check_deflines("18S", &sets, &expected);
}

/// `BeginContainer`-free BER of one Seq-id as NCBI's writer gives it (`MSerial_AsnBinary`).
fn ber(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    let mut enc = AsnBinaryEncoder::new(&mut out);
    SeqId::Db(id).write(&mut enc);
    assert_eq!(enc.error(), None);
    out
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

#[test]
fn every_seq_id_encodes_as_makeblastdb_stored_it() {
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for (name, seqtype) in DATABASES {
        let db = open(name, seqtype);
        for (at, volume) in db.volumes().iter().enumerate() {
            assert_eq!(at, 0);
            for oid in 0..volume.num_oids() {
                records.push((format!("{name} {oid}"), volume.header_bytes(oid).unwrap()));
            }
        }
    }
    for (oid, record) in sample_18s() {
        records.push((format!("18S {oid}"), record));
    }
    let mut ids = 0;
    for (what, record) in &records {
        let set = decode_full_set(record).unwrap();
        for id in set.ids() {
            assert!(contains(record, &ber(id)), "{what}: {id:?}");
            ids += 1;
        }
    }
    assert!(ids > 150, "{ids}");
}

/// The first word of the subject Bioseq's generated title (`GenerateDefline` of its title;
/// no fixture record with a `lcl|Subject_` id has a Seq-id that changes the title).
fn title_word(bioseq: &BlastDefline, protein: bool) -> Vec<u8> {
    let title = bioseq.title.as_deref().unwrap_or(b"");
    let defline = generate_defline(title, protein, false).expect("a known encoding");
    defline.text.split(|&c| c == b' ').next().unwrap().to_vec()
}

/// The ids and titles that `tabular_id` and the joins give for one sequence, in the column
/// order `sseqid sallseqid sacc saccver sallacc sgi sallgi stitle salltitles staxid staxids`.
fn tabular_row(set: &DeflineSet, protein: bool) -> Vec<String> {
    let aln = set.alignment_id().expect("an id");
    let bioseq = &set.deflines[set.bioseq_defline(aln).expect("a target")];
    let word = title_word(bioseq, protein);
    let target = seq_id_list(&bioseq.ids, &word);
    let staxid = set.first_taxid();
    let staxids = set.tax_ids();
    vec![
        with(|out| tabular_id(&target, TabularId::Full, out)),
        with(|out| set.tabular_all(TabularId::Full, &word, out)),
        with(|out| tabular_id(&target, TabularId::Accession, out)),
        with(|out| tabular_id(&target, TabularId::AccessionVersion, out)),
        with(|out| set.tabular_all(TabularId::Accession, &word, out)),
        with(|out| tabular_id(&target, TabularId::Gi, out)),
        with(|out| set.tabular_all(TabularId::Gi, &word, out)),
        with(|out| set.first_title(out)),
        with(|out| set.all_titles(out)),
        if staxid == 0 {
            "N/A".to_string()
        } else {
            staxid.to_string()
        },
        if staxids.is_empty() {
            "N/A".to_string()
        } else {
            staxids
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(";")
        },
    ]
}

const CASES: [(&str, &str, SeqType); 4] = [
    ("blastp_prot_ids", "prot_ids", SeqType::Protein),
    ("blastn_nucl_ids", "nucl_ids", SeqType::Nucleotide),
    ("blastp_prot_extra", "prot_extra", SeqType::Protein),
    ("blastp_variants", "variants", SeqType::Protein),
];

/// Matches each distinct NCBI row (one per subject: the reports have one HSP per subject,
/// and the queries repeat subjects) with the sequence whose computed row `computed` gives
/// it; a sequence matches at most one row.
fn match_rows<T: PartialEq + std::fmt::Debug>(
    case: &str,
    expected: Vec<T>,
    computed: impl Fn(usize) -> Option<T>,
    sequences: usize,
) -> Vec<usize> {
    let mut distinct: Vec<T> = Vec::new();
    for row in expected {
        if !distinct.contains(&row) {
            distinct.push(row);
        }
    }
    assert!(!distinct.is_empty(), "{case}");
    let rows: Vec<Option<T>> = (0..sequences).map(computed).collect();
    let mut matched = HashSet::new();
    distinct
        .iter()
        .map(|row| {
            let at = rows
                .iter()
                .position(|candidate| candidate.as_ref() == Some(row))
                .unwrap_or_else(|| panic!("{case}: no sequence gives {row:?}"));
            assert!(matched.insert(at), "{case}: sequence {at} gives two rows");
            at
        })
        .collect()
}

#[test]
fn tabular_id_and_title_fields_are_ncbis() {
    for (case, name, seqtype) in CASES {
        let sets = sets(&open(name, seqtype));
        let mut expected = read_tsv(&format!("{case}.fmt6.tsv"));
        for row in &mut expected {
            // Approved difference (line DBID, DECISIONS 5, Q2): NCBI prints the type of the
            // first id registered under the same accession key in the process (`dbj|`, from
            // another sequence); LOSAT prints the definition line's own `gb|`.
            if row[0] == "gi|1515|dbj|AB000015.1|" {
                row[0] = "gi|1515|gb|AB000015.1|".to_string();
            }
        }
        let protein = seqtype == SeqType::Protein;
        match_rows(
            case,
            expected,
            |oid| Some(tabular_row(&sets[oid], protein)),
            sets.len(),
        );
    }
}

/// Outfmt 5 `Hit_id`, `Hit_def` and `Hit_accession` of one sequence.
fn xml_row(set: &DeflineSet) -> Vec<String> {
    let (id, mut def) = set.hit_id_and_def(true);
    if def.is_empty() {
        def = b"No definition line".to_vec();
    }
    let aln = set.alignment_id().unwrap();
    let target = &set.deflines[set.bioseq_defline(aln).unwrap()].ids;
    let best = find_best(target, text_score).unwrap();
    vec![
        text(&id),
        text(&def),
        with(|out| format_label(best, false, out)),
    ]
}

/// A Seq-id in ASN.1 text (outfmt 8, blastdbcmd `%d`) with white space collapsed.
fn asn_text(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    enc.begin_document("Seq-id");
    SeqId::Db(id).write(&mut enc);
    enc.end_document();
    assert_eq!(enc.error(), None);
    let text = text(&out);
    let text = text.strip_prefix("Seq-id ::= ").expect("a document");
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// A Seq-id in JSON (outfmt 12) without the white space outside strings.
fn json_text(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    let mut enc = JsonEncoder::new(&mut out);
    enc.begin_document("Seq-id");
    SeqId::Db(id).write(&mut enc);
    enc.end_document();
    assert_eq!(enc.error(), None);
    let (mut compact, mut in_string, mut escaped) = (String::new(), false, false);
    for c in text(&out).chars() {
        if in_string {
            in_string = escaped || c != '"';
            escaped = !escaped && c == '\\';
        } else if c == '"' {
            in_string = true;
        } else if c.is_whitespace() {
            continue;
        }
        compact.push(c);
    }
    // `{"Seq_id":<value>}`: the value.
    let value = &compact[compact.find(':').expect("a key") + 1..];
    value.strip_suffix('}').expect("a document").to_string()
}

#[test]
fn xml_hit_fields_and_alignment_ids_are_ncbis() {
    for (case, name, seqtype) in CASES {
        let sets = sets(&open(name, seqtype));
        let hits = read_tsv(&format!("{case}.fmt5.tsv"));
        match_rows(case, hits, |oid| Some(xml_row(&sets[oid])), sets.len());
        // The Seq-align subject id of each hit in outfmt 8 and 12. NCBI's XML and outfmt 0
        // put consecutive Seq-aligns whose subject ids `Match` into one hit (`pdb|1ABC|` and
        // `pdb|1ABC|AA` of `variants`), so the hits are matched apart from the Seq-aligns.
        let aln_text = |oid: usize| sets[oid].alignment_id().map(|id| vec![asn_text(id)]);
        let aln_json = |oid: usize| sets[oid].alignment_id().map(|id| vec![json_text(id)]);
        let by_text = match_rows(
            case,
            read_tsv(&format!("{case}.fmt8.tsv")),
            aln_text,
            sets.len(),
        );
        let by_json = match_rows(
            case,
            read_tsv(&format!("{case}.fmt12.tsv")),
            aln_json,
            sets.len(),
        );
        assert_eq!(
            by_text, by_json,
            "{case}: outfmt 8 and 12 name the same sequences"
        );
    }
}

#[test]
fn every_seq_id_writes_as_ncbis_blastdbcmd_writes_it() {
    for (name, seqtype) in DATABASES {
        let sets = sets(&open(name, seqtype));
        let rows = read_tsv(&format!("{name}.asn.tsv"));
        let stored: usize = sets.iter().map(|set| set.ids().count()).sum();
        assert_eq!(rows.len(), stored, "{name}: Seq-ids");
        for row in rows {
            let [oid, line, at, expected] = &row[..] else {
                panic!("{name}: {row:?}");
            };
            let (oid, line, at): (usize, usize, usize) = (
                oid.parse().unwrap(),
                line.parse().unwrap(),
                at.parse().unwrap(),
            );
            let id = &sets[oid].deflines[line].ids[at];
            assert_eq!(&asn_text(id), expected, "{name} {oid} {line} {at}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:340-359
// (`s_WrapOutputLine`: a title of more than 60 bytes gets an end of line after the first
// white space at or after each multiple of 60).
fn wrap60(title: &str) -> String {
    let bytes = title.as_bytes();
    if bytes.len() <= 60 {
        return title.to_string();
    }
    let mut out = String::new();
    let mut do_wrap = false;
    for (i, &c) in bytes.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        out.push(c as char);
        if do_wrap && c.is_ascii_whitespace() {
            out.push('\n');
            do_wrap = false;
        }
    }
    out
}

/// The definition-line block of the outfmt 0 alignment of one sequence: every definition
/// line with its id (`alignment_line_id`) and title; an empty title is the subject Bioseq's
/// (its target definition line's), `unnamed protein product` for an empty protein one.
fn block(set: &DeflineSet, protein: bool, show_gis: bool) -> String {
    let aln = set.alignment_id().unwrap();
    let bioseq_title = set.deflines[set.bioseq_defline(aln).unwrap()]
        .title
        .clone()
        .unwrap_or_default();
    let mut lines = Vec::new();
    for (i, defline) in set.deflines.iter().enumerate() {
        let mut title = defline.title.clone().unwrap_or_default();
        if title.is_empty() {
            title = bioseq_title.clone();
        }
        if title.is_empty() && protein {
            title = b"unnamed protein product".to_vec();
        }
        let id = with(|out| alignment_line_id(defline, show_gis, out));
        let prefix = if i == 0 { ">" } else { " " };
        lines.push(format!("{prefix}{id} {}", wrap60(&text(&title))));
    }
    lines.join("\n")
}

#[test]
fn outfmt0_id_parts_are_ncbis() {
    for (case, name, seqtype) in CASES {
        let protein = seqtype == SeqType::Protein;
        let sets = sets(&open(name, seqtype));
        let mut checked = 0;
        for row in read_tsv(&format!("{case}.fmt0.tsv")) {
            let (kind, line) = (&row[0], &row[1]);
            let show_gis = kind.starts_with("show_gis");
            if kind.ends_with("_block") {
                let matching = sets
                    .iter()
                    .filter(|set| block(set, protein, show_gis) == *line)
                    .count();
                assert_eq!(matching, 1, "{case}: sequences that give the block\n{line}");
            } else {
                // The description: the id of `GetDisplayIds`, two spaces after `gi|N|`, then
                // the title of the subject Bioseq (truncated, then the scores).
                let found = sets.iter().any(|set| {
                    let aln = set.alignment_id().unwrap();
                    let shown = set.display_defline(aln).unwrap();
                    let defline = &set.deflines[shown.index];
                    let gi = if shown.found { defline.gi() } else { 0 };
                    let mut start = Vec::new();
                    let wrote_gi =
                        description_id(defline.display_id().unwrap(), gi, show_gis, &mut start);
                    if wrote_gi {
                        start.extend_from_slice(b"  ");
                    }
                    let mut title = set.deflines[set.bioseq_defline(aln).unwrap()]
                        .title
                        .clone()
                        .unwrap_or_default();
                    if title.is_empty() && protein {
                        title = b"unnamed protein product".to_vec();
                    }
                    start.extend_from_slice(&title[..title.len().min(10)]);
                    line.starts_with(&text(&start))
                });
                assert!(
                    found,
                    "{case}: no sequence gives the description line\n{line}"
                );
            }
            checked += 1;
        }
        assert!(checked > 10, "{case}: {checked}");
    }
}

fn textseq(accession: Option<&str>, name: Option<&str>, version: Option<i32>) -> TextseqId {
    TextseqId {
        name: name.map(|s| s.as_bytes().to_vec()),
        accession: accession.map(|s| s.as_bytes().to_vec()),
        release: None,
        version,
    }
}

/// Forms that no fixture holds (makeblastdb does not write them from FASTA), by NCBI's rules
/// quoted in `label` and `rank`.
#[test]
fn printings_of_forms_without_a_fixture() {
    let giim = ParsedSeqId::Giim(GiimportId {
        id: 77,
        db: Some(b"embl".to_vec()),
        release: None,
    });
    assert_eq!(with(|o| fasta(&giim, o)), "gim|77");
    assert_eq!(with(|o| content(&giim, true, o)), "77");
    let nat = ParsedSeqId::NamedAnnotTrack(textseq(Some("na1"), None, Some(2)));
    assert_eq!(with(|o| fasta(&nat, o)), "nat|na1.2|");
    assert_eq!(with(|o| content(&nat, true, o)), "NA1.2");
    // A version 0: not in the FASTA form, `.0` in the content label.
    let zero = ParsedSeqId::Genbank(textseq(Some("ab1"), Some("LOC"), Some(0)));
    assert_eq!(with(|o| fasta(&zero, o)), "gb|ab1|LOC");
    assert_eq!(with(|o| content(&zero, true, o)), "AB1.0");
    assert_eq!(with(|o| content(&zero, false, o)), "AB1");
    // No accession: the name, not upper-cased, with the version when asked.
    let named = ParsedSeqId::Pir(textseq(None, Some("s1"), Some(3)));
    assert_eq!(with(|o| content(&named, true, o)), "s1.3");
    assert_eq!(with(|o| content(&named, false, o)), "s1");
    assert_eq!(with(|o| bare(&named, o)), "pir||s1");
    // A patent with an application number and a document type.
    let patent = ParsedSeqId::Patent(PatentSeqId {
        seqid: 4,
        cit: IdPat {
            country: b"JP".to_vec(),
            id: PatentNumber::AppNumber(b"2001-1".to_vec()),
            doc_type: Some(b"A".to_vec()),
        },
    });
    assert_eq!(with(|o| fasta(&patent, o)), "pat|JP|2001-1A|4");
    assert_eq!(with(|o| content(&patent, true, o)), "JP2001-1_4");
    // PDB chains: '|' and 0 in the FASTA form; a chain at or below ' ' leaves the label
    // without one; an unset chain reads as 32.
    let pdb = |chain| {
        ParsedSeqId::Pdb(PdbSeqId {
            mol: b"1XYZ".to_vec(),
            chain,
            rel: None,
            chain_id: None,
        })
    };
    assert_eq!(with(|o| fasta(&pdb(Some(i32::from(b'|'))), o)), "pdb|1XYZ|");
    assert_eq!(with(|o| fasta(&pdb(Some(0)), o)), "pdb|1XYZ| ");
    assert_eq!(with(|o| fasta(&pdb(None), o)), "pdb|1XYZ| ");
    assert_eq!(with(|o| content(&pdb(Some(32)), true, o)), "1XYZ");
    assert_eq!(with(|o| content(&pdb(Some(66)), true, o)), "1XYZ_B");
    // An unset chain: written only by JSON, which writes the defaults (DEFAULT 32).
    assert_eq!(asn_text(&pdb(None)), r#"pdb { mol "1XYZ" }"#);
    assert_eq!(
        json_text(&pdb(None)),
        r#"{"pdb":{"mol":"1XYZ","chain":32}}"#
    );
    assert_eq!(asn_text(&pdb(Some(32))), r#"pdb { mol "1XYZ", chain 32 }"#);
    // The formatter's label of a general id is its tag, of one with an empty tag `db:`.
    let general = |tag: &str| {
        ParsedSeqId::General(Dbtag {
            db: b"db1".to_vec(),
            tag: ObjectId::Str(tag.as_bytes().to_vec()),
        })
    };
    assert_eq!(with(|o| format_label(&general("x7"), true, o)), "x7");
    assert_eq!(with(|o| format_label(&general(""), true, o)), "db1:");
    // Ranks: TRACE after ti; a Textseq-id without version, accession or name ranks lower.
    let trace = ParsedSeqId::General(Dbtag {
        db: b"TRACE".to_vec(),
        tag: ObjectId::Id(1),
    });
    let ti = ParsedSeqId::General(Dbtag {
        db: b"ti".to_vec(),
        tag: ObjectId::Id(1),
    });
    assert_eq!((text_score(&ti), text_score(&trace)), (600, 605));
    assert_eq!(
        text_score(&ParsedSeqId::Other(textseq(None, None, None))),
        109
    );
    assert_eq!(best_rank_score(&ParsedSeqId::Gi(1)), 510);
    let ids = [trace.clone(), ti.clone()];
    assert_eq!(find_best(&ids, text_score), Some(&ti));
    assert_eq!(find_best(&ids[..1], text_score), Some(&trace));
    assert_eq!(find_best(&[], text_score), None);
}

#[test]
fn seq_id_comparisons_follow_ncbi() {
    let gb = ParsedSeqId::Genbank(textseq(Some("AB000015"), None, Some(1)));
    let emb = ParsedSeqId::Embl(textseq(Some("ab000015"), None, None));
    let emb2 = ParsedSeqId::Embl(textseq(Some("AB000016"), None, None));
    let other = ParsedSeqId::Other(textseq(Some("AB000015"), None, Some(1)));
    assert_eq!(compare(&gb, &emb), Sic::Yes);
    assert_eq!(compare(&gb, &emb2), Sic::No);
    assert_eq!(compare(&gb, &other), Sic::Diff);
    let local = |tag| ParsedSeqId::Local(tag);
    assert_eq!(
        compare(
            &local(ObjectId::Str(b"5".to_vec())),
            &local(ObjectId::Id(5))
        ),
        Sic::Yes
    );
    assert_eq!(
        compare(
            &local(ObjectId::Str(b"05".to_vec())),
            &local(ObjectId::Id(5))
        ),
        Sic::No
    );
    assert_eq!(
        compare(
            &local(ObjectId::Str(b"Ab".to_vec())),
            &local(ObjectId::Str(b"aB".to_vec()))
        ),
        Sic::Yes
    );
    let general = |db: &str, tag: i32| {
        ParsedSeqId::General(Dbtag {
            db: db.as_bytes().to_vec(),
            tag: ObjectId::Id(tag),
        })
    };
    assert_eq!(compare(&general("db", 1), &general("DB", 1)), Sic::Yes);
    assert_eq!(compare(&general("db", 1), &general("db", 2)), Sic::No);
    assert_eq!(compare(&general("db", 1), &general("xx", 1)), Sic::Diff);
    // The Bioseq's definition line: the first id comparable with the target decides.
    let set = DeflineSet {
        deflines: vec![
            BlastDefline {
                ids: vec![ParsedSeqId::Gi(1), ParsedSeqId::Gi(2)],
                title: Some(Vec::new()),
                ..BlastDefline::default()
            },
            BlastDefline {
                ids: vec![ParsedSeqId::Gi(2)],
                title: Some(Vec::new()),
                ..BlastDefline::default()
            },
        ],
    };
    assert_eq!(set.bioseq_defline(&ParsedSeqId::Gi(2)), Some(1));
    assert_eq!(set.alignment_id(), Some(&ParsedSeqId::Gi(1)));
    let shown = set.display_defline(&ParsedSeqId::Gi(3)).unwrap();
    assert_eq!((shown.index, shown.found), (1, false));
}

/// Rewrites a BER encoding with indefinite lengths into definite lengths (NCBI's reader
/// reads both).
fn to_definite(data: &[u8]) -> Vec<u8> {
    fn length(out: &mut Vec<u8>, len: usize) {
        if len < 0x80 {
            out.push(len as u8);
        } else {
            let bytes = len.to_be_bytes();
            let skip = bytes.iter().take_while(|&&b| b == 0).count();
            out.push(0x80 | (bytes.len() - skip) as u8);
            out.extend_from_slice(&bytes[skip..]);
        }
    }
    fn value(data: &[u8], pos: &mut usize, out: &mut Vec<u8>) {
        let tag = data[*pos];
        *pos += 1;
        if tag & 0x20 != 0 && data[*pos] == 0x80 {
            *pos += 1;
            let mut inner = Vec::new();
            while data[*pos..*pos + 2] != [0, 0] {
                value(data, pos, &mut inner);
            }
            *pos += 2;
            out.push(tag);
            length(out, inner.len());
            out.extend_from_slice(&inner);
        } else {
            let first = data[*pos];
            *pos += 1;
            let len = if first < 0x80 {
                usize::from(first)
            } else {
                let n = usize::from(first - 0x80);
                let len = data[*pos..*pos + n]
                    .iter()
                    .fold(0, |a, &b| a << 8 | usize::from(b));
                *pos += n;
                len
            };
            out.push(tag);
            length(out, len);
            out.extend_from_slice(&data[*pos..*pos + len]);
            *pos += len;
        }
    }
    let mut out = Vec::new();
    let mut pos = 0;
    value(data, &mut pos, &mut out);
    assert_eq!(pos, data.len());
    out
}

#[test]
fn definite_lengths_read_as_indefinite_ones() {
    for (oid, record) in sample_18s() {
        let definite = to_definite(&record);
        assert_ne!(definite, record);
        assert_eq!(
            decode_full_set(&definite).unwrap(),
            decode_full_set(&record).unwrap(),
            "18S {oid}"
        );
    }
}

#[test]
fn damaged_records_are_errors() {
    let (_, record) = sample_18s().remove(0);
    let decoded = decode_full_set(&record).unwrap();
    assert_eq!(decode_full_set(&[]).unwrap(), DeflineSet::default());
    // Every record cut short is an error; a changed byte gives an error or a value, never
    // a panic.
    for (oid, record) in sample_18s() {
        for cut in 1..record.len() {
            let error = decode_full_set(&record[..cut]).expect_err("cut short");
            assert!(
                error.to_string().starts_with("header record, byte "),
                "18S {oid} cut at {cut}: {error}"
            );
        }
        let mut changed = record.clone();
        for at in 0..record.len() {
            for byte in [0x00, 0x80, 0xff, record[at] ^ 0x01, record[at] ^ 0x20] {
                changed[at] = byte;
                let _ = decode_full_set(&changed);
            }
            changed[at] = record[at];
        }
    }
    let mut longer = record.clone();
    longer.push(0);
    assert_eq!(
        decode_full_set(&longer).unwrap_err().message,
        "unexpected data after Blast-def-line-set"
    );
    // A definition line without its Seq-id list: `30 80 30 80 a0 80 1a 01 x 00 00 00 00 00 00`.
    let no_seqid = [
        0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x01, b'x', 0, 0, 0, 0, 0, 0,
    ];
    assert_eq!(
        decode_full_set(&no_seqid).unwrap_err().message,
        "a Blast-def-line has no seqid"
    );
    // An alternative 20 of Seq-id does not exist.
    let id = [
        0x30, 0x80, 0x30, 0x80, 0xa1, 0x80, 0x30, 0x80, 0xb4, 0x80, 0x02, 0x01, 0x01,
    ];
    assert_eq!(
        decode_full_set(&id).unwrap_err().message,
        "tag 0xb4 is not an alternative of Seq-id"
    );
    // A title with bytes outside ' ' ..= '~' reads with '#' for each, and counts them.
    let title = [
        0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x03, b'a', 0xe9, 0x7f, 0, 0, 0xa1, 0x80, 0x30,
        0x80, 0xa0, 0x80, 0xa1, 0x80, 0x1a, 0x02, b'x', 0x01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let set = decode_full_set(&title).unwrap();
    let defline = &set.deflines[0];
    assert_eq!(defline.title.as_deref(), Some(&b"a##"[..]));
    assert_eq!(
        (defline.replaced_title_bytes, defline.replaced_id_bytes),
        (2, 1)
    );
    assert_eq!(
        defline.ids,
        [ParsedSeqId::Local(ObjectId::Str(b"x#".to_vec()))]
    );
    assert!(decoded.deflines[0].ids.len() == 2);
}

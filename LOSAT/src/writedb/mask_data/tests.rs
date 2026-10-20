//! Tests of the mask-column writer.
//!
//! Byte equality with NCBI: the databases of `LOSAT/tests/seqdb/db` that
//! `make_mask_fixtures.py` built in format version 4 (dust, segmasker and
//! windowmasker masks, `-mask_id`, several volumes) and the cases of
//! `make_mask_writer_fixtures.py` (`tests/seqdb/mask_writer`: Seq-locs that are
//! combined, algorithm ids, masks that match nothing, ranges past the end,
//! volumes without a masked sequence, ...) are rebuilt through the library
//! with the mask files of the same runs and compared file by file; the text
//! makeblastdb prints for a failed run is compared as well. The tests run no
//! NCBI program.
//!
//! The mask files are ASN.1 text. `writedb::mask_info` will read the files
//! (text and binary) for makeblastdb; until it lands, `parse_mask_text` below
//! reads the text files of the fixtures (a small subset of ASN.1 text: the
//! shapes dustmasker, segmasker and the fixture script write). The dust and
//! segmasker masks are `dust.asnt` / `seg.asnt` of `tests/seqdb/mask_writer`,
//! written by `dustmasker` / `segmasker -outfmt maskinfo_asn1_text` on
//! `mask_nucl.fna` / `mask_prot.faa`: the same masks as the binary files
//! `mask_dust.asnb` / `mask_seg.asnb` of the fixtures (a byte-equal database
//! shows it).

use std::path::{Path, PathBuf};

use super::registry::{
    FILTER_PROGRAM_DUST, FILTER_PROGRAM_OTHER, FILTER_PROGRAM_REPEAT, FILTER_PROGRAM_SEG,
    FILTER_PROGRAM_WINDOWMASKER,
};
use super::{ranges_are_empty, MaskInfoRegistry, MaskedRangeSet, MaskedRanges};
use crate::seqdb::{MaskRanges, SeqDb, SeqType};
use crate::writedb::build_db::BuildDatabase;
use crate::writedb::fasta_source::Warn;
use crate::writedb::mask_info::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskSeqId};
use crate::writedb::tests::{files_named, seqdb_fixtures, stored_date, time_of, TempDir};
use crate::writedb::{WriteDb, WriteDbError};

// ---------------------------------------------------------------------
// ASN.1 text: the part of it that the mask files use
// ---------------------------------------------------------------------

#[derive(Debug, PartialEq)]
enum Tok {
    Open,
    Close,
    Comma,
    Assign,
    Str(String),
    Word(String),
}

fn tokenize(text: &str) -> Vec<Tok> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = text.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        match c {
            _ if c.is_whitespace() => i += 1,
            '{' => {
                tokens.push(Tok::Open);
                i += 1;
            }
            '}' => {
                tokens.push(Tok::Close);
                i += 1;
            }
            ',' => {
                tokens.push(Tok::Comma);
                i += 1;
            }
            ':' => {
                assert_eq!(&chars[i..i + 3], [':', ':', '=']);
                tokens.push(Tok::Assign);
                i += 3;
            }
            '"' => {
                let end = chars[i + 1..].iter().position(|&c| c == '"').unwrap();
                tokens.push(Tok::Str(chars[i + 1..i + 1 + end].iter().collect()));
                i += end + 2;
            }
            _ => {
                let end = chars[i..]
                    .iter()
                    .position(|c| c.is_whitespace() || "{},\"".contains(*c))
                    .unwrap_or(chars.len() - i);
                tokens.push(Tok::Word(chars[i..i + end].iter().collect()));
                i += end;
            }
        }
    }
    tokens
}

#[derive(Debug, Clone)]
enum Node {
    Int(i64),
    Str(String),
    Seq(Vec<Node>),
    /// `name value`, the value being absent for `name` alone.
    Named(String, Option<Box<Node>>),
}

struct Parser {
    tokens: Vec<Tok>,
    pos: usize,
}

impl Parser {
    fn next(&mut self) -> Tok {
        let token = std::mem::replace(&mut self.tokens[self.pos], Tok::Comma);
        self.pos += 1;
        token
    }

    fn peek(&self) -> Option<&Tok> {
        self.tokens.get(self.pos)
    }

    fn value(&mut self) -> Node {
        match self.next() {
            Tok::Open => {
                let mut items = Vec::new();
                if self.peek() == Some(&Tok::Close) {
                    self.next();
                    return Node::Seq(items);
                }
                loop {
                    items.push(self.value());
                    match self.next() {
                        Tok::Comma => {}
                        Tok::Close => return Node::Seq(items),
                        other => panic!("unexpected {other:?}"),
                    }
                }
            }
            Tok::Str(s) => Node::Str(s),
            Tok::Word(w) => {
                if let Ok(n) = w.parse::<i64>() {
                    return Node::Int(n);
                }
                let inner = match self.peek() {
                    Some(Tok::Open | Tok::Str(_) | Tok::Word(_)) => Some(Box::new(self.value())),
                    _ => None,
                };
                Node::Named(w, inner)
            }
            other => panic!("unexpected {other:?}"),
        }
    }
}

fn field<'a>(items: &'a [Node], name: &str) -> Option<&'a Node> {
    items.iter().find_map(|node| match node {
        Node::Named(key, Some(value)) if key == name => Some(&**value),
        _ => None,
    })
}

fn seq(node: &Node) -> &[Node] {
    match node {
        Node::Seq(items) => items,
        other => panic!("not a sequence: {other:?}"),
    }
}

fn int(node: &Node) -> i64 {
    match node {
        Node::Int(n) => *n,
        other => panic!("not an integer: {other:?}"),
    }
}

fn word(node: &Node) -> &str {
    match node {
        Node::Named(w, None) => w,
        other => panic!("not a word: {other:?}"),
    }
}

/// `Blast-filter-program` and `Na-strand` names.
fn program_number(node: &Node) -> i32 {
    match node {
        Node::Int(n) => *n as i32,
        _ => match word(node) {
            "not-set" => 0,
            "dust" => 10,
            "seg" => 20,
            "windowmasker" => 30,
            "repeat" => 40,
            "other" => 100,
            "max" => 255,
            other => panic!("unknown program {other}"),
        },
    }
}

fn strand_number(node: &Node) -> u8 {
    match word(node) {
        "unknown" => 0,
        "plus" => 1,
        "minus" => 2,
        "both" => 3,
        "both-rev" => 4,
        "other" => 255,
        other => panic!("unknown strand {other}"),
    }
}

fn interval_of(fields: &[Node]) -> (MaskInterval, MaskSeqId) {
    // `id local id N` names a sequence of the database; `id local str "x"`
    // is a Seq-id that no sequence of a database without -parse_seqids has.
    let id = match field(fields, "id") {
        Some(Node::Named(kind, Some(inner))) if kind == "local" => {
            if let Some(n) = field(std::slice::from_ref(&**inner), "id") {
                MaskSeqId::LocalId(int(n) as i32)
            } else if let Some(Node::Str(text)) = field(std::slice::from_ref(&**inner), "str") {
                MaskSeqId::Other(text.clone().into_bytes())
            } else {
                panic!("local id: {inner:?}");
            }
        }
        other => panic!("id: {other:?}"),
    };
    let interval = MaskInterval {
        from: int(field(fields, "from").unwrap()) as u32,
        to: int(field(fields, "to").unwrap()) as u32,
        strand: field(fields, "strand").map(strand_number),
    };
    (interval, id)
}

fn loc_of(node: &Node) -> MaskLoc {
    let Node::Named(kind, Some(body)) = node else {
        panic!("not a Seq-loc: {node:?}");
    };
    let (intervals, id) = match kind.as_str() {
        "int" => {
            let (interval, id) = interval_of(seq(body));
            (vec![interval], id)
        }
        "packed-int" => {
            let pairs: Vec<_> = seq(body).iter().map(|n| interval_of(seq(n))).collect();
            let id = pairs[0].1.clone();
            (pairs.into_iter().map(|(iv, _)| iv).collect(), id)
        }
        other => panic!("unknown Seq-loc {other}"),
    };
    MaskLoc {
        id: Some(id),
        intervals,
    }
}

/// The `Blast-db-mask-info` objects of an ASN.1 text file.
fn parse_mask_text(text: &str) -> Vec<MaskInfo> {
    let mut parser = Parser {
        tokens: tokenize(text),
        pos: 0,
    };
    let mut objects = Vec::new();
    while parser.peek().is_some() {
        assert_eq!(parser.next(), Tok::Word("Blast-db-mask-info".to_string()));
        assert_eq!(parser.next(), Tok::Assign);
        let value = parser.value();
        let items = seq(&value);
        let Node::Str(options) = field(items, "algo-options").unwrap() else {
            panic!("algo-options");
        };
        let masks = seq(field(items, "masks").unwrap());
        let locs = seq(field(masks, "masks").unwrap());
        objects.push(MaskInfo {
            algo_id: int(field(items, "algo-id").unwrap()) as i32,
            algo_program: program_number(field(items, "algo-program").unwrap()),
            algo_options: options.clone().into_bytes(),
            masks: MaskList {
                masks: locs.iter().map(loc_of).collect(),
                more: word(field(masks, "more").unwrap()) == "TRUE",
            },
        });
    }
    objects
}

#[test]
fn the_text_reader_reads_the_fixture_masks() {
    let dir = seqdb_fixtures().join("mask_writer");
    let dust = parse_mask_text(&std::fs::read_to_string(dir.join("dust.asnt")).unwrap());
    assert_eq!(dust.len(), 1);
    assert_eq!(dust[0].algo_program, 10);
    assert_eq!(dust[0].algo_options, b"window=64; level=20; linker=1");
    let intervals: Vec<(Option<MaskSeqId>, Vec<(u32, u32)>)> = dust[0]
        .masks
        .masks
        .iter()
        .map(|m| {
            (
                m.id.clone(),
                m.intervals.iter().map(|i| (i.from, i.to)).collect(),
            )
        })
        .collect();
    assert_eq!(
        intervals[0],
        (Some(MaskSeqId::LocalId(1)), vec![(199, 280)])
    );
    assert_eq!(
        intervals[3],
        (Some(MaskSeqId::LocalId(6)), vec![(150, 219), (620, 689)])
    );
    let strand = parse_mask_text(&std::fs::read_to_string(dir.join("strand.asnt")).unwrap());
    assert_eq!(strand[0].masks.masks[0].intervals[0].strand, Some(2));
    assert_eq!(strand[0].masks.masks[0].intervals[1].strand, Some(1));
    let objects = parse_mask_text(&std::fs::read_to_string(dir.join("objs.asnt")).unwrap());
    assert_eq!(objects.len(), 2);
    assert_eq!(objects[1].algo_options, b"-ignored 2");
    let none = parse_mask_text(&std::fs::read_to_string(dir.join("none.asnt")).unwrap());
    assert!(none[0].masks.masks.is_empty());
}

// ---------------------------------------------------------------------
// The combination of Seq-locs and the registry on their own
// ---------------------------------------------------------------------

fn loc(id: i32, ranges: &[(u32, u32)]) -> MaskLoc {
    MaskLoc {
        id: Some(MaskSeqId::LocalId(id)),
        intervals: ranges
            .iter()
            .map(|&(from, to)| MaskInterval {
                from,
                to,
                strand: None,
            })
            .collect(),
    }
}

fn stranded(id: i32, ranges: &[(u32, u32, u8)]) -> MaskLoc {
    MaskLoc {
        id: Some(MaskSeqId::LocalId(id)),
        intervals: ranges
            .iter()
            .map(|&(from, to, strand)| MaskInterval {
                from,
                to,
                strand: Some(strand),
            })
            .collect(),
    }
}

/// The ranges of sequence `id` for the algorithm `algo` after the Seq-locs
/// are inserted in order.
fn combined(algo: i32, id: i32, locs: &[MaskLoc]) -> Result<Vec<(u32, u32)>, WriteDbError> {
    let mut set = MaskedRangeSet::new();
    for l in locs {
        set.insert(algo, &MaskSeqId::LocalId(id), l)?;
    }
    let ranges = set.get_ranges(&[MaskSeqId::LocalId(id)])?;
    assert_eq!(ranges.len(), 1);
    Ok(ranges[0].offsets.clone())
}

/// What NCBI's makeblastdb stored (blastdbcmd `%m`) for the same Seq-locs:
/// the oracle runs behind `fixtures` in `make_mask_writer_fixtures.py` and
/// the experiments of the M4 writer brief.
#[test]
fn seq_locs_of_one_sequence_are_merged_as_ncbi_does() {
    // One Seq-loc is kept as given: unsorted, overlapping, repeated.
    assert_eq!(
        combined(30, 1, &[loc(1, &[(300, 309), (0, 9), (5, 20), (8, 8)])]).unwrap(),
        [(300, 310), (0, 10), (5, 21), (8, 9)]
    );
    assert_eq!(
        combined(30, 1, &[loc(1, &[(100, 150), (100, 150)])]).unwrap(),
        [(100, 151), (100, 151)]
    );
    // A single-base range, a range past 31 bits, and a reversed one.
    assert_eq!(combined(30, 1, &[loc(1, &[(7, 7)])]).unwrap(), [(7, 8)]);
    assert_eq!(combined(30, 1, &[loc(1, &[(6, 5)])]).unwrap(), [(6, 6)]);
    assert_eq!(
        combined(30, 1, &[loc(1, &[(0, u32::MAX)])]).unwrap(),
        [(0, 0)]
    );
    // Two Seq-locs: sorted, overlapping and abutting ranges joined.
    assert_eq!(
        combined(
            30,
            1,
            &[
                loc(1, &[(300, 309), (0, 9)]),
                loc(1, &[(5, 20), (310, 320), (100, 100)])
            ]
        )
        .unwrap(),
        [(0, 21), (100, 101), (300, 321)]
    );
    // Contained ranges and abutting single bases.
    assert_eq!(
        combined(
            30,
            1,
            &[
                loc(1, &[(50, 60)]),
                loc(1, &[(0, 100)]),
                loc(1, &[(20, 30)]),
                loc(1, &[(101, 101)]),
                loc(1, &[(500, 510)]),
                loc(1, &[(505, 505)])
            ]
        )
        .unwrap(),
        [(0, 102), (500, 511)]
    );
    // Gaps of one base stay apart.
    assert_eq!(
        combined(
            30,
            1,
            &[loc(1, &[(0, 0)]), loc(1, &[(1, 1)]), loc(1, &[(3, 3)])]
        )
        .unwrap(),
        [(0, 2), (3, 4)]
    );
    // An empty first Seq-loc still makes the second one sorted and merged.
    assert_eq!(
        combined(30, 1, &[loc(1, &[]), loc(1, &[(30, 40), (0, 5), (3, 8)])]).unwrap(),
        [(0, 9), (30, 41)]
    );
}

#[test]
fn reverse_strand_ranges_are_merged_apart_and_come_last() {
    const PLUS: u8 = 1;
    const MINUS: u8 = 2;
    const BOTH: u8 = 3;
    const BOTH_REV: u8 = 4;
    // One minus-strand Seq-loc alone stays in the order it was given.
    assert_eq!(
        combined(30, 1, &[stranded(1, &[(300, 309, MINUS), (0, 9, MINUS)])]).unwrap(),
        [(300, 310), (0, 10)]
    );
    // Two minus-strand Seq-locs: by decreasing end.
    assert_eq!(
        combined(
            30,
            1,
            &[
                stranded(1, &[(300, 309, MINUS), (0, 9, MINUS)]),
                stranded(1, &[(5, 20, MINUS), (310, 320, MINUS)])
            ]
        )
        .unwrap(),
        [(300, 321), (0, 21)]
    );
    // Plus and minus ranges are not joined with each other.
    assert_eq!(
        combined(
            30,
            1,
            &[
                stranded(1, &[(300, 309, MINUS), (0, 9, PLUS)]),
                stranded(1, &[(5, 20, MINUS), (310, 320, PLUS), (7, 8, MINUS)])
            ]
        )
        .unwrap(),
        [(0, 10), (310, 321), (300, 310), (5, 21)]
    );
    // `both` and unknown strands count as plus, `both-rev` as minus.
    assert_eq!(
        combined(
            30,
            1,
            &[
                stranded(1, &[(0, 9, BOTH)]),
                stranded(1, &[(10, 19, PLUS)]),
                loc(1, &[(20, 29)]),
                stranded(1, &[(40, 49, BOTH_REV)]),
                stranded(1, &[(30, 39, BOTH_REV)])
            ]
        )
        .unwrap(),
        [(0, 30), (30, 50)]
    );
    assert_eq!(
        combined(
            30,
            1,
            &[
                stranded(1, &[(0, 9, MINUS)]),
                stranded(1, &[(5, 19, MINUS)])
            ]
        )
        .unwrap(),
        [(0, 20)]
    );
}

#[test]
fn combinations_with_empty_or_whole_ranges_are_refused() {
    // NCBI drops an empty range (end before start) from a merged Seq-loc;
    // LOSAT does not port that.
    let error = combined(30, 2, &[loc(2, &[(6, 5)]), loc(2, &[(20, 25)])]).unwrap_err();
    assert!(
        error.to_string().contains("LOSAT does not support"),
        "{error}"
    );
    assert!(combined(30, 2, &[loc(2, &[(0, u32::MAX)]), loc(2, &[(1, 2)])]).is_err());
    // A single Seq-loc passes whatever it holds.
    assert_eq!(
        combined(30, 2, &[loc(2, &[(6, 5), (20, 25)])]).unwrap(),
        [(6, 6), (20, 26)]
    );
    // A negative algorithm id cannot index NCBI's vector.
    let mut set = MaskedRangeSet::new();
    assert!(set
        .insert(-1, &MaskSeqId::LocalId(1), &loc(1, &[(1, 2)]))
        .is_err());
}

/// NCBI combines the Seq-locs of any Seq-id and only looks up those of the
/// sequences it writes, so a combination LOSAT refuses is an error only for a
/// sequence that is looked up.
#[test]
fn a_refused_combination_is_an_error_only_when_it_is_looked_up() {
    let mut set = MaskedRangeSet::new();
    let other = MaskSeqId::Other(b"x".to_vec());
    let far = MaskSeqId::LocalId(1_000_000);
    let near = MaskSeqId::LocalId(2);
    for id in [&other, &far, &near] {
        set.insert(30, id, &loc(2, &[(6, 5)])).unwrap();
        set.insert(30, id, &loc(2, &[(20, 25)])).unwrap();
        // More Seq-locs after the refused one change nothing.
        set.insert(30, id, &loc(2, &[(1, 2)])).unwrap();
    }
    set.insert(31, &MaskSeqId::LocalId(1), &loc(1, &[(0, 9)]))
        .unwrap();
    // Sequence 1 has no refused combination (algorithm 30 has no Seq-loc of it).
    let ranges = set.get_ranges(&[MaskSeqId::LocalId(1)]).unwrap();
    assert_eq!(ranges[1].offsets, [(0, 10)]);
    // A sequence with one does.
    let error = set.get_ranges(&[near]).unwrap_err();
    assert!(
        error.to_string().contains("LOSAT does not support"),
        "{error}"
    );
}

#[test]
fn ranges_follow_the_first_insert_of_each_algorithm() {
    let mut set = MaskedRangeSet::new();
    let id1 = MaskSeqId::LocalId(1);
    let id2 = MaskSeqId::LocalId(2);
    set.insert(31, &id1, &loc(1, &[(1, 2)])).unwrap();
    set.insert(30, &id2, &loc(2, &[(3, 4)])).unwrap();
    set.insert(31, &id2, &loc(2, &[(5, 6)])).unwrap();
    let ranges = set.get_ranges(&[id1.clone()]).unwrap().to_vec();
    assert_eq!(
        ranges,
        [
            MaskedRanges {
                algorithm_id: 31,
                offsets: vec![(1, 3)]
            },
            MaskedRanges {
                algorithm_id: 30,
                offsets: vec![]
            }
        ]
    );
    assert!(!ranges_are_empty(&ranges));
    // No id of the sequence, no ranges; a different kind of Seq-id never matches.
    let other = MaskSeqId::Other(vec![1, 2, 3]);
    assert!(ranges_are_empty(set.get_ranges(&[other]).unwrap()));
    assert!(ranges_are_empty(set.get_ranges(&[]).unwrap()));
    // Several ids of one sequence: their Seq-locs are added.
    let both = set.get_ranges(&[id1, id2]).unwrap();
    assert_eq!(both[0].offsets, [(1, 3), (5, 7)]);
    assert_eq!(both[1].offsets, [(3, 5)]);
}

fn algo(registry: &mut MaskInfoRegistry, program: i32, options: &str, name: &str) -> String {
    match registry.add_program(program, options.as_bytes(), name.as_bytes()) {
        Ok(id) => id.to_string(),
        Err(error) => error.to_string(),
    }
}

/// The ids and the messages of makeblastdb 2.17.0 (`reg/` experiments of the M4
/// writer brief, the cases `ids_*`, `dup*` and `invalid*` of
/// `make_mask_writer_fixtures.py`).
#[test]
fn the_registry_assigns_ids_as_ncbi_does() {
    let mut r = MaskInfoRegistry::new();
    assert!(!r.is_registered(10));
    assert_eq!(algo(&mut r, FILTER_PROGRAM_DUST, "", "d0"), "10");
    assert_eq!(algo(&mut r, FILTER_PROGRAM_DUST, "-a 1", "d1"), "11");
    assert_eq!(algo(&mut r, FILTER_PROGRAM_DUST, "-b 2", "d2"), "12");
    assert!(r.is_registered(10) && r.is_registered(12) && !r.is_registered(13));
    // A second set of default options is refused whatever the file is called.
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_DUST, "", "other name"),
        "BLAST Database creation error: Masking algorithm with default arguments already provided"
    );
    // The same program, name and options is a duplicate.
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_DUST, "-a 1", "d1"),
        "BLAST Database creation error: Duplicate masking algorithm found."
    );
    assert_eq!(algo(&mut r, FILTER_PROGRAM_SEG, "", "s"), "20");
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_WINDOWMASKER, "-x 1", "w1"),
        "31"
    );
    assert_eq!(algo(&mut r, FILTER_PROGRAM_WINDOWMASKER, "", "w0"), "30");
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_WINDOWMASKER, "-x 2", "w2"),
        "32"
    );
    // `repeat` and `other` take the first free id from their own value on;
    // their options do not matter.
    assert_eq!(algo(&mut r, FILTER_PROGRAM_REPEAT, "", "r0"), "40");
    assert_eq!(algo(&mut r, FILTER_PROGRAM_REPEAT, "", "r1"), "41");
    assert_eq!(algo(&mut r, FILTER_PROGRAM_OTHER, "", "o0"), "100");
    assert_eq!(algo(&mut r, FILTER_PROGRAM_OTHER, "-x:1", "o1"), "101");
    // The ids of `-mask_id` come from the same range.
    assert_eq!(r.add_id(b"alpha").unwrap(), 102);
    assert_eq!(
        r.add_id(b"alpha").unwrap_err().to_string(),
        "BLAST Database creation error: Duplicate masking algorithm found."
    );
    // Programs that are not filtering programs.
    for (program, shown) in [(0, "0"), (55, "55"), (255, "255"), (-1, "-1")] {
        assert_eq!(
            algo(&mut r, program, "", "x"),
            format!("BLAST Database creation error: Invalid filtering program: {shown}")
        );
    }
}

#[test]
fn the_registry_runs_out_of_ids_in_each_range() {
    // dust: 10 for the defaults and 11 to 19 for others.
    let mut r = MaskInfoRegistry::new();
    for i in 0..10 {
        let options = if i == 0 {
            String::new()
        } else {
            format!("-x {i}")
        };
        assert_eq!(
            r.add_program(FILTER_PROGRAM_DUST, options.as_bytes(), b"f")
                .unwrap(),
            10 + i
        );
    }
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_DUST, "-x 10", "f"),
        "BLAST Database creation error: Too many IDs in range 11-20"
    );
    // repeat: 40 to 99.
    let mut r = MaskInfoRegistry::new();
    for i in 0..60 {
        assert_eq!(
            r.add_program(FILTER_PROGRAM_REPEAT, b"", format!("r{i}").as_bytes())
                .unwrap(),
            40 + i
        );
    }
    assert_eq!(
        algo(&mut r, FILTER_PROGRAM_REPEAT, "", "r60"),
        "BLAST Database creation error: Too many IDs in range 40-100"
    );
    // -mask_id: 100 to 254.
    let mut r = MaskInfoRegistry::new();
    for i in 0..155 {
        assert_eq!(r.add_id(format!("m{i}").as_bytes()).unwrap(), 100 + i);
    }
    assert_eq!(
        r.add_id(b"m155").unwrap_err().to_string(),
        "BLAST Database creation error: Too many IDs in range 100-255"
    );
}

// ---------------------------------------------------------------------
// WriteDb alone
// ---------------------------------------------------------------------

fn bioseq(local_id: usize, length: usize) -> crate::writedb::bioseq::Bioseq {
    crate::writedb::bioseq::Bioseq {
        local_id: local_id.to_string(),
        title: format!("s{local_id}").into_bytes(),
        is_aa: true,
        length: length as u64,
        seq_data: Some(crate::writedb::bioseq::SeqData::Ncbieaa(vec![b'A'; length])),
    }
}

#[test]
fn set_mask_data_checks_what_ncbi_checks() {
    let tmp = TempDir::new("mask-checks");
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 0);
    // No sequence: "Need sequence data."
    let one = [MaskedRanges {
        algorithm_id: 20,
        offsets: vec![(0, 5)],
    }];
    assert_eq!(
        db.set_mask_data(&one).unwrap_err().to_string(),
        "BLAST Database creation error: Need sequence data."
    );
    let id = db
        .register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"seg")
        .unwrap();
    assert_eq!(id, 20);
    db.add_sequence(bioseq(1, 10)).unwrap();
    // An id nobody registered.
    let unknown = [MaskedRanges {
        algorithm_id: 21,
        offsets: vec![(0, 5)],
    }];
    assert_eq!(
        db.set_mask_data(&unknown).unwrap_err().to_string(),
        "BLAST Database creation error: Error: Algorithm IDs must be registered before use. \
         Unknown algorithm ID = 21"
    );
    // A range with end before start, or past the end of the sequence.
    for offsets in [vec![(6, 5)], vec![(0, 11)], vec![(0, 5), (3, 12)]] {
        let bad = [MaskedRanges {
            algorithm_id: 20,
            offsets,
        }];
        assert_eq!(
            db.set_mask_data(&bad).unwrap_err().to_string(),
            "BLAST Database creation error: Error: Masked data offsets out of bounds."
        );
    }
    // Ranges of an unregistered id with no offsets are skipped.
    let skipped = [MaskedRanges {
        algorithm_id: 77,
        offsets: vec![],
    }];
    db.set_mask_data(&skipped).unwrap();
    // The blobs of a sequence are set once.
    assert_eq!(
        db.set_mask_data(&one).unwrap_err().to_string(),
        "BLAST Database creation error: Error: Already have blob for this sequence and column"
    );
    db.abandon();
}

#[test]
fn a_column_is_not_created_in_an_open_volume() {
    let tmp = TempDir::new("mask-late");
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 0);
    db.add_sequence(bioseq(1, 10)).unwrap();
    db.add_sequence(bioseq(2, 10)).unwrap();
    // The first sequence is published, so a volume exists.
    let error = db
        .register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"seg")
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("before the first sequence is published"),
        "{error}"
    );
    db.abandon();
}

/// Registering after the first sequence was added, before it is published, is
/// what NCBI allows.
#[test]
fn algorithms_are_registered_until_the_first_volume_is_written() {
    let tmp = TempDir::new("mask-early");
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 0);
    db.set_column_date(0, b"10/20/2026 12:00:00");
    db.add_sequence(bioseq(1, 10)).unwrap();
    assert_eq!(
        db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"seg")
            .unwrap(),
        20
    );
    db.set_mask_data(&[MaskedRanges {
        algorithm_id: 20,
        offsets: vec![(2, 4)],
    }])
    .unwrap();
    db.close().unwrap();
    let written = files_named(&tmp.0, "db.");
    assert_eq!(
        written,
        ["db.paa", "db.pab", "db.pac", "db.phr", "db.pin", "db.psq"]
    );
    let reader = SeqDb::open_with_masks(&tmp.base("db"), SeqType::Protein).unwrap();
    assert_eq!(
        reader.mask_ranges(0, 20).unwrap(),
        MaskRanges::Ranges(vec![(2, 4)])
    );
}

/// The date of each column index file is the one a test sets, or the time at
/// which the volume was created, as `CTime::AsString()` writes it.
#[test]
fn column_dates_come_from_the_hook_or_the_clock() {
    let tmp = TempDir::new("mask-dates");
    // 100 bytes per file: every sequence is a volume.
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 100);
    db.set_column_date(1, b"01/02/2003 04:05:06");
    db.set_column_date_fallback(b"09/09/2009 09:09:09");
    let id = db
        .register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"seg")
        .unwrap();
    for i in 1..=3 {
        db.add_sequence(bioseq(i, 40)).unwrap();
        db.set_mask_data(&[MaskedRanges {
            algorithm_id: id,
            offsets: vec![(0, 5)],
        }])
        .unwrap();
    }
    db.close().unwrap();
    let dates: Vec<String> = (0..3)
        .map(|i| column_date(&tmp.0.join(format!("db.0{i}.paa"))).unwrap())
        .collect();
    assert_eq!(dates[1], "01/02/2003 04:05:06");
    for date in [&dates[0], &dates[2]] {
        let b = date.as_bytes();
        assert_eq!(b.len(), 19, "{date}");
        assert!(
            b.iter().enumerate().all(|(i, c)| match i {
                2 | 5 => *c == b'/',
                10 => *c == b' ',
                13 | 16 => *c == b':',
                _ => c.is_ascii_digit(),
            }),
            "{date}"
        );
    }
}

// ---------------------------------------------------------------------
// Databases rebuilt through BuildDatabase
// ---------------------------------------------------------------------

/// What makeblastdb's `x_ProcessMaskData` and `x_ProcessInputData` are given.
struct Case {
    name: String,
    protein: bool,
    /// The FASTA files in the order of `-in`.
    inputs: Vec<PathBuf>,
    /// The mask files in the order of `-mask_data`: the name as given on the
    /// command line and the ASN.1 text.
    masks: Vec<(String, String)>,
    /// `-mask_id` and `-mask_desc`, split at the commas.
    ids: Option<Vec<String>>,
    descs: Option<Vec<String>>,
    max_file_size: u64,
    title: String,
}

/// `x_ProcessMaskData` as far as the files are concerned: registers the
/// algorithm of each file and adds its masks to the set.
fn process_mask_data(db: &mut BuildDatabase<'_>, case: &Case) -> Result<(), WriteDbError> {
    for (i, (name, text)) in case.masks.iter().enumerate() {
        let mut algo_id: Option<i32> = None;
        for info in parse_mask_text(text) {
            let id = match algo_id {
                Some(id) => id,
                None => {
                    let id = if let Some(ids) = &case.ids {
                        let desc = case.descs.as_ref().map_or("", |d| d[i].as_str());
                        db.register_masking_algorithm_id(
                            ids[i].as_bytes(),
                            desc.as_bytes(),
                            &info.algo_options,
                        )?
                    } else {
                        db.register_masking_algorithm(
                            info.algo_program,
                            &info.algo_options,
                            name.as_bytes(),
                        )?
                    };
                    algo_id = Some(id);
                    id
                }
            };
            assert!(
                !info.masks.more,
                "continuations are not used by the fixtures"
            );
            if db.mask_data_source_mut().is_none() && !info.masks.masks.is_empty() {
                db.set_mask_data_source(MaskedRangeSet::new());
            }
            for mask in &info.masks.masks {
                let seq_id = mask
                    .id
                    .as_ref()
                    .expect("the fixtures give every mask an id");
                if let Some(set) = db.mask_data_source_mut() {
                    set.insert(id, seq_id, mask)?;
                }
            }
        }
    }
    Ok(())
}

/// `x_BuildDatabase` from the mask files on: the result, and the text
/// makeblastdb posts to the standard error stream (the diagnostics of the
/// build, the message of the exception that ends it, and the message of the
/// destructor of `CBuildDatabase`, in that order).
fn run_case(
    case: &Case,
    out_base: &str,
    index_date: &str,
    column_dates: &[Option<String>],
) -> (Result<(), WriteDbError>, Vec<u8>) {
    let mut log = Vec::new();
    let mut messages = Vec::new();
    let result = (|| {
        let mut warn = |line: &[u8]| {
            messages.extend_from_slice(line);
            messages.push(b'\n');
        };
        let warn_ref: Warn<'_> = &mut warn;
        let mut db = BuildDatabase::new(
            out_base,
            &case.title,
            case.protein,
            &mut log,
            &time_of(index_date),
        )?;
        for (i, date) in column_dates.iter().enumerate() {
            if let Some(date) = date {
                db.set_column_date(i, date.as_bytes());
            }
        }
        db.set_max_file_size(case.max_file_size);
        let mut result = process_mask_data(&mut db, case);
        if result.is_ok() {
            for input in &case.inputs {
                let file = std::fs::File::open(input).unwrap();
                result = db.add_fasta(file, &mut *warn_ref);
                if result.is_err() {
                    break;
                }
            }
        }
        if result.is_ok() {
            result = db.end_build(false, &mut *warn_ref).map(|_| ());
        }
        if let Err(error) = &result {
            // CATCH_ALL
            warn_ref(error.to_string().as_bytes());
        }
        db.post_destruction_diagnostics(&mut *warn_ref);
        result
    })();
    (result, messages)
}

const DB_EXTENSIONS: [&str; 14] = [
    "nin", "nsq", "nhr", "naa", "nab", "nac", "nal", "pin", "psq", "phr", "paa", "pab", "pac",
    "pal",
];

fn is_database_file(name: &str) -> bool {
    name.rsplit_once('.')
        .is_some_and(|(_, ext)| DB_EXTENSIONS.contains(&ext))
}

/// The names of the volumes of a database in `dir`: those the alias file
/// lists, or the database itself.
fn volume_bases(dir: &Path, name: &str, protein: bool) -> Vec<String> {
    let alias = dir.join(format!("{name}.{}al", if protein { 'p' } else { 'n' }));
    match std::fs::read_to_string(alias) {
        Ok(text) => text
            .lines()
            .find_map(|line| line.strip_prefix("DBLIST "))
            .expect("an alias file lists its volumes")
            .split(' ')
            .map(str::to_string)
            .collect(),
        Err(_) => vec![name.to_string()],
    }
}

/// The date in a column index file: after the 32 bytes of fixed fields, the
/// title and the date, each a variable-length integer and the bytes.
fn column_date(path: &Path) -> Option<String> {
    let bytes = std::fs::read(path).ok()?;
    let mut at = 32;
    let mut strings = Vec::new();
    for _ in 0..2 {
        // Lengths of 63 or less are one byte.
        let len = usize::from(bytes[at]);
        assert!(len < 64, "{path:?}");
        strings.push(String::from_utf8(bytes[at + 1..at + 1 + len].to_vec()).unwrap());
        at += 1 + len;
    }
    assert_eq!(strings[0], "BlastDb/MaskData");
    Some(strings.remove(1))
}

fn stored_title(index_file: &Path, seqtype: SeqType) -> String {
    let bytes = std::fs::read(index_file).unwrap();
    let index = crate::seqdb::index::IndexFile::parse(&bytes, seqtype, index_file).unwrap();
    String::from_utf8(index.title).unwrap()
}

/// Every file `out_dir/name.*` equals the fixture's file of the same
/// extension and the fixture has no other file.
fn assert_same_files(out_dir: &Path, name: &str, fixture_dir: &Path, fixture: &str) {
    let suffix = |file: &str, stem: &str| file[stem.len()..].to_string();
    let written = files_named(out_dir, &format!("{name}."));
    let expected: Vec<String> = files_named(fixture_dir, &format!("{fixture}."))
        .iter()
        .filter(|file| is_database_file(file))
        .map(|file| suffix(file, fixture))
        .collect();
    let actual: Vec<String> = written.iter().map(|file| suffix(file, name)).collect();
    assert_eq!(actual, expected, "the files of {fixture}");
    for file in &written {
        let rest = suffix(file, name);
        let expected = std::fs::read(fixture_dir.join(format!("{fixture}{rest}"))).unwrap();
        let actual = std::fs::read(out_dir.join(file)).unwrap();
        assert!(
            actual == expected,
            "{fixture}{rest}: {} bytes written, {} expected; first difference at byte {:?}",
            actual.len(),
            expected.len(),
            actual.iter().zip(&expected).position(|(a, b)| a != b)
        );
    }
}

fn ranges_text(result: &MaskRanges) -> String {
    match result {
        MaskRanges::Ranges(ranges) if ranges.is_empty() => "none".to_string(),
        MaskRanges::Ranges(ranges) => ranges.iter().map(|(a, b)| format!("{a}-{b};")).collect(),
        other => panic!("{other:?}"),
    }
}

/// The database `out_base` reads back as NCBI's blastdbcmd reads the
/// fixture: the table of algorithms and the ranges of every sequence.
fn assert_reads_back(out_base: &str, protein: bool, expected_dir: &Path, fixture: &str) {
    let seqtype = if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let db = SeqDb::open_with_masks(out_base, seqtype).unwrap();
    let algorithms_file = expected_dir.join(format!("{fixture}.algorithms.txt"));
    let table = db.mask_algorithms().unwrap();
    if algorithms_file.exists() {
        let expected = std::fs::read_to_string(&algorithms_file).unwrap();
        assert_eq!(table.descriptions().unwrap(), expected, "{fixture}");
    } else {
        assert!(table.ids().is_empty(), "{fixture}");
    }
    let rows = std::fs::read_to_string(expected_dir.join(format!("{fixture}.masks.tsv")))
        .unwrap_or_default();
    for line in rows.lines() {
        let mut fields = line.split('\t');
        let algo: i32 = fields.next().unwrap().parse().unwrap();
        let oid: u32 = fields.next().unwrap().parse().unwrap();
        let expected = fields.next().unwrap();
        let result = db.mask_ranges(oid, algo).unwrap();
        assert_eq!(
            ranges_text(&result),
            expected,
            "{fixture}: algorithm {algo}, sequence {oid}"
        );
    }
}

/// Rebuilds the NCBI database `fixture` of `fixture_dir` from `case` and
/// compares every file. The index date, the title and the dates of the
/// column files are the ones in NCBI's files.
fn assert_rebuilds(case: &Case, fixture_dir: &Path, fixture: &str, expected_stderr: &[u8]) {
    let letter = if case.protein { 'p' } else { 'n' };
    let seqtype = if case.protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let bases = volume_bases(fixture_dir, fixture, case.protein);
    let first_index = fixture_dir.join(format!("{}.{letter}in", bases[0]));
    let index_date = stored_date(&first_index, seqtype);
    let column_dates: Vec<Option<String>> = bases
        .iter()
        .map(|base| column_date(&fixture_dir.join(format!("{base}.{letter}aa"))))
        .collect();
    let tmp = TempDir::new(&case.name);
    let (result, messages) = run_case(case, &tmp.base(fixture), &index_date, &column_dates);
    result.unwrap_or_else(|error| panic!("{fixture}: {error}"));
    assert_eq!(
        String::from_utf8_lossy(&messages),
        String::from_utf8_lossy(expected_stderr),
        "{fixture}: standard error"
    );
    assert_same_files(&tmp.0, fixture, fixture_dir, fixture);
    assert_eq!(
        stored_title(&first_index, seqtype),
        case.title,
        "{fixture}: the title is the one of the fixture"
    );
    let expected_dir = fixture_dir.parent().unwrap().join("expected");
    if expected_dir.join(format!("{fixture}.masks.tsv")).exists() {
        assert_reads_back(&tmp.base(fixture), case.protein, &expected_dir, fixture);
    }
}

fn read(dir: &Path, name: &str) -> String {
    std::fs::read_to_string(dir.join(name)).unwrap_or_else(|error| panic!("{name}: {error}"))
}

/// The databases of `make_mask_fixtures.py` in format version 4. The mask
/// files are named as that script gave them to makeblastdb.
fn fixture_case(name: &str) -> Case {
    let seqdb = seqdb_fixtures();
    let writer = seqdb.join("mask_writer");
    let dust = ("dust.asnb", read(&writer, "dust.asnt"));
    let seg = ("seg.asnb", read(&writer, "seg.asnt"));
    let wm = ("wm.asnt", read(&seqdb, "mask_wm.asnt"));
    let wm_opts = ("wm_opts.asnt", read(&seqdb, "mask_wm_opts.asnt"));
    let wm_colon = ("wm_colon.asnt", read(&seqdb, "mask_wm_colon.asnt"));
    let (protein, masks, ids, descs, max_file_size): (
        bool,
        Vec<(&str, String)>,
        Option<&[&str]>,
        Option<&[&str]>,
        u64,
    ) = match name {
        "mask_nucl_v4" => (false, vec![dust, wm], None, None, 3_000_000_000),
        "mask_nucl_mv" => (false, vec![dust, wm], None, None, 400),
        "mask_nucl_ids" => (
            false,
            vec![dust, wm],
            Some(&["dustX", "wmY"]),
            Some(&["dust: from dustmasker", "hand written"]),
            3_000_000_000,
        ),
        "mask_nucl_z" => (false, vec![wm_opts], Some(&["wmZ"]), None, 3_000_000_000),
        "mask_nucl_y" => (false, vec![wm_opts], Some(&["wmY"]), None, 3_000_000_000),
        "mask_nucl_opts" => (false, vec![wm_opts, wm], None, None, 3_000_000_000),
        "mask_nucl_colon" => (false, vec![wm_colon], None, None, 3_000_000_000),
        "mask_nucl_colon_id" => (false, vec![wm_colon], Some(&["c:id"]), None, 3_000_000_000),
        "mask_prot_v4" => (true, vec![seg], None, None, 3_000_000_000),
        "mask_nucl_plain" => (false, vec![], None, None, 3_000_000_000),
        other => panic!("no fixture {other}"),
    };
    let strings = |list: &[&str]| list.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    // makeblastdb's default title is the name of the input file; the two
    // fixtures with a title of their own state it in their index file.
    let fasta = if protein {
        "mask_prot.faa"
    } else {
        "mask_nucl.fna"
    };
    let title = match name {
        "mask_nucl_v4" => "LOSAT mask fixture, nucleotide, version 4",
        "mask_nucl_mv" => "masked nucleotide in volumes",
        _ => fasta,
    };
    Case {
        name: name.to_string(),
        protein,
        inputs: vec![seqdb.join(fasta)],
        masks: masks
            .into_iter()
            .map(|(file, text)| (file.to_string(), text))
            .collect(),
        ids: ids.map(strings),
        descs: descs.map(strings),
        max_file_size,
        title: title.to_string(),
    }
}

/// Every database that `make_mask_fixtures.py` built in format version 4,
/// from its FASTA file and mask files.
#[test]
fn the_writer_rebuilds_makeblastdb_s_masked_databases_byte_for_byte() {
    let fixtures = seqdb_fixtures().join("db");
    for name in [
        "mask_nucl_v4",
        "mask_nucl_mv",
        "mask_nucl_ids",
        "mask_nucl_z",
        "mask_nucl_y",
        "mask_nucl_opts",
        "mask_nucl_colon",
        "mask_nucl_colon_id",
        "mask_prot_v4",
        "mask_nucl_plain",
    ] {
        assert_rebuilds(&fixture_case(name), &fixtures, name, b"");
    }
}

/// The cases of `make_mask_writer_fixtures.py`.
fn writer_cases() -> Vec<Case> {
    let seqdb = seqdb_fixtures();
    let dir = seqdb.join("mask_writer");
    let table = read(&dir, "cases.tsv");
    let resolve = |file: &str| {
        if dir.join(file).exists() {
            dir.join(file)
        } else {
            seqdb.join(file)
        }
    };
    let list = |field: &str| {
        (field != "-").then(|| field.split(',').map(str::to_string).collect::<Vec<_>>())
    };
    table
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 8, "{line}");
            Case {
                name: f[0].to_string(),
                protein: f[1] == "prot",
                inputs: f[2].split(' ').map(resolve).collect(),
                masks: list(f[3])
                    .unwrap_or_default()
                    .iter()
                    .map(|file| (file.clone(), read(&dir, file)))
                    .collect(),
                ids: list(f[4]),
                descs: list(f[5]),
                max_file_size: f[6]
                    .strip_suffix('B')
                    .map_or(3_000_000_000, |n| n.parse().unwrap()),
                title: f[7].to_string(),
            }
        })
        .collect()
}

#[test]
fn the_writer_matches_makeblastdb_on_the_cases_of_the_writer_fixtures() {
    let dir = seqdb_fixtures().join("mask_writer");
    let cases = writer_cases();
    assert!(cases.len() >= 30);
    let mut built = 0;
    let mut failed = 0;
    for case in &cases {
        let name = &case.name;
        let expected_status = read(&dir.join("expected"), &format!("{name}.rc"));
        let expected_stderr =
            std::fs::read(dir.join("expected").join(format!("{name}.stderr"))).unwrap();
        if expected_status.trim() == "0" {
            built += 1;
            assert_rebuilds(case, &dir.join("db"), name, &expected_stderr);
        } else {
            failed += 1;
            assert_eq!(expected_status.trim(), "1", "{name}");
            assert!(files_named(&dir.join("db"), &format!("{name}.")).is_empty());
            let tmp = TempDir::new(name);
            let (result, messages) = run_case(case, &tmp.base(name), "Oct 20, 2026  1:00 PM", &[]);
            assert!(result.is_err(), "{name}: the run should fail");
            assert_eq!(
                String::from_utf8_lossy(&messages),
                String::from_utf8_lossy(&expected_stderr),
                "{name}: standard error"
            );
            // The files of a failed run are removed.
            assert_eq!(
                files_named(&tmp.0, &format!("{name}.")),
                Vec::<String>::new(),
                "{name}: files left behind"
            );
        }
    }
    assert!(
        built >= 20 && failed >= 10,
        "{built} built, {failed} failed"
    );
}

/// A database built without mask data lists and writes no column files.
#[test]
fn a_database_without_mask_data_has_no_column_files() {
    let tmp = TempDir::new("mask-plain");
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 100);
    for i in 1..=5 {
        db.add_sequence(bioseq(i, 40)).unwrap();
    }
    db.close().unwrap();
    for file in db.list_files() {
        assert!(
            !file.ends_with("aa") && !file.ends_with("ab") && !file.ends_with("ac"),
            "{file}"
        );
    }
    assert!(files_named(&tmp.0, "db.")
        .iter()
        .all(|f| is_database_file(f)));
    assert!(files_named(&tmp.0, "db.").iter().all(|f| {
        let ext = f.rsplit('.').next().unwrap();
        ["pin", "psq", "phr", "pal"].contains(&ext)
    }));
}

/// The message of the destructor of `CBuildDatabase`: posted only when mask data
/// were read and no sequence matched.
#[test]
fn the_message_of_a_run_without_matching_masks_is_ncbi_s() {
    let mut log = Vec::new();
    let tmp = TempDir::new("mask-destructor");
    let mut db = BuildDatabase::new(
        &tmp.base("db"),
        "t",
        true,
        &mut log,
        &time_of("Oct 12, 2026  3:22 PM"),
    )
    .unwrap();
    let mut messages: Vec<u8> = Vec::new();
    let mut warn = |line: &[u8]| {
        messages.extend_from_slice(line);
        messages.push(b'\n');
    };
    // No mask data source: no message.
    db.post_destruction_diagnostics(&mut warn);
    assert!(!db.found_matching_masks());
    db.set_mask_data_source(MaskedRangeSet::new());
    db.post_destruction_diagnostics(&mut warn);
    assert_eq!(
        String::from_utf8_lossy(&messages),
        "Error: [makeblastdb] No sequences matched any of the masks provided.\nPlease ensure that \
         the -parse_seqids option is used in the\nfiltering program as well as makeblastdb.\n"
    );
}

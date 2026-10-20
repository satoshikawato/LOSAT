//! Unit tests of the mask-column writer: the combination of Seq-locs, the registry, `WriteDb`
//! alone. The databases rebuilt from the real mask files (`process_mask_data` and
//! `BuildDatabase`, compared with NCBI's files byte by byte) are in `tests_e2e.rs`.

use std::path::Path;

use super::registry::{
    FILTER_PROGRAM_DUST, FILTER_PROGRAM_OTHER, FILTER_PROGRAM_REPEAT, FILTER_PROGRAM_SEG,
    FILTER_PROGRAM_WINDOWMASKER,
};
use super::{ranges_are_empty, MaskInfoRegistry, MaskedRangeSet, MaskedRanges};
use crate::seqdb::{MaskRanges, SeqDb, SeqType};
use crate::writedb::mask_info::{MaskInterval, MaskLoc, MaskSeqId};
use crate::writedb::tests::{files_named, TempDir};
use crate::writedb::{WriteDb, WriteDbError};

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

pub(super) fn bioseq(local_id: usize, length: usize) -> crate::writedb::bioseq::Bioseq {
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

/// The date in a column index file: after the 32 bytes of fixed fields, the
/// title and the date, each a variable-length integer and the bytes.
pub(super) fn column_date(path: &Path) -> Option<String> {
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

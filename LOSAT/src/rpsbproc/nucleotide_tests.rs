//! The refusals of a nucleotide query's SITES and MOTIFS lines: NCBI's `TranslateAll` of fewer
//! than 3 bases or reverse-complement letters, a site residue outside its frame's translation,
//! and the codon and motif positions whose `int` arithmetic overflows.

use super::*;
use crate::rpsbproc::archive::{Align, Segs, StdLoc, StdPair};
use crate::rpsbproc::domains::parse_std_align;
use crate::rpsbproc::segset::Seg;

/// The value of `result`, or a panic naming its refusal.
fn value<T>(result: Result<T, PrintError>) -> T {
    match result {
        Ok(value) => value,
        Err(PrintError::Unsupported(message) | PrintError::Exception(message)) => {
            panic!("{}", String::from_utf8_lossy(&message))
        }
    }
}

/// The message of a refusal.
fn refusal<T>(result: Result<T, PrintError>) -> String {
    match result {
        Err(PrintError::Unsupported(message)) => String::from_utf8(message).unwrap(),
        Err(PrintError::Exception(what)) => panic!("{}", String::from_utf8_lossy(&what)),
        Ok(_) => panic!("not refused"),
    }
}

fn frame(rfidx: u32, translation: &[u8]) -> Frame<'_> {
    Frame {
        rfidx,
        label: RF_IDS[rfidx as usize].to_string(),
        translation,
    }
}

/// A hit of a translated query on one strand (`parse_std_align`).
fn hit(from: u32, to: u32, strand: i32) -> DomAlign {
    let align = Align {
        segs: Segs::Other("std"),
        score: Some(10),
        e_value: Some(1e-5),
        bit_score: Some(40.5),
        num_ident: None,
        score_of_other_type: None,
    };
    let seg = StdPair {
        query: StdLoc::Int {
            from,
            to,
            strand: Some(strand),
        },
        hit: StdLoc::Int {
            from: 0,
            to: 731,
            strand: None,
        },
    };
    parse_std_align(&[seg], &align, 7).unwrap()
}

// biodata_core.cpp:113-145 (quoted at `frames::translate_all`): NCBI's `TranslateAll` erases
// the first byte of a translation made from fewer than 3 bases, or from fewer than 3 letters of
// the reverse complement (which drops what is not an IUPACna letter): an empty string, which
// crashes NCBI (PS16c oracle `short2`, 2 bases: segmentation fault). Refused.
#[test]
fn translations_of_fewer_than_3_bases_or_letters() {
    assert_eq!(
        refusal(translations(b"AC")),
        "a nucleotide query of 2 bases whose reverse complement has 2 IUPACna letters, fewer than 3 in either (NCBI's TranslateAll erases from an empty string, undefined), is not supported by LOSAT's rpsbproc"
    );
    assert_eq!(
        refusal(translations(b"AC-")),
        "a nucleotide query of 3 bases whose reverse complement has 2 IUPACna letters, fewer than 3 in either (NCBI's TranslateAll erases from an empty string, undefined), is not supported by LOSAT's rpsbproc"
    );
    assert_eq!(value(translations(b"ATGAAA"))[0], b"MK");
}

// rpsbproc_main.cpp:1457 (quoted at `codon`): `(char)toupper(dimTranslated[rfidx][iterRP->curr])`.
// At the end of the translation `operator[]` gives its terminating 0; elsewhere outside it NCBI
// reads outside the string (PS16c oracle `natoprot`, whose two 0.5.1 builds print different
// bytes): refused, the residue numbered from 1.
#[test]
fn residues_outside_the_translation() {
    let translation = frame(4, b"mKx");
    assert_eq!(value(translated_residue(&translation, 0)), b'M');
    assert_eq!(value(translated_residue(&translation, 2)), b'X');
    assert_eq!(value(translated_residue(&translation, 3)), 0);
    assert_eq!(
        refusal(translated_residue(&translation, 4)),
        "a site at residue 5 of the translation of reading frame -2 (3 residues; NCBI reads outside it) is not supported by LOSAT's rpsbproc"
    );
    assert_eq!(
        refusal(translated_residue(&translation, -1)),
        "a site at residue 0 of the translation of reading frame -2 (3 residues; NCBI reads outside it) is not supported by LOSAT's rpsbproc"
    );
}

// rpsbproc_main.cpp:949,1456-1457,1631-1633 (quoted at `codon`): the plus-strand codon of a
// residue, `xres + 1` in `int` and `xres +/- 2 + 1` in `unsigned int` (wraps); `xres + 1` at
// `INT_MAX` and `MinusPr2PlusNA`'s `pr - 1` at `INT_MIN` overflow `int`: refused.
#[test]
fn codons_and_their_overflow() {
    // Frame +2: 10 * 3 + 1 = 31. Frame -1 of 5,000 bases: 5000 - 9 * 3 - 3 - 1 = 4969.
    assert_eq!(value(codon(&frame(1, b""), 10, 5000)), (32, 34));
    assert_eq!(value(codon(&frame(3, b""), 10, 5000)), (4970, 4968));
    // Frame +1, residue -1: xres = -3; the second number is (2^32 - 3) + 3 modulo 2^32.
    assert_eq!(value(codon(&frame(0, b""), -1, 5000)), (-2, 0));
    // Frame +2, residue 715,827,882: 715827882 * 3 + 1 = 2^31 - 1.
    assert_eq!(
        refusal(codon(&frame(1, b""), 715_827_882, 5000)),
        "a site at residue 715827883 of reading frame 2 whose position on the query overflows NCBI's int (undefined) is not supported by LOSAT's rpsbproc"
    );
    assert_eq!(
        refusal(codon(&frame(3, b""), i32::MIN, 5000)),
        "a site at residue -2147483647 of reading frame -1 whose position on the query overflows NCBI's int (undefined) is not supported by LOSAT's rpsbproc"
    );
}

// rpsbproc_main.cpp:1746-1747,1788-1789 (quoted at `motif_ends`): `GetLeft() + 1` in `int` and
// `GetRight() +/- 2 + 1` in `unsigned int` after `Pr2NaConvert`; a set outside `int`, a
// minus-strand residue at `INT_MIN` (`domains::pr2na_ends`) and `GetLeft() + 1` at `INT_MAX`
// overflow `int`: refused.
#[test]
fn motif_ends_and_their_overflow() {
    let set = |from: i64, to: i64| SegSet::from_segs(vec![Seg::new(from, to)]);
    // Frame +2 (388 % 3 = 1): 134 * 3 + 1 = 403, 138 * 3 + 1 + 2 = 417, then + 2 + 1.
    let plus = hit(388, 2583, 1);
    assert_eq!(value(motif_ends(&set(134, 138), &plus, 1, "2")), (404, 420));
    let refused = "a motif of reading frame 2 whose positions overflow NCBI's int (undefined) is not supported by LOSAT's rpsbproc";
    assert_eq!(
        refusal(motif_ends(&set(715_827_882, 715_827_882), &plus, 1, "2")),
        refused
    );
    assert_eq!(
        refusal(motif_ends(&set(5, i64::from(i32::MAX) + 1), &plus, 1, "2")),
        refused
    );
    // Frame -1 of a minus-strand hit (alignedLen (599 + 1) / 3 = 200, na length 601).
    let minus = hit(301, 599, 2);
    assert_eq!(value(motif_ends(&set(10, 12), &minus, 3, "-1")), (563, 569));
    assert_eq!(
        refusal(motif_ends(&set(i64::from(i32::MIN), 12), &minus, 3, "-1")),
        "a motif of reading frame -1 whose positions overflow NCBI's int (undefined) is not supported by LOSAT's rpsbproc"
    );
}

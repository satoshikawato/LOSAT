//! Tests of `seq_vector`: the reads of a row of a request's Bioseq (`RowVector`) and their
//! exceptions, against the rows of NCBI BLAST+ 2.17.0 (probes `probe3.py`, `probe4.py` of agent E3,
//! 2026-11-27; their rows are frozen in `tests/blast_formatter`). One case for each mechanism: the
//! query's rows (`of6_lq_*`, `of6_dq*`) read as the subject's do.

use super::*;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::{read_document, table::named};

/// `n` letters of IUPACna.
fn letters(n: usize) -> String {
    "ACGT".repeat(n.div_ceil(4))[..n].to_string()
}

/// A raw nucleotide `inst` of `length` letters that holds `held` of them.
fn raw(length: u32, held: usize) -> String {
    format!(
        "repr raw, mol na, length {length}, seq-data iupacna \"{}\"",
        letters(held)
    )
}

/// A delta nucleotide `inst` of `length` letters: its literals' lengths and the letters each
/// holds (`None`: a gap).
fn delta(length: u32, literals: &[(u32, Option<usize>)]) -> String {
    let literals: Vec<String> = literals
        .iter()
        .map(|&(n, held)| match held {
            Some(held) => format!(
                "literal {{ length {n}, seq-data iupacna \"{}\" }}",
                letters(held)
            ),
            None => format!("literal {{ length {n} }}"),
        })
        .collect();
    format!(
        "repr delta, mol na, length {length}, ext delta {{ {} }}",
        literals.join(", ")
    )
}

/// The reads of the row over the letters `from..=to` of the plus strand of the Bioseq whose
/// `inst` is `inst`, on the plus strand (`plus`) or the minus strand: `ok` and the letters the
/// map gives there (of the plus strand), or the message of NCBI's exception.
fn row_read(inst: &str, plus: bool, from: u32, to: u32) -> (String, Vec<u8>) {
    let text = format!(
        "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"s\" }}, inst {{ {inst} }} }} }} }}"
    );
    let (value, posts) = read_document(text.as_bytes(), named::BIOSEQ_SET, Names::STRATEGY);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.unwrap_or_else(|error| panic!("{error:?}"));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let map = SeqMap::of(Bioseq::of(entry.choice().1)).unwrap_or_else(|e| panic!("{e:?}"));
    match RowVector::new(&map, plus).read(from, to) {
        Ok(()) => {
            let mut out = Vec::new();
            map.letters(from, to + 1, &mut out);
            ("ok".into(), out)
        }
        Err(RenderError::Exception(e)) => {
            (String::from_utf8_lossy(e.get_msg()).into_owned(), vec![])
        }
        Err(other) => panic!("{other:?}"),
    }
}

/// The message of the read of `row_read`.
fn outcome(inst: &str, plus: bool, from: u32, to: u32) -> String {
    row_read(inst, plus, from, to).0
}

fn out_of_range(end: u32) -> String {
    format!("reference out of range of Seq-inst data: {end}")
}

fn cannot_get(from: u32, to: u32) -> String {
    format!("CSeqVector_CI::GetSeqData: cannot get seq-data in range: {from}-{to}")
}

/// A plus-strand row of 2048 letters from the first letter of a Bioseq that declares 3000 and
/// holds fewer: each fill of the cache takes up to 1024 letters from where the read is, and the
/// read fills once more after its last chunk (clipped at the segment's end), so the error names
/// the end of the first fill that passes the letters held, not the row's end (probe3
/// `of6_lp_1000_3000`: 1024; `lp_2000_3000`, `lq_2000_3000`: 2048; `lp_2048_3000`,
/// `lp_2100_3000`, `lq_2048_3000`: 3000). A Bioseq that holds every letter reads
/// (`lp_3000_3000`), as does one whose length the read ends at, with no fill after it
/// (`lp_2048_2048`, `lq_2048_2048`). The fills count from the row's first letter: a row from
/// letter 100 throws at 1124 (`of6_slen_short_legacy_sseq`).
#[test]
fn a_plus_row_throws_at_the_end_of_the_fill_that_passes_the_data() {
    for (held, want) in [
        (1000, out_of_range(1024)),
        (2000, out_of_range(2048)),
        (2048, out_of_range(3000)),
        (2100, out_of_range(3000)),
        (3000, "ok".into()),
    ] {
        assert_eq!(
            outcome(&raw(3000, held), true, 0, 2047),
            want,
            "held {held}"
        );
    }
    let (ok, letters_read) = row_read(&raw(2048, 2048), true, 0, 2047);
    assert_eq!(ok, "ok");
    assert_eq!(letters_read, letters(2048).into_bytes());
    assert_eq!(outcome(&raw(1500, 400), true, 100, 219), out_of_range(1124));
}

/// The minus-strand row of the same alignment (probe3 `of6_lm_*`): the vector reads the map back
/// from its end, so its fills of the vector's positions 952-1976 and 1976-3000 read the Seq-data
/// at 1024-2048 and 0-1024: a Bioseq that holds 2048 letters reads (`lm_2048_3000`,
/// `lm_2100_3000`), and one that holds fewer throws at 2048 whatever it holds
/// (`lm_2000_3000`, `lm_1000_3000`).
#[test]
fn a_minus_row_reads_the_data_back_from_the_row_s_end() {
    for (held, want) in [
        (1000, out_of_range(2048)),
        (2000, out_of_range(2048)),
        (2048, "ok".into()),
        (2100, "ok".into()),
    ] {
        assert_eq!(
            outcome(&raw(3000, held), false, 0, 2047),
            want,
            "held {held}"
        );
    }
    assert_eq!(outcome(&raw(2048, 2048), false, 0, 2047), "ok");
}

/// Rows of delta Bioseqs (probe4 `of6_d*`, `of6_dq*`): a literal that holds fewer letters than it
/// declares throws at the end of its fill, counted from the literal's start (`d2`: the second
/// literal of 250 holds 100, the row from letter 100 throws at 250); segments that end before
/// the length: the range check of the read (`CanGetRange`, `d5` plus: `0-2048`; `d4` minus:
/// `952-3000`), or where no segment holds the vector's first position of the read, the
/// iterator's `cannot locate segment at` (`d5m`: 952; `dq5` plus, after the last segment: 120);
/// segments past the length are cut at it (`d10`); a gap reads as `N` (`d11m`).
#[test]
fn a_delta_row_throws_at_its_literal_or_where_no_segment_holds_the_read() {
    assert_eq!(
        outcome(
            &delta(400, &[(150, Some(150)), (250, Some(100))]),
            true,
            100,
            219
        ),
        out_of_range(250)
    );
    let short = delta(3000, &[(1000, Some(1000)), (1000, Some(1000))]);
    assert_eq!(outcome(&short, true, 0, 2047), cannot_get(0, 2048));
    assert_eq!(
        outcome(&short, false, 0, 2047),
        "CSeqVector_CI: cannot locate segment at 952"
    );
    let longer = delta(3000, &[(1000, Some(1000)), (1500, Some(1500))]);
    assert_eq!(outcome(&longer, false, 0, 2047), cannot_get(952, 3000));
    assert_eq!(outcome(&longer, true, 0, 2047), "ok");
    assert_eq!(
        outcome(&delta(200, &[(60, Some(60)), (60, Some(60))]), true, 0, 119),
        "CSeqVector_CI: cannot locate segment at 120"
    );
    assert_eq!(
        outcome(
            &delta(2000, &[(1000, Some(1000)), (2000, Some(2000))]),
            true,
            0,
            2047
        ),
        "ok"
    );
    let (ok, letters_read) = row_read(
        &delta(3000, &[(1000, Some(1000)), (500, None), (1500, Some(1500))]),
        false,
        0,
        2047,
    );
    assert_eq!(ok, "ok");
    assert_eq!(
        letters_read,
        [
            letters(1000).into_bytes(),
            vec![b'N'; 500],
            letters(548).into_bytes()
        ]
        .concat()
    );
}

/// `kInvalidSeqPos` (4294967295): a Bioseq that declares it has no vector size, `Invalid
/// sequence length` at the row's first read (`of6_qlen_huge`); a delta literal of that length
/// after a first segment: the plus read's range check moves to it (probe4 `d6`: `cannot get
/// seq-data in range: 0-2048`), the minus read's first position is found past it (`d6m`:
/// `Invalid sequence length`).
#[test]
fn a_length_of_k_invalid_seq_pos_throws_where_the_read_resolves_it() {
    assert_eq!(
        outcome(&raw(u32::MAX, 8), true, 0, 7),
        "Invalid sequence length"
    );
    let invalid = delta(3000, &[(1000, Some(1000)), (u32::MAX, Some(2000))]);
    assert_eq!(outcome(&invalid, true, 0, 2047), cannot_get(0, 2048));
    assert_eq!(outcome(&invalid, false, 0, 2047), "Invalid sequence length");
}

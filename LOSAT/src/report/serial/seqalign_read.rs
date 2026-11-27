//! The Seq-aligns of a BLAST archive as NCBI's tabular formats read them (`-outfmt 6`, 7, 10 and
//! 20 of blast_formatter, `report::tabular_seqalign`): the split of an ungapped search's
//! Seq-aligns (`PrepareBlastUngappedSeqalign`), the scores (`GetAlnScores`, `GetNamedScore`), each
//! row's Seq-id and range, and the Dense-seg of two rows that the formatter builds
//! (`CreateDensegFromDendiag`, or the Dense-seg as it is). The mirror of `seqalign`, which writes
//! them.
//!
//! What NCBI reads of a form that BLAST's archives do not hold (a Dense-seg of other than two
//! rows, vectors of other sizes than its segments', a row without letters) reads memory the
//! archive does not decide or is not modelled: it is LOSAT's refusal by name
//! ([`ReadStop::Unsupported`]). A Std-seg (the translated programs) is refused the same way by
//! the formatter before any row (`blast_formatter::refusals`).

use crate::report::archive::read::{invalid_row, seq_align_id, Stop};
use crate::report::serial::asn_types::view::V;
use crate::report::serial::seqalign::Strand;
use crate::seqid::ParsedSeqId;

/// Where reading a Seq-align stops.
#[derive(Debug)]
pub(crate) enum ReadStop {
    /// Where NCBI's code stops (`report::archive::read::Stop`).
    Stop(Stop),
    /// A form that LOSAT does not read as NCBI: what it is.
    Unsupported(&'static str),
}

impl From<Stop> for ReadStop {
    fn from(stop: Stop) -> Self {
        ReadStop::Stop(stop)
    }
}

/// One Seq-align of a result set as the formatter reads it: the archive's, or one Dense-diag of
/// an ungapped search's Seq-align that `PrepareBlastUngappedSeqalign` made a Seq-align of its own
/// (`piece`), whose scores are the Dense-diag's (the Dense-diag's `GetScores` set as the Seq-align's scores).
#[derive(Clone, Copy, Debug)]
pub(crate) struct ReadAlign<'a> {
    /// The archive's Seq-align (for a piece, the one it was cut from).
    pub(crate) align: V<'a>,
    piece: Option<V<'a>>,
}

/// What `GetAlnScores` reads (align_format_util.cpp:811-848): each -1 where the list has none
/// (`comp_adjustment_method` 0).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AlnScores {
    pub(crate) score: i64,
    pub(crate) bit_score: f64,
    pub(crate) evalue: f64,
    pub(crate) num_ident: i64,
    /// The "sum_n" score: the size of the HSP's linked set (outfmt 0's `Expect(N)`).
    pub(crate) sum_n: i64,
    /// The "comp_adjustment_method" score (outfmt 0's `Method:`).
    pub(crate) comp_adjustment_method: i64,
    /// Whether the list holds a "use_this_gi" score (the GIs the description lines show).
    pub(crate) use_this_gi: bool,
}

/// The Dense-seg of two rows that the tabular formats build of a Seq-align
/// (tabular.cpp:925-941): its segments' starts (two per segment, -1 for a gap) and lengths, and
/// the strands of its first segment (`GetSeqStrand`: plus where the Dense-seg has none;
/// `CDense_seg::Reverse`, Dense_seg.cpp:644-681, flips the strands vector in place without
/// reversing segment order, so the reversed Dense-seg's first-segment strands are the flipped
/// `strands`).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PairDenseg {
    pub(crate) starts: Vec<i32>,
    pub(crate) lens: Vec<u32>,
    pub(crate) strands: [Strand; 2],
    /// Whether the Dense-seg has strands and its first is minus: the tabular formats reverse it
    /// then (tabular.cpp:943-952).
    pub(crate) first_minus: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3160-3212
/// The Seq-aligns of a result set as the formatter reads them: where the search is ungapped
/// (`ungapped`), a Seq-align of more than one Dense-diag or Std-seg becomes one per segment
/// (`PrepareBlastUngappedSeqalign`); every other one is kept.
pub(crate) fn read_set<'a>(set: &[V<'a>], ungapped: bool) -> Vec<ReadAlign<'a>> {
    let mut out = Vec::with_capacity(set.len());
    for &align in set {
        let (variant, segs) = align.field("segs").choice();
        if ungapped && matches!(variant, "dendiag" | "std") && segs.len() > 1 {
            out.extend(segs.items().map(|piece| ReadAlign {
                align,
                piece: Some(piece),
            }));
        } else {
            out.push(ReadAlign { align, piece: None });
        }
    }
    out
}

/// The `Na-strand` of an ENUMERATED value: plus (1) and minus (2); every other value is read as
/// neither (`CAlnMap::IsPositiveStrand` takes it for plus, `CDense_seg::Reverse` keeps it).
fn strand(value: i64) -> Strand {
    match value {
        1 => Strand::Plus,
        2 => Strand::Minus,
        _ => Strand::Unknown,
    }
}

impl<'a> ReadAlign<'a> {
    /// The variant of the segments (`segs.Which()`): a piece's is `dendiag` or `std`.
    pub(crate) fn variant(&self) -> &'static str {
        self.align.field("segs").choice().0
    }

    /// The Dense-diags of a `dendiag` Seq-align: the piece's, else every one.
    fn diags(&self) -> Vec<V<'a>> {
        match self.piece {
            Some(piece) => vec![piece],
            None => self.align.field("segs").choice().1.items().collect(),
        }
    }

    /// `GetSeq_id` of `row`.
    pub(crate) fn seq_id(&self, row: usize) -> Result<&'a ParsedSeqId, Stop> {
        match self.piece {
            Some(piece) => match piece.get("ids").and_then(|ids| ids.item(row)) {
                Some(id) => Ok(id.seq_id()),
                None => Err(invalid_row()),
            },
            None => seq_align_id(self.align, row),
        }
    }

    /// The Seq-align's own score list (`GetScore()`, `None` where it has none): a piece's is its
    /// segment's.
    fn own_scores(&self) -> Option<V<'a>> {
        match self.piece {
            Some(piece) => piece.get("scores"),
            None => self.align.get("score"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:811-848
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:229-270
    /// `GetAlnScores`: the scores of the Seq-align's list, and where it has no `e_value` or
    /// `sum_e` those of its first segment (a Dense-diag's, a Std-seg's or the Dense-seg's) read
    /// on top of the values already found. A score without an id is NCBI's unset member (`GetId()`), a value of the
    /// other kind its `CInvalidChoiceSelection`; the `use_this_gi` scores are not read by the
    /// tabular formats, but a `gi:` score that is no number throws there and is refused.
    pub(crate) fn scores(&self) -> Result<AlnScores, ReadStop> {
        let mut scores = AlnScores {
            score: -1,
            bit_score: -1.0,
            evalue: -1.0,
            num_ident: -1,
            sum_n: -1,
            comp_adjustment_method: 0,
            use_this_gi: false,
        };
        if let Some(list) = self.own_scores() {
            if read_scores(list, &mut scores)? {
                return Ok(scores);
            }
        }
        let (variant, segs) = self.align.field("segs").choice();
        let first = match (self.piece, variant) {
            (Some(piece), _) => Some(piece),
            (None, "std" | "dendiag") => match segs.item(0) {
                Some(first) => Some(first),
                None => {
                    return Err(
                        Stop::Undefined("front() of a Seq-align's empty segment list").into(),
                    )
                }
            },
            (None, "denseg") => Some(segs),
            _ => None,
        };
        if let Some(list) = first.and_then(|first| first.get("scores")) {
            read_scores(list, &mut scores)?;
        }
        Ok(scores)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:520-592
    /// `GetNamedScore` of an integer: the first score of the Seq-align's own list with the string
    /// id `name`, a real value cut to an `int` (one past the `int` range is refused: the cast is
    /// undefined).
    pub(crate) fn named_int(&self, name: &str) -> Result<Option<i64>, ReadStop> {
        match self.named(name) {
            None => Ok(None),
            Some(value) => match value.choice() {
                ("int", n) => Ok(Some(n.int())),
                (_, real) => {
                    let real = real.real();
                    if real.is_nan() || real <= f64::from(i32::MIN) - 1.0 || real >= 2147483648.0 {
                        return Err(
                            Stop::Undefined("a real score past the int range cast to int").into(),
                        );
                    }
                    Ok(Some(real.trunc() as i64))
                }
            },
        }
    }

    /// `GetNamedScore` of a real number.
    pub(crate) fn named_real(&self, name: &str) -> Option<f64> {
        self.named(name).map(|value| match value.choice() {
            ("int", n) => n.int() as f64,
            (_, real) => real.real(),
        })
    }

    /// The value that `x_GetNamedScore` finds for `name`.
    fn named(&self, name: &str) -> Option<V<'a>> {
        self.own_scores()?
            .items()
            .find_map(|score| match score.get("id").map(|id| id.choice()) {
                Some(("str", text)) if text.str() == name.as_bytes() => Some(score.field("value")),
                _ => None,
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1150-1210
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:925-941
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:241-270
    /// The Dense-seg that the tabular formats read: a Dense-diag Seq-align's
    /// (`CreateDensegFromDendiag`: the segments of its Dense-diags in order), or the Dense-seg as
    /// it is. The Dense-seg must have two rows, two ids, each vector one entry per row and
    /// segment (strands also none), at least one segment and letters in each row: NCBI's
    /// `CAlnMap` reads its vectors without a check.
    pub(crate) fn dense_seg(&self) -> Result<PairDenseg, ReadStop> {
        let mut starts = Vec::new();
        let mut lens = Vec::new();
        let mut strands: Vec<i64> = Vec::new();
        match self.variant() {
            "denseg" => {
                let ds = self.align.field("segs").choice().1;
                check_rows(ds)?;
                let numseg = ds.field("numseg").int();
                for len in ds.field("lens").ints() {
                    lens.push(len);
                }
                if usize::try_from(numseg).ok() != Some(lens.len()) {
                    return Err(ReadStop::Unsupported(
                        "a Dense-seg whose number of segments is not that of its lengths",
                    ));
                }
                starts.extend(ds.field("starts").ints());
                if let Some(list) = ds.get("strands") {
                    strands.extend(list.ints());
                }
            }
            "dendiag" => {
                for diag in self.diags() {
                    check_rows(diag)?;
                    let diag_starts: Vec<i64> = diag.field("starts").ints().collect();
                    if diag_starts.len() != 2 {
                        return Err(ReadStop::Unsupported(
                            "a Dense-diag whose starts are not one per row",
                        ));
                    }
                    starts.extend(diag_starts);
                    lens.push(diag.field("len").int());
                    if let Some(list) = diag.get("strands") {
                        let diag_strands: Vec<i64> = list.ints().collect();
                        if diag_strands.len() != 2 {
                            return Err(ReadStop::Unsupported(
                                "a Dense-diag whose strands are not one per row",
                            ));
                        }
                        strands.extend(diag_strands);
                    } else if !strands.is_empty() {
                        return Err(ReadStop::Unsupported(
                            "Dense-diags of which only some have strands",
                        ));
                    }
                }
                if !strands.is_empty() && strands.len() != starts.len() {
                    return Err(ReadStop::Unsupported(
                        "Dense-diags of which only some have strands",
                    ));
                }
            }
            _ => {
                return Err(ReadStop::Unsupported(
                    "segments other than a Dense-seg or Dense-diags",
                ))
            }
        }
        let numseg = lens.len();
        if numseg == 0 {
            return Err(Stop::Undefined("the last segment of an alignment without one").into());
        }
        if starts.len() != 2 * numseg || !(strands.is_empty() || strands.len() == 2 * numseg) {
            return Err(ReadStop::Unsupported(
                "a Dense-seg whose starts or strands are not one per row and segment",
            ));
        }
        let lens: Vec<u32> = lens
            .into_iter()
            .map(|len| u32::try_from(len).ok().filter(|&len| len < i32::MAX as u32))
            .collect::<Option<_>>()
            .ok_or(ReadStop::Unsupported(
                "a segment length outside 0 to 2^31 - 2",
            ))?;
        let starts: Vec<i32> = starts
            .into_iter()
            .map(|start| i32::try_from(start).ok().filter(|&start| start >= -1))
            .collect::<Option<_>>()
            .ok_or(ReadStop::Unsupported(
                "a segment start outside -1 to 2^31 - 1",
            ))?;
        for row in 0..2 {
            if !(0..numseg).any(|seg| starts[seg * 2 + row] >= 0) {
                return Err(ReadStop::Unsupported("an alignment row without letters"));
            }
        }
        let row_strand = |at: usize| strands.get(at).map_or(Strand::Plus, |&s| strand(s));
        Ok(PairDenseg {
            first_minus: strands.first() == Some(&2),
            strands: [row_strand(0), row_strand(1)],
            starts,
            lens,
        })
    }
}

/// A Dense-seg's or Dense-diag's rows: two (`dim`, 2 by default) and two ids.
fn check_rows(value: V<'_>) -> Result<(), ReadStop> {
    if value.get("dim").is_some_and(|dim| dim.int() != 2) || value.field("ids").len() != 2 {
        return Err(ReadStop::Unsupported("an alignment of other than two rows"));
    }
    Ok(())
}

/// `s_GetBlastScore` of `list` into `scores`: whether it has an `e_value` or `sum_e`.
fn read_scores(list: V<'_>, scores: &mut AlnScores) -> Result<bool, ReadStop> {
    let mut has_score = false;
    for score in list.items() {
        let Some(id) = score.get("id") else {
            return Err(Stop::Unassigned {
                type_name: "Score",
                member: "id",
            }
            .into());
        };
        let ("str", name) = id.choice() else {
            continue;
        };
        let value = score.field("value");
        let int = || match value.choice() {
            ("int", n) => Ok(n.int()),
            (current, _) => Err(ReadStop::Stop(Stop::InvalidSelection {
                module: "NCBI-Seqalign",
                type_name: "Score.value",
                required: "int",
                current,
            })),
        };
        let real = || match value.choice() {
            ("real", x) => Ok(x.real()),
            (current, _) => Err(ReadStop::Stop(Stop::InvalidSelection {
                module: "NCBI-Seqalign",
                type_name: "Score.value",
                required: "real",
                current,
            })),
        };
        match name.str() {
            b"score" => scores.score = int()?,
            b"bit_score" => scores.bit_score = real()?,
            b"e_value" | b"sum_e" => {
                scores.evalue = real()?;
                has_score = true;
            }
            b"use_this_gi" => {
                int()?;
                scores.use_this_gi = true;
            }
            b"sum_n" => scores.sum_n = int()?,
            b"comp_adjustment_method" => scores.comp_adjustment_method = int()?,
            b"num_ident" => scores.num_ident = int()?,
            text if text.starts_with(b"gi:") => {
                if gi_of_name(text).is_none() {
                    return Err(ReadStop::Unsupported(
                        "a gi: score whose name NCBI BLAST+ cannot convert to a GI, where it \
                         throws",
                    ));
                }
            }
            _ => {}
        }
    }
    Ok(has_score)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:262-266
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3370-3384
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:868-940
/// The GI of a `gi:` score's name as `s_GetBlastScore` converts it: every `gi:` removed
/// (`NStr::Replace`, which looks again at the place of each removal), then `StringToNumeric<TGi>`
/// (an `Int8` in base 10 without flags: one optional sign, at least one digit and nothing else,
/// within the `Int8` range). `None` where the conversion throws.
pub(crate) fn gi_of_name(name: &[u8]) -> Option<i64> {
    let mut text = name.to_vec();
    let mut at = 0;
    while let Some(found) = text[at..].windows(3).position(|window| window == b"gi:") {
        at += found;
        text.drain(at..at + 3);
    }
    let (negative, digits) = match text.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, &text[..]),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let mut value: i64 = 0;
    for &digit in digits {
        let digit = i64::from(digit - b'0');
        value = value.checked_mul(10)?;
        value = if negative {
            value.checked_sub(digit)?
        } else {
            value.checked_add(digit)?
        };
    }
    Some(value)
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:147-170
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:165-240
/// `GetSeqRange` of `row` of a Seq-align whose Dense-seg is `ds` (`dense_seg`), from and to as
/// `TSeqPos`: on a minus row (the first segment's strand) the start of its last segment with
/// letters and the end of its first, else the start of its first and the end of its last.
/// `ds` comes from a Dense-seg Seq-align; a Dense-diag Seq-align combines its diagonals' ranges
/// (`ReadAlign::seq_range`).
fn denseg_range(ds: &PairDenseg, row: usize) -> (u32, u32) {
    let numseg = ds.lens.len();
    let with_letters = |seg: &usize| ds.starts[seg * 2 + row] >= 0;
    let first = (0..numseg).find(with_letters).expect("a row with letters");
    let last = (0..numseg)
        .rev()
        .find(with_letters)
        .expect("a row with letters");
    let (start_seg, stop_seg) = if ds.strands[row] == Strand::Minus {
        (last, first)
    } else {
        (first, last)
    };
    let start = ds.starts[start_seg * 2 + row] as u32;
    let stop = (ds.starts[stop_seg * 2 + row] as u32)
        .wrapping_add(ds.lens[stop_seg])
        .wrapping_sub(1);
    (start, stop)
}

impl ReadAlign<'_> {
    /// `GetSeqRange` of `row` of the Seq-align whose Dense-seg is `ds` (`dense_seg` of this one):
    /// a Dense-seg's (`denseg_range`), or the combined ranges of its Dense-diags (each from its
    /// start to its start plus its length less 1, `CRange::CombineWith`).
    pub(crate) fn seq_range(&self, ds: &PairDenseg, row: usize) -> (u32, u32) {
        if self.variant() == "denseg" {
            return denseg_range(ds, row);
        }
        let mut range: Option<(u32, u32)> = None;
        for (seg, &len) in ds.lens.iter().enumerate() {
            let from = ds.starts[seg * 2 + row] as u32;
            let to_open = from.wrapping_add(len);
            if to_open <= from {
                continue;
            }
            range = Some(match range {
                Some((low, high_open)) => (low.min(from), high_open.max(to_open)),
                None => (from, to_open),
            });
        }
        match range {
            Some((from, to_open)) => (from, to_open.wrapping_sub(1)),
            None => (u32::MAX, u32::MAX.wrapping_sub(1)),
        }
    }
}

#[cfg(test)]
#[path = "seqalign_read_tests.rs"]
mod tests;

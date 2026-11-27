//! The Seq-aligns of a BLAST archive as NCBI's tabular formats read them (`-outfmt 6`, 7, 10 and
//! 20 of blast_formatter, `report::tabular_seqalign`): the split of an ungapped search's
//! Seq-aligns (`PrepareBlastUngappedSeqalign`), the scores (`GetAlnScores`, `GetNamedScore`), each
//! row's Seq-id and range (`GetSeqRange`, a Std-seg's its own), and the Dense-seg of two rows
//! that the formatter builds (`CreateDensegFromDendiag`, `CreateDensegFromStdseg` and for two
//! nucleotide rows `CreateTranslatedDensegFromNADenseg`, or the Dense-seg as it is). The mirror of
//! `seqalign`, which writes them; the one reader of the archive path (line RG item F3).
//!
//! What NCBI reads of a form that BLAST's archives do not hold (a Dense-seg of other than two
//! rows, vectors of other sizes than its segments', a row without letters, a Std-seg location
//! that is neither an interval nor empty, a Std-seg of more rows than the first) reads memory the
//! archive does not decide or is not modelled: it is LOSAT's refusal by name
//! ([`ReadStop::Unsupported`]). NCBI's own exceptions of the conversions and ranges are
//! [`Stop::Exception`]s with NCBI's texts.

use crate::report::archive::read::{invalid_row, seq_align_id, seq_loc_id, Stop};
use crate::report::serial::asn_types::view::V;
use crate::report::serial::seqalign::{
    row_range_exception, seq_align_exception, RowRangeError, Strand,
};
use crate::seqid::defline_set::{compare, Sic};
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

/// What `GetAlnScores` reads (align_format_util.cpp:811-848): each -1 where the list has none.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AlnScores {
    pub(crate) score: i64,
    pub(crate) bit_score: f64,
    pub(crate) evalue: f64,
    pub(crate) num_ident: i64,
}

/// The Dense-seg of two rows that the tabular formats build of a Seq-align
/// (tabular.cpp:925-941): its segments' starts (two per segment, -1 for a gap) and lengths (in
/// the rows' residues where it has widths), the strands of its first segment (`GetSeqStrand`:
/// plus where the Dense-seg has none; `CDense_seg::Reverse`, Dense_seg.cpp:644-681, flips the
/// strands vector in place without reversing segment order, so the reversed Dense-seg's
/// first-segment strands are the flipped `strands`) and the rows' widths.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PairDenseg {
    pub(crate) starts: Vec<i32>,
    pub(crate) lens: Vec<u32>,
    pub(crate) strands: [Strand; 2],
    /// The rows' widths (`CDense_seg::GetWidths`): 3 on a translated row of a Std-seg's Dense-seg
    /// (both rows of TBLASTX's), `[1, 1]` where the Dense-seg has none.
    pub(crate) widths: [u32; 2],
    /// Whether it was made of Std-segs (`kTranslated`, tabular.cpp:918): its rows' coordinates
    /// and frames follow each row's strand, and it is never reversed.
    pub(crate) translated: bool,
    /// Whether the tabular formats reverse it (tabular.cpp:943-952): not translated, with
    /// strands, the first minus.
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

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:320-357
    /// `GetSeq_id` of `row`: a piece's (the Seq-align of one Dense-diag or Std-seg that
    /// `PrepareBlastUngappedSeqalign` makes) from its ids, a Std-seg's else from its location of
    /// the row (`seq_loc_id`, as `seq_align_id` reads a Seq-align of Std-segs).
    pub(crate) fn seq_id(&self, row: usize) -> Result<&'a ParsedSeqId, Stop> {
        let Some(piece) = self.piece else {
            return seq_align_id(self.align, row);
        };
        if let Some(id) = piece.get("ids").and_then(|ids| ids.item(row)) {
            return Ok(id.seq_id());
        }
        if self.variant() == "std" {
            if let Some(id) = piece.field("loc").item(row).and_then(seq_loc_id) {
                return Ok(id);
            }
        }
        Err(invalid_row())
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:913-941
    /// The Dense-seg that the tabular formats read of this Seq-align: a Std-seg Seq-align's
    /// (`std_dense_seg`; of two nucleotide rows, `both_nucleotide`, TBLASTX's, translated by
    /// `CreateTranslatedDensegFromNADenseg`), a Dense-diag Seq-align's or the Dense-seg
    /// (`pair_dense_seg`).
    pub(crate) fn dense_seg(&self, both_nucleotide: bool) -> Result<PairDenseg, ReadStop> {
        match self.variant() {
            "std" => self.std_dense_seg(both_nucleotide),
            _ => self.pair_dense_seg(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1150-1210
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:241-270
    /// The Dense-seg of a Dense-diag Seq-align (`CreateDensegFromDendiag`: the segments of its
    /// Dense-diags in order), or the Dense-seg as it is. The Dense-seg must have two rows, two
    /// ids, each vector one entry per row and segment (strands also none), at least one segment
    /// and letters in each row: NCBI's `CAlnMap` reads its vectors without a check.
    fn pair_dense_seg(&self) -> Result<PairDenseg, ReadStop> {
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
                    "segments other than a Dense-seg, Dense-diags or Std-segs",
                ))
            }
        }
        if !(strands.is_empty() || strands.len() == starts.len()) {
            return Err(ReadStop::Unsupported(
                "a Dense-seg whose starts or strands are not one per row and segment",
            ));
        }
        let first = |at: usize| strands.get(at).map_or(Strand::Plus, |&s| strand(s));
        let first_minus = strands.first() == Some(&2);
        pair_of(
            starts,
            lens,
            [first(0), first(1)],
            [1, 1],
            false,
            first_minus,
        )
    }

    /// The Std-segs of a Std-seg Seq-align: the piece's, else every one.
    fn std_segs(&self) -> Vec<V<'a>> {
        match self.piece {
            Some(piece) => vec![piece],
            None => self.align.field("segs").choice().1.items().collect(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:726-948
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:953-1000
    /// `CreateDensegFromStdseg` of the Std-segs (without a `SSeqIdChooser`), then where both
    /// rows are nucleotide (`both_nucleotide`) `CreateTranslatedDensegFromNADenseg`. Each Std-seg
    /// is a segment, each location a row (an interval its start and length, an empty location a
    /// gap); the first Std-seg's rows decide the dimension, the ids (each later location's id
    /// must equal its row's, `SerialEquals`) and each row's strand (the first interval's that has
    /// one, each later one equal); a segment's length is its first row's with letters, divided by
    /// 3 where another row's is a third of it (that row of width 1, the first of width 3) or
    /// where its first row with letters has the width 3 of another segment; widths that a
    /// segment does not decide are none where no segment decides any (the rows of equal
    /// lengths, TBLASTX's). NCBI's exceptions (`CSeqalignException::eInvalidInputAlignment`) are
    /// its own; a location of another kind, which NCBI reads through a null Seq-id, and a
    /// Std-seg of more rows than the first, read past NCBI's vectors, are refused.
    fn std_dense_seg(&self, both_nucleotide: bool) -> Result<PairDenseg, ReadStop> {
        let thrown = |line: u32, message: String| -> ReadStop {
            Stop::Exception(seq_align_exception(
                line,
                "CSeqalignException::eInvalidInputAlignment",
                message,
            ))
            .into()
        };
        let lengths_not_accurate = |line| {
            thrown(
                line,
                "CreateDensegFromStdseg(): Std-seg segment lengths not accurate!".into(),
            )
        };
        let segs = self.std_segs();
        let mut lens: Vec<u32> = Vec::with_capacity(segs.len());
        let mut starts: Vec<i64> = Vec::new();
        let mut strands: Vec<i64> = Vec::new();
        let mut widths: Vec<u32> = Vec::new();
        let mut determined: Vec<bool> = Vec::with_capacity(segs.len());
        let mut ids: Vec<&ParsedSeqId> = Vec::new();
        let mut row_lens: Vec<u32> = Vec::new();
        let mut dim = 0usize;
        for (seg, ss) in segs.iter().enumerate() {
            lens.push(0);
            row_lens.clear();
            determined.push(false);
            let mut row = 0usize;
            for loc in ss.field("loc").items() {
                if seg == 0 {
                    widths.push(0);
                    strands.push(0);
                } else if row >= dim {
                    return Err(ReadStop::Unsupported(
                        "a Std-seg of more rows than the first, which NCBI BLAST+ reads past its \
                         vectors",
                    ));
                }
                let (variant, value) = loc.choice();
                let seq_id = match variant {
                    "int" => {
                        let from = seq_pos(value.field("from"))?;
                        let to = seq_pos(value.field("to"))?;
                        if from > i32::MAX as u32 {
                            return Err(ReadStop::Unsupported(
                                "a Std-seg interval that starts past 2^31 - 1",
                            ));
                        }
                        starts.push(i64::from(from));
                        let row_len = to.wrapping_sub(from).wrapping_add(1);
                        row_lens.push(row_len);
                        let seg_len = &mut lens[seg];
                        let mut width = 0;
                        if *seg_len == 0 {
                            *seg_len = row_len;
                        } else if row_len.wrapping_mul(3) == *seg_len {
                            *seg_len /= 3;
                            width = 1;
                        } else if row_len / 3 == *seg_len {
                            width = 3;
                        } else if row_len != *seg_len {
                            return Err(lengths_not_accurate(801));
                        }
                        if width > 0 {
                            determined[seg] = true;
                            if widths[row] > 0 && widths[row] != width {
                                return Err(lengths_not_accurate(809));
                            }
                            widths[row] = width;
                        }
                        if let Some(strand) = value.get("strand") {
                            let strand = strand.int();
                            if seg == 0 || strands[row] == 0 {
                                strands[row] = strand;
                            } else if strands[row] != strand {
                                return Err(thrown(
                                    827,
                                    "CreateDensegFromStdseg(): Inconsistent strands!".into(),
                                ));
                            }
                        }
                        value.field("id").seq_id()
                    }
                    "empty" => {
                        starts.push(-1);
                        if seg == 0 {
                            strands[row] = 0;
                        }
                        row_lens.push(0);
                        value.seq_id()
                    }
                    _ => {
                        return Err(ReadStop::Unsupported(
                            "a Std-seg location other than an interval or an empty location, \
                             whose Seq-id NCBI BLAST+ reads through a null pointer",
                        ))
                    }
                };
                if seg == 0 {
                    ids.push(seq_id);
                } else if ids[row] != seq_id {
                    let mut text = b"CreateDensegFromStdseg(): Seq-ids: ".to_vec();
                    crate::seqid::label::fasta(ids[row], &mut text);
                    text.extend_from_slice(b" and ");
                    crate::seqid::label::fasta(seq_id, &mut text);
                    text.extend_from_slice(
                        b" are not identical! Without the OM it cannot be determined if they \
                          belong to the same sequence. Define and pass ChooseSeqId to resolve \
                          seq-ids.",
                    );
                    return Err(Stop::Exception(seq_align_exception(
                        861,
                        "CSeqalignException::eInvalidInputAlignment",
                        text,
                    ))
                    .into());
                }
                row += 1;
                if seg == 0 {
                    dim += 1;
                }
            }
            let ss_dim = ss.get("dim").map_or(2, |dim| dim.int());
            if i64::try_from(dim).ok() != Some(ss_dim) || row != dim {
                return Err(thrown(
                    873,
                    "CreateDensegFromStdseg(): Inconsistent dimentions!".into(),
                ));
            }
            if determined[seg] {
                let seg_len = lens[seg];
                for row in 0..dim {
                    let row_len = row_lens[row];
                    if row_len > 0 {
                        let width = if row_len == seg_len.wrapping_mul(3) {
                            3
                        } else if row_len == seg_len {
                            1
                        } else {
                            0
                        };
                        if widths[row] > 0 && widths[row] != width {
                            return Err(lengths_not_accurate(889));
                        }
                        widths[row] = width;
                    }
                }
            }
        }
        // The lengths of the segments that no width decided, from their first row with letters.
        let (mut failure, mut success) = (false, false);
        for seg in 0..lens.len() {
            if determined[seg] {
                success = true;
                continue;
            }
            if let Some(row) = (0..dim).find(|&row| starts[seg * dim + row] >= 0) {
                match widths[row] {
                    3 => lens[seg] /= 3,
                    0 => failure = true,
                    _ => {}
                }
            }
        }
        let has_widths = match (failure, success) {
            (true, true) => {
                return Err(thrown(
                    929,
                    "CreateDensegFromStdseg(): Some widths cannot be determined!".into(),
                ))
            }
            (true, false) => false,
            _ => true,
        };
        if dim != 2 {
            return Err(ReadStop::Unsupported("an alignment of other than two rows"));
        }
        let mut pair_widths = [1, 1];
        if both_nucleotide {
            if has_widths {
                return Err(thrown(
                    965,
                    "CSeq_align::CreateTranslatedDensegFromNADenseg(): Widths already exist for \
                     the original alignment"
                        .into(),
                ));
            }
            for (seg, len) in lens.iter_mut().enumerate() {
                if *len % 3 != 0 {
                    return Err(thrown(
                        985,
                        format!(
                            "CSeq_align::CreateTranslatedDensegFromNADenseg(): Length of segment \
                             {seg} is not divisible by 3."
                        ),
                    ));
                }
                *len /= 3;
            }
            pair_widths = [3, 3];
        } else if has_widths {
            if widths.contains(&0) {
                return Err(ReadStop::Unsupported(
                    "a Std-seg Seq-align whose rows have no letters in the segments that decide \
                     the widths",
                ));
            }
            pair_widths = [widths[0], widths[1]];
        }
        pair_of(
            starts,
            lens.into_iter().map(i64::from).collect(),
            [strand(strands[0]), strand(strands[1])],
            pair_widths,
            true,
            false,
        )
    }
}

/// A `TSeqPos` of a Seq-loc (an interval's `from` and `to`, a point's `point`): an INTEGER of
/// the `Uint4` range (NCBI's reader takes no other).
fn seq_pos(value: V<'_>) -> Result<u32, ReadStop> {
    u32::try_from(value.int())
        .map_err(|_| ReadStop::Unsupported("a Seq-loc position outside the Uint4 range"))
}

/// The pair Dense-seg of `starts` (two per segment) and `lens`, the first segment's `strands`,
/// the rows' `widths`, whether it was `translated` and whether the tabular formats reverse it
/// (`first_minus`): at least one segment, the starts from -1 to 2^31 - 1, the lengths to 2^31 - 2
/// and letters in each row (`CAlnMap`'s rows of gaps only throw in NCBI's report: not modelled).
fn pair_of(
    starts: Vec<i64>,
    lens: Vec<i64>,
    strands: [Strand; 2],
    widths: [u32; 2],
    translated: bool,
    first_minus: bool,
) -> Result<PairDenseg, ReadStop> {
    let numseg = lens.len();
    if numseg == 0 {
        return Err(Stop::Undefined("the last segment of an alignment without one").into());
    }
    if starts.len() != 2 * numseg {
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
    Ok(PairDenseg {
        starts,
        lens,
        strands,
        widths,
        translated,
        first_minus,
    })
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
            b"use_this_gi" | b"sum_n" | b"comp_adjustment_method" => {
                int()?;
            }
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
    /// `GetSeqRange` of `row` (`TSeqPos` from and to): a Std-seg Seq-align's (`std_range`), else
    /// that of its Dense-seg (`seq_range`).
    pub(crate) fn row_range(&self, row: usize) -> Result<(u32, u32), ReadStop> {
        if self.variant() == "std" {
            return self.std_range(row);
        }
        let ds = self.dense_seg(false)?;
        Ok(self.seq_range(&ds, row))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-230
    /// `GetSeqRange` of `row` of a Std-seg Seq-align (its Std-segs: the piece's, else every one):
    /// the row's intervals and points combined (`CRange::CombineWith`, which skips an empty
    /// range). NCBI's exceptions: a Seq-id of the row other than the first Std-seg's (by
    /// `CSeq_id::Compare`; where the first Std-seg's location of the row is no interval or point,
    /// the row's Seq-id stays unset and the next interval's differs), a row past a Std-seg's
    /// locations, a Std-seg of other locations than the Seq-align's `dim`, no letters.
    pub(crate) fn std_range(&self, row: usize) -> Result<(u32, u32), ReadStop> {
        let exception = |line: u32, kind: &'static str, message: &'static str| -> ReadStop {
            Stop::Exception(seq_align_exception(line, kind, message)).into()
        };
        // The Seq-align's `dim` (a piece's Seq-align, made by `PrepareBlastUngappedSeqalign`, has
        // none: showalign.cpp:3169-3178).
        let dim = match self.piece {
            Some(_) => None,
            None => self.align.get("dim").map(|dim| dim.int()),
        };
        // `from` and `to_open` (`TSeqPos`) of the combined range, none while it is empty.
        let mut range: Option<(u32, u32)> = None;
        let mut seq_id: Option<&ParsedSeqId> = None;
        for (seg_i, seg) in self.std_segs().into_iter().enumerate() {
            let mut row_i = 0usize;
            for loc in seg.field("loc").items() {
                if row_i == row {
                    let (variant, value) = loc.choice();
                    let place = match variant {
                        "int" => Some((
                            value.field("id").seq_id(),
                            seq_pos(value.field("from"))?,
                            seq_pos(value.field("to"))?,
                            186,
                        )),
                        "pnt" => {
                            let point = seq_pos(value.field("point"))?;
                            Some((value.field("id").seq_id(), point, point, 201))
                        }
                        _ => None,
                    };
                    if let Some((id, from, to, line)) = place {
                        if seg_i == 0 {
                            seq_id = Some(id);
                        } else if seq_id.is_none_or(|kept| compare(kept, id) != Sic::Yes) {
                            return Err(exception(
                                line,
                                "CSeqalignException::eInvalidRowNumber",
                                "CSeq_align::GetSeqRange(): Row seqids not consistent. Cannot \
                                 determine range.",
                            ));
                        }
                        let to_open = to.wrapping_add(1);
                        if to_open > from {
                            range = Some(match range {
                                Some((low, high)) => (low.min(from), high.max(to_open)),
                                None => (from, to_open),
                            });
                        }
                    }
                }
                row_i += 1;
            }
            if row >= row_i {
                return Err(exception(
                    215,
                    "CSeqalignException::eInvalidRowNumber",
                    "CSeq_align::GetSeqRange(): Invalid row number",
                ));
            }
            if dim.is_some_and(|dim| i64::try_from(row_i).ok() != Some(dim)) {
                return Err(exception(
                    220,
                    "CSeqalignException::eInvalidAlignment",
                    "CSeq_align::GetSeqRange(): loc.size is inconsistent with dim",
                ));
            }
        }
        match range {
            Some((from, to_open)) => Ok((from, to_open - 1)),
            None => Err(Stop::Exception(row_range_exception(RowRangeError::Empty)).into()),
        }
    }

    /// `GetSeqRange` of `row` of the Dense-seg or Dense-diag Seq-align whose Dense-seg is `ds`
    /// (`dense_seg` of this one): a Dense-seg's (`denseg_range`), or the combined ranges of its
    /// Dense-diags (each from its start to its start plus its length less 1,
    /// `CRange::CombineWith`).
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

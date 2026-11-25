//! NCBI's sequence ranges, `-query_loc` and `-subject_loc` (`ParseSequenceRange` and the
//! range checks of `CBlastFastaInputSource`): one range, given once, for every record that
//! the input of its role reads.
//!
//! The search sees only the letters of each record's interval; NCBI keeps the record
//! itself for the report, so LOSAT cuts each record to its interval and keeps, per record,
//! where the interval lies in the record (`RecordPlacement`): the reports shift the
//! coordinates by the interval's start (`RemapToQueryLoc`, `s_RemapToSubjectLoc`) and
//! print the record's length.

use anyhow::Result;

use crate::blastinput::app;
use crate::blastinput::fasta_reader::InputRecord;

/// The role whose input a range applies to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RangeRole {
    Query,
    Subject,
}

impl RangeRole {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1996-1998
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2541-2545
    fn error_prefix(self) -> &'static str {
        match self {
            Self::Query => "Invalid specification of query location",
            Self::Subject => "Invalid specification of subject location",
        }
    }
}

/// A range as given on the command line, 0-based and closed (`TSeqRange` after
/// one is taken from each of `from` and `to`); it is not limited to the length of any record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SequenceRange {
    pub from: usize,
    pub to: usize,
}

impl SequenceRange {
    /// The length of the range as given (`TSeqRange::GetLength`), which the tabular
    /// formats print as `qlen` even where the range ends past a record.
    pub fn length(&self) -> usize {
        self.to - self.from + 1
    }
}

/// NCBI's `ParseSequenceRange`: `start-stop`, 1-based and closed.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:145-179
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-227
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
/// `NStr::Split` keeps empty tokens (`1--5` has three, `-10` an empty first one).
/// `NStr::StringToInt` (flags 0, `ncbi_str::string_to_int`) converts the start, then the stop,
/// before the other checks; a part that it cannot convert (`a`, `1.5`, ` 5`) throws its
/// `CStringException`, which ends the program with `CATCH_ALL`'s report (`NativeError::exception`,
/// exit status 255) where the range is read, before anything is written.
pub fn parse_sequence_range(value: &str, role: RangeRole) -> Result<SequenceRange> {
    let prefix = role.error_prefix();
    let tokens: Vec<&str> = value.split('-').collect();
    if tokens.len() != 2 || tokens[0].is_empty() || tokens[1].is_empty() {
        return Err(app::engine_error(format!("{prefix} (Format: start-stop)")));
    }
    let mut numbers = [0i32; 2];
    for (number, token) in numbers.iter_mut().zip(&tokens) {
        *number = crate::ncbi_str::string_to_int(token.as_bytes(), 0)
            .map_err(|e| crate::cli::NativeError::exception(&e))?;
    }
    let [from, to] = numbers;
    if from <= 0 || to <= 0 {
        return Err(app::engine_error(format!(
            "{prefix} (range elements cannot be less than or equal to 0)"
        )));
    }
    if from == to {
        return Err(app::engine_error(format!(
            "{prefix} (range cannot be empty)"
        )));
    }
    if from > to {
        return Err(app::engine_error(format!(
            "{prefix} (start cannot be larger than stop)"
        )));
    }
    Ok(SequenceRange {
        from: from as usize - 1,
        to: to as usize - 1,
    })
}

/// Parses the range of a role when it was given (`args.Exist` is true for the argument and `args`
/// holds a value for it).
pub fn parse_optional_range(value: Option<&str>, role: RangeRole) -> Result<Option<SequenceRange>> {
    value
        .map(|value| parse_sequence_range(value, role))
        .transpose()
}

/// What a range makes of one record (`x_FastaToSeqLoc`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordInterval {
    /// The 0-based closed interval `[from, to]` of the record; `to` is the record's last
    /// letter where the range ends past it.
    Letters { from: usize, to: usize },
    /// The range starts just past the record's end: an interval without letters, which
    /// NCBI searches as a sequence without data.
    Empty,
    /// The range starts more than one letter past the record's end: NCBI's
    /// `Invalid from coordinate (greater than sequence length)`.
    PastEnd,
}

/// The interval of a record of `length` letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-460
/// The parser leaves `to > from`, so `Invalid sequence range` cannot occur. A start equal
/// to the record's length passes the check and gives the interval `[seqlen, seqlen-1]`.
pub fn record_interval(range: &SequenceRange, length: usize) -> RecordInterval {
    if range.from > length {
        RecordInterval::PastEnd
    } else if range.from == length {
        RecordInterval::Empty
    } else {
        RecordInterval::Letters {
            from: range.from,
            to: if range.to < length {
                range.to
            } else {
                length - 1
            },
        }
    }
}

/// Where the searched letters of a record lie in the record: the start of its interval
/// and the record's length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RecordPlacement {
    pub offset: usize,
    pub length: usize,
}

/// The placements of the records of one role; empty without a range (every record is
/// searched whole).
#[derive(Clone, Debug, Default)]
pub struct Placements(pub Vec<RecordPlacement>);

impl Placements {
    /// The start of record `index`'s interval in the record (0 without a range).
    pub fn offset(&self, index: usize) -> usize {
        self.0.get(index).map_or(0, |placement| placement.offset)
    }

    /// The length of record `index` itself, or `searched` (the length of the letters that
    /// were searched) without a range.
    pub fn length(&self, index: usize, searched: usize) -> usize {
        self.0
            .get(index)
            .map_or(searched, |placement| placement.length)
    }

    /// Whether a range applies (some record may then be searched in part).
    pub fn ranged(&self) -> bool {
        !self.0.is_empty()
    }

    /// Where record `index`'s letters lie in the record, `None` without a range.
    pub fn get(&self, index: usize) -> Option<RecordPlacement> {
        self.0.get(index).copied()
    }
}

/// The subjects as NCBI searches them with `-subject_loc`: every record cut to its interval
/// (an empty one for an interval without letters), and where each lies in its record.
/// A record whose interval starts more than one letter past its end stops NCBI, in file
/// order, while the subjects are read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// The range's `CInputException` is not caught here.
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
pub fn cut_subjects<R: InputRecord>(
    records: &[R],
    range: Option<&SequenceRange>,
) -> Result<Option<(Vec<R>, Placements)>> {
    let Some(range) = range else {
        return Ok(None);
    };
    let mut cut_records = Vec::with_capacity(records.len());
    let mut placements = Vec::with_capacity(records.len());
    for record in records {
        let (from, to_exclusive) = record_cut(range, record.seq().len())?;
        cut_records.push(record.cut(from, to_exclusive));
        placements.push(RecordPlacement {
            offset: from,
            length: record.seq().len(),
        });
    }
    Ok(Some((cut_records, Placements(placements))))
}

/// The letters `[from, to_exclusive)` of a record of `length` letters that `cut_subjects` keeps
/// (`from == to_exclusive == length` for an interval without letters), or the error of an interval
/// that starts past the end.
pub fn record_cut(range: &SequenceRange, length: usize) -> Result<(usize, usize)> {
    match record_interval(range, length) {
        RecordInterval::Letters { from, to } => Ok((from, to + 1)),
        RecordInterval::Empty => Ok((length, length)),
        RecordInterval::PastEnd => Err(app::options_error(
            "Invalid from coordinate (greater than sequence length)",
        )),
    }
}

/// How many subject records NCBI's reader reads (and writes the title warnings of) before
/// the range stops it: every record, or the records up to the first one whose interval
/// starts more than one letter past its end (`GetAllSeqs` reads a record whole, then
/// checks its range).
pub fn subjects_read<R: InputRecord>(records: &[R], range: Option<&SequenceRange>) -> usize {
    let Some(range) = range else {
        return records.len();
    };
    records
        .iter()
        .position(|record| record_interval(range, record.seq().len()) == RecordInterval::PastEnd)
        .map_or(records.len(), |index| index + 1)
}

/// The query input as NCBI reads it, one record at a time: where each searched query lies
/// in its record, its position in the input, and the records that are not searched.
#[derive(Clone, Debug, Default)]
pub struct QueryInput {
    /// Where each searched query lies in its record (empty without `-query_loc`).
    pub placements: Placements,
    /// The 0-based position in the input of each searched query: NCBI's reader numbers
    /// the records it reads (`Query_<n>`), skipped ones too.
    pub ordinals: Vec<usize>,
    /// The length of each record of the input, and whether it is skipped (its interval
    /// starts more than one letter past its end).
    pub input_lengths: Vec<usize>,
    pub skipped: Vec<bool>,
    /// The length of the range as given (the tabular `qlen`), with `-query_loc`.
    pub range_length: Option<usize>,
}

impl QueryInput {
    /// Every record searched whole (no `-query_loc`).
    pub fn whole<R: InputRecord>(records: &[R]) -> Self {
        Self {
            placements: Placements::default(),
            ordinals: (0..records.len()).collect(),
            input_lengths: records.iter().map(|record| record.seq().len()).collect(),
            skipped: vec![false; records.len()],
            range_length: None,
        }
    }

    /// The input position of searched query `index` (`Query_<n>` is this plus one).
    pub fn ordinal(&self, index: usize) -> usize {
        self.ordinals.get(index).copied().unwrap_or(index)
    }

    /// The searched queries of the input records `input`.
    pub fn searched_in(&self, input: std::ops::Range<usize>) -> std::ops::Range<usize> {
        let first = self
            .ordinals
            .partition_point(|&ordinal| ordinal < input.start);
        let last = self
            .ordinals
            .partition_point(|&ordinal| ordinal < input.end);
        first..last
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-146
    // No NCBI counterpart for the slice itself: a batch is counted from its first record alone, so
    // the input records from the start of a batch have the batches of the whole input from there.
    /// The input records `input` as an input of their own: their lengths and skips, and the
    /// searched queries among them (`searched_in`) with their placements, numbered from
    /// `input.start`. EXPERIMENT (LOSAT_X_BPQBLOCK): a block of whole query batches of the search
    /// of all queries (`blast_engine.rs`: `x_prelim_block_plan`, and `run_preliminary`'s loop over
    /// its blocks, `blast_engine/prelim_run.rs`), whose `batches` are the whole input's batches of those records.
    pub fn slice(&self, input: std::ops::Range<usize>) -> Self {
        let searched = self.searched_in(input.clone());
        Self {
            placements: Placements(if self.placements.0.is_empty() {
                Vec::new()
            } else {
                self.placements.0[searched.clone()].to_vec()
            }),
            ordinals: self.ordinals[searched]
                .iter()
                .map(|&ordinal| ordinal - input.start)
                .collect(),
            input_lengths: self.input_lengths[input.clone()].to_vec(),
            skipped: self.skipped[input].to_vec(),
            range_length: self.range_length,
        }
    }

    /// The query batches of the input with a fixed batch size: NCBI reads batches while
    /// records remain (`End` of the input is false), so a batch can hold only skipped records (NCBI then
    /// stops with `Empty CBlastQueryVector`, and no batch follows). A batch size of 0
    /// reads no record: one empty batch. Each record counts its own letters: under
    /// `-parse_deflines` NCBI's reader counts a record otherwise only where it repeats an id of
    /// its scope (`query_batch::BatchCount`), and the run stops at the batch of the first such
    /// record (`query_batch::QueryScopes::add_batch`), so the batches before it are NCBI's and the
    /// program's loop (`query_batch::read_query_batches`) cuts that batch as NCBI does.
    pub fn batches(&self, batch_size: u32) -> Vec<RangedBatch> {
        let mut batches = Vec::new();
        let mut start = 0;
        while start < self.input_lengths.len() {
            let end = next_ranged_batch_end(&self.input_lengths, &self.skipped, start, batch_size);
            let searched = self.searched_in(start..end);
            let empty = searched.is_empty();
            batches.push(RangedBatch {
                input: start..end,
                searched,
            });
            if empty {
                break;
            }
            start = end;
        }
        batches
    }
}

/// The queries as NCBI searches them with `-query_loc`: the searched queries, each cut to
/// its interval, in input order, and the input they come from.
#[derive(Clone, Debug)]
pub struct RangedQueries<R = crate::blastinput::fasta_reader::FastaRecord> {
    pub records: Vec<R>,
    pub input: QueryInput,
}

/// Cuts every query record to its interval. NCBI's batch reader ignores the exception of a
/// record whose interval starts more than one letter past its end, so that record is not
/// searched and nothing is reported for it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:144-155
pub fn cut_queries<R: InputRecord>(records: &[R], range: &SequenceRange) -> RangedQueries<R> {
    let mut ranged = RangedQueries {
        records: Vec::with_capacity(records.len()),
        input: QueryInput {
            placements: Placements(Vec::with_capacity(records.len())),
            ordinals: Vec::with_capacity(records.len()),
            input_lengths: Vec::with_capacity(records.len()),
            skipped: Vec::with_capacity(records.len()),
            range_length: Some(range.length()),
        },
    };
    for (ordinal, record) in records.iter().enumerate() {
        let length = record.seq().len();
        ranged.input.input_lengths.push(length);
        let (from, to_exclusive) = match record_interval(range, length) {
            RecordInterval::Letters { from, to } => (from, to + 1),
            RecordInterval::Empty => (length, length),
            RecordInterval::PastEnd => {
                ranged.input.skipped.push(true);
                continue;
            }
        };
        ranged.input.skipped.push(false);
        ranged.records.push(record.cut(from, to_exclusive));
        ranged.input.placements.0.push(RecordPlacement {
            offset: from,
            length,
        });
        ranged.input.ordinals.push(ordinal);
    }
    ranged
}

/// One query batch of the input: the input records `input` (skipped ones included), and
/// the searched queries among them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RangedBatch {
    pub input: std::ops::Range<usize>,
    pub searched: std::ops::Range<usize>,
}

/// The end (an input record index) of the query batch that starts at input record `start`.
/// A skipped record adds nothing to the residues read; the others add the length of the
/// whole record, not of its interval.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
/// `size_read` and the batch size are `TSeqPos` (32-bit unsigned): an `Int4` batch size
/// converts to it (a negative one becomes 2^32 plus the value).
pub fn next_ranged_batch_end(
    input_lengths: &[usize],
    skipped: &[bool],
    start: usize,
    batch_size: u32,
) -> usize {
    next_ranged_batch(input_lengths, skipped, start, batch_size).0
}

/// `next_ranged_batch_end`, and whether the batch reached its size (`size_read` is no
/// longer below the batch size). A batch that did not reach it ends where the records end,
/// and NCBI's reader reads on while lines are left (`End()` is false), so what ends the
/// reading after the last record (an `eEOF`, a reader error) comes in that batch; after a
/// batch that reached its size it comes in the next one.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-146
pub fn next_ranged_batch(
    input_lengths: &[usize],
    skipped: &[bool],
    start: usize,
    batch_size: u32,
) -> (usize, bool) {
    next_counted_batch(input_lengths.len(), start, batch_size, &mut |record| {
        if skipped[record] {
            0
        } else {
            input_lengths[record]
        }
    })
}

/// `next_ranged_batch` of `records` input records, each counted by `count` (of its index, every
/// record read in turn): 0 for a skipped record, else NCBI's `sequence::GetLength` of the
/// record's id in the scope, which is the record's own length but for a repeated id under
/// `-parse_deflines` (`query_batch::BatchCount`).
pub fn next_counted_batch(
    records: usize,
    start: usize,
    batch_size: u32,
    count: &mut dyn FnMut(usize) -> usize,
) -> (usize, bool) {
    let mut end = start;
    let mut size_read: u32 = 0;
    while end < records && size_read < batch_size {
        size_read = size_read.wrapping_add(count(end) as u32);
        end += 1;
    }
    (end, size_read >= batch_size)
}

/// NCBI's frame of an aligned nucleotide row as its reports print it: from the record
/// coordinate of the row's start (the low end on the plus strand, the high end on the
/// minus strand) and the length of the whole record.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:411-422
/// `start` is the 0-based record coordinate; the tabular `qframe`/`sframe` use the same
/// rule (`CAlignFormatUtil::GetFrame`, align_format_util.cpp:1233-1245). Without a range
/// it equals the frame of the search.
pub fn record_frame(plus: bool, start: usize, record_length: usize) -> i32 {
    if plus {
        (start % 3) as i32 + 1
    } else {
        -(((record_length - start - 1) % 3) as i32 + 1)
    }
}

/// A mask of the searched letters `[from, to]` (0-based) as the reports show it: moved into the
/// record by the start of the interval (`offset`, 0 without a range) and cut at the interval's end
/// (`length` letters searched); `None` for an empty mask, one past the interval, or one that
/// covers the whole interval.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:825-843
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:936-938
pub fn record_mask(from: i32, to: i32, length: usize, offset: usize) -> Option<(i32, i32)> {
    let target = (offset as i32, (offset + length) as i32 - 1);
    if from > to || from > target.1 || from + target.0 > target.1 {
        return None;
    }
    let mapped = (target.0 + from, (target.0 + to).min(target.1));
    (mapped != target).then_some(mapped)
}

#[cfg(test)]
#[path = "seq_range_tests.rs"]
mod tests;

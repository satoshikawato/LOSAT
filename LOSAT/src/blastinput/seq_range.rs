//! NCBI's sequence ranges, `-query_loc` and `-subject_loc` (`ParseSequenceRange` and the
//! range checks of `CBlastFastaInputSource`): one range, given once, for every record that
//! the input of its role reads.
//!
//! The search sees only the letters of each record's interval; NCBI keeps the record
//! itself for the report, so LOSAT cuts each record to its interval and keeps, per record,
//! where the interval lies in the record (`RecordPlacement`): the reports shift the
//! coordinates by the interval's start (`RemapToQueryLoc`, `s_RemapToSubjectLoc`) and
//! print the record's length.

use anyhow::{bail, Result};

use crate::blastinput::app;
use crate::blastinput::fasta_reader::InputRecord;

/// The role whose input a range applies to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RangeRole {
    Query,
    Subject,
}

impl RangeRole {
    fn option(self) -> &'static str {
        match self {
            Self::Query => "-query_loc",
            Self::Subject => "-subject_loc",
        }
    }

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1996-1998
    // ```c++
    //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
    //                                      "Invalid specification of query location");
    // ```
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2541-2545
    // ```c++
    //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
    //             subj_range =
    //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
    //                             "Invalid specification of subject location");
    //         }
    // ```
    fn error_prefix(self) -> &'static str {
        match self {
            Self::Query => "Invalid specification of query location",
            Self::Subject => "Invalid specification of subject location",
        }
    }
}

/// A range as given on the command line, 0-based and closed (`TSeqRange` after
/// `from--, to--`); it is not limited to the length of any record.
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input_aux.cpp:145-179
/// ```c++
/// TSeqRange
/// ParseSequenceRange(const string& range_str,
///                    const char* error_prefix /* = NULL */)
/// {
///     static const char* kDfltErrorPrefix = "Failed to parse sequence range";
///     static const string kDelimiters("-");
///     string error_msg(error_prefix ? error_prefix : kDfltErrorPrefix);
///
///     vector<string> tokens;
///     NStr::Split(range_str, kDelimiters, tokens);
///     if (tokens.size() != 2 || tokens.front().empty() || tokens.back().empty()) {
///         error_msg += " (Format: start-stop)";
///         NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
///     }
///     int from = NStr::StringToInt(tokens.front());
///     int to = NStr::StringToInt(tokens.back());
///     if (from <= 0 || to <= 0) {
///         error_msg += " (range elements cannot be less than or equal to 0)";
///         NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
///     }
///     if (from == to) {
///         error_msg += " (range cannot be empty)";
///         NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
///     }
///     if (from > to) {
///         error_msg += " (start cannot be larger than stop)";
///         NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
///     }
///     from--, to--;   // decrement to make range 0-based
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:216-227
/// ```c++
///     catch (const blast::CBlastException& e) {                               \
///         ...
///         } else {                                                            \
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
/// `NStr::Split` keeps empty tokens (`1--5` has three, `-10` an empty first one).
/// `NStr::StringToInt` converts the start, then the stop, before the other checks; a part
/// that it cannot convert stops NCBI with a `CStringException`, whose text names the
/// source files of the build (exit 255): LOSAT rejects that value explicitly, as it does
/// the environment variables that NCBI converts the same way (plan TD-15).
pub fn parse_sequence_range(value: &str, role: RangeRole, program: &str) -> Result<SequenceRange> {
    let prefix = role.error_prefix();
    let tokens: Vec<&str> = value.split('-').collect();
    if tokens.len() != 2 || tokens[0].is_empty() || tokens[1].is_empty() {
        return Err(app::engine_error(&format!("{prefix} (Format: start-stop)")));
    }
    let mut numbers = [0i32; 2];
    for (number, token) in numbers.iter_mut().zip(&tokens) {
        match crate::blastinput::query_batch::ncbi_string_to_int(std::ffi::OsStr::new(token)) {
            Some(value) => *number = value,
            None => bail!(
                "the {} value {value:?} has a part {token:?} that NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's {program}",
                role.option()
            ),
        }
    }
    let [from, to] = numbers;
    if from <= 0 || to <= 0 {
        return Err(app::engine_error(&format!(
            "{prefix} (range elements cannot be less than or equal to 0)"
        )));
    }
    if from == to {
        return Err(app::engine_error(&format!(
            "{prefix} (range cannot be empty)"
        )));
    }
    if from > to {
        return Err(app::engine_error(&format!(
            "{prefix} (start cannot be larger than stop)"
        )));
    }
    Ok(SequenceRange {
        from: from as usize - 1,
        to: to as usize - 1,
    })
}

/// Parses the range of a role when it was given (`args.Exist(...) && args[...]`).
pub fn parse_optional_range(
    value: Option<&str>,
    role: RangeRole,
    program: &str,
) -> Result<Option<SequenceRange>> {
    value
        .map(|value| parse_sequence_range(value, role, program))
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-460
/// ```c++
///     // sanity checks for the range
///     const TSeqPos from = m_Config.GetRange().GetFrom() == kEmptyRange.GetFrom()
///         ? 0 : m_Config.GetRange().GetFrom();
///     const TSeqPos to = m_Config.GetRange().GetTo() == kEmptyRange.GetTo()
///         ? 0 : m_Config.GetRange().GetTo();
///
///     // Get the sequence length
///     const TSeqPos seqlen = seq_entry->GetSeq().GetInst().GetLength();
///     ...
///     if (to > 0 && to < from) {
///         NCBI_THROW(CInputException, eInvalidRange,
///                    "Invalid sequence range");
///     }
///     if (from > seqlen) {
///         NCBI_THROW(CInputException, eInvalidRange,
///                    "Invalid from coordinate (greater than sequence length)");
///     }
///     // N.B.: if the to coordinate is greater than or equal to the sequence
///     // length, we fix that silently
///
///
///     // set sequence range
///     retval->SetInt().SetFrom(from);
///     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
/// ```
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
}

/// The subjects as NCBI searches them with `-subject_loc`: every record cut to its interval
/// (an empty one for an interval without letters), and where each lies in its record.
/// A record whose interval starts more than one letter past its end stops NCBI, in file
/// order, while the subjects are read.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// ```c++
/// CRef<CBlastQueryVector>
/// CBlastInput::GetAllSeqs(CScope& scope)
/// {
///     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
///
///     while (!End()) {
///         try { retval->AddQuery(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             auto err = e.GetErrCode();
///             if (err == CObjReaderParseException::eEOF) {
///                 break;
///             } else if (err == CObjReaderParseException::eNoDefline) {
///                 ...
///             }
///             throw;
///         }
///     }
///
///     return retval;
/// ```
/// The range's `CInputException` is not caught here.
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:172-176
/// ```c++
///     catch (const blast::CInputException& e) {                               \
///         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
///         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
///         exit_code = BLAST_INPUT_ERROR;                                      \
/// ```
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
        let length = record.seq().len();
        let (from, to_exclusive) = match record_interval(range, length) {
            RecordInterval::Letters { from, to } => (from, to + 1),
            RecordInterval::Empty => (length, length),
            RecordInterval::PastEnd => {
                return Err(app::options_error(
                    "Invalid from coordinate (greater than sequence length)",
                ))
            }
        };
        cut_records.push(record.cut(from, to_exclusive));
        placements.push(RecordPlacement {
            offset: from,
            length,
        });
    }
    Ok(Some((cut_records, Placements(placements))))
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

/// The explicit rejection of an interval without letters (a range that starts just past
/// a record's end), which NCBI searches as a sequence without data ("Sequence contains no
/// data"); LOSAT does not search records without letters either.
/// `ordinals` are the input positions of the records (empty: their indexes).
pub fn check_no_empty_interval<R: InputRecord>(
    records: &[R],
    placements: &Placements,
    ordinals: &[usize],
    role: RangeRole,
    program: &str,
) -> Result<()> {
    for (index, (record, placement)) in records.iter().zip(&placements.0).enumerate() {
        let index = ordinals.get(index).copied().unwrap_or(index);
        if record.seq().is_empty() && placement.length > 0 {
            let (role_name, option) = match role {
                RangeRole::Query => ("query", "-query_loc"),
                RangeRole::Subject => ("subject", "-subject_loc"),
            };
            bail!(
                "{role_name} record {} ({}) has no letters in the {option} range (the range starts just past its end); NCBI BLAST+ searches it as a sequence without data, which is not supported by LOSAT's {program}",
                index + 1,
                first_word(&record.title_bytes())
            );
        }
    }
    Ok(())
}

/// The title of a record up to its first space, as text (the `bio` ID of a `bio` record),
/// for LOSAT's messages.
fn first_word(title: &[u8]) -> String {
    let end = title
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(title.len());
    String::from_utf8_lossy(&title[..end]).into_owned()
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

    /// The query batches of the input with a fixed batch size: NCBI reads batches while
    /// records remain (`!input.End()`), so a batch can hold only skipped records (NCBI then
    /// stops with `Empty CBlastQueryVector`, and no batch follows). A batch size of 0
    /// reads no record: one empty batch.
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:144-155
/// ```c++
///         CRef<CBlastSearchQuery> q;
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
///                 break;
///             }
///             throw;
///         }
///         catch (const exception&) {
///             continue; //SB-2307. ignore well formed, not found accession
///         }
/// ```
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:138-170
/// ```c++
///     TSeqPos size_read = 0;
///
///     while (size_read < GetBatchSize()) {
///
///         if (End())
///             break;
///         ...
///         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
///
///         if (loc->IsInt()) {
///             size_read += sequence::GetLength(loc->GetInt().GetId(),
///                                              q->GetScope());
///         ...
///         retval->AddQuery(q);
///     }
/// ```
/// `size_read` and the batch size are `TSeqPos` (32-bit unsigned), as in
/// `query_batch::next_query_batch_end`; without skipped records the two agree.
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:140-146
/// ```c++
///     while (size_read < GetBatchSize()) {
///
///         if (End())
///             break;
///
///         CRef<CBlastSearchQuery> q;
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
/// ```
pub fn next_ranged_batch(
    input_lengths: &[usize],
    skipped: &[bool],
    start: usize,
    batch_size: u32,
) -> (usize, bool) {
    let mut end = start;
    let mut size_read: u32 = 0;
    while end < input_lengths.len() && size_read < batch_size {
        if !skipped[end] {
            size_read = size_read.wrapping_add(input_lengths[end] as u32);
        }
        end += 1;
    }
    (end, size_read >= batch_size)
}

/// NCBI's frame of an aligned nucleotide row as its reports print it: from the record
/// coordinate of the row's start (the low end on the plus strand, the high end on the
/// minus strand) and the length of the whole record.
///
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:411-422
/// ```c++
/// static int s_GetFrame (int start, ENa_strand strand, const CSeq_id& id,
///                        CScope& sp)
/// {
///     int frame = 0;
///     if (strand == eNa_strand_plus) {
///         frame = (start % 3) + 1;
///     } else if (strand == eNa_strand_minus) {
///         frame = -(((int)sp.GetBioseqHandle(id).GetBioseqLength() - start - 1)
///                   % 3 + 1);
///
///     }
///     return frame;
/// }
/// ```
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

#[cfg(test)]
mod tests {
    use super::*;

    fn message(error: anyhow::Error) -> String {
        format!("{error:#}")
    }

    #[test]
    fn spellings_follow_ncbis_parser() {
        let parse = |text: &str| parse_sequence_range(text, RangeRole::Query, "BLASTN");
        assert_eq!(parse("10-20").unwrap(), SequenceRange { from: 9, to: 19 });
        assert_eq!(parse("+1-5").unwrap(), SequenceRange { from: 0, to: 4 });
        assert_eq!(parse("01-05").unwrap(), SequenceRange { from: 0, to: 4 });
        assert_eq!(parse("1-+5").unwrap(), SequenceRange { from: 0, to: 4 });
        assert_eq!(
            parse("1-2147483647").unwrap(),
            SequenceRange {
                from: 0,
                to: 2147483646
            }
        );
        for text in ["", "-", "--", "-10", "10-", "10", "1--5", "1-2-3", "-5-10"] {
            assert_eq!(
                message(parse(text).unwrap_err()),
                "BLAST engine error: Invalid specification of query location (Format: start-stop)\n",
                "{text:?}"
            );
        }
        for text in ["0-10", "10-0"] {
            assert!(message(parse(text).unwrap_err())
                .contains("(range elements cannot be less than or equal to 0)"));
        }
        for text in ["1-1", "5-5"] {
            assert!(message(parse(text).unwrap_err()).contains("(range cannot be empty)"));
        }
        assert!(message(parse("20-10").unwrap_err()).contains("(start cannot be larger than stop)"));
        for text in [
            " 1-5",
            "1-5 ",
            "a-5",
            "1-5a",
            "1.0-5",
            "0x10-0x20",
            "1-2147483648",
            "0-a",
        ] {
            let error = message(parse(text).unwrap_err());
            assert!(
                error.contains("cannot convert to an int") && error.contains("LOSAT's BLASTN"),
                "{text:?}: {error}"
            );
        }
        assert!(
            message(parse_sequence_range("5-3", RangeRole::Subject, "TBLASTN").unwrap_err())
                .starts_with("BLAST engine error: Invalid specification of subject location")
        );
    }

    #[test]
    fn intervals_follow_the_record_length() {
        let range = SequenceRange { from: 9, to: 19 };
        assert_eq!(
            record_interval(&range, 30),
            RecordInterval::Letters { from: 9, to: 19 }
        );
        assert_eq!(
            record_interval(&range, 15),
            RecordInterval::Letters { from: 9, to: 14 }
        );
        assert_eq!(
            record_interval(&range, 10),
            RecordInterval::Letters { from: 9, to: 9 }
        );
        assert_eq!(record_interval(&range, 9), RecordInterval::Empty);
        assert_eq!(record_interval(&range, 8), RecordInterval::PastEnd);
    }

    #[test]
    fn skipped_queries_keep_their_ordinals_and_batches_count_whole_records() {
        let records: Vec<crate::blastinput::fasta_reader::FastaRecord> =
            [8000usize, 400, 3000, 200]
                .iter()
                .enumerate()
                .map(|(index, &length)| {
                    crate::blastinput::fasta_reader::FastaRecord::new(
                        format!("Query_{}", index + 1),
                        format!("q{index}").as_bytes(),
                        &vec![b'A'; length],
                    )
                })
                .collect();
        let ranged = cut_queries(&records, &SequenceRange { from: 500, to: 999 });
        let input = &ranged.input;
        assert_eq!(input.ordinals, vec![0, 2]);
        assert_eq!(input.skipped, vec![false, true, false, true]);
        assert_eq!(input.range_length, Some(500));
        assert_eq!(ranged.records[0].seq().len(), 500);
        assert_eq!(input.placements.offset(1), 500);
        assert_eq!(input.placements.length(1, 500), 3000);
        // 8000 letters reach a batch of 5000: the next batch holds q1 (skipped), q2, q3.
        assert_eq!(
            input.batches(5000),
            vec![
                RangedBatch {
                    input: 0..1,
                    searched: 0..1
                },
                RangedBatch {
                    input: 1..4,
                    searched: 1..2
                }
            ]
        );
        // A batch of 3000 letters closes after q2; q3 alone is skipped: an empty batch.
        assert_eq!(
            input.batches(3000),
            vec![
                RangedBatch {
                    input: 0..1,
                    searched: 0..1
                },
                RangedBatch {
                    input: 1..3,
                    searched: 1..2
                },
                RangedBatch {
                    input: 3..4,
                    searched: 2..2
                }
            ]
        );
    }

    #[test]
    fn frames_are_those_of_the_whole_record() {
        assert_eq!(record_frame(true, 0, 100), 1);
        assert_eq!(record_frame(true, 4, 100), 2);
        assert_eq!(record_frame(false, 99, 100), -1);
        assert_eq!(record_frame(false, 97, 100), -3);
    }
}

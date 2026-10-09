//! Query batches of a search (NCBI `CBlastInput::GetNextSeqBatch`).

use std::ops::Range;

use crate::algorithm::blastn::blast_engine::{QueryReading, NO_DATA_MESSAGE};
use crate::blastinput::fasta_reader::FastaRecord;

/// Splits queries of the given lengths into NCBI's query batches.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:138-170
/// ```c++
///     TSeqPos size_read = 0;
///     ...
///     while (size_read < GetBatchSize()) {
///         ...
///             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
///         ...
///         retval->AddQuery(q);
/// ```
/// The query that reaches the batch size stays in the current batch. The batch size
/// depends on the program (`blast_input_aux.cpp:104-119`).
pub fn query_batches(lengths: &[usize], batch_size: usize) -> Vec<Range<usize>> {
    let mut batches = Vec::new();
    let mut start = 0;
    while start < lengths.len() {
        let mut end = start;
        let mut residues = 0usize;
        while end < lengths.len() && residues < batch_size {
            residues += lengths[end];
            end += 1;
        }
        batches.push(start..end);
        start = end;
    }
    batches
}

/// The end (a query index) of the query batch that starts at `start`: queries are added
/// while the residues read are fewer than the batch size.
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
///             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
///         ...
///         retval->AddQuery(q);
///     }
/// ```
/// The batch size and `size_read` are `TSeqPos` (32-bit unsigned): an `Int4` batch size
/// converts to it (a negative one becomes 2^32 plus the value).
pub fn next_query_batch_end(lengths: &[usize], start: usize, batch_size: u32) -> usize {
    let mut end = start;
    let mut size_read: u32 = 0;
    while end < lengths.len() && size_read < batch_size {
        size_read = size_read.wrapping_add(lengths[end] as u32);
        end += 1;
    }
    end
}

/// `NStr::StringToInt` with the default flags: an optional sign and decimal digits, within
/// the range of an `int`, and nothing else (no spaces); `None` where NCBI throws its
/// `CStringException`.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbistr.cpp:635-643
/// ```c
/// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
/// {
///     S2N_CONVERT_GUARD_EX(flags);
///     Int8 value = StringToInt8(str, flags, base);
///     if ( value < kMin_Int  ||  value > kMax_Int ) {
///         S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
///     }
///     return (int) value;
/// }
/// ```
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbistr.cpp:810-862
/// ```c
///     // Determine sign
///     bool sign = false;
///     switch (str[pos]) {
///     case '-':
///         sign = true;
///         /*FALLTHRU*/
///     case '+':
///         pos++;
///         break;
///     ...
///     // Last checks
///     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
///         S2N_CONVERT_ERROR_INVAL(Int8);
///     }
/// ```
/// `i32::from_str` reads the same strings.
pub fn ncbi_string_to_int(value: &std::ffi::OsStr) -> Option<i32> {
    value.to_str()?.parse::<i32>().ok()
}

/// The query batch size of blastn from the environment variable `BATCH_SIZE` (0 when it
/// is not set: NCBI then uses `CBatchSizeMixer`). `Err` holds the value that NCBI cannot
/// convert (an empty value too).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
/// ```c
///     // used for experimentation purposes
///     char* batch_sz_str = getenv("BATCH_SIZE");
///     if (batch_sz_str) {
///         retval = NStr::StringToInt(batch_sz_str);
///         _TRACE("DEBUG: Using query batch size " << retval);
///         return retval;
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:99
/// ```c
///     if (! use_default) return 0;
/// ```
pub fn get_query_batch_size(value: Option<&std::ffi::OsStr>) -> Result<i32, String> {
    match value {
        None => Ok(0),
        Some(value) => {
            ncbi_string_to_int(value).ok_or_else(|| value.to_string_lossy().into_owned())
        }
    }
}

/// NCBI's adaptive query batch size of blastn: the first batch aims at 1/200 of the target
/// hits, and each later batch at the target over the (mixed) ratio of successful ungapped
/// extensions per residue of the batches before.
///
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:57-83
/// ```c++
/// class CBatchSizeMixer
/// {
/// private:
///     const double k_MixIn;        // mixing factor between batches
///     Int4 m_TargetHits;           // the target hits per batch
///     double m_Ratio;              // the hits to batch size ratio
///     Int4 m_BatchSize;            // the batch size for next run
///     const Int4 k_MinBatchSize;   // the minimum allowable batch size
///     const Int4 k_MaxBatchSize;   // the maximum allowable batch size
///     const Int4 k_MinTargetHits;  // the minimum target hits per batch
///
/// public:
///     CBatchSizeMixer(Int4 max_batch_size)
///           : k_MixIn        (0.3),
///             m_TargetHits   (1000000),
///             m_Ratio        (-1.0),
///             m_BatchSize    (5000),
///             k_MinBatchSize (100),
///             k_MaxBatchSize (max_batch_size),
/// 	    k_MinTargetHits(1000000) { }
///
///     void SetTargetHits(Int4 target) {
///         m_TargetHits = (target < k_MinTargetHits) ? k_MinTargetHits : target;
///         m_BatchSize = m_TargetHits / 200;
///     }
/// ```
pub struct BatchSizeMixer {
    target_hits: i32,
    ratio: f64,
    batch_size: i32,
    max_batch_size: i32,
}

impl BatchSizeMixer {
    const MIX_IN: f64 = 0.3;
    const MIN_BATCH_SIZE: i32 = 100;
    const MIN_TARGET_HITS: i32 = 1_000_000;

    pub fn new(max_batch_size: i32) -> Self {
        Self {
            target_hits: Self::MIN_TARGET_HITS,
            ratio: -1.0,
            batch_size: 5000,
            max_batch_size,
        }
    }

    pub fn set_target_hits(&mut self, target: i32) {
        self.target_hits = target.max(Self::MIN_TARGET_HITS);
        self.batch_size = self.target_hits / 200;
    }

    /// The next batch size; `hits` is the number of successful ungapped extensions of the
    /// batch just searched (`None` before the first batch).
    ///
    /// NCBI reference: c++/src/app/blast/blast_app_util.cpp:65-81
    /// ```c++
    /// Int4 CBatchSizeMixer::GetBatchSize(Int4 hits)
    /// {
    ///      if (hits >= 0) {
    ///          double ratio = 1.0 * (hits+1) / m_BatchSize;
    ///          m_Ratio = (m_Ratio < 0) ? ratio
    ///                  : k_MixIn * ratio + (1.0 - k_MixIn) * m_Ratio;
    ///          m_BatchSize = (Int4) (1.0 * m_TargetHits / m_Ratio);
    ///      }
    ///      if (m_BatchSize > k_MaxBatchSize) {
    ///          m_BatchSize = k_MaxBatchSize;
    ///          m_Ratio = -1.0;
    ///      } else if (m_BatchSize < k_MinBatchSize) {
    ///          m_BatchSize = k_MinBatchSize;
    ///          m_Ratio = -1.0;
    ///      }
    ///      return m_BatchSize;
    /// }
    /// ```
    /// The `(Int4)` conversion is the x86-64 one of the oracle (`INT_MIN` beyond `Int4`), so a
    /// batch with fewer than about 2.3 extensions per million target hits and residue is
    /// followed by the smallest batch.
    pub fn batch_size(&mut self, hits: Option<i32>) -> i32 {
        if let Some(hits) = hits.filter(|&hits| hits >= 0) {
            let ratio = f64::from(hits.wrapping_add(1)) / f64::from(self.batch_size);
            self.ratio = if self.ratio < 0.0 {
                ratio
            } else {
                Self::MIX_IN * ratio + (1.0 - Self::MIX_IN) * self.ratio
            };
            self.batch_size = crate::core::blast_util::ncbi_int4_from_double(
                f64::from(self.target_hits) / self.ratio,
            );
        }
        if self.batch_size > self.max_batch_size {
            self.batch_size = self.max_batch_size;
            self.ratio = -1.0;
        } else if self.batch_size < Self::MIN_BATCH_SIZE {
            self.batch_size = Self::MIN_BATCH_SIZE;
            self.ratio = -1.0;
        }
        self.batch_size
    }
}

/// The query batches as NCBI's batch loop reads them, up to the batch that stops the run.
pub(crate) struct QueryBatches {
    /// The reader's messages of each batch's records (skipped ones too), before the
    /// report of the batch's first searched query, by searched query.
    pub(crate) warnings: Vec<Vec<u8>>,
    /// The input records and the searched queries of the batches that NCBI searches.
    pub(crate) input_done: usize,
    pub(crate) searched_done: usize,
    /// What stops the run at the batch after them: the messages that NCBI writes as it
    /// reads that batch, and its error.
    pub(crate) stop: Option<(Vec<u8>, anyhow::Error)>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:134-171
// ```c
// CRef<CBlastQueryVector>
// CBlastInput::GetNextSeqBatch(CScope& scope)
// {
//     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
//     TSeqPos size_read = 0;
//
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
// ```c
//     if (queries.Empty()) {
//         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
// ```c
//         } catch (const CException& e) {
//             ...
//             CRef<CSearchMessage> m
//                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
//             messages[index].push_back(m);
//             s_InvalidateQueryContexts(qinfo, index);
//         }
//     ...
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
// ```c
// TSearchMessages::ToString() const
// {
//     string retval;
//     ITERATE(vector<TQueryMessages>, qm, *this) {
//         if (qm->empty()) {
//             continue;
//         }
//         ITERATE(TQueryMessages, msg, *qm) {
//             retval += (*msg)->GetMessage() + " ";
//         }
//     }
//     return retval;
// }
// ```
/// Reads the query batches as NCBI's loop reads them (`input` and `searched` are the
/// records read and the queries searched, `seq_range::cut_queries`): a batch counts the
/// whole records up to the batch size, and the reader's messages of its records come
/// before the report of its first searched query. The run stops at the first batch that
/// NCBI fails, after the reports of the batches before (without the epilog) and the
/// messages of that batch's records:
/// - the reading ends there with a reader error (or LOSAT's rejection of a Seq-id line);
/// - the batch has no query (a batch size of 0, a batch of records that `-query_loc`
///   skips, or a batch that an `eEOF` ends before its first record): `Empty
///   CBlastQueryVector`;
/// - every query of the batch has no letters: the context of a protein query is valid at
///   the set-up exactly when it has letters, so NCBI's set-up fails with one `Sequence
///   contains no data` per query (`GetMessage` puts the severity first) and `CATCH_ALL`
///   writes `BLAST engine error: ` (exit 3).
///
/// A batch whose queries have letters but no Karlin-Altschul parameters is set up and not
/// searched (`run_resolved_in_pool`).
pub(crate) fn read_query_batches(
    input_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &crate::blastinput::seq_range::QueryInput,
    mut reading: QueryReading,
    batch_size: u32,
) -> QueryBatches {
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); searched.len()];
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = crate::blastinput::seq_range::next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        let reading_ends_here =
            input_end == input.input_lengths.len() && !reached_size && reading.reads_past_records();
        // `CFastaReader` writes its messages about a batch's records (their lines and
        // titles) when it reads them, after the report of the batch before.
        let mut batch_warnings: Vec<u8> = input_records[input_start..input_end]
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        let mut stop_error = None;
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                batch_warnings.extend_from_slice(&error_warnings);
                stop_error = Some(error.into_app_error());
            }
        }
        if stop_error.is_none() && start == end {
            stop_error = Some(
                crate::cli::NativeError {
                    exit: 3,
                    message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
                }
                .into(),
            );
        }
        if stop_error.is_none()
            && searched[start..end]
                .iter()
                .all(|record| record.seq().is_empty())
        {
            let mut message = String::from("BLAST engine error: ");
            for _ in start..end {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            stop_error = Some(crate::cli::NativeError { exit: 3, message }.into());
        }
        if let Some(error) = stop_error {
            return QueryBatches {
                warnings,
                input_done: input_start,
                searched_done: start,
                stop: Some((batch_warnings, error)),
            };
        }
        warnings[start].extend_from_slice(&batch_warnings);
        input_start = input_end;
    }
    QueryBatches {
        warnings,
        input_done: input_start,
        searched_done: searched.len(),
        stop: None,
    }
}

#[cfg(test)]
mod tests {
    use super::{
        get_query_batch_size, ncbi_string_to_int, next_query_batch_end, query_batches,
        BatchSizeMixer,
    };

    #[test]
    fn the_query_that_reaches_the_batch_size_stays_in_the_batch() {
        assert_eq!(query_batches(&[19_999, 2, 1], 20_000), vec![0..2, 2..3]);
        assert_eq!(
            query_batches(&[], 20_000),
            Vec::<std::ops::Range<usize>>::new()
        );
    }

    #[test]
    fn batch_size_strings_are_read_as_ncbis_string_to_int() {
        use std::ffi::OsStr;
        for (text, value) in [
            ("100", Some(100)),
            ("+100", Some(100)),
            ("-1", Some(-1)),
            ("007", Some(7)),
            ("2147483647", Some(i32::MAX)),
            ("-2147483648", Some(i32::MIN)),
            ("2147483648", None),
            ("", None),
            (" 100", None),
            ("100 ", None),
            ("1e5", None),
            ("abc", None),
            ("+", None),
            ("1,000", None),
            ("0x10", None),
        ] {
            assert_eq!(ncbi_string_to_int(OsStr::new(text)), value, "{text:?}");
        }
        assert_eq!(get_query_batch_size(None), Ok(0));
        assert_eq!(get_query_batch_size(Some(OsStr::new("5000"))), Ok(5000));
        assert_eq!(
            get_query_batch_size(Some(OsStr::new(""))),
            Err(String::new())
        );
    }

    #[test]
    fn a_negative_batch_size_is_a_huge_tseqpos() {
        assert_eq!(
            next_query_batch_end(&[3000, 1999, 2, 7], 0, -1i32 as u32),
            4
        );
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 1), 1);
    }

    #[test]
    fn the_next_batch_takes_queries_until_the_batch_size_is_reached() {
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 5000), 3);
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 3, 5000), 4);
        assert_eq!(next_query_batch_end(&[3000, 2000, 2], 0, 5000), 2);
    }

    // The values follow NCBI's CBatchSizeMixer (blast_app_util.cpp:65-81).
    #[test]
    fn the_mixer_is_ncbis() {
        let mut mixer = BatchSizeMixer::new(1_000_000 - 1000);
        assert_eq!(mixer.batch_size(None), 5000);
        // 4999 hits: ratio 1.0, 1e6 residues, above the maximum.
        assert_eq!(mixer.batch_size(Some(4999)), 999_000);
        // 999 hits in 999000 residues: ratio 0.001 (reset after the clamp), 1e9, clamped.
        assert_eq!(mixer.batch_size(Some(998)), 999_000);
        // 9999 hits in 5000 residues: ratio 2.0, 500000 residues.
        let mut mixer = BatchSizeMixer::new(4_999_000);
        mixer.batch_size(None);
        assert_eq!(mixer.batch_size(Some(9999)), 500_000);
        // Mixed: 0.3 * (20000 / 500000) + 0.7 * 2.0 = 1.412 -> 708215.
        assert_eq!(mixer.batch_size(Some(19_999)), 708_215);
        // No hits: 1e6 / (1 / 5000) overflows Int4, INT_MIN on x86-64, so the smallest batch.
        let mut mixer = BatchSizeMixer::new(999_000);
        mixer.batch_size(None);
        assert_eq!(mixer.batch_size(Some(0)), 100);
        // A large subject set raises the target and the first batch, within the maximum.
        let mut mixer = BatchSizeMixer::new(4_999_000);
        mixer.set_target_hits(3_000_000_000_i64.min(i64::from(i32::MAX)) as i32);
        assert_eq!(mixer.batch_size(None), 4_999_000);
        let mut mixer = BatchSizeMixer::new(999_000);
        mixer.set_target_hits(60_000_000);
        assert_eq!(mixer.batch_size(None), 300_000);
    }
}

//! Query batches of a search (NCBI `CBlastInput::GetNextSeqBatch`).

use std::ops::Range;

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
pub fn next_query_batch_end(lengths: &[usize], start: usize, batch_size: i32) -> usize {
    let mut end = start;
    let mut residues = 0usize;
    while end < lengths.len() && (residues as i64) < i64::from(batch_size) {
        residues += lengths[end];
        end += 1;
    }
    end
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

#[cfg(test)]
mod tests {
    use super::{next_query_batch_end, query_batches, BatchSizeMixer};

    #[test]
    fn the_query_that_reaches_the_batch_size_stays_in_the_batch() {
        assert_eq!(query_batches(&[19_999, 2, 1], 20_000), vec![0..2, 2..3]);
        assert_eq!(
            query_batches(&[], 20_000),
            Vec::<std::ops::Range<usize>>::new()
        );
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

//! Query batches of a search (NCBI `CBlastInput::GetNextSeqBatch`).

use std::ops::Range;

/// Splits queries of the given lengths into NCBI's query batches.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:137-165
/// ```c++
///     while (size_read < GetBatchSize()) {
///         ...
///             size_read += sequence::GetLength(*q->GetQuerySeqLoc(), q->GetScope());
///             retval->AddQuery(q);
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

#[cfg(test)]
mod tests {
    use super::query_batches;

    #[test]
    fn the_query_that_reaches_the_batch_size_stays_in_the_batch() {
        assert_eq!(query_batches(&[19_999, 2, 1], 20_000), vec![0..2, 2..3]);
        assert_eq!(
            query_batches(&[], 20_000),
            Vec::<std::ops::Range<usize>>::new()
        );
    }
}

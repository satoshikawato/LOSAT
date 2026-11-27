#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::super::SeqRanges;
use super::*;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn chunk(length: usize) -> SubjectChunk {
    SubjectChunk {
        offset: 0,
        length,
        seq_ranges: SeqRanges::Inline {
            range: (0, length as i32),
        },
        masked: true,
        overlap: 0,
    }
}

/// The batches of `next_chunk_batch` of the chunks of `lengths` with at most `max_jobs` jobs
/// each: per batch the lengths of each job's chunks. Fails if a batch is empty or the chunks
/// come out of order, or if `next_chunk` is called after it ended.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn batches(lengths: &[usize], max_jobs: usize) -> Vec<Vec<Vec<usize>>> {
    let mut source = lengths.iter().enumerate();
    let mut ended = false;
    let (mut chunks, mut job_ends, mut kept) = (Vec::new(), Vec::new(), None);
    let mut found = Vec::new();
    let mut order = Vec::new();
    loop {
        let done = next_chunk_batch(&mut chunks, &mut job_ends, &mut kept, max_jobs, || {
            assert!(!ended, "a chunk asked for after the last");
            let next = source.next().map(|(index, &length)| SubjectChunk {
                offset: index,
                ..chunk(length)
            });
            ended = next.is_none();
            next
        });
        if chunks.is_empty() {
            assert!(done && kept.is_none());
            break;
        }
        assert!(job_ends.len() <= max_jobs && job_ends.last() == Some(&chunks.len()));
        order.extend(chunks.iter().map(|chunk| chunk.offset));
        found.push(
            batch_jobs(&chunks, &job_ends)
                .iter()
                .map(|job| job.iter().map(|chunk| chunk.length).collect())
                .collect(),
        );
        if done {
            assert!(kept.is_none());
            break;
        }
    }
    assert_eq!(order, (0..lengths.len()).collect::<Vec<_>>());
    found
}

// The batches of the pool's chunk jobs (review F4 of M4a-2b, review 2 S2 and T2): jobs are
// runs of consecutive chunks of at most MAX_DBSEQ_LEN letters in all, at least one chunk
// each, and a batch holds up to `max_jobs` of them; the chunk that would start one more job
// is kept for the next batch. A hard-masked subject's short chunks (chunks of no letters
// too) are one job, searched serially; full chunks are a job each, as are the chunks of a
// range over MAX_DBSEQ_LEN letters after a short one (`hard1` of
// `tests/blastn_db_hard_long.rs`: 1,001, 5,000,000 and 299,100 letters).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn chunk_batches_hold_jobs_of_at_most_ncbi_s_chunk_length() {
    const MAX: usize = MAX_DBSEQ_LEN;
    assert_eq!(batches(&[1000, 0, 2000], 4), [vec![vec![1000, 0, 2000]]]);
    assert_eq!(
        batches(&[MAX, MAX, 500_100], 8),
        [vec![vec![MAX], vec![MAX], vec![500_100]]]
    );
    assert_eq!(
        batches(&[MAX - 10, 10, 0, 1], 4),
        [vec![vec![MAX - 10, 10, 0], vec![1]]]
    );
    assert_eq!(
        batches(&[1001, MAX, 299_100], 4),
        [vec![vec![1001], vec![MAX], vec![299_100]]]
    );
    assert!(batches(&[], 4).is_empty());
    // The kept chunk: five jobs at two jobs per batch (two threads' worth of the probe `jb`
    // of review 2 at one thread), the chunk that starts the third job of a batch opens the
    // next one, with the chunks that join it.
    assert_eq!(
        batches(&[MAX, MAX, MAX - 5, 5, 3, MAX, 7], 2),
        [
            vec![vec![MAX], vec![MAX]],
            vec![vec![MAX - 5, 5], vec![3]],
            vec![vec![MAX], vec![7]],
        ]
    );
    assert_eq!(
        batches(&[2, 3, MAX, 4], 1),
        [vec![vec![2, 3]], vec![vec![MAX]], vec![vec![4]]]
    );
    // Random lengths: the jobs are the greedy runs (a chunk opens a new job only when the job
    // before cannot take it), every batch but the last is full.
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % bound) as usize
    };
    for _ in 0..2000 {
        let lengths: Vec<usize> = (0..next(30))
            .map(|_| match next(4) {
                0 => MAX,
                1 => next(MAX as u64 + 1),
                _ => next(5_000),
            })
            .collect();
        let max_jobs = 1 + next(6);
        let found = batches(&lengths, max_jobs);
        let jobs: Vec<&Vec<usize>> = found.iter().flatten().collect();
        for (index, batch) in found.iter().enumerate() {
            assert!(!batch.is_empty() && batch.len() <= max_jobs);
            assert!(index + 1 == found.len() || batch.len() == max_jobs);
        }
        for (index, job) in jobs.iter().enumerate() {
            let letters: usize = job.iter().sum();
            assert!(job.len() == 1 || letters <= MAX, "{lengths:?}");
            if let Some(following) = jobs.get(index + 1) {
                assert!(letters + following[0] > MAX, "{lengths:?} {max_jobs}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1526
// PD-LOSAT-NCBI-DEFECTS 1.21 (rule 3): `record_chunk_merge_error` orders cell entries by ordinal
// id while naming the sequence by its mapped database OID (`OidMap::oid`), so filtered databases
// (`-gilist`, `-seqidlist`, `-taxids`, alias ranges) report the real OID rather than the ordinal.
#[test]
fn record_chunk_merge_error_names_mapped_oid_and_orders_by_ordinal() {
    let cell: std::sync::Mutex<Option<(u32, anyhow::Error)>> = std::sync::Mutex::new(None);
    let oid_of = |ordinal: u32| ordinal + 100;

    record_chunk_merge_error(&cell, 42, oid_of, "BLASTN");
    {
        let guard = cell.lock().unwrap();
        let (key, error) = guard.as_ref().expect("recorded ordinal 42");
        assert_eq!(*key, 42);
        assert!(
            error.is::<crate::algorithm::common::hard_chunk_scan::ShortChunkSeedRefusal>(),
            "expected ShortChunkSeedRefusal, got {error:?}"
        );
        let text = error.to_string();
        assert!(text.contains("BLASTN: sequence of OID 142"), "{text}");
        assert!(!text.contains("OID 42"), "{text}");
        assert!(text.contains("undefined in NCBI BLAST+"), "{text}");
    }

    // A later subject (higher ordinal) does not replace an earlier ordinal.
    record_chunk_merge_error(&cell, 50, oid_of, "BLASTN");
    {
        let guard = cell.lock().unwrap();
        let (key, error) = guard.as_ref().expect("kept ordinal 42");
        assert_eq!(*key, 42);
        assert!(error.to_string().contains("BLASTN: sequence of OID 142"));
    }

    // An earlier subject (lower ordinal) replaces the higher ordinal and names its mapped OID.
    record_chunk_merge_error(&cell, 5, oid_of, "BLASTN");
    {
        let guard = cell.lock().unwrap();
        let (key, error) = guard.as_ref().expect("replaced with ordinal 5");
        assert_eq!(*key, 5);
        let text = error.to_string();
        assert!(text.contains("BLASTN: sequence of OID 105"), "{text}");
        assert!(!text.contains("OID 5"), "{text}");
    }
}

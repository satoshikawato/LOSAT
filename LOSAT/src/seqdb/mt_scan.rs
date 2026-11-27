//! NCBI's multithreaded database scan: the chunks `CSeqDBImpl::x_FillSeqBuffer` cuts and how far
//! its reads go past the end of the iteration range.
//!
//! A search with two or more threads in DB-split mode takes the subjects in chunks of raw OIDs.
//! A chunk reads its first OID, then reads the OID after each one it takes, and takes the OID read
//! while its length is not negative, the bytes taken so far leave room for it (the budget:
//! [`mt_chunk_budget`]) and the volume-relative OID after it is below the range's end.
//! NCBI compares that volume-relative OID with the end as a database-wide OID, so in a volume
//! after the first it reads on past the end, as far as the budget, the volume and the first
//! negative length allow. Those reads cannot change what a search reports (the chunk is cut at the
//! end afterwards), but a read of a nucleotide's last byte past the mapping crashes NCBI.

use super::volume::ScanProblem;
use super::{alias, SeqDb, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:496-501,597
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697
/// The bases one chunk of `x_FillSeqBuffer` may take after its first subject: the slice size (the
/// largest file the atlas measured so far, at most 1 GiB: `largest_file`) over four times the
/// threads, plus one. The one owner of the formula (the scheduling of tblastn's `-subject` cuts
/// its chunks by it too).
pub(crate) fn mt_chunk_budget(threads: usize, largest_file: u64) -> u64 {
    // NCBI's `Uint8` division by the `int` product of 4 and `m_NumThreads`; usize fits u64.
    let divisor = u64::try_from(threads.max(1))
        .unwrap_or(u64::MAX)
        .saturating_mul(4);
    largest_file.min(1 << 30) / divisor + 1
}

/// The inputs of [`last_chunk_reach`]: the scan of one volume.
pub(super) struct ChunkWalk {
    pub threads: usize,
    /// The largest file measured when the scan cuts the volume's first chunk, and when it cuts
    /// a later one (the volume's sequence file is measured by the first read).
    pub first_file: u64,
    pub later_file: u64,
    /// The volume-relative OID the scan's first chunk of the volume starts at.
    pub start: u32,
    pub num_oids: u32,
    /// The database-wide OID of the volume's first sequence, and the range's end.
    pub volume_start: u32,
    pub range_end: u32,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:324-373,684-710
/// The highest volume-relative OID that the chunks of `x_FillSeqBuffer` read in the volume of
/// `walk`, from the chunk at `walk.start` to the one that reaches the range's end, with
/// `length(oid)` the length `x_GetSequence` returns (negative: the chunk stops, and no more is
/// read), or the error of a read that crashes.
pub(super) fn last_chunk_reach<E>(
    walk: &ChunkWalk,
    mut length: impl FnMut(u32) -> Result<i32, E>,
) -> Result<u32, E> {
    let range_end = u64::from(walk.range_end);
    let mut begin = walk.start;
    let mut reach = begin;
    let mut first = true;
    loop {
        let file = if first {
            walk.first_file
        } else {
            walk.later_file
        };
        // The budget is at most 2^27 + 1, as is every length taken from it (`Int8`).
        let mut room = i64::try_from(mt_chunk_budget(walk.threads, file)).unwrap_or(i64::MAX);
        let mut at = begin;
        let mut len = length(at)?;
        at += 1;
        if len < 0 {
            return Ok(reach);
        }
        let mut taken = 0u32;
        loop {
            room -= i64::from(len);
            taken += 1;
            len = if at < walk.num_oids {
                reach = reach.max(at);
                length(at)?
            } else {
                // A read past the volume's last OID returns -1 without reading.
                -1
            };
            at += 1;
            if !(len >= 0 && room >= i64::from(len) && u64::from(at) < range_end) {
                break;
            }
        }
        // The chunk ends where it stopped taking, at most at the range's end.
        let next = u64::from(walk.volume_start) + u64::from(begin) + u64::from(taken);
        if next >= range_end {
            return Ok(reach);
        }
        begin = (next - u64::from(walk.volume_start)) as u32;
        if begin >= walk.num_oids {
            return Ok(reach);
        }
        first = false;
    }
}

impl SeqDb {
    /// The volume after the first whose OIDs include `end`, the iteration range's end, with an
    /// OID of the volume after it: where NCBI's volume-relative comparison with the end lets a
    /// chunk read past it.
    fn volume_past_range_end(&self, end: u32) -> Option<usize> {
        if end >= self.num_oids() {
            return None;
        }
        let at = self.volume_starts.partition_point(|&start| start <= end) - 1;
        (at > 0 && self.volume_starts[at] != end).then_some(at)
    }

    /// The first OID of volume `at` from its OID `from` on (volume-relative) that stops or
    /// crashes NCBI's scan, `None` where there is none.
    fn first_scan_problem_from(
        &self,
        at: usize,
        from: u32,
    ) -> Result<Option<ScanProblem>, SeqDbError> {
        const CHUNK: u32 = 131_072;
        let volume = &self.volumes[at];
        let start = self.volume_starts[at];
        let mut cur = from;
        while cur < volume.num_oids() {
            let hi = volume.num_oids().min(cur.saturating_add(CHUNK));
            if let Some(problem) = volume.x_first_scan_problem_in(cur, hi, start + cur)? {
                return Ok(Some(problem));
            }
            cur = hi;
        }
        Ok(None)
    }

    /// The largest file that NCBI's atlas has measured when the scan cuts its first chunk of
    /// volume `at`, for a scan that starts at the OID `begin`: the alias files and index files
    /// opened with the database, and the sequence files of the volumes the scan has read so far.
    fn largest_file_before_chunks(&self, begin: u32, at: usize) -> Result<u64, SeqDbError> {
        let mut largest = match &self.search_path {
            Some(path) => alias::largest_resolution_file(&self.db_name, self.seqtype, path),
            None => 0,
        };
        for volume in &self.volumes {
            let index = super::volume::with_extension(volume.base(), self.seqtype, "in");
            if let Ok(meta) = std::fs::metadata(index) {
                largest = largest.max(meta.len());
            }
        }
        let first = self.volume_starts.partition_point(|&start| start <= begin) - 1;
        for volume in &self.volumes[first..at] {
            if volume.num_oids() > 0 {
                largest = largest.max(volume.sequence_file_len()?);
            }
        }
        Ok(largest)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:324-373,684-710
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1785-1838
    /// The errors of NCBI's reads past the end `end` of the iteration range `begin..end` in a
    /// volume after the first (`last_chunk_reach`), when the scan has read the OIDs up to
    /// `end` without meeting a negative length (the caller's check of the OIDs before it).
    /// NCBI's reads there are `last_chunk_reach`'s; of the OIDs after `end`, the first that
    /// stops the scan decides: a negative length ends it with nothing left to fail, a read that
    /// crashes fails the search when the chunks reach it. A volume whose OIDs past `end` all read
    /// cleanly takes no walk. `Ok` in every other case.
    pub(super) fn mt_check_past_range_end(&self, begin: u32, end: u32) -> Result<(), SeqDbError> {
        if self.mt_threads == 0 || self.seqtype != SeqType::Nucleotide || begin >= end {
            return Ok(());
        }
        let Some(at) = self.volume_past_range_end(end) else {
            return Ok(());
        };
        let volume = &self.volumes[at];
        let start = self.volume_starts[at];
        if end + 1 - start >= volume.num_oids() {
            return Ok(());
        }
        // The walk needs the lengths of the OIDs from the volume's first chunk on: only a crash
        // past the range's end can change the outcome.
        let Some(ScanProblem::Crash(_)) = self.first_scan_problem_from(at, end + 1 - start)?
        else {
            return Ok(());
        };
        let first_file = self.largest_file_before_chunks(begin, at)?;
        let walk = ChunkWalk {
            threads: self.mt_threads,
            first_file,
            later_file: first_file.max(volume.sequence_file_len()?),
            start: begin.max(start) - start,
            num_oids: volume.num_oids(),
            volume_start: start,
            range_end: end,
        };
        last_chunk_reach(&walk, |oid| volume.x_scan_length_i32(oid, start + oid)).map(|_| ())
    }
}

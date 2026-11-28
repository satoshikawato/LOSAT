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

use super::volume::with_extension;
use super::{Atlas, SeqDb, SeqDbError, SeqType};

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

/// LOSAT's refusal of a scan whose reads past the range's end crash or not by whether NCBI's data
/// loader made a lookup that LOSAT does not know it made (`SeqDb::mt_check_past_range_end`).
const MAYBE_LOOKUP_DECIDES: &str = "whether NCBI BLAST+ 2.17.0's multithreaded scan crashes reading past the end of the alias file's OID range depends on whether its reports looked a subject's Seq-align id up in the database, which they do or not in an order of the loads that LOSAT does not know (another subject's definition lines also hold the id)";

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

    /// Whether the files that a query batch reads after its scan can change the scan of a later
    /// batch: NCBI's multithreaded scan of a nucleotide database whose iteration range ends
    /// inside a volume after the first, whose chunks' budgets the largest file measured sets
    /// (`mt_check_past_range_end`). The readers note a batch's reads (`note_header_reads`,
    /// `note_lookup_reads`) only then.
    pub(crate) fn notes_batch_reads(&self) -> bool {
        self.mt_threads > 0
            && self.seqtype == SeqType::Nucleotide
            && self.filter.as_deref().is_some_and(|filter| {
                let (begin, end) = filter.range();
                begin < end && self.volume_past_range_end(end).is_some()
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:126-133
    /// Notes in the atlas the header file of the volume of each sequence of `ordinals` (ordinal
    /// ids: an alias file's range or list numbers its sequences from 0, [`SeqDb::oid`]) whose
    /// header NCBI reads (its first read opens the file): the subjects of a batch's Seq-aligns.
    pub(crate) fn note_header_reads(&self, ordinals: impl IntoIterator<Item = u32>) {
        let count = self
            .filter
            .as_deref()
            .map_or(self.num_oids(), |filter| filter.count());
        let mut last = None;
        for ordinal in ordinals {
            if ordinal >= count {
                continue;
            }
            let oid = self.oid(ordinal).get();
            let at = self.volume_starts.partition_point(|&start| start <= oid) - 1;
            if last != Some(at) {
                self.atlas
                    .size(with_extension(self.volumes[at].base(), self.seqtype, "hr"));
                last = Some(at);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:389-406
    /// Notes in the atlas the files that the data loader's lookups of the Seq-align ids `ids` read
    /// (`lookup::SeqidLookup::note_sized_files`): the loads of a query batch's reports.
    pub(crate) fn note_lookup_reads<'i>(
        &self,
        ids: impl IntoIterator<Item = &'i crate::seqid::ParsedSeqId>,
    ) -> Result<(), SeqDbError> {
        let lookup = self.seqid_lookup()?;
        for id in ids {
            lookup.note_sized_files(id, &self.atlas);
        }
        Ok(())
    }

    /// [`SeqDb::note_lookup_reads`] of lookups that NCBI's data loader may or may not make (a
    /// place before the lookup may answer, in an order LOSAT does not know): the files of each
    /// are noted apart ([`Atlas::maybe`]).
    pub(crate) fn note_maybe_lookup_reads<'i>(
        &self,
        ids: impl IntoIterator<Item = &'i crate::seqid::ParsedSeqId>,
    ) -> Result<(), SeqDbError> {
        let lookup = self.seqid_lookup()?;
        for id in ids {
            let files = Atlas::default();
            lookup.note_sized_files(id, &files);
            self.atlas.maybe(&files);
        }
        Ok(())
    }

    /// Takes the files that the readers of another database of the run (the repeat database of
    /// `-filtering_db`, `note_single_thread_search`) have measured, noted in `other`: NCBI's
    /// atlas is the process's.
    pub(crate) fn absorb_atlas(&self, other: &Atlas) {
        self.atlas.absorb(other);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:374-399
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:358-385
    // NCBI's search of the repeats runs without threads (`CLocalBlast` of a new `CSeqDB` that
    // shares the process's atlas): its scan reads each OID that the list and the iteration range
    // keep, so the sequence file of each volume with one, after the set-up's read of OID 0 where
    // the statistics read it; its Seq-aligns read the headers of the subjects of its hits.
    /// Notes in `atlas` the files that this database's reader has measured, and those that
    /// NCBI's search of it without threads measures, whose Seq-aligns are on the sequences
    /// `subjects`.
    pub(crate) fn note_single_thread_search(
        &self,
        subjects: impl IntoIterator<Item = u32>,
        atlas: &Atlas,
    ) {
        if self.first_oid_length.is_some() {
            if let Some(volume) = self.volumes.iter().find(|volume| volume.num_oids() > 0) {
                self.atlas
                    .size(with_extension(volume.base(), self.seqtype, "sq"));
            }
        }
        let (begin, end) = self
            .filter
            .as_deref()
            .map_or((0, self.num_oids()), |filter| filter.range());
        for (at, volume) in self.volumes.iter().enumerate() {
            let (start, stop) = (self.volume_starts[at], self.volume_starts[at + 1]);
            let mut next = start.max(begin);
            if next >= stop.min(end) {
                continue;
            }
            if self.check_or_find_oid(&mut next).unwrap_or(false) && next < stop.min(end) {
                self.atlas
                    .size(with_extension(volume.base(), self.seqtype, "sq"));
            }
        }
        self.note_header_reads(subjects);
        atlas.absorb(&self.atlas);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-938
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1783-1785
    /// Notes in the atlas the files that NCBI measures before its scan cuts the first chunk of
    /// the volume `at`, in a scan that starts at the OID `begin`: the sequence file of OID 0's
    /// volume where the statistics read OID 0 (each batch's `BLAST_GapAlignSetUp`, before its
    /// scan), and the sequence file of each volume with OIDs that the scan read before (its
    /// first read opens it); then gives the largest file measured.
    fn largest_file_before_chunks(&self, begin: u32, at: usize) -> u64 {
        if self.first_oid_length.is_some() {
            if let Some(volume) = self.volumes.iter().find(|volume| volume.num_oids() > 0) {
                self.atlas
                    .size(with_extension(volume.base(), self.seqtype, "sq"));
            }
        }
        let first = self.volume_starts.partition_point(|&start| start <= begin) - 1;
        for volume in &self.volumes[first..at] {
            if volume.num_oids() > 0 {
                self.atlas
                    .size(with_extension(volume.base(), self.seqtype, "sq"));
            }
        }
        self.atlas.largest()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:684-710
    /// Whether a read past the end `end` of the iteration range that NCBI's chunks can make in the
    /// volume `at` crashes. A chunk starts before `end` and takes an OID while the bases it took
    /// leave room for it, so it reads an OID `j` past `end` only where the lengths of the OIDs
    /// from `end` to `j - 1` fit `budget` (the largest budget of the volume's chunks); a negative
    /// length ends every chunk that reads it. The OID at `end` was read and checked before
    /// (`SeqDb::mt_check_trailing_excluded`).
    fn crash_within_reach(&self, at: usize, end: u32, budget: u64) -> Result<bool, SeqDbError> {
        let volume = &self.volumes[at];
        let start = self.volume_starts[at];
        let budget = i64::try_from(budget).unwrap_or(i64::MAX);
        let mut oid = end - start;
        let mut taken = 0i64;
        while oid < volume.num_oids() {
            let len = match volume.x_scan_length_i32(oid, start + oid) {
                Ok(len) => len,
                Err(_) => return Ok(true),
            };
            if len < 0 {
                return Ok(false);
            }
            taken += i64::from(len);
            if taken > budget {
                return Ok(false);
            }
            oid += 1;
        }
        Ok(false)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:324-373,684-710
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1785-1838
    /// The errors of NCBI's reads past the end `end` of the iteration range `begin..end` in a
    /// volume after the first (`last_chunk_reach`), when the scan has read the OIDs up to
    /// `end` without meeting a negative length (the caller's check of the OIDs before it).
    /// NCBI's reads there are `last_chunk_reach`'s, with the budgets of the largest file that
    /// NCBI's atlas has measured when each chunk is cut (`super::atlas::Atlas`: before the
    /// volume's first chunk, then with the volume's sequence file, which the first read
    /// measures); of the OIDs after `end`, the first that stops the scan decides: a negative
    /// length ends it with nothing left to fail, a read that crashes fails the search when the
    /// chunks reach it. A volume whose OIDs past `end` that a chunk can reach all read cleanly
    /// takes no walk (`crash_within_reach`). `Ok` in every other case. Where NCBI may or may not
    /// have made a data loader's lookup that would measure a larger file ([`Atlas::maybe_largest`]),
    /// the reads are known only where both largest files give the same: else LOSAT refuses the
    /// database.
    pub(super) fn mt_check_past_range_end(&self, begin: u32, end: u32) -> Result<(), SeqDbError> {
        if self.mt_threads == 0 || self.seqtype != SeqType::Nucleotide || begin >= end {
            return Ok(());
        }
        let Some(at) = self.volume_past_range_end(end) else {
            return Ok(());
        };
        let volume = &self.volumes[at];
        let start = self.volume_starts[at];
        let first_file = self.largest_file_before_chunks(begin, at);
        let sequences = with_extension(volume.base(), self.seqtype, "sq");
        self.atlas.size(&sequences);
        let later_file = self.atlas.largest();
        let reads = |first_file: u64, later_file: u64| -> Result<(), SeqDbError> {
            if !self.crash_within_reach(at, end, mt_chunk_budget(self.mt_threads, later_file))? {
                return Ok(());
            }
            let walk = ChunkWalk {
                threads: self.mt_threads,
                first_file,
                later_file,
                start: begin.max(start) - start,
                num_oids: volume.num_oids(),
                volume_start: start,
                range_end: end,
            };
            last_chunk_reach(&walk, |oid| volume.x_scan_length_i32(oid, start + oid)).map(|_| ())
        };
        let known = reads(first_file, later_file);
        for maybe in self.atlas.maybe_largest() {
            if maybe <= first_file {
                continue;
            }
            if reads(first_file.max(maybe), later_file.max(maybe)).is_ok() != known.is_ok() {
                return Err(SeqDbError::Unsupported {
                    path: sequences,
                    message: MAYBE_LOOKUP_DECIDES.to_string(),
                });
            }
        }
        known
    }
}

#[cfg(test)]
mod tests {
    use super::{last_chunk_reach, mt_chunk_budget, ChunkWalk, MAYBE_LOOKUP_DECIDES};

    // Line RG D2 tail R3 (fixture `tests/fixtures/blastdbcheck/nonasc_xprog/s3db`, frozen from NCBI
    // BLAST+ 2.17.0 by the probe `D/d2/tail/work_r2/probe/freeze_s3.py`): the 2-thread scan of the
    // alias `r` (OIDs 150..154) with OID 157 of `v2` made to crash completes where NCBI has
    // measured no string ISAM, as when the subject's id `gb|AB900001.1|` is answered by another
    // sequence's blob without a lookup. That lookup would measure `v0.nsd` (15,380 bytes), whose
    // budget lets the chunks reach OID 157: made, the scan crashes; made or not, LOSAT does not
    // know whether NCBI crashes and refuses the database.
    #[test]
    fn a_lookup_that_may_not_be_made_decides_or_refuses() {
        use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/blastdbcheck/nonasc_xprog/s3db");
        let dir = crate::utils::test_dir::TempDir::unique("mt-maybe-lookup");
        for entry in std::fs::read_dir(&fixture).unwrap() {
            let entry = entry.unwrap();
            std::fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
        }
        // `v2.nin`: the ambiguity offset of OID 6 of `v2` (OID 157) 2^28, past the sequence file.
        let nin = dir.join("v2.nin");
        let mut bytes = std::fs::read(&nin).unwrap();
        let word = |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
        let mut at = 8;
        at += 4 + word(at);
        at += 4 + word(at);
        let count = word(at);
        let amb = at + 16 + 8 * (count + 1);
        bytes[amb + 24..amb + 28].copy_from_slice(&0x1000_0000u32.to_be_bytes());
        std::fs::write(&nin, &bytes).unwrap();
        let held = crate::seqid::ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
            accession: Some(b"AB900001".to_vec()),
            version: Some(1),
            ..Default::default()
        });
        let scan = |note: fn(&SeqDb, &crate::seqid::ParsedSeqId)| {
            let path = SearchPath::new(&dir, None, None);
            let mut db = SeqDb::open_search_in(b"r", SeqType::Nucleotide, &path).unwrap();
            db.set_mt_stop_at_negative(true);
            db.set_mt_threads(2);
            note(&db, &held);
            let mut first = 0;
            while let Some(range) = db.next_block(first, 1 << 20) {
                first = range.end_ordinal();
                db.read_nucleotide_search_block(range)?;
            }
            Ok::<(), SeqDbError>(())
        };
        scan(|_, _| {}).unwrap();
        let made = scan(|db, id| db.note_lookup_reads([id]).unwrap()).unwrap_err();
        assert!(
            matches!(&made, SeqDbError::Format { message, .. }
                if message.contains("passes the end of the file")),
            "{made:?}"
        );
        let unknown = scan(|db, id| db.note_maybe_lookup_reads([id]).unwrap()).unwrap_err();
        assert!(
            matches!(&unknown, SeqDbError::Unsupported { path, message }
                if message == MAYBE_LOOKUP_DECIDES && path.ends_with("v2.nsq")),
            "{unknown:?}"
        );
    }

    /// The lengths of the 30 OIDs of the volume after the first (3 OIDs) in the cases below.
    fn lens() -> Vec<i32> {
        (0..30).map(|i| 10 + (i * 7) % 13).collect()
    }

    /// `last_chunk_reach` on the volume after a 3-OID volume, as a database-wide OID.
    fn reach(range_end: u32, threads: usize, files: (u64, u64), begin: u32) -> u32 {
        let lens = lens();
        let walk = ChunkWalk {
            threads,
            first_file: files.0,
            later_file: files.1,
            start: begin.max(3) - 3,
            num_oids: 30,
            volume_start: 3,
            range_end,
        };
        3 + last_chunk_reach(&walk, |oid| Ok::<_, ()>(lens[oid as usize])).unwrap()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:684-710
    /// The highest OID read, with the cases of the probes of line RG D2 (tail): the expected
    /// numbers are those of NCBI's loop, which compares the volume-relative OID with the range's
    /// end (the model in the probes' `sim.py` agreed with NCBI BLAST+ 2.17.0 on 150 databases).
    #[test]
    fn chunks_read_past_the_range_end_in_a_volume_after_the_first() {
        // A budget of 251 bases takes about 14 of the OIDs of 10..22 bases, but the volume-relative
        // OID stops each chunk at 5 (end 6): the OID of rank 5 is the last read, database OID 8.
        assert_eq!(reach(6, 2, (2000, 2000), 0), 8);
        assert_eq!(reach(20, 2, (2000, 2000), 0), 22);
        assert_eq!(reach(20, 4, (2000, 2000), 0), 22);
        assert_eq!(reach(20, 8, (2000, 2000), 0), 21);
        // The volume's sequence file counts from the second chunk on: 3000 bytes raise the
        // budget of the chunks after the first.
        assert_eq!(reach(12, 2, (500, 3000), 0), 14);
        assert_eq!(reach(12, 2, (500, 500), 0), 13);
        // Chunks too small to take a second OID read the OID after the first one only.
        assert_eq!(reach(8, 2, (160, 160), 0), 8);
        assert_eq!(reach(25, 2, (60, 60), 0), 25);
        // A scan that starts inside the volume (range begin 5).
        assert_eq!(reach(8, 2, (300, 300), 5), 9);
    }

    #[test]
    fn a_negative_length_ends_the_reads_and_a_crash_is_the_error() {
        let lens = lens();
        let walk = ChunkWalk {
            threads: 2,
            first_file: 2000,
            later_file: 2000,
            start: 0,
            num_oids: 30,
            volume_start: 3,
            range_end: 20,
        };
        // OID 4 of the volume is negative: nothing after it is read.
        let at = last_chunk_reach(&walk, |oid| {
            Ok::<_, ()>(if oid == 4 { -1 } else { lens[oid as usize] })
        });
        assert_eq!(at, Ok(4));
        // A read that crashes is reached by the chunk (range end 20, reach 22: volume OID 19).
        let at = last_chunk_reach(&walk, |oid| {
            if oid == 19 {
                Err(oid)
            } else {
                Ok(lens[oid as usize])
            }
        });
        assert_eq!(at, Err(19));
        // And one past the reach is not.
        let at = last_chunk_reach(&walk, |oid| {
            if oid == 20 {
                Err(oid)
            } else {
                Ok(lens[oid as usize])
            }
        });
        assert_eq!(at, Ok(19));
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:496-501,597
    #[test]
    fn the_budget_is_the_slice_over_four_times_the_threads_plus_one() {
        assert_eq!(mt_chunk_budget(2, 4000), 501);
        assert_eq!(mt_chunk_budget(8, 163), 6);
        assert_eq!(mt_chunk_budget(4, 0), 1);
        assert_eq!(mt_chunk_budget(4, u64::MAX), (1 << 30) / 16 + 1);
    }
}

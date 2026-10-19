//! One database volume: its index, sequence and header files, and blocks of
//! sequences read by position.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use super::file::DbFile;
use super::header::{decode_defline_set, Defline};
use super::index::{self, IndexFile};
use super::{SeqDbError, SeqType, XDbIdxRead};

/// One volume of a database (`NAME.pin` / `.psq` / `.phr`, or the `n` files).
#[derive(Debug)]
pub struct Volume {
    base: PathBuf,
    index: IndexFile,
    seq: OnceLock<DbFile>,
    hdr: OnceLock<DbFile>,
}

/// File name of `base` with the three-letter extension `{p|n}{ext}`.
pub(super) fn with_extension(base: &Path, seqtype: SeqType, ext: &str) -> PathBuf {
    let mut name = base.as_os_str().to_os_string();
    name.push(".");
    name.push(seqtype.letter());
    name.push(ext);
    PathBuf::from(name)
}

impl Volume {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:104-109
    // ```c++
    //     m_Idx.Reset(new CSeqDBIdxFile(atlas, name, prot_nucl));
    //
    //     m_VolEnd = m_VolStart + m_Idx->GetNumOIDs();
    //
    //     // To allow for empty volumes, we must tolerate the absence of all
    //     // files other than the index file.
    // ```
    /// Opens the volume `base` (a path without extension) by reading its
    /// index file. The sequence and header files are opened on first use.
    /// EXPERIMENT (LOSAT_X_DBIDXREAD): the index file stays open and its arrays
    /// are read where they are needed (`IndexFile::open`); without the switch
    /// the whole file is read and its arrays are held (`IndexFile::parse`).
    pub fn open(base: &Path, seqtype: SeqType) -> Result<Self, SeqDbError> {
        Self::open_as(base, seqtype, super::x_dbidxread())
    }

    /// `open` with the index read as `mode` says.
    pub(super) fn open_as(
        base: &Path,
        seqtype: SeqType,
        mode: XDbIdxRead,
    ) -> Result<Self, SeqDbError> {
        let index_path = with_extension(base, seqtype, "in");
        let index_file = DbFile::open(&index_path)?;
        let index = match mode {
            XDbIdxRead::Off => IndexFile::parse(&index_file.read_all()?, seqtype, &index_path)?,
            XDbIdxRead::On => IndexFile::open(index_file, seqtype, false)?,
            XDbIdxRead::Shadow => IndexFile::open(index_file, seqtype, true)?,
        };
        let volume = Self {
            base: base.to_path_buf(),
            index,
            seq: OnceLock::new(),
            hdr: OnceLock::new(),
        };
        volume.reject_column_files()?;
        Ok(volume)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3333-3343
    // ```c++
    //     string alpha("abcdefghijklmnopqrstuvwxyz");
    //     string ei("??a"), ed("??b"), ed2("??c");
    //
    //     ei[0] = ed[0] = ed2[0] = (m_IsAA ? 'p' : 'n');
    //
    //     map<string,int> unique_titles;
    //
    //     for(size_t i = 0; i < alpha.size(); i++) {
    //         ei[1] = ed[1] = ed2[1] = alpha[i];
    //
    //         if (CSeqDBColumn::ColumnExists(m_VolName, ei, m_Atlas)) {
    // ```
    // Column files hold per-sequence data such as the masks of
    // `makeblastdb -mask_data`. LOSAT does not read them, so a volume that has
    // any is rejected instead of being searched without them.
    fn reject_column_files(&self) -> Result<(), SeqDbError> {
        for letter in b'a'..=b'z' {
            let ext = [letter, b'a'];
            let ext = std::str::from_utf8(&ext).unwrap_or("aa");
            let path = with_extension(&self.base, self.index.seqtype, ext);
            if path.is_file() {
                return Err(SeqDbError::Unsupported {
                    path,
                    message: "the volume has column data (for example the masks of makeblastdb \
                              -mask_data); LOSAT does not read column files"
                        .to_string(),
                });
            }
        }
        Ok(())
    }

    /// Path of the volume without extension.
    pub fn base(&self) -> &Path {
        &self.base
    }

    /// The parsed index file.
    pub fn index(&self) -> &IndexFile {
        &self.index
    }

    /// Number of sequences in the volume.
    pub fn num_oids(&self) -> u32 {
        self.index.num_oids
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:117-134
    // ```c++
    // void
    // CSeqDBVol::x_OpenSeqFile(void) const {
    //     CFastMutexGuard mtx_gurad(m_MtxSeq);
    //     if (!m_SeqFileOpened && m_Idx->GetNumOIDs() != 0) {
    //         m_Seq.Reset(new CSeqDBSeqFile(m_Atlas, m_VolName, (m_IsAA?'p':'n')));
    //     }
    //     m_SeqFileOpened = true;
    // }
    //
    // void
    // CSeqDBVol::x_OpenHdrFile(void) const {
    //     CFastMutexGuard mtx_gurad(m_MtxHdr);
    //     if (!m_HdrFileOpened && m_Idx->GetNumOIDs() != 0) {
    //         m_Hdr.Reset(new CSeqDBHdrFile(m_Atlas, m_VolName, (m_IsAA?'p':'n')));
    //     }
    //     m_HdrFileOpened = true;
    // }
    // ```
    fn lazy_file<'a>(
        &self,
        slot: &'a OnceLock<DbFile>,
        ext: &str,
    ) -> Result<&'a DbFile, SeqDbError> {
        if let Some(file) = slot.get() {
            return Ok(file);
        }
        let file = DbFile::open(&with_extension(&self.base, self.index.seqtype, ext))?;
        // Two threads may both open the file; the first one stored is used.
        Ok(slot.get_or_init(|| file))
    }

    fn seq_file(&self) -> Result<&DbFile, SeqDbError> {
        self.lazy_file(&self.seq, "sq")
    }

    fn hdr_file(&self) -> Result<&DbFile, SeqDbError> {
        self.lazy_file(&self.hdr, "hr")
    }

    fn check_range(&self, first: u32, end: u32) -> Result<(), SeqDbError> {
        if first >= end || end > self.index.num_oids {
            return Err(SeqDbError::Range(format!(
                "sequences {first}..{end} are not a non-empty range of the {} sequences of {}",
                self.index.num_oids,
                self.base.display()
            )));
        }
        Ok(())
    }

    /// Bytes before the first sequence that a block reads with it: the
    /// protein sentinel byte.
    fn lead(&self) -> u64 {
        match self.index.seqtype {
            SeqType::Protein => 1,
            SeqType::Nucleotide => 0,
        }
    }

    /// Bytes of the sequence file that a block of the sequences `first..end`
    /// occupies (protein: with the sentinel bytes around them; nucleotide:
    /// with their ambiguity data).
    pub fn block_bytes(&self, first: u32, end: u32) -> Result<u64, SeqDbError> {
        Ok(self.index.sequence_start(end)? - self.index.sequence_start(first)? + self.lead())
    }

    /// The end of the block that starts at `first` (`SeqDb::next_block`): the
    /// last end in `first + 1 ..= num_oids` whose `block_bytes` fit
    /// `max_bytes`, at least `first + 1`. The offsets ascend, so the ends that
    /// fit form a prefix.
    pub(super) fn block_end(&self, first: u32, max_bytes: u64) -> Result<u32, SeqDbError> {
        if let Some(offsets) = self.index.held_sequence_offsets() {
            return Ok(held_block_end(offsets, first, max_bytes, self.lead()));
        }
        let end = self.x_read_block_end(first, max_bytes, X_BLOCK_END_WINDOW)?;
        if let Some(offsets) = self.index.x_shadow_sequence_offsets() {
            index::X_SHADOW_BLOCK_ENDS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let expected = held_block_end(offsets, first, max_bytes, self.lead());
            if end != expected {
                index::x_shadow_difference(format!(
                    "{}: the block from {first} within {max_bytes} bytes ends at {end}, the held \
                     array's at {expected}",
                    self.base.display()
                ));
            }
        }
        Ok(end)
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): `block_end` over the offsets in the
    /// file: the end of the volume when it fits (one entry read), else the
    /// offsets after `first` read forward in windows of `window` entries,
    /// doubled up to `X_BLOCK_END_MAX_WINDOW`, to the first end that does not
    /// fit.
    pub(super) fn x_read_block_end(
        &self,
        first: u32,
        max_bytes: u64,
        window: usize,
    ) -> Result<u32, SeqDbError> {
        let n = self.num_oids();
        let lead = self.lead();
        if max_bytes < lead {
            return Ok(first + 1);
        }
        // `block_bytes(first, end) <= max_bytes` is `sequence_start(end) <= limit`.
        let limit = self
            .index
            .sequence_start(first)?
            .saturating_add(max_bytes - lead);
        if self.index.sequence_start(n)? <= limit {
            return Ok(n);
        }
        let mut window = window.max(1);
        let mut next = first + 1;
        loop {
            // The end `n` does not fit, so the scan stops at or before it.
            let count = window.min((n + 1 - next) as usize);
            let offsets = self.index.sequence_offsets(next, next + count as u32 - 1)?;
            if let Some(at) = offsets.iter().position(|&offset| u64::from(offset) > limit) {
                return Ok((next + at as u32).saturating_sub(1).max(first + 1));
            }
            next += count as u32;
            window = (window * 2).min(X_BLOCK_END_MAX_WINDOW);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1838
    // ```c++
    //     if ('p' == seqtype) {
    //         // Subtract one, for the inter-sequence null.
    //
    //         end_offset --;
    //
    //         length = int(end_offset - start_offset);
    //
    //         // Although we subtracted one above to get the correct length,
    //         // we expand the range here by one byte in both directions.
    //         // The normal consumer of this data relies on them, and can
    //         // walk off memory if a sequence ends on a slice boundary.
    //         *buffer = m_Seq->GetFileDataPtr(start_offset-1) + 1;
    //         ...
    //     } else if ('n' == seqtype) {
    //         ...
    //         *buffer = m_Seq->GetFileDataPtr(start_offset);
    //         ...
    //         int whole_bytes = int(end_offset - start_offset - 1);
    //
    //         char last_char = (*buffer)[whole_bytes];
    //
    //         int remainder = last_char & 3;
    //         length = (whole_bytes * 4) + remainder;
    //     }
    // ```
    // NCBI points into the mapped file one sequence at a time; LOSAT copies
    // the bytes of a contiguous range of sequences with one positioned read.
    // The protein range starts one byte early, so every sequence of the block
    // has its sentinel byte on both sides as in NCBI's buffer.
    /// Reads the sequences `first..end` (volume-relative, non-empty) as one
    /// block. `first_ordinal` is the database ordinal id of `first`.
    pub fn read_block(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.check_range(first, end)?;
        let seqtype = self.index.seqtype;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        let data_start = u64::from(seq_offsets[0]) - self.lead();
        let data_end = u64::from(seq_offsets[seq_offsets.len() - 1]);
        let file = self.seq_file()?;
        let data = file.read_range(data_start, data_end)?;
        let block = SequenceBlock {
            seqtype,
            first_ordinal,
            path: file.path().to_path_buf(),
            data,
            data_start,
            seq_offsets,
            amb_offsets,
        };
        Ok(block)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2316
    // ```c++
    // CTempString
    // CSeqDBVol::x_GetHdrAsn1Binary(int oid) const
    // {
    //     TIndx hdr_start = 0;
    //     TIndx hdr_end   = 0;
    //
    //     //m_Atlas.Lock(locked);
    //
    //     if (!m_HdrFileOpened) x_OpenHdrFile();
    //
    //     m_Idx->GetHdrStartEnd(oid, hdr_start, hdr_end);
    //     const char * asn_region = m_Hdr->GetFileDataPtr(hdr_start);
    //
    //     return CTempString(asn_region, hdr_end - hdr_start);
    // }
    // ```
    /// The encoded `Blast-def-line-set` of the sequence `oid` (volume-relative).
    pub fn header_bytes(&self, oid: u32) -> Result<Vec<u8>, SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let (start, end) = self.index.header_range(oid)?;
        self.hdr_file()?.read_range(start, end)
    }

    /// The definition lines of the sequence `oid` (volume-relative), with
    /// their ordinal ids still relative to the volume.
    pub fn deflines(&self, oid: u32) -> Result<Vec<Defline>, SeqDbError> {
        let record = self.header_bytes(oid)?;
        decode_defline_set(
            &record,
            &with_extension(&self.base, self.index.seqtype, "hr"),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2364-2384
    // ```c++
    //     int length = int(end_offset - start_offset);
    //
    //     if (length) {
    //         int total = length / 4;
    //         ...
    //         ambchars.resize(total);
    //
    //         for(int i = 0; i<total; i++) {
    //             ambchars[i] = SeqDB_GetStdOrd((const int *)(& buffer[i]));
    //         }
    // ```
    /// Whether the ambiguity data of nucleotide `oid` (volume-relative;
    /// `ordinal` its database ordinal id, for the message) decode as
    /// `SequenceBlock::ambiguity_runs` decodes them: the same error, read
    /// without the sequence's letters (only the bytes between its ambiguity
    /// offset and the next sequence; nothing when it has none). A protein
    /// volume has no ambiguity data.
    pub fn check_ambiguity(&self, oid: u32, ordinal: u32) -> Result<(), SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let Some((amb_start, next)) = self.index.ambiguity_words(oid)? else {
            return Ok(());
        };
        self.check_ambiguity_data(amb_start, next, ordinal)
    }

    /// `check_ambiguity` of the sequence whose ambiguity data are
    /// `amb_start..next` of the sequence file.
    fn check_ambiguity_data(
        &self,
        amb_start: u32,
        next: u32,
        ordinal: u32,
    ) -> Result<(), SeqDbError> {
        if u64::from(amb_start) >= u64::from(next) {
            return Ok(());
        }
        let file = self.seq_file()?;
        let block = SequenceBlock {
            seqtype: SeqType::Nucleotide,
            first_ordinal: ordinal,
            path: file.path().to_path_buf(),
            data: file.read_range(u64::from(amb_start), u64::from(next))?,
            data_start: u64::from(amb_start),
            // The block holds the ambiguity data only: its letters' range is
            // empty, and nothing but `ambiguity_runs` reads it.
            seq_offsets: vec![amb_start, next],
            amb_offsets: vec![amb_start],
        };
        block.ambiguity_runs(0).map(|_| ())
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): calls `visit(oid, self.deflines(oid))`
    /// for the sequences `first..end` (volume-relative, non-empty) in order
    /// until `visit` returns false; returns whether every call returned true.
    /// The header offsets of the range are read at once and the records in
    /// pieces of at most `piece_bytes` (one longer record alone), each with one
    /// read; a record is decoded from its piece as `deflines` decodes it. A
    /// sequence whose record passes the end of the header file, and every
    /// sequence of a range or piece that cannot be read, goes through
    /// `deflines`, which gives its own result.
    pub(super) fn visit_deflines(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        visit: &mut dyn FnMut(u32, Result<Vec<Defline>, SeqDbError>) -> bool,
    ) -> bool {
        let one_by_one =
            |oids: std::ops::Range<u32>,
             visit: &mut dyn FnMut(u32, Result<Vec<Defline>, SeqDbError>) -> bool| {
                oids.into_iter().all(|oid| visit(oid, self.deflines(oid)))
            };
        if self.check_range(first, end).is_err() {
            return one_by_one(first..end, visit);
        }
        let (Ok(offsets), Ok(file)) = (self.index.header_offsets(first, end), self.hdr_file())
        else {
            return one_by_one(first..end, visit);
        };
        let hr_path = with_extension(&self.base, self.index.seqtype, "hr");
        let mut i = 0usize;
        let count = (end - first) as usize;
        while i < count {
            let start = u64::from(offsets[i]);
            if u64::from(offsets[i + 1]) > file.len() {
                // `read_range`'s error, from `deflines`.
                if !one_by_one(first + i as u32..first + i as u32 + 1, visit) {
                    return false;
                }
                i += 1;
                continue;
            }
            let mut stop = i + 1;
            while stop < count
                && u64::from(offsets[stop + 1]) <= file.len()
                && u64::from(offsets[stop + 1]) - start <= piece_bytes
            {
                stop += 1;
            }
            let piece_end = u64::from(offsets[stop]);
            let Ok(data) = file.read_range(start, piece_end) else {
                if !one_by_one(first + i as u32..first + stop as u32, visit) {
                    return false;
                }
                i = stop;
                continue;
            };
            for j in i..stop {
                let lo = (u64::from(offsets[j]) - start) as usize;
                let hi = (u64::from(offsets[j + 1]) - start) as usize;
                if !visit(
                    first + j as u32,
                    decode_defline_set(&data[lo..hi], &hr_path),
                ) {
                    return false;
                }
            }
            index::X_SHADOW_RUN_RECORDS
                .fetch_add((stop - i) as u64, std::sync::atomic::Ordering::Relaxed);
            i = stop;
        }
        true
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): calls `visit(oid, self.check_ambiguity(oid, ordinal))`
    /// for the sequences `first..end` (volume-relative, non-empty; `first_ordinal` the
    /// database ordinal id of `first`) in order until `visit` returns false; returns whether
    /// every call returned true. The sequence and ambiguity offsets of the range are read at
    /// once; the ambiguity data are read as `check_ambiguity` reads them, only for the
    /// sequences that have any. A range whose offsets cannot be read goes through
    /// `check_ambiguity` one sequence at a time.
    pub(super) fn visit_ambiguity(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        visit: &mut dyn FnMut(u32, Result<(), SeqDbError>) -> bool,
    ) -> bool {
        let offsets = match self.check_range(first, end) {
            Ok(()) if self.index.seqtype == SeqType::Nucleotide => self
                .index
                .sequence_offsets(first, end)
                .and_then(|seq| Ok((seq, self.index.ambiguity_offsets(first, end)?))),
            _ => Err(SeqDbError::Range(String::new())),
        };
        let Ok((seq, amb)) = offsets else {
            return (first..end).all(|oid| {
                visit(
                    oid,
                    self.check_ambiguity(oid, first_ordinal + (oid - first)),
                )
            });
        };
        (0..(end - first) as usize).all(|i| {
            visit(
                first + i as u32,
                self.check_ambiguity_data(amb[i], seq[i + 1], first_ordinal + i as u32),
            )
        })
    }
}

/// Entries of the first window that `Volume::x_read_block_end` reads.
const X_BLOCK_END_WINDOW: usize = 1024;
/// The largest window of `Volume::x_read_block_end` (1 MiB of offsets).
const X_BLOCK_END_MAX_WINDOW: usize = 1 << 18;

/// `Volume::block_end` over the held sequence offsets (`lead`: `Volume::lead`).
pub(super) fn held_block_end(offsets: &[u32], first: u32, max_bytes: u64, lead: u64) -> u32 {
    let n = (offsets.len() - 1) as u32;
    let block_bytes =
        |end: u32| u64::from(offsets[end as usize]) - u64::from(offsets[first as usize]) + lead;
    let fitting = super::partition_point_u32(first + 1..n + 1, |end| block_bytes(end) <= max_bytes);
    first + (fitting.max(1) as u32)
}

/// A contiguous range of sequences of one volume, held in memory.
///
/// Sequences are addressed by their position `i` in the block
/// (`0..self.len()`); the database ordinal id of position `i` is
/// `self.first_ordinal() + i`.
#[derive(Debug, Clone)]
pub struct SequenceBlock {
    seqtype: SeqType,
    first_ordinal: u32,
    path: PathBuf,
    data: Vec<u8>,
    data_start: u64,
    /// Start offsets of the block's sequences in the file, plus the end.
    seq_offsets: Vec<u32>,
    /// Ambiguity start offsets of the block's sequences (nucleotide).
    amb_offsets: Vec<u32>,
}

/// A run of one ambiguity code in a nucleotide sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AmbiguityRun {
    /// Offset of the first base of the run.
    pub offset: u32,
    /// Number of bases.
    pub length: u32,
    /// The ambiguity code (NCBI4na, 0-15).
    pub ncbi4na: u8,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:590-607
// ```c++
// unsigned SeqDB_ncbina8_to_blastna8[] = {
//     15, /* Gap, 0 */
//     0,  /* A,   1 */
//     1,  /* C,   2 */
//     6,  /* M,   3 */
//     2,  /* G,   4 */
//     4,  /* R,   5 */
//     9,  /* S,   6 */
//     13, /* V,   7 */
//     3,  /* T,   8 */
//     8,  /* W,   9 */
//     5,  /* Y,  10 */
//     12, /* H,  11 */
//     7,  /* K,  12 */
//     11, /* D,  13 */
//     10, /* B,  14 */
//     14  /* N,  15 */
// };
// ```
const NCBINA8_TO_BLASTNA8: [u8; 16] = [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];

/// Sentinel byte NCBI puts before and after a BLASTNA sequence.
pub const BLASTNA_SENTINEL: u8 = 15;

impl SequenceBlock {
    /// Sequence type of the block.
    pub fn seqtype(&self) -> SeqType {
        self.seqtype
    }

    /// Database ordinal id of the first sequence of the block.
    pub fn first_ordinal(&self) -> u32 {
        self.first_ordinal
    }

    /// Number of sequences in the block.
    pub fn len(&self) -> usize {
        self.seq_offsets.len() - 1
    }

    /// Whether the block holds no sequence (never true for a block that was read).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bytes of the sequence file held by the block.
    pub fn byte_len(&self) -> usize {
        self.data.len()
    }

    fn slice(&self, start: u64, end: u64) -> &[u8] {
        let lo = (start - self.data_start) as usize;
        let hi = (end - self.data_start) as usize;
        &self.data[lo..hi]
    }

    fn format(&self, i: usize, message: &str) -> SeqDbError {
        SeqDbError::Format {
            path: self.path.clone(),
            message: format!(
                "sequence {}: {message}",
                u64::from(self.first_ordinal) + i as u64
            ),
        }
    }

    fn expect(&self, seqtype: SeqType) {
        assert!(
            self.seqtype == seqtype,
            "sequence block of the other sequence type"
        );
    }

    /// The residues of protein `i` in NCBIstdaa, one byte each.
    pub fn protein(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Protein);
        self.slice(
            u64::from(self.seq_offsets[i]),
            u64::from(self.seq_offsets[i + 1]) - 1,
        )
    }

    /// Protein `i` with the sentinel byte (NUL in a well-formed volume)
    /// before and after its residues, as NCBI's buffer has them.
    pub fn protein_with_sentinels(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Protein);
        self.slice(
            u64::from(self.seq_offsets[i]) - 1,
            u64::from(self.seq_offsets[i + 1]),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1803
    // ```c++
    //     if ('p' == seqtype) {
    //         // Subtract one, for the inter-sequence null.
    //
    //         end_offset --;
    //
    //         length = int(end_offset - start_offset);
    //
    //         // Although we subtracted one above to get the correct length,
    //         // we expand the range here by one byte in both directions.
    //         // The normal consumer of this data relies on them, and can
    //         // walk off memory if a sequence ends on a slice boundary.
    //         *buffer = m_Seq->GetFileDataPtr(start_offset-1) + 1;
    //         if (! (*buffer - 1)) return -1;
    // ```
    // NCBI's search reads these bytes in place. A LOSAT block is a private
    // copy of them, so a search may rewrite it before it scans it (BLASTP reads
    // O as X); the layout stays the one of the file.
    /// The bytes of a protein block as read: the sentinel byte before the
    /// first sequence, then every sequence followed by its sentinel byte
    /// (`protein_with_sentinels(i)` is a window of them).
    pub fn protein_bytes_mut(&mut self) -> &mut [u8] {
        self.expect(SeqType::Protein);
        &mut self.data
    }

    /// The 2-bit data of nucleotide `i` as stored: four bases per byte from
    /// the high bits down, the last byte holding the remaining bases and,
    /// in its two low bits, their count.
    pub fn nucleotide_packed(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Nucleotide);
        self.slice(
            u64::from(self.seq_offsets[i]),
            u64::from(self.amb_offsets[i]),
        )
    }

    /// Length of nucleotide `i` in bases. A length above NCBI's signed
    /// 32-bit `int` is a format error.
    pub fn nucleotide_len(&self, i: usize) -> Result<usize, SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let whole_bytes = packed.len() - 1;
        let remainder = u64::from(packed[whole_bytes] & 3);
        let length = whole_bytes as u64 * 4 + remainder;
        if length > i32::MAX as u64 {
            return Err(self.format(i, "the sequence is longer than a signed 32-bit length"));
        }
        Ok(length as usize)
    }

    /// The ambiguity data of nucleotide `i` as stored (empty when the
    /// sequence has no ambiguous base).
    pub fn ambiguity_bytes(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Nucleotide);
        self.slice(
            u64::from(self.amb_offsets[i]),
            u64::from(self.seq_offsets[i + 1]),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2364-2384
    // ```c++
    //     int length = int(end_offset - start_offset);
    //
    //     if (length) {
    //         int total = length / 4;
    //         ...
    //         // This is probably unnecessary
    //         total &= 0x7FFFFFFF;
    //
    //         ambchars.resize(total);
    //
    //         for(int i = 0; i<total; i++) {
    //             ambchars[i] = SeqDB_GetStdOrd((const int *)(& buffer[i]));
    //         }
    //     } else {
    //         ambchars.clear();
    //     }
    // ```
    fn ambiguity_words(&self, i: usize) -> Vec<u32> {
        self.ambiguity_bytes(i)
            .chunks_exact(4)
            .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:642-713
    // ```c++
    // inline Uint4 s_ResLenNew(const vector<Int4> & ambchars, Uint4 i)
    // {
    //     return (ambchars[i] >> 16) & 0xFFF;
    // }
    // ...
    // inline Uint4 s_ResPosNew(const vector<Int4> & ambchars, Uint4 i)
    // {
    //     return ambchars[i+1];
    // }
    // ...
    // inline Uint4 s_ResVal(const vector<Int4> & ambchars, Uint4 i)
    // {
    //     return (ambchars[i] >> 28) & 0xF;
    // }
    // ...
    // inline Uint4 s_ResLenOld(const vector<Int4> & ambchars, Uint4 i)
    // {
    //     return (ambchars[i] >> 24) & 0xF;
    // }
    // ...
    // inline Uint4 s_ResPosOld(const vector<Int4> & ambchars, Uint4 i)
    // {
    //     return ambchars[i] & 0xFFFFFF; // RES_OFFSET
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:810-832
    // ```c++
    //     Uint4 amb_num = amb_chars[0];
    //
    //     /* Check if highest order bit set. */
    //     bool new_format = (amb_num & 0x80000000) != 0;
    //
    //     if (new_format)  amb_num &= 0x7FFFFFFF;
    //
    //     for(Uint4 i = 1; i < amb_num+1; i++) {
    //         Int4  row_len  = 0;
    //         Int4  position = 0;
    //         Uint1 trans_ch = 0;
    //
    // 	if (new_format) {
    //             trans_ch  = s_ResVal   (amb_chars, i);
    //             row_len   = s_ResLenNew(amb_chars, i) + 1;
    //             position  = s_ResPosNew(amb_chars, i);
    // 	} else {
    //             trans_ch  = s_ResVal   (amb_chars, i);
    //             row_len   = s_ResLenOld(amb_chars, i) + 1;
    //             position  = s_ResPosOld(amb_chars, i);
    // 	}
    //
    //         if (new_format) ++i;
    // ```
    // The entry count comes from the data; a count that passes the end of the
    // ambiguity data, where NCBI reads past its vector, is a format error.
    /// The ambiguity runs of nucleotide `i`, in stored order.
    pub fn ambiguity_runs(&self, i: usize) -> Result<Vec<AmbiguityRun>, SeqDbError> {
        let words = self.ambiguity_words(i);
        let mut runs = Vec::new();
        let Some(&first) = words.first() else {
            return Ok(runs);
        };
        let new_format = first & 0x8000_0000 != 0;
        let amb_num = (first & 0x7FFF_FFFF) as usize;
        let words_per_entry = if new_format { 2 } else { 1 };
        runs.reserve(amb_num.min(words.len()) / words_per_entry);
        let mut index = 1;
        // In the new format the loop above advances `i` by two per entry, so
        // there `amb_num` counts 32-bit words, not entries.
        while index < amb_num + 1 {
            if index + words_per_entry > words.len() {
                return Err(self.format(i, "the ambiguity data are shorter than their count"));
            }
            let word = words[index];
            let ncbi4na = ((word >> 28) & 0xF) as u8;
            let (length, offset) = if new_format {
                (((word >> 16) & 0xFFF) + 1, words[index + 1])
            } else {
                (((word >> 24) & 0xF) + 1, word & 0x00FF_FFFF)
            };
            runs.push(AmbiguityRun {
                offset,
                length,
                ncbi4na,
            });
            index += words_per_entry;
        }
        Ok(runs)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1619-1707
    // ```c++
    //     int base_length = x_GetSequence(oid, &tmp);
    //     ...
    //     if (base_length < 1)
    //         NCBI_THROW(CSeqDBException, eFileErr, "Error: could not get sequence or range.");
    //     ...
    //         bool sentinel = (nucl_code == kSeqDBNuclBlastNA8);
    //         *buffer = x_AllocType(base_length + (sentinel ? 2 : 0), alloc_type);
    //         char *seq = *buffer - range.begin + (sentinel ? 1 : 0);
    //         ...
    //             s_SeqDBMapNA2ToNA8(tmp, seq, range);
    //             s_SeqDBRebuildDNA_NA8(seq, ambchars, range);
    //             s_SeqDBMaskSequence(seq, masks, (char)14, range);
    //             if (sentinel) s_SeqDBMapNcbiNA8ToBlastNA8(seq, range);
    //         ...
    //         if (sentinel) {
    //             (*buffer)[0] = (char)15;
    //             (*buffer)[base_length+1] = (char)15;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:834-842
    // ```c++
    //         if (position + row_len <= region.begin)
    //             continue;
    //
    //         if(position >= region.end)
    //         	break;
    //
    //         for (int j = 0; j < row_len; ++j, ++position)
    //             if ( position >= region.begin  && position < region.end)
    //                 seq[position] = trans_ch;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:619-625
    // ```c++
    // static void
    // s_SeqDBMapNcbiNA8ToBlastNA8(char              * buf,
    //                             const SSeqDBSlice & range)
    // {
    //     for(int i = range.begin; i < range.end; i++)
    //         buf[i] = SeqDB_ncbina8_to_blastna8[ buf[i] & 0xF ];
    // }
    // ```
    // The whole sequence is the range; no mask is applied (LOSAT reads no
    // mask data). `s_SeqDBMapNA2ToNA8` (seqdbvol.cpp:455-471) turns each
    // 2-bit code `p` into the NCBI8na value `1 << p`.
    /// Nucleotide `i` in BLASTNA, one byte per base with the ambiguity codes
    /// restored, and the sentinel byte 15 before and after the bases (the
    /// buffer NCBI's traceback gets): `len + 2` bytes.
    pub fn nucleotide_blastna(&self, i: usize) -> Result<Vec<u8>, SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        let mut buffer = vec![0u8; base_length + 2];
        {
            let seq = &mut buffer[1..=base_length];
            for (position, base) in seq.iter_mut().enumerate() {
                let code = (packed[position >> 2] >> (6 - 2 * (position & 3))) & 3;
                *base = 1 << code;
            }
            let region_end = base_length as i64;
            for run in self.ambiguity_runs(i)? {
                let mut position = i64::from(run.offset as i32);
                let row_len = i64::from(run.length);
                if position + row_len <= 0 {
                    continue;
                }
                if position >= region_end {
                    break;
                }
                for _ in 0..row_len {
                    if position >= 0 && position < region_end {
                        seq[position as usize] = run.ncbi4na;
                    }
                    position += 1;
                }
            }
            for base in seq.iter_mut() {
                *base = NCBINA8_TO_BLASTNA8[usize::from(*base & 0xF)];
            }
        }
        buffer[0] = BLASTNA_SENTINEL;
        buffer[base_length + 1] = BLASTNA_SENTINEL;
        Ok(buffer)
    }
}

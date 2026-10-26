//! One database volume: its index, sequence and header files, and blocks of
//! sequences read by position.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use super::file::DbFile;
use super::header::{decode_defline_taxids, decode_deflines_as, decode_full, Defline};
use super::index::IndexFile;
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
        let mode = super::x_dbidxread();
        if mode != XDbIdxRead::Off {
            crate::utils::xreach::x_reach!("LOSAT_X_DBIDXREAD", true);
        }
        Self::open_as(base, seqtype, mode)
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
            XDbIdxRead::On => IndexFile::open(index_file, seqtype)?,
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

    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the number of residues of protein `oid` (volume-relative), read
    /// from the index alone: `read_block(oid, oid + 1, _)`'s `protein(0).len()`.
    pub fn protein_length(&self, oid: u32) -> Result<u32, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein length of a nucleotide volume"
        );
        self.check_range(oid, oid.saturating_add(1))?;
        // One read of the two offsets, as `read_block` reads them; the residues end before the
        // sentinel byte at the start of the next sequence's offset.
        let offsets = self.index.sequence_offsets(oid, oid + 1)?;
        offsets[1]
            .checked_sub(offsets[0])
            .and_then(|bytes| bytes.checked_sub(1))
            .ok_or_else(|| self.index.sequence_changed(oid, 2))
    }

    /// Bytes of the sequence file that a block of the sequences `first..end`
    /// occupies (protein: with the sentinel bytes around them; nucleotide:
    /// with their ambiguity data).
    pub fn block_bytes(&self, first: u32, end: u32) -> Result<u64, SeqDbError> {
        let (start, stop) = (
            self.index.sequence_start(first)?,
            self.index.sequence_start(end)?,
        );
        // Ascending when the arrays are held (checked at open); read from the file, the two
        // entries can disagree only if the file changed after the open.
        let bytes = stop.checked_sub(start).ok_or_else(|| {
            self.index
                .sequence_changed(first, (end - first) as usize + 1)
        })?;
        Ok(bytes + self.lead())
    }

    /// The end of the block that starts at `first` (`SeqDb::next_block`): the
    /// last end in `first + 1 ..= num_oids` whose `block_bytes` fit
    /// `max_bytes`, at least `first + 1`. The offsets ascend, so the ends that
    /// fit form a prefix.
    pub(super) fn block_end(&self, first: u32, max_bytes: u64) -> Result<u32, SeqDbError> {
        if let Some(offsets) = self.index.held_sequence_offsets() {
            return Ok(held_block_end(offsets, first, max_bytes, self.lead()));
        }
        self.x_read_block_end(first, max_bytes, X_BLOCK_END_WINDOW)
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
            // The end `n` does not fit, so the scan stops at or before it, unless the file
            // changed after the open and the entry of `n` reads differently.
            if next > n {
                return Err(self.index.sequence_changed(first, (n - first) as usize + 1));
            }
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
        self.read_block_in(first, end, first_ordinal, BlockBuffer::New)
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_block` with the bytes in a buffer of
    /// `File::x_read_range_mapped`.
    pub fn x_read_block_mapped(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_block_in(first, end, first_ordinal, BlockBuffer::Mapped)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): `read_block` with the bytes in `buf` (`File::x_read_range_into`;
    /// `SequenceBlock::x_into_data` gives the buffer back).
    pub fn x_read_block_into(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        buf: Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_block_in(first, end, first_ordinal, BlockBuffer::Given(buf))
    }

    fn read_block_in(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        buffer: BlockBuffer,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.check_range(first, end)?;
        let seqtype = self.index.seqtype;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        // Each sequence's ambiguity data inside its data (checked at open): read from the file,
        // the entries can break it only if the file changed after the open.
        if !self.index.holds_offsets()
            && !amb_offsets
                .iter()
                .enumerate()
                .all(|(i, &amb)| seq_offsets[i] < amb && amb <= seq_offsets[i + 1])
        {
            return Err(self.index.ambiguity_changed(first, amb_offsets.len()));
        }
        let data_start = u64::from(seq_offsets[0]) - self.lead();
        let data_end = u64::from(seq_offsets[seq_offsets.len() - 1]);
        let file = self.seq_file()?;
        let data = match buffer {
            BlockBuffer::New => file.read_range(data_start, data_end)?,
            BlockBuffer::Mapped => file.x_read_range_mapped(data_start, data_end)?,
            BlockBuffer::Given(buf) => file.x_read_range_into(data_start, data_end, buf)?,
        };
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

    /// EXPERIMENT (LOSAT_X_TBNDBSTREAM): what `read_block(first, end, first_ordinal)` of one
    /// nucleotide (`end == first + 1`) gives `nucleotide_len(0)` and `ambiguity_runs(0)`, without
    /// its 2-bit data: the index reads and checks of `read_block` in its order and the checks of
    /// the file range it reads (`DbFile::x_check_range`), then the last byte of the 2-bit data and
    /// the ambiguity data only. Returns the length, with `nucleotide_len`'s error, and a block of
    /// the ambiguity data whose `ambiguity_runs(0)` is the one of `read_block`'s block (as in
    /// `check_ambiguity_data`). An I/O error in the bytes between is not seen here; the search
    /// reads them later.
    pub fn x_nucleotide_check(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
    ) -> Result<(usize, SequenceBlock), SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Nucleotide && end == first.wrapping_add(1),
            "the check of one nucleotide"
        );
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        let (start, amb_start, next) = (seq_offsets[0], amb_offsets[0], seq_offsets[1]);
        if !self.index.holds_offsets() && !(start < amb_start && amb_start <= next) {
            return Err(self.index.ambiguity_changed(first, amb_offsets.len()));
        }
        let file = self.seq_file()?;
        file.x_check_range(u64::from(start), u64::from(next))?;
        let last = file.read_range(u64::from(amb_start) - 1, u64::from(amb_start))?;
        let length = u64::from(amb_start - start - 1) * 4 + u64::from(last[0] & 3);
        let ambiguity = SequenceBlock {
            seqtype: SeqType::Nucleotide,
            first_ordinal,
            path: file.path().to_path_buf(),
            data: file.read_range(u64::from(amb_start), u64::from(next))?,
            data_start: u64::from(amb_start),
            // The block holds the ambiguity data only: its letters' range is
            // empty, and nothing but `ambiguity_runs` reads it.
            seq_offsets: vec![amb_start, next],
            amb_offsets: vec![amb_start],
        };
        if length > i32::MAX as u64 {
            return Err(ambiguity.format(0, "the sequence is longer than a signed 32-bit length"));
        }
        Ok((length as usize, ambiguity))
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the numbers of residues of the proteins `first..end`
    /// (volume-relative, non-empty), read from the index alone with one read of their offsets:
    /// `read_block(first, end, _)`'s `protein(i).len()`.
    pub fn x_protein_lengths(&self, first: u32, end: u32) -> Result<Vec<u32>, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein lengths of a nucleotide volume"
        );
        self.check_range(first, end)?;
        let offsets = self.index.sequence_offsets(first, end)?;
        // The residues end before the sentinel byte at the start of the next sequence's offset;
        // read from the file, the entries can break that only if the file changed after the open.
        offsets
            .windows(2)
            .map(|pair| {
                pair[1]
                    .checked_sub(pair[0])
                    .and_then(|bytes| bytes.checked_sub(1))
            })
            .collect::<Option<Vec<u32>>>()
            .ok_or_else(|| self.index.sequence_changed(first, offsets.len()))
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

    /// The `taxid` of each definition line of the sequence `oid` (volume-relative), 0 for none.
    pub fn defline_taxids(&self, oid: u32) -> Result<Vec<i32>, SeqDbError> {
        let record = self.header_bytes(oid)?;
        decode_defline_taxids(
            &record,
            &with_extension(&self.base, self.index.seqtype, "hr"),
        )
    }

    /// The definition lines of the sequence `oid` (volume-relative) with every member and
    /// every Seq-id (`header::decode_full`); a `general` id `BL_ORD_ID` keeps its
    /// volume-relative tag.
    pub fn defline_set(&self, oid: u32) -> Result<crate::seqid::DeflineSet, SeqDbError> {
        let record = self.header_bytes(oid)?;
        decode_full(
            &record,
            &with_extension(&self.base, self.index.seqtype, "hr"),
        )
    }

    /// The definition lines of the sequence `oid` (volume-relative), with
    /// their ordinal ids still relative to the volume.
    pub fn deflines(&self, oid: u32) -> Result<Vec<Defline>, SeqDbError> {
        self.deflines_as(oid, false)
    }

    /// `deflines`, with the `taxid` of each definition line decoded when `want_taxid` (the
    /// decode of `defline_taxids`, which also gives the titles).
    pub fn deflines_as(&self, oid: u32, want_taxid: bool) -> Result<Vec<Defline>, SeqDbError> {
        let record = self.header_bytes(oid)?;
        decode_deflines_as(
            &record,
            &with_extension(&self.base, self.index.seqtype, "hr"),
            want_taxid,
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1838
    // ```c++
    //         int whole_bytes = int(end_offset - start_offset - 1);
    //
    //         char last_char = (*buffer)[whole_bytes];
    //
    //         int remainder = last_char & 3;
    //         length = (whole_bytes * 4) + remainder;
    // ```
    /// The packed bytes of nucleotide `oid` (volume-relative) from `skip` bytes in, `count`
    /// of them, with the offsets of its sequence and its ambiguity data: `(bytes, packed
    /// start, ambiguity start, next sequence's start)`.
    fn packed_part(
        &self,
        oid: u32,
        ordinal: u32,
        skip: u64,
        count: Option<u64>,
    ) -> Result<(Vec<u8>, u32, u32, u32), SeqDbError> {
        let (start, amb_start, next) = self.part_offsets(oid, ordinal)?;
        let file = self.seq_file()?;
        let packed = u64::from(amb_start - start);
        let from = u64::from(start) + skip.min(packed);
        let to = match count {
            Some(count) => (from + count).min(u64::from(amb_start)),
            None => u64::from(amb_start),
        };
        Ok((file.read_range(from, to)?, start, amb_start, next))
    }

    /// The offsets in the sequence file of nucleotide `oid` (volume-relative): `(start of its
    /// packed bases, start of its ambiguity data, start of the next sequence)`.
    fn part_offsets(&self, oid: u32, ordinal: u32) -> Result<(u32, u32, u32), SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let seq = self.index.sequence_offsets(oid, oid + 1)?;
        let amb = self.index.ambiguity_offsets(oid, oid + 1)?;
        let (start, amb_start, next) = (seq[0], amb[0], seq[1]);
        if !(start < amb_start && amb_start <= next) {
            return Err(SeqDbError::Format {
                path: self.seq_file()?.path().to_path_buf(),
                message: format!("sequence {ordinal}: the ambiguity offset lies outside its data"),
            });
        }
        Ok((start, amb_start, next))
    }

    /// The length in bases of nucleotide `oid` (volume-relative; `ordinal` its database
    /// ordinal id, for the message), from the last byte of its packed bases only.
    pub fn nucleotide_length(&self, oid: u32, ordinal: u32) -> Result<usize, SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let amb = self.index.ambiguity_offsets(oid, oid + 1)?;
        let seq = self.index.sequence_offsets(oid, oid + 1)?;
        let (start, amb_start) = (u64::from(seq[0]), u64::from(amb[0]));
        let file = self.seq_file()?;
        if amb_start <= start {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!("sequence {ordinal}: the ambiguity offset lies outside its data"),
            });
        }
        let last = file.read_range(amb_start - 1, amb_start)?;
        let length = (amb_start - start - 1) * 4 + u64::from(last[0] & 3);
        if length > i32::MAX as u64 {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!(
                    "sequence {ordinal}: the sequence is longer than a signed 32-bit length"
                ),
            });
        }
        Ok(length as usize)
    }

    /// The IUPAC letters `[from, to)` of nucleotide `oid` (volume-relative), as `ordinal`'s
    /// stored data give them: the 2-bit bases of the range and the ambiguity runs that overlap
    /// it. Only the bytes of the range and the sequence's ambiguity data are read.
    pub fn nucleotide_iupac_range(
        &self,
        oid: u32,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        let ambiguity = self.nucleotide_ambiguity(oid, ordinal)?;
        self.nucleotide_iupac_range_with(oid, ordinal, from, to, &ambiguity)
    }

    /// The length of nucleotide `oid` (volume-relative; `ordinal` its database ordinal id) and
    /// its decoded ambiguity runs, read once for any number of
    /// [`Volume::nucleotide_iupac_range_with`] calls on the same sequence.
    pub fn nucleotide_ambiguity(
        &self,
        oid: u32,
        ordinal: u32,
    ) -> Result<NucleotideAmbiguity, SeqDbError> {
        let length = self.nucleotide_length(oid, ordinal)?;
        let (packed_start, amb_start, next) = self.part_offsets(oid, ordinal)?;
        let mut runs = Vec::new();
        if u64::from(amb_start) < u64::from(next) {
            let file = self.seq_file()?;
            let words: Vec<u32> = file
                .read_range(u64::from(amb_start), u64::from(next))?
                .chunks_exact(4)
                .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
                .collect();
            runs = decode_ambiguity_runs(&words).map_err(|message| SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!("sequence {ordinal}: {message}"),
            })?;
        }
        Ok(NucleotideAmbiguity::new(
            length,
            u64::from(packed_start),
            runs,
        ))
    }

    /// [`Volume::nucleotide_iupac_range`] with the length and the ambiguity runs that
    /// [`Volume::nucleotide_ambiguity`] read for the same `oid`: one read of the range's packed
    /// bytes, and a binary search for the first run that can overlap the range.
    pub fn nucleotide_iupac_range_with(
        &self,
        oid: u32,
        ordinal: u32,
        from: usize,
        to: usize,
        ambiguity: &NucleotideAmbiguity,
    ) -> Result<Vec<u8>, SeqDbError> {
        let length = ambiguity.length;
        if from > to || to > length {
            return Err(SeqDbError::Range(format!(
                "letters {from}..{to} are outside sequence {ordinal}, which has {length}"
            )));
        }
        self.check_range(oid, oid.saturating_add(1))?;
        let first_byte = (from / 4) as u64;
        let end_byte = (to as u64).div_ceil(4);
        let packed = self.seq_file()?.read_range(
            ambiguity.packed_start + first_byte,
            ambiguity.packed_start + end_byte,
        )?;
        let mut letters: Vec<u8> = (from..to)
            .map(|at| {
                let byte = packed[at / 4 - first_byte as usize];
                b"ACGT"[usize::from((byte >> (6 - 2 * (at % 4))) & 3)]
            })
            .collect();
        ambiguity.overlay(from, to, &mut letters);
        Ok(letters)
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
        self.visit_deflines_as(first, end, piece_bytes, false, visit)
    }

    /// `visit_deflines`, with the `taxid` of each definition line decoded when `want_taxid`
    /// (the results of `deflines_as`): one pass over the headers gives the titles and the taxids.
    pub(super) fn visit_deflines_as(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        want_taxid: bool,
        visit: &mut dyn FnMut(u32, Result<Vec<Defline>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_headers(
            first,
            end,
            piece_bytes,
            &|oid| self.deflines_as(oid, want_taxid),
            &|record, path| decode_deflines_as(record, path, want_taxid),
            visit,
        )
    }

    /// `visit_deflines` for `defline_taxids` (OF5b: the taxdb pre-pass of BLASTX 15/16 reads the
    /// taxid of every sequence's first definition line): the same pieces and fallbacks, with
    /// the decoder of `defline_taxids`.
    pub(super) fn visit_defline_taxids(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        visit: &mut dyn FnMut(u32, Result<Vec<i32>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_headers(
            first,
            end,
            piece_bytes,
            &|oid| self.defline_taxids(oid),
            &decode_defline_taxids,
            visit,
        )
    }

    /// The shared walk of `visit_deflines` and `visit_defline_taxids`: `single(oid)` is the
    /// one-sequence read that gives the result for a sequence that is not read in a piece,
    /// `decode(record, path)` decodes a record read in a piece.
    fn visit_headers<T>(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        single: &dyn Fn(u32) -> Result<T, SeqDbError>,
        decode: &dyn Fn(&[u8], &Path) -> Result<T, SeqDbError>,
        visit: &mut dyn FnMut(u32, Result<T, SeqDbError>) -> bool,
    ) -> bool {
        let one_by_one =
            |oids: std::ops::Range<u32>,
             visit: &mut dyn FnMut(u32, Result<T, SeqDbError>) -> bool| {
                oids.into_iter().all(|oid| visit(oid, single(oid)))
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
                if !visit(first + j as u32, decode(&data[lo..hi], &hr_path)) {
                    return false;
                }
            }
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

/// Where `Volume::read_block_in` puts the bytes of a block.
enum BlockBuffer {
    /// A new buffer (`read_block`).
    New,
    /// EXPERIMENT (LOSAT_X_BNMAPPED): `File::x_read_range_mapped`'s buffer.
    Mapped,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP): the caller's buffer.
    Given(Vec<u8>),
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

/// The length of one nucleotide sequence and its ambiguity runs, read once by
/// [`Volume::nucleotide_ambiguity`] for any number of range reads
/// ([`Volume::nucleotide_iupac_range_with`]).
#[derive(Debug, Clone)]
pub struct NucleotideAmbiguity {
    length: usize,
    /// Offset in the sequence file of the first packed byte.
    packed_start: u64,
    /// The runs in stored order: a later run overwrites an earlier one.
    runs: Vec<AmbiguityRun>,
    /// Whether the runs are in increasing order of offset and do not overlap, as NCBI writes
    /// them: then the runs that overlap a range are found by binary search. Otherwise every run
    /// is applied in stored order.
    ordered: bool,
}

impl NucleotideAmbiguity {
    fn new(length: usize, packed_start: u64, runs: Vec<AmbiguityRun>) -> Self {
        let ordered = runs.windows(2).all(|pair| {
            u64::from(pair[0].offset) + u64::from(pair[0].length) <= u64::from(pair[1].offset)
        });
        Self {
            length,
            packed_start,
            runs,
            ordered,
        }
    }

    /// The length of the sequence in bases.
    pub fn length(&self) -> usize {
        self.length
    }

    /// Writes the IUPAC letters of the runs that overlap `[from, to)` into `letters`, which
    /// holds the bases of that range.
    fn overlay(&self, from: usize, to: usize, letters: &mut [u8]) {
        overlay_ambiguity_runs(&self.runs, self.ordered, from, to, letters);
    }
}

/// Writes the IUPAC letters of `runs` over `letters` (the bases `[from, to)`), a later run over
/// an earlier one. With `ordered` runs (increasing offsets, no overlap) only those from the
/// first run that ends after `from` to the first that starts at or after `to` are visited.
fn overlay_ambiguity_runs(
    runs: &[AmbiguityRun],
    ordered: bool,
    from: usize,
    to: usize,
    letters: &mut [u8],
) {
    let first = if ordered {
        runs.partition_point(|run| (run.offset as usize + run.length as usize) <= from)
    } else {
        0
    };
    for run in &runs[first..] {
        if ordered && run.offset as usize >= to {
            break;
        }
        let run_from = (run.offset as usize).max(from);
        let run_to = (run.offset as usize + run.length as usize).min(to);
        for at in run_from..run_to.max(run_from) {
            letters[at - from] = NCBI4NA_TO_IUPAC[usize::from(run.ncbi4na & 0xF)];
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNTBRANGE): the ambiguity codes of `runs` over the bases of `ranges`
/// (sorted, disjoint) in `seq` (the whole sequence's bases): the runs in stored order, with the
/// whole decode's tests (`nucleotide_blastna_in`: a run that ends before the sequence is skipped,
/// the first one that starts at or after its end stops the loop), each run written over the
/// ranges it overlaps. So a range's bytes are the whole decode's for ordered and unordered runs.
fn x_overlay_runs_in_ranges(runs: &[AmbiguityRun], ranges: &[(usize, usize)], seq: &mut [u8]) {
    let region_end = seq.len() as i64;
    for run in runs {
        let position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        let run_from = position.max(0) as usize;
        let run_to = (position + row_len).min(region_end) as usize;
        let first = ranges.partition_point(|r| r.1 <= run_from);
        for &(from, to) in &ranges[first..] {
            if from >= run_to {
                break;
            }
            for at in run_from.max(from)..run_to.min(to) {
                seq[at] = run.ncbi4na;
            }
        }
    }
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

/// `NCBI4NA_TO_IUPACNA` (c++/src/algo/blast/core/blast_encoding.c:80-83).
const NCBI4NA_TO_IUPAC: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

/// Sentinel byte NCBI puts before and after a BLASTNA sequence.
pub const BLASTNA_SENTINEL: u8 = 15;

/// EXPERIMENT (LOSAT_X_BNTBRANGE): NCBI's `FENCE_SENTRY` (seqdbgeneral.hpp:1229, blast_util.h:364),
/// the byte around each range of a partly decoded sequence (`x_nucleotide_blastna_ranges`).
pub const FENCE_SENTRY: u8 = 201;

/// The ambiguity runs of the ambiguity words of one sequence, in stored order
/// (`SequenceBlock::ambiguity_runs`); the error is its message.
fn decode_ambiguity_runs(words: &[u32]) -> Result<Vec<AmbiguityRun>, &'static str> {
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
            return Err("the ambiguity data are shorter than their count");
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

impl SequenceBlock {
    /// EXPERIMENT (LOSAT_X_BNMAPPED): gives the block's bytes to `utils::x_mapped::X_BYTES` (a
    /// buffer of `File::x_read_range_mapped` of 1 MiB or more; smaller ones are dropped). No NCBI
    /// counterpart: only where the bytes go.
    pub fn x_recycle(mut self) {
        crate::utils::x_mapped::X_BYTES.give(&mut self.data);
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): the block's buffer, for the next `x_read_block_into`. No
    /// NCBI counterpart: only where the bytes are kept.
    pub fn x_into_data(self) -> Vec<u8> {
        self.data
    }

    /// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the bytes of a block of one nucleotide cut to its 2-bit
    /// data (`nucleotide_packed(0)`, which starts the block), in the block's buffer with its
    /// capacity. No NCBI counterpart: only where the bytes are kept.
    pub fn x_into_nucleotide_packed(mut self) -> Vec<u8> {
        self.expect(SeqType::Nucleotide);
        assert_eq!(self.len(), 1, "the packed bytes of a block of one sequence");
        assert_eq!(
            u64::from(self.seq_offsets[0]),
            self.data_start,
            "a nucleotide block starts with its first sequence"
        );
        let end = self.nucleotide_packed(0).len();
        self.data.truncate(end);
        self.data
    }

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
        decode_ambiguity_runs(&self.ambiguity_words(i)).map_err(|message| self.format(i, message))
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
        self.nucleotide_blastna_in(i, false)
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `nucleotide_blastna` in a buffer taken from or sized
    /// by `utils::x_mapped::X_BYTES.reserve` (a mapped buffer from 1 MiB on, given back with
    /// `x_mapped::X_BYTES.give`). No NCBI counterpart: only where the bytes are kept.
    pub fn x_nucleotide_blastna_mapped(&self, i: usize) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_blastna_in(i, true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1565-1595
    // ```c++
    //     *buffer = x_AllocType(base_length + (sentinel ? 2 : 0), alloc_type);
    //     char *seq = *buffer + (sentinel ? 1 : 0);
    //
    //     vector<Int4> ambchars;
    // 	x_GetAmbChar(oid, ambchars);
    //     ITERATE(CSeqDB::TSequenceRanges, riter, *(partial_ranges)) {
    //     	int begin(riter->first);
    //         int end(riter->second);
    //
    // 	    if (begin) seq[begin - 1] = (char) FENCE_SENTRY;
    // 	    if (end < base_length) seq[end] = (char) FENCE_SENTRY;
    // 	}
    //     ITERATE(CSeqDB::TSequenceRanges, riter, *(partial_ranges)) {
    //         SSeqDBSlice slice(max(0, (int)riter->first), min(base_length, (int)riter->second));
    // 	    s_SeqDBMapNA2ToNA8(tmp, seq, slice);
    // 	    s_SeqDBRebuildDNA_NA8(seq, ambchars, slice);
    // 	    s_SeqDBMaskSequence(seq, masks, (char)14, slice);
    // 	    if (sentinel){
    // 	    	s_SeqDBMapNcbiNA8ToBlastNA8(seq, slice);
    // 	    }
    // 	}
    // ```
    // The ambiguity runs are applied in stored order, each to the ranges it overlaps, with the
    // whole decode's tests (`nucleotide_blastna_in`), so every byte of a range is the whole
    // decode's (NCBI's slices stop at the first run past each slice).
    /// EXPERIMENT (LOSAT_X_BNTBRANGE): `nucleotide_blastna` of the bases of `ranges` only
    /// (`[from, to)`, sorted and disjoint, as NCBI's partial-fetch ranges are; a range past the
    /// end is cut there), with `FENCE_SENTRY` before and after
    /// each range that does not reach an end of the sequence (`GetAmbigPartialSeq`): `len + 2`
    /// bytes, zero outside the ranges. The buffer is zeroed by the allocator; from
    /// `x_mapped::X_MAPPED_FROM` on its capacity is above glibc's mmap ceiling, so glibc maps it
    /// (unless an arena has a free chunk that large) and writes none of it: the pages outside the
    /// ranges are not resident. Other allocators may zero the whole buffer (Wasm's dlmalloc does).
    pub fn x_nucleotide_blastna_ranges(
        &self,
        i: usize,
        ranges: &[(usize, usize)],
    ) -> Result<Vec<u8>, SeqDbError> {
        debug_assert!(ranges.windows(2).all(|pair| pair[0].1 <= pair[1].0));
        let packed = self.nucleotide_packed(i);
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        let mut buffer = vec![0u8; crate::utils::x_mapped::x_mapped_capacity(base_length + 2, 1)];
        buffer.truncate(base_length + 2);
        {
            let seq = &mut buffer[1..=base_length];
            let ranges = ranges
                .iter()
                .map(|&(from, to)| (from.min(base_length), to.min(base_length)));
            for (from, to) in ranges.clone() {
                if from > 0 {
                    seq[from - 1] = FENCE_SENTRY;
                }
                if to < base_length {
                    seq[to] = FENCE_SENTRY;
                }
            }
            for (from, to) in ranges.clone() {
                for (position, base) in seq[from..to.max(from)].iter_mut().enumerate() {
                    let position = from + position;
                    let code = (packed[position >> 2] >> (6 - 2 * (position & 3))) & 3;
                    *base = 1 << code;
                }
            }
            let ranges: Vec<(usize, usize)> = ranges.filter(|r| r.0 < r.1).collect();
            x_overlay_runs_in_ranges(&self.ambiguity_runs(i)?, &ranges, seq);
            for &(from, to) in &ranges {
                for base in &mut seq[from..to] {
                    *base = NCBINA8_TO_BLASTNA8[usize::from(*base & 0xF)];
                }
            }
        }
        buffer[0] = BLASTNA_SENTINEL;
        buffer[base_length + 1] = BLASTNA_SENTINEL;
        Ok(buffer)
    }

    fn nucleotide_blastna_in(&self, i: usize, mapped: bool) -> Result<Vec<u8>, SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        let mut buffer = if mapped {
            let mut buffer = Vec::new();
            crate::utils::x_mapped::X_BYTES.reserve(&mut buffer, base_length + 2);
            buffer.resize(base_length + 2, 0u8);
            buffer
        } else {
            vec![0u8; base_length + 2]
        };
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

#[cfg(test)]
mod tests {
    use super::{overlay_ambiguity_runs, AmbiguityRun, NucleotideAmbiguity};

    /// A small deterministic generator (xorshift), so the test needs no fixture and no crate.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self, below: u64) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 % below
        }
    }

    /// The letters of the whole sequence with every run applied in stored order (the decode
    /// that `SequenceBlock` does), the reference for any range.
    fn whole_decode(length: usize, runs: &[AmbiguityRun]) -> Vec<u8> {
        let mut letters = vec![b'A'; length];
        overlay_ambiguity_runs(runs, false, 0, length, &mut letters);
        letters
    }

    fn check(length: usize, runs: Vec<AmbiguityRun>, rng: &mut Rng, ranges: usize) -> bool {
        let whole = whole_decode(length, &runs);
        let ambiguity = NucleotideAmbiguity::new(length, 0, runs);
        for _ in 0..ranges {
            let from = rng.next(length as u64 + 1) as usize;
            let to = from + rng.next((length - from) as u64 + 1) as usize;
            let mut letters = vec![b'A'; to - from];
            ambiguity.overlay(from, to, &mut letters);
            assert_eq!(letters, whole[from..to], "range {from}..{to}");
        }
        ambiguity.ordered
    }

    /// The ambiguity loop of `nucleotide_blastna_in` over the whole sequence.
    fn whole_runs(runs: &[AmbiguityRun], seq: &mut [u8]) {
        let region_end = seq.len() as i64;
        for run in runs {
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
    }

    // LOSAT_X_BNTBRANGE: the runs over ranges give the whole decode's bytes inside the ranges,
    // for ordered runs, shuffled and overlapping ones, offsets of 2^31 and more (negative as
    // NCBI's Int4), and runs past the end (which stop the whole decode's loop).
    #[test]
    fn runs_over_ranges_equal_the_whole_decode_inside_them() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        for case in 0..400 {
            let length = 1 + rng.next(3000) as usize;
            let mut runs = Vec::new();
            for _ in 0..rng.next(40) {
                let offset = match rng.next(20) {
                    0 => 0x8000_0000u32 + rng.next(1 << 20) as u32,
                    1 => length as u32 + rng.next(100) as u32,
                    2 => u32::MAX - rng.next(50) as u32,
                    _ => rng.next(length as u64) as u32,
                };
                let longest = if rng.next(10) == 0 { 5000 } else { 20 };
                runs.push(AmbiguityRun {
                    offset,
                    length: 1 + rng.next(longest) as u32,
                    ncbi4na: rng.next(16) as u8,
                });
            }
            if case % 2 == 0 {
                runs.sort_by_key(|run| run.offset);
            }
            let mut ranges = Vec::new();
            let mut at = 0usize;
            while at < length {
                let from = at + rng.next(200) as usize;
                let to = (from + rng.next(300) as usize).min(length);
                if from < to {
                    ranges.push((from, to));
                }
                at = to + 1 + rng.next(100) as usize;
            }
            let mut whole = vec![1u8; length];
            whole_runs(&runs, &mut whole);
            let mut ranged = vec![1u8; length];
            super::x_overlay_runs_in_ranges(&runs, &ranges, &mut ranged);
            for &(from, to) in &ranges {
                assert_eq!(
                    ranged[from..to],
                    whole[from..to],
                    "case {case}, range {from}..{to}"
                );
            }
        }
    }

    #[test]
    fn many_runs_and_many_ranges_equal_the_whole_decode() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        // 3000 ordered runs of 1..=20 bases with gaps of 0..=30 bases (NCBI's layout).
        let mut runs = Vec::new();
        let mut at = 5u32;
        for _ in 0..3000 {
            let length = 1 + rng.next(20) as u32;
            runs.push(AmbiguityRun {
                offset: at,
                length,
                ncbi4na: 2 + rng.next(14) as u8,
            });
            at += length + rng.next(31) as u32;
        }
        let length = at as usize + 10;
        assert!(check(length, runs.clone(), &mut rng, 5000));
        // Ranges that start inside a run, end inside one, or lie between two.
        for run in runs.iter().step_by(97) {
            let mut one = vec![b'A'; 2];
            let ambiguity = NucleotideAmbiguity::new(length, 0, runs.clone());
            ambiguity.overlay(run.offset as usize, run.offset as usize + 1, &mut one[..1]);
            assert_ne!(one[0], b'A');
        }
        // Stored runs that overlap or are out of order: every run applies in stored order, so
        // a later run still overwrites an earlier one.
        let mut messy = runs.clone();
        messy.swap(10, 2900);
        messy.push(AmbiguityRun {
            offset: 7,
            length: 5000,
            ncbi4na: 3,
        });
        messy.push(AmbiguityRun {
            offset: 100,
            length: 3,
            ncbi4na: 9,
        });
        assert!(!check(length, messy, &mut rng, 5000));
    }

    #[test]
    fn a_range_outside_every_run_stays_unchanged() {
        let runs = vec![AmbiguityRun {
            offset: 10,
            length: 4,
            ncbi4na: 15,
        }];
        let ambiguity = NucleotideAmbiguity::new(40, 0, runs);
        for (from, to) in [(0, 10), (14, 40), (0, 0), (40, 40)] {
            let mut letters = vec![b'C'; to - from];
            ambiguity.overlay(from, to, &mut letters);
            assert!(letters.iter().all(|&letter| letter == b'C'), "{from}..{to}");
        }
        let mut letters = vec![b'C'; 6];
        ambiguity.overlay(8, 14, &mut letters);
        assert_eq!(letters, b"CCNNNN");
    }
}

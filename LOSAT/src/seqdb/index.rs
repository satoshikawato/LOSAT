//! Index file of a database volume (`.pin` / `.nin`), format versions 4 and 5.
//!
//! EXPERIMENT (LOSAT_X_DBIDXREAD, in the strict set): the offset arrays stay in the file and
//! are read by position when a range of them is needed (`IndexFile::open`); without the switch
//! the whole file is read and the arrays are held (`IndexFile::parse`). Both forms check the
//! arrays when the volume is opened and give the same values and the same errors.

use std::path::{Path, PathBuf};

use super::file::{mapped_length, DbFile};
use super::{SeqDbError, SeqType};

/// Parsed index file: the volume's summary values and its offset arrays.
#[derive(Debug)]
pub struct IndexFile {
    /// Format version, 4 or 5.
    #[cfg(test)]
    pub version: u32,
    /// Sequence type stored in the file.
    pub seqtype: SeqType,
    /// Volume number (version 5 only).
    #[cfg(test)]
    pub volume_number: Option<u32>,
    /// Volume title, as stored.
    pub title: Vec<u8>,
    /// Name of the LMDB file (version 5 only). LOSAT does not read that file.
    pub lmdb_file: Option<Vec<u8>>,
    /// Creation date, as stored (may end in NUL padding).
    pub date: Vec<u8>,
    /// Number of sequences (ordinal ids) in the volume.
    pub num_oids: u32,
    /// Total number of residues or bases in the volume.
    pub volume_length: u64,
    /// Length of the longest sequence.
    pub max_length: u32,
    offsets: Offsets,
}

/// The three offset arrays (`num_oids + 1` entries each; the ambiguity array is empty for
/// protein).
#[derive(Debug)]
struct HeldOffsets {
    hdr: Vec<u32>,
    seq: Vec<u32>,
    amb: Vec<u32>,
}

/// EXPERIMENT (LOSAT_X_DBIDXREAD): where the offset arrays are read from.
#[derive(Debug)]
enum Offsets {
    /// Read with the file and held (`IndexFile::parse`).
    Held(HeldOffsets),
    /// Read from the open index file by position (`IndexFile::open`).
    Read(ReadOffsets),
}

/// EXPERIMENT (LOSAT_X_DBIDXREAD): the index file and the positions of its arrays.
#[derive(Debug)]
struct ReadOffsets {
    file: DbFile,
    hdr_at: u64,
    seq_at: u64,
    amb_at: u64,
}

/// One of the three offset arrays.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Region {
    Header,
    Sequence,
    Ambiguity,
}

impl HeldOffsets {
    fn get(&self, region: Region) -> &[u32] {
        match region {
            Region::Header => &self.hdr,
            Region::Sequence => &self.seq,
            Region::Ambiguity => &self.amb,
        }
    }
}

/// Why a read of the index file's fixed part fails.
#[derive(Debug)]
enum ReadFault {
    /// The read passes the end of the bytes: `wanted` bytes at offset `at`, of `len`.
    PastEnd {
        path: PathBuf,
        wanted: usize,
        at: u64,
        len: u64,
    },
    /// Any other error.
    Error(SeqDbError),
}

impl From<SeqDbError> for ReadFault {
    fn from(error: SeqDbError) -> Self {
        Self::Error(error)
    }
}

impl From<ReadFault> for SeqDbError {
    fn from(fault: ReadFault) -> Self {
        match fault {
            ReadFault::PastEnd {
                path,
                wanted,
                at,
                len,
            } => SeqDbError::Format {
                path,
                message: format!(
                    "index file ends inside a field ({wanted} bytes wanted at offset {at}, \
                     file has {len})"
                ),
            },
            ReadFault::Error(error) => error,
        }
    }
}

/// Where the bytes of an index file come from while its fixed part is parsed.
trait Source {
    fn len(&self) -> u64;
    /// Fills `buf` from `offset`; the range lies inside `len()`.
    fn read_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError>;
}

impl Source for [u8] {
    fn len(&self) -> u64 {
        <[u8]>::len(self) as u64
    }

    fn read_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError> {
        let at = offset as usize;
        buf.copy_from_slice(&self[at..at + buf.len()]);
        Ok(())
    }
}

impl Source for DbFile {
    fn len(&self) -> u64 {
        DbFile::len(self)
    }

    fn read_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError> {
        self.read_exact_at(buf, offset)
    }
}

struct Cursor<'a, S: Source + ?Sized> {
    source: &'a S,
    pos: u64,
    path: &'a Path,
}

impl<'a, S: Source + ?Sized> Cursor<'a, S> {
    /// The offset of the `len` bytes at the cursor, which it moves past them.
    fn skip(&mut self, len: usize) -> Result<u64, ReadFault> {
        let end = self
            .pos
            .checked_add(len as u64)
            .filter(|&end| end <= self.source.len());
        let Some(end) = end else {
            return Err(ReadFault::PastEnd {
                path: self.path.to_path_buf(),
                wanted: len,
                at: self.pos,
                len: self.source.len(),
            });
        };
        let start = self.pos;
        self.pos = end;
        Ok(start)
    }

    fn fill(&mut self, buf: &mut [u8]) -> Result<(), ReadFault> {
        let at = self.skip(buf.len())?;
        Ok(self.source.read_at(buf, at)?)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:71-79
    fn be_u32(&mut self) -> Result<u32, ReadFault> {
        let mut bytes = [0u8; 4];
        self.fill(&mut bytes)?;
        Ok(u32::from_be_bytes(bytes))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:81-89
    // `SeqDB_GetBroken` reads the one little-endian field of the format.
    fn broken_u64(&mut self) -> Result<u64, ReadFault> {
        let mut raw = [0u8; 8];
        self.fill(&mut raw)?;
        Ok(u64::from_le_bytes(raw))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:91-104
    fn string(&mut self) -> Result<Vec<u8>, ReadFault> {
        let len = self.be_u32()? as usize;
        // The length is checked before the buffer is allocated.
        let at = self.skip(len)?;
        let mut value = vec![0u8; len];
        self.source.read_at(&mut value, at)?;
        Ok(value)
    }

    /// The offset of an array of `count` offsets at the cursor, which it moves past them.
    fn offsets(&mut self, count: usize) -> Result<u64, ReadFault> {
        self.skip(count.checked_mul(4).ok_or_else(|| {
            ReadFault::Error(SeqDbError::Format {
                path: self.path.to_path_buf(),
                message: "index file offset array size overflows".to_string(),
            })
        })?)
    }
}

/// The values of an index file before its arrays, and where the arrays start.
struct Fixed {
    #[cfg(test)]
    version: u32,
    #[cfg(test)]
    volume_number: Option<u32>,
    title: Vec<u8>,
    lmdb_file: Option<Vec<u8>>,
    date: Vec<u8>,
    num_oids: u32,
    volume_length: u64,
    max_length: u32,
    hdr_at: u64,
    seq_at: u64,
    amb_at: u64,
}

/// The index file as NCBI's reader sees it through its memory mapping: the bytes of the file,
/// then zeros up to the end of the file's last 4 KiB page (`file::mapped_length`). A read past
/// the last page is outside the mapping: not data of the file.
struct Mapped<'a>(&'a DbFile);

impl Source for Mapped<'_> {
    fn len(&self) -> u64 {
        mapped_length(self.0.len())
    }

    fn read_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError> {
        let held = self.0.len().saturating_sub(offset).min(buf.len() as u64) as usize;
        if held > 0 {
            self.0.read_exact_at(&mut buf[..held], offset)?;
        }
        buf[held..].fill(0);
        Ok(())
    }
}

/// A byte slice padded with zeros to its last 4 KiB mapped page (`file::mapped_length`), matching
/// [`Mapped`] for `IndexFile::parse` (`LOSAT_X_DBIDXREAD=0`).
struct MappedSlice<'a>(&'a [u8]);

impl Source for MappedSlice<'_> {
    fn len(&self) -> u64 {
        mapped_length(self.0.len() as u64)
    }

    fn read_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError> {
        let at = usize::try_from(offset).unwrap_or(usize::MAX);
        let held = self.0.len().saturating_sub(at).min(buf.len());
        if held > 0 {
            buf[..held].copy_from_slice(&self.0[at..at + held]);
        }
        buf[held..].fill(0);
        Ok(())
    }
}

/// What a volume's index file holds before its offset arrays, as NCBI reads it when it opens
/// the database (`CSeqDBIdxFile::CSeqDBIdxFile`).
#[derive(Debug)]
pub(crate) struct VolumeHeader {
    /// The title, as stored.
    pub title: Vec<u8>,
    /// The name of the LMDB file (version 5 only), as stored.
    pub lmdb_file: Option<Vec<u8>>,
    /// The date, as stored.
    pub date: Vec<u8>,
    /// The number of sequences.
    pub num_oids: u32,
    /// The total length of the sequences.
    pub volume_length: u64,
}

/// The fields of an index file before its offset arrays.
struct Header {
    #[cfg_attr(not(test), allow(dead_code))]
    version: u32,
    db_seqtype: u32,
    #[cfg_attr(not(test), allow(dead_code))]
    volume_number: Option<u32>,
    title: Vec<u8>,
    lmdb_file: Option<Vec<u8>>,
    date: Vec<u8>,
    num_oids: u32,
    volume_length: u64,
    max_length: u32,
}

/// The big-endian offsets in `bytes`.
fn decode_words(bytes: &[u8], out: &mut [u32]) {
    for (word, raw) in out.iter_mut().zip(bytes.chunks_exact(4)) {
        *word = u32::from_be_bytes([raw[0], raw[1], raw[2], raw[3]]);
    }
}

impl IndexFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:134-236
    // The constructor reads the values through the mapped file and checks nothing but the
    // version and the type. A file that ends before its header does (a write that failed)
    // reads as zeros up to the end of the last page: the version of a file of fewer than 4 bytes
    // is not 4 or 5, a missing string is empty, a missing count is 0. A file of no bytes has
    // no page: NCBI maps it as a null pointer and crashes. A read beyond the last page is outside
    // the mapping: it reads whatever memory follows, which LOSAT cannot know.
    /// The header of the index file `file` of a database of `seqtype`, as NCBI's constructor
    /// reads it, whatever follows it (the arrays are not looked at).
    ///
    /// The error is NCBI's own where NCBI throws, and where NCBI does not give a result LOSAT can
    /// reproduce, LOSAT's refusal: [`SeqDbError::NcbiCrashes`] for a file of no bytes and
    /// [`SeqDbError::Unsupported`] for a read that leaves the mapping of the file, where what NCBI
    /// reads, or whether it faults, depends on its memory.
    pub(crate) fn read_mapped_header(
        file: &DbFile,
        seqtype: SeqType,
    ) -> Result<VolumeHeader, SeqDbError> {
        let path = file.path();
        let mut cursor = Cursor {
            source: &Mapped(file),
            pos: 0,
            path,
        };
        let header = Self::read_header(&mut cursor, path).map_err(|fault| match fault {
            ReadFault::PastEnd { .. } if file.is_empty() => SeqDbError::NcbiCrashes {
                path: path.to_path_buf(),
                message: "the index file has no bytes, which NCBI maps as a null pointer"
                    .to_string(),
            },
            ReadFault::PastEnd {
                wanted, at, len, ..
            } => SeqDbError::Unsupported {
                path: path.to_path_buf(),
                message: format!(
                    "the index file ends inside its header ({wanted} bytes wanted at offset \
                     {at}, the file has {} bytes) and NCBI BLAST+'s mapping of it ends at \
                     {len}: NCBI reads memory past its mapping there, whose bytes, or \
                     whether NCBI faults, depend on its memory map",
                    file.len()
                ),
            },
            ReadFault::Error(error) => error,
        })?;
        let file_seqtype = if header.db_seqtype == 1 {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        if file_seqtype != seqtype {
            return Err(SeqDbError::Format {
                path: path.to_path_buf(),
                message: "Error: requested sequence type does not match DB.".to_string(),
            });
        }
        Ok(VolumeHeader {
            title: header.title,
            lmdb_file: header.lmdb_file,
            date: header.date,
            num_oids: header.num_oids,
            volume_length: header.volume_length,
        })
    }

    /// The values before the arrays: what the constructor of `CSeqDBIdxFile` reads, in its
    /// order, up to the longest sequence's length.
    fn read_header<S: Source + ?Sized>(
        cursor: &mut Cursor<'_, S>,
        path: &Path,
    ) -> Result<Header, ReadFault> {
        let version = cursor.be_u32()?;
        let dbv5 = version == 5;
        if !dbv5 && version != 4 {
            return Err(ReadFault::Error(SeqDbError::Format {
                path: path.to_path_buf(),
                message: "Error: Not a valid version 4 or 5 database.".to_string(),
            }));
        }
        let db_seqtype = cursor.be_u32()?;
        let volume_number = if dbv5 { Some(cursor.be_u32()?) } else { None };
        let title = cursor.string()?;
        let lmdb_file = if dbv5 { Some(cursor.string()?) } else { None };
        let date = cursor.string()?;
        let num_oids = cursor.be_u32()?;
        let volume_length = cursor.broken_u64()?;
        let max_length = cursor.be_u32()?;
        Ok(Header {
            version,
            db_seqtype,
            volume_number,
            title,
            lmdb_file,
            date,
            num_oids,
            volume_length,
            max_length,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:175-235
    /// The values before the arrays and where the arrays start; the arrays' lengths are
    /// checked against the file's.
    fn parse_fixed<S: Source + ?Sized>(
        source: &S,
        seqtype: SeqType,
        path: &Path,
    ) -> Result<Fixed, SeqDbError> {
        let format = |message: &str| SeqDbError::Format {
            path: path.to_path_buf(),
            message: message.to_string(),
        };
        let mut cursor = Cursor {
            source,
            pos: 0,
            path,
        };
        let Header {
            version,
            db_seqtype,
            volume_number,
            title,
            lmdb_file,
            date,
            num_oids,
            volume_length,
            max_length,
        } = Self::read_header(&mut cursor, path)?;
        // Only the tests look at the version and the volume number.
        #[cfg(not(test))]
        let _ = (version, volume_number);

        let file_seqtype = if db_seqtype == 1 {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        if file_seqtype != seqtype {
            return Err(format("Error: requested sequence type does not match DB."));
        }
        // NCBI holds the count in a signed 32-bit `int`.
        if num_oids > i32::MAX as u32 {
            return Err(format(
                "index file holds more sequences than a signed 32-bit count",
            ));
        }

        let region = num_oids as usize + 1;
        let hdr_at = cursor.offsets(region)?;
        let seq_at = cursor.offsets(region)?;
        let amb_at = match seqtype {
            SeqType::Nucleotide => cursor.offsets(region)?,
            SeqType::Protein => 0,
        };
        Ok(Fixed {
            #[cfg(test)]
            version,
            #[cfg(test)]
            volume_number,
            title,
            lmdb_file,
            date,
            num_oids,
            volume_length,
            max_length,
            hdr_at,
            seq_at,
            amb_at,
        })
    }

    fn with_offsets(fixed: Fixed, seqtype: SeqType, offsets: Offsets) -> Self {
        Self {
            #[cfg(test)]
            version: fixed.version,
            seqtype,
            #[cfg(test)]
            volume_number: fixed.volume_number,
            title: fixed.title,
            lmdb_file: fixed.lmdb_file,
            date: fixed.date,
            num_oids: fixed.num_oids,
            volume_length: fixed.volume_length,
            max_length: fixed.max_length,
            offsets,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:134-236
    // NCBI maps the index file (zeros to the end of the last 4 KiB page) and does not validate
    // offset monotonicity at open: non-ascending offsets are handled at sequence/header read time.
    /// Parses the bytes of an index file opened as `seqtype`, and holds its arrays.
    pub fn parse(bytes: &[u8], seqtype: SeqType, path: &Path) -> Result<Self, SeqDbError> {
        let mapped = MappedSlice(bytes);
        let fixed = Self::parse_fixed(&mapped, seqtype, path)?;
        let region = fixed.num_oids as usize + 1;
        let array = |at: u64| -> Result<Vec<u32>, SeqDbError> {
            let mut raw = vec![0u8; 4 * region];
            mapped.read_at(&mut raw, at)?;
            let mut words = vec![0u32; region];
            decode_words(&raw, &mut words);
            Ok(words)
        };
        let held = HeldOffsets {
            hdr: array(fixed.hdr_at)?,
            seq: array(fixed.seq_at)?,
            amb: match seqtype {
                SeqType::Nucleotide => array(fixed.amb_at)?,
                SeqType::Protein => Vec::new(),
            },
        };
        Ok(Self::with_offsets(fixed, seqtype, Offsets::Held(held)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:134-236
    // NCBI maps the index file and reads an offset where it needs one (`x_GetHdr()[oid]`, cited
    // at `header_range`). LOSAT reads the fixed part when it opens the file (zero-filling the
    // last 4 KiB mapped page) and then the ranges of the arrays it needs with positioned reads.
    /// EXPERIMENT (LOSAT_X_DBIDXREAD): opens the index file `file` as `seqtype`, the arrays left
    /// in the file. The errors are those of `parse`.
    pub fn open(file: DbFile, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let file = file.zero_pad_last_page();
        let path = file.path().to_path_buf();
        let fixed = Self::parse_fixed(&Mapped(&file), seqtype, &path)?;
        let (hdr_at, seq_at, amb_at) = (fixed.hdr_at, fixed.seq_at, fixed.amb_at);
        Ok(Self::with_offsets(
            fixed,
            seqtype,
            Offsets::Read(ReadOffsets {
                file,
                hdr_at,
                seq_at,
                amb_at,
            }),
        ))
    }

    #[cfg(test)]
    pub(super) fn open_checked(
        file: DbFile,
        seqtype: SeqType,
        _chunk_entries: usize,
    ) -> Result<Self, SeqDbError> {
        Self::open(file, seqtype)
    }

    /// Fills `out` with the entries `first..first + out.len()` of `region`.
    fn read_words(
        &self,
        region: Region,
        first: usize,
        out: &mut [u32],
        scratch: &mut Vec<u8>,
    ) -> Result<(), SeqDbError> {
        let entries = match (region, self.seqtype) {
            (Region::Ambiguity, SeqType::Protein) => 0,
            _ => self.num_oids as usize + 1,
        };
        assert!(
            first <= entries && out.len() <= entries - first,
            "index entries {first}..{} outside the {entries} of the array",
            first + out.len()
        );
        match &self.offsets {
            Offsets::Held(held) => out.copy_from_slice(&held.get(region)[first..first + out.len()]),
            Offsets::Read(read) => {
                let at = match region {
                    Region::Header => read.hdr_at,
                    Region::Sequence => read.seq_at,
                    Region::Ambiguity => read.amb_at,
                } + 4 * first as u64;
                let mut stack = [0u8; 64];
                let bytes = if out.len() * 4 <= stack.len() {
                    &mut stack[..out.len() * 4]
                } else {
                    scratch.resize(out.len() * 4, 0);
                    &mut scratch[..]
                };
                read.file.read_exact_at(bytes, at)?;
                decode_words(bytes, out);
            }
        }
        Ok(())
    }

    /// The error of an index offset range that cannot be read as ascending sequence or ambiguity
    /// bounds (`Volume::block_bytes`, `Volume::x_read_block_end`, `Volume::read_block`).
    fn changed(&self, region: Region, first: usize, count: usize) -> SeqDbError {
        let path = match &self.offsets {
            Offsets::Read(read) => read.file.path().to_path_buf(),
            Offsets::Held(_) => std::path::PathBuf::new(),
        };
        SeqDbError::Format {
            path,
            message: format!(
                "{region:?} offsets {first}..{} are not in ascending order",
                first + count
            ),
        }
    }

    /// The error of `changed` for the sequence offsets `first..first + count`.
    pub(super) fn sequence_changed(&self, first: u32, count: usize) -> SeqDbError {
        self.changed(Region::Sequence, first as usize, count)
    }

    /// The error of `changed` for the ambiguity offsets `first..first + count`.
    pub(super) fn ambiguity_changed(&self, first: u32, count: usize) -> SeqDbError {
        self.changed(Region::Ambiguity, first as usize, count)
    }

    /// The entries `first..first + count` of `region`.
    fn words(&self, region: Region, first: u32, count: usize) -> Result<Vec<u32>, SeqDbError> {
        let mut out = vec![0u32; count];
        self.read_words(region, first as usize, &mut out, &mut Vec::new())?;
        Ok(out)
    }

    fn pair(&self, region: Region, first: u32) -> Result<(u64, u64), SeqDbError> {
        let mut out = [0u32; 2];
        self.read_words(region, first as usize, &mut out, &mut Vec::new())?;
        Ok((u64::from(out[0]), u64::from(out[1])))
    }

    fn word(&self, region: Region, at: u32) -> Result<u32, SeqDbError> {
        let mut out = [0u32; 1];
        self.read_words(region, at as usize, &mut out, &mut Vec::new())?;
        Ok(out[0])
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:664-669
    /// Byte range of the header (`Blast-def-line-set`) of `oid` in the header file.
    pub fn header_range(&self, oid: u32) -> Result<(u64, u64), SeqDbError> {
        self.pair(Region::Header, oid)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:671-681
    /// Byte range of the sequence data of `oid` in the sequence file. For a
    /// protein the range includes the NUL byte that follows the residues; for
    /// a nucleotide it ends where the ambiguity data start.
    pub fn sequence_range(&self, oid: u32) -> Result<(u64, u64), SeqDbError> {
        match self.seqtype {
            SeqType::Protein => self.pair(Region::Sequence, oid),
            SeqType::Nucleotide => Ok((
                u64::from(self.word(Region::Sequence, oid)?),
                u64::from(self.word(Region::Ambiguity, oid)?),
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:651-662
    /// Byte range of the ambiguity data of `oid` (nucleotide volumes only).
    #[cfg(test)]
    pub fn ambiguity_range(&self, oid: u32) -> Result<Option<(u64, u64)>, SeqDbError> {
        Ok(self
            .ambiguity_words(oid)?
            .map(|(start, end)| (u64::from(start), u64::from(end))))
    }

    /// The ambiguity range of `oid` as stored (`ambiguity_range`).
    #[cfg(test)]
    pub(super) fn ambiguity_words(&self, oid: u32) -> Result<Option<(u32, u32)>, SeqDbError> {
        match self.seqtype {
            SeqType::Protein => Ok(None),
            SeqType::Nucleotide => Ok(Some((
                self.word(Region::Ambiguity, oid)?,
                self.word(Region::Sequence, oid + 1)?,
            ))),
        }
    }

    /// Offset in the sequence file where the data of `oid` start; `oid` may be
    /// `num_oids`, which gives the end of the last sequence's data.
    pub fn sequence_start(&self, oid: u32) -> Result<u64, SeqDbError> {
        Ok(u64::from(self.word(Region::Sequence, oid)?))
    }

    /// The sequence start offsets of `first..=end` (`end` may be `num_oids`): the starts of
    /// the sequences `first..end` and the end of the last one's data.
    pub fn sequence_offsets(&self, first: u32, end: u32) -> Result<Vec<u32>, SeqDbError> {
        self.words(Region::Sequence, first, (end - first) as usize + 1)
    }

    /// The ambiguity start offsets of the sequences `first..end` (empty for protein).
    pub fn ambiguity_offsets(&self, first: u32, end: u32) -> Result<Vec<u32>, SeqDbError> {
        match self.seqtype {
            SeqType::Protein => Ok(Vec::new()),
            SeqType::Nucleotide => self.words(Region::Ambiguity, first, (end - first) as usize),
        }
    }

    /// The header start offsets of `first..=end` (`end` may be `num_oids`): the starts of the
    /// headers of the sequences `first..end` and the end of the last one's.
    pub fn header_offsets(&self, first: u32, end: u32) -> Result<Vec<u32>, SeqDbError> {
        self.words(Region::Header, first, (end - first) as usize + 1)
    }

    /// The sequence start offsets, `num_oids + 1` entries, when the arrays are held
    /// (`parse`, `LOSAT_X_DBIDXREAD=0`).
    pub(super) fn held_sequence_offsets(&self) -> Option<&[u32]> {
        match &self.offsets {
            Offsets::Held(held) => Some(&held.seq),
            Offsets::Read(_) => None,
        }
    }

    /// Whether the arrays are held in memory (`parse`) rather than read from the file.
    pub fn holds_offsets(&self) -> bool {
        matches!(self.offsets, Offsets::Held(_))
    }
}

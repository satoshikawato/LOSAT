//! The GI-indexed masks of a database (`makeblastdb -gi_mask`): the files `<name>.gni` (index),
//! `<name>.gno` (offsets) and `<name>.gnd` or `<name>.NN.gnd` (data) of each name of the alias
//! file's `MASKLIST` (`alias::AliasTree::gi_masks`), and the GI of each sequence that keys them.
//!
//! NCBI reads them with `CSeqDBGiMask` (`seqdbgimask.cpp`, `seqdbgimask.hpp`), and the GI with
//! `CSeqDBImpl::x_GetSeqGI` and the volume's `CSeqDBGiIndex` (the `.nog` file). This module is
//! the port; `mask.rs` dispatches to it where the database has GI masks, and
//! `blastinput::db_mask` and `losat blastdbcmd` read them through `mask.rs`.
//!
//! NCBI maps the files and reads its arrays in place. LOSAT reads the index, offset and `.nog`
//! files, whose reads are binary searches and one GI per sequence, in chunks of 64 KiB on first
//! use (`MappedFile`), kept for the life of the database as NCBI's mapping keeps the pages it
//! touched. A data file is read through the reader's window (`DataWindow`, one per
//! `mask::MaskBlock`): makeblastdb writes the records of the sequences in their order, so a
//! reader that asks for one sequence after another reads each part of the file once and keeps
//! one window, not the file.
//!
//! Errors: NCBI's `CSeqDBException` texts are NCBI's ([`SeqDbError::NotFound`],
//! [`SeqDbError::Format`], [`SeqDbError::Range`], [`SeqDbError::MemoryMap`]). A failed
//! `SEQDB_FILE_ASSERT` (its text names a source line) is [`SeqDbError::Integrity`] in LOSAT's
//! words; a read outside a file, where NCBI reads past its mapping, is [`SeqDbError::Format`] in
//! LOSAT's words. The `.nog`, which NCBI reads without a check, reads as NCBI's mapping of it
//! does (zeros to the end of its last page); NCBI's crash on an empty one is
//! [`SeqDbError::NcbiCrashes`], and a GI past the mapping, where NCBI reads memory it never
//! mapped, [`SeqDbError::Unsupported`]: LOSAT's refusals.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use super::column::Blob;
use super::file::{mapped_length, DbFile};
use super::mask::push_padded;
use super::volume::{with_extension, Volume};
use super::{SearchPath, SeqDbError, SeqType};
use crate::utils::os_bytes::path_of_bytes;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:124-133
/// The extensions of the index, offset and data files that a host of this byte order reads; their
/// arrays are in the host's order.
const EXTENSIONS: [&str; 3] = if cfg!(target_endian = "little") {
    ["gni", "gno", "gnd"]
} else {
    ["gmi", "gmo", "gmd"]
};

/// The bytes of the fixed fields at the start of an index file (`kFixedFieldBytes`).
const FIXED_FIELD_BYTES: u64 = 32;

/// The size of the chunks in which `MappedFile` reads a file.
const CHUNK: u64 = 1 << 16;

/// The least number of bytes of a data file that a [`DataWindow`] reads at once.
const WINDOW: u64 = 1 << 16;

/// What a [`MappedFile`] is, in LOSAT's errors of its reads.
const INDEX_FILE: &str = "GI-mask index file";
const OFFSET_FILE: &str = "GI-mask offset file";
const NOG_FILE: &str = "GI file of the volume (.nog)";

/// A file read in chunks on first use and kept, as NCBI's memory map of it.
#[derive(Debug)]
struct MappedFile {
    file: DbFile,
    /// What the file is ([`INDEX_FILE`], ...), for the errors.
    kind: &'static str,
    chunks: Box<[OnceLock<Box<[u8]>>]>,
}

impl MappedFile {
    fn open(path: &Path, kind: &'static str) -> Result<Self, SeqDbError> {
        let file = DbFile::open(path)?;
        let count = file.len().div_ceil(CHUNK);
        let chunks = (0..count).map(|_| OnceLock::new()).collect();
        Ok(Self { file, kind, chunks })
    }

    fn len(&self) -> u64 {
        self.file.len()
    }

    fn path(&self) -> &Path {
        self.file.path()
    }

    fn chunk(&self, index: usize) -> Result<&[u8], SeqDbError> {
        if let Some(bytes) = self.chunks[index].get() {
            return Ok(bytes);
        }
        let start = index as u64 * CHUNK;
        let bytes = self
            .file
            .read_range(start, (start + CHUNK).min(self.len()))?;
        Ok(self.chunks[index].get_or_init(|| bytes.into_boxed_slice()))
    }

    /// The bytes at `start..start + out.len()` into `out`. A range outside the file, which NCBI
    /// reads outside its mapping, is a format error in LOSAT's words.
    fn read_at(&self, start: i64, out: &mut [u8]) -> Result<(), SeqDbError> {
        let end = start.saturating_add(out.len() as i64);
        if start < 0 || end as u64 > self.len() {
            return Err(SeqDbError::Format {
                path: self.path().to_path_buf(),
                message: format!(
                    "a read of bytes {start}..{end} of the {} is outside the file ({} bytes)",
                    self.kind,
                    self.len()
                ),
            });
        }
        let mut at = start as u64;
        let mut filled = 0;
        while filled < out.len() {
            let chunk = self.chunk((at / CHUNK) as usize)?;
            let from = (at % CHUNK) as usize;
            let take = (chunk.len() - from).min(out.len() - filled);
            out[filled..filled + take].copy_from_slice(&chunk[from..from + take]);
            filled += take;
            at += take as u64;
        }
        Ok(())
    }

    /// The `N` bytes at `at`: from their chunk, a copy of a constant length, where one chunk
    /// holds them all (every read of the arrays but one that crosses a chunk's end); otherwise
    /// [`MappedFile::read_at`], with its error outside the file. The binary searches read the
    /// arrays a number at a time, which a copy of a length known only at run time made a call
    /// of `memmove` each (review round of M4c-2, 2026-11-24: 50 % of a masked search's cycles on
    /// a database of 1,000,000 sequences).
    fn array_at<const N: usize>(&self, at: i64) -> Result<[u8; N], SeqDbError> {
        let mut bytes = [0u8; N];
        if let Ok(start) = u64::try_from(at) {
            let from = (start % CHUNK) as usize;
            if from + N <= CHUNK as usize && start + N as u64 <= self.len() {
                let chunk = self.chunk((start / CHUNK) as usize)?;
                bytes.copy_from_slice(&chunk[from..from + N]);
                return Ok(bytes);
            }
        }
        self.read_at(at, &mut bytes)?;
        Ok(bytes)
    }

    /// The 4 bytes at `at` as a number in the host's order (the arrays of the GI-mask files).
    fn u32_at(&self, at: i64) -> Result<u32, SeqDbError> {
        self.array_at::<4>(at).map(u32::from_ne_bytes)
    }

    /// [`check_range`] of this file.
    fn check_range(&self, start: i64, end: i64, what: &str) -> Result<(), SeqDbError> {
        check_range(&self.file, start, end, what)
    }

    /// The `N` bytes at `at` as NCBI reads them through its mapping of the file, which it does
    /// not check against the file's length (the `.nog`, `CSeqDBGiIndex::GetSeqGI`): the file's
    /// bytes, then zeros to the end of its last page ([`mapped_length`]). `None` where they are
    /// not all in that mapping (before the file, or past its last page).
    fn mapped<const N: usize>(&self, at: i64) -> Result<Option<[u8; N]>, SeqDbError> {
        let end = at.saturating_add(N as i64);
        if at < 0 || end as u64 > mapped_length(self.len()) {
            return Ok(None);
        }
        if end as u64 <= self.len() {
            return self.array_at::<N>(at).map(Some);
        }
        let mut bytes = [0u8; N];
        let held = self.len().saturating_sub(at as u64).min(N as u64) as usize;
        if held > 0 {
            self.read_at(at, &mut bytes[..held])?;
        }
        Ok(Some(bytes))
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:123-130
/// The two `SEQDB_FILE_ASSERT`s of `GetFileDataPtr` (with a lease and the range's start and
/// end) on `file`, `what` naming the range in LOSAT's error.
fn check_range(file: &DbFile, start: i64, end: i64, what: &str) -> Result<(), SeqDbError> {
    if start < end && end as u64 <= file.len() && end >= 0 {
        return Ok(());
    }
    Err(SeqDbError::Integrity {
        path: file.path().to_path_buf(),
        message: format!(
            "{what} (bytes {start}..{end}) is empty or passes the end of the file ({} bytes)",
            file.len()
        ),
    })
}

/// The bytes of a GI-mask data file that one reader keeps from one record to the next (a
/// `mask::MaskBlock` holds one): at least [`WINDOW`] bytes read at once from the start of a
/// record, so that the records of the sequences one after the other, which makeblastdb writes in
/// that order, take one read a window, and the reader holds a window, not the file. No NCBI
/// counterpart: NCBI maps the data files.
#[derive(Debug, Default)]
pub(super) struct DataWindow {
    /// The files ([`GiMaskFiles`]'s `id`) and the data volume the bytes are of; `None`: none.
    of: Option<(u64, usize)>,
    /// Where the bytes start in that file.
    start: u64,
    bytes: Vec<u8>,
}

impl DataWindow {
    /// The `len` bytes at `start` of `file`, the data volume that `of` names, which holds them
    /// (the caller checked the range).
    fn get(
        &mut self,
        of: (u64, usize),
        file: &DbFile,
        start: u64,
        len: u64,
    ) -> Result<&[u8], SeqDbError> {
        let held = self.of == Some(of)
            && start >= self.start
            && start + len <= self.start + self.bytes.len() as u64;
        if !held {
            let end = start.saturating_add(len.max(WINDOW)).min(file.len());
            self.of = None;
            self.bytes.resize((end - start) as usize, 0);
            file.read_exact_at(&mut self.bytes, start)?;
            self.of = Some(of);
            self.start = start;
        }
        let from = (start - self.start) as usize;
        Ok(&self.bytes[from..from + len as usize])
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:204-245
// The port reads each key where NCBI does (`key_at`), so a list with equal GIs (a version 4
// database, survey D7) gives the entry NCBI's search lands on. The indexes are 64-bit, where
// NCBI's `int` sum could overflow.
/// `s_BinarySearch` over the `n` keys that `key_at` reads: whether `key` was found, and the index
/// (the entry found; otherwise -1 for a key outside the keys, else the entry before it).
fn binary_search(
    n: i64,
    key: u32,
    key_at: impl Fn(i64) -> Result<u32, SeqDbError>,
) -> Result<(bool, i64), SeqDbError> {
    let (mut lower, mut upper) = (0i64, n - 1);
    if key > key_at(upper)? || key < key_at(lower)? {
        return Ok((false, -1));
    }
    if key == key_at(upper)? {
        return Ok((true, upper));
    }
    if key == key_at(lower)? {
        return Ok((true, lower));
    }
    let mut idx = (lower + upper) / 2;
    while idx != lower {
        let at = key_at(idx)?;
        if key > at {
            lower = idx;
            idx = (lower + upper) / 2;
        } else if key < at {
            upper = idx;
            idx = (lower + upper) / 2;
        } else {
            return Ok((true, idx));
        }
    }
    Ok((false, idx))
}

/// The numbers of the opened [`GiMaskFiles`], each its own (the key of a [`DataWindow`]).
static OPENED: AtomicU64 = AtomicU64::new(0);

/// The opened files of one algorithm (`x_Open`, `x_ReadFields`).
#[derive(Debug)]
pub(super) struct GiMaskFiles {
    /// This opening's own number (`OPENED`).
    id: u64,
    index: MappedFile,
    offsets: MappedFile,
    data: Vec<DbFile>,
    gi_size: i32,
    offset_size: i32,
    page_size: i32,
    num_index: i32,
    num_gis: i32,
    index_start: i32,
    desc: Vec<u8>,
}

/// The fields of an index file (`x_ReadFields`) and the number of data volumes.
struct Fields {
    num_vols: i32,
    gi_size: i32,
    offset_size: i32,
    page_size: i32,
    num_index: i32,
    num_gis: i32,
    index_start: i32,
    desc: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:200-251
// `ReferTo` keeps the blob's read offset, so the strings are read from byte 32 of the header.
// The sums are 64-bit, where NCBI's `int` product could overflow.
/// `x_ReadFields` of the index file `index`.
fn read_fields(index: &MappedFile) -> Result<Fields, SeqDbError> {
    let format = |message: String| SeqDbError::Format {
        path: index.path().to_path_buf(),
        message,
    };
    index.check_range(
        0,
        FIXED_FIELD_BYTES as i64,
        "the fixed fields of the GI-mask index",
    )?;
    let mut header = [0u8; FIXED_FIELD_BYTES as usize];
    index.read_at(0, &mut header)?;
    let mut blob = Blob::new(&header, 0);
    let mut field = || blob.read_i32().map_err(format);
    if field()? != 1 {
        return Err(format(
            "Gi-mask file uses unknown format_version.".to_string(),
        ));
    }
    let num_vols = field()?;
    let gi_size = field()?;
    let offset_size = field()?;
    let page_size = field()?;
    let num_index = field()?;
    let num_gis = field()?;
    let index_start = field()?;
    let integrity = |message: String| SeqDbError::Integrity {
        path: index.path().to_path_buf(),
        message,
    };
    if index_start < 0 {
        return Err(integrity(format!(
            "the GI-mask index states a negative start of its arrays ({index_start})"
        )));
    }
    if index.len() < index_start as u64 {
        return Err(integrity(format!(
            "the GI-mask index has {} bytes, fewer than the start of its arrays ({index_start})",
            index.len()
        )));
    }
    index.check_range(0, i64::from(index_start), "the header of the GI-mask index")?;
    let mut head = vec![0u8; index_start as usize];
    index.read_at(0, &mut head)?;
    let mut blob = Blob::new(&head, FIXED_FIELD_BYTES as i64);
    let desc = blob.read_string_var().map_err(format)?.to_vec();
    let date = blob.read_string_var().map_err(format)?;
    if desc.is_empty() || date.is_empty() {
        return Err(integrity(
            "the GI-mask index has an empty description or date".to_string(),
        ));
    }
    let begin = i64::from(index_start);
    let end = begin + i64::from(num_index) * (i64::from(gi_size) + i64::from(offset_size));
    index.check_range(begin, end, "the arrays of the GI-mask index")?;
    Ok(Fields {
        num_vols,
        gi_size,
        offset_size,
        page_size,
        num_index,
        num_gis,
        index_start,
        desc,
    })
}

impl GiMaskFiles {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:137-178
    // The files are looked for along the search path (the current directory, `BLASTDB`, the
    // registry), never beside the alias file (survey D5). A data file that is not found maps the
    // empty name before `found` is tested (`Init` of the volume's data lease with an empty
    // name): the corelib's "To be memory mapped the file must exist: ''"
    // (`CSeqDBException::eMemoryMappingFailure`, seqdbatlas.hpp:709-719), so "Could not open
    // gi-mask data file(s)." is never thrown (survey D8). The data files of an earlier algorithm
    // stay in NCBI's vectors, but `[vol]` is opened again, so each algorithm reads its own files.
    /// `x_Open` of the mask `name`: its index, offset and data files.
    fn open(name: &str, search: &SearchPath) -> Result<Self, SeqDbError> {
        let index_path = search.resolve_file(format!("{name}.{}", EXTENSIONS[0]).as_bytes());
        let offset_path = search.resolve_file(format!("{name}.{}", EXTENSIONS[1]).as_bytes());
        let (Some(index_path), Some(offset_path)) = (index_path, offset_path) else {
            return Err(SeqDbError::NotFound(
                b"Could not open gi-mask index files.".to_vec(),
            ));
        };
        let index = MappedFile::open(&path_of_bytes(&index_path), INDEX_FILE)?;
        let offsets = MappedFile::open(&path_of_bytes(&offset_path), OFFSET_FILE)?;
        let fields = read_fields(&index)?;
        let mut data = Vec::new();
        for vol in 0..fields.num_vols.max(0) {
            let file = if fields.num_vols == 1 {
                format!("{name}.{}", EXTENSIONS[2])
            } else {
                format!("{name}.{}{}.{}", vol / 10, vol % 10, EXTENSIONS[2])
            };
            let Some(path) = search.resolve_file(file.as_bytes()) else {
                return Err(SeqDbError::MemoryMap(super::oid_list::not_a_file(b"")));
            };
            data.push(DbFile::open(&path_of_bytes(&path))?);
        }
        Ok(Self {
            id: OPENED.fetch_add(1, Ordering::Relaxed),
            index,
            offsets,
            data,
            gi_size: fields.gi_size,
            offset_size: fields.offset_size,
            page_size: fields.page_size,
            num_index: fields.num_index,
            num_gis: fields.num_gis,
            index_start: fields.index_start,
            desc: fields.desc,
        })
    }

    /// The `(vol, off)` of the record of `key`, or `None` without one: the index, then the page
    /// of the offset file (`GetMaskData`, cited at [`GiMaskFiles::ranges_into`]).
    fn locate(&self, key: u32) -> Result<Option<(i32, i32)>, SeqDbError> {
        let index_start = i64::from(self.index_start);
        let num_index = i64::from(self.num_index);
        let (found, page) =
            binary_search(num_index, key, |i| self.index.u32_at(index_start + 4 * i))?;
        let (file, at) = if found {
            (&self.index, index_start + 4 * (num_index + page * 2))
        } else {
            if page == -1 {
                return Ok(None);
            }
            let page_size = i64::from(self.page_size);
            let num_gis = i64::from(self.num_gis);
            let mut pagesize = page_size;
            if page * page_size + pagesize > num_gis {
                pagesize = num_gis - page * page_size;
            }
            let entry = i64::from(self.gi_size) + i64::from(self.offset_size);
            let begin = page * page_size * entry;
            let end = begin + pagesize * entry;
            self.offsets
                .check_range(begin, end, "a page of the GI-mask offsets")?;
            let (found, idx) =
                binary_search(pagesize, key, |i| self.offsets.u32_at(begin + 4 * i))?;
            if !found {
                return Ok(None);
            }
            (&self.offsets, begin + 4 * (pagesize + idx * 2))
        };
        // NCBI assigns the `Uint4` values to `int`s.
        Ok(Some((file.u32_at(at)? as i32, file.u32_at(at + 4)? as i32)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:73-112
    // The arrays are indexed as 4-byte GIs and 4-byte (vol, off) pairs whatever the header's
    // sizes say; each read is checked against the file (NCBI reads its mapping). A `vol` outside
    // the data files (NCBI's debug-only `_ASSERT`, then an index past its vector) is a format
    // error in LOSAT's words; the sums are 64-bit, where NCBI's `int` ones could overflow.
    /// `GetMaskData` of the key `key` (the sequence's GI as a `Uint4`, `GiMasks::seq_gi_key`):
    /// its ranges appended to `out` in the order stored, none without a record; the data file is
    /// read through the reader's `window`.
    pub(super) fn ranges_into(
        &self,
        key: u32,
        out: &mut Vec<(u32, u32)>,
        window: &mut DataWindow,
    ) -> Result<(), SeqDbError> {
        let Some((vol, off)) = self.locate(key)? else {
            return Ok(());
        };
        let Some(at) = usize::try_from(vol).ok().filter(|&at| at < self.data.len()) else {
            return Err(SeqDbError::Format {
                path: self.index.path().to_path_buf(),
                message: format!(
                    "the GI-mask record of GI {key} is in data volume {vol}, of {} volumes",
                    self.data.len()
                ),
            });
        };
        let off = i64::from(off);
        let count = self.record_bytes(at, off, off + 4, "the count of a GI-mask record", window)?;
        let n = i64::from(i32::from_ne_bytes([count[0], count[1], count[2], count[3]]));
        let what = "the ranges of a GI-mask record";
        let pairs = self.record_bytes(at, off + 4, off + 8 * n + 4, what, window)?;
        out.extend(pairs.chunks_exact(8).map(|pair| {
            (
                u32::from_ne_bytes([pair[0], pair[1], pair[2], pair[3]]),
                u32::from_ne_bytes([pair[4], pair[5], pair[6], pair[7]]),
            )
        }));
        Ok(())
    }

    /// The bytes `start..end` of the data volume `at` through `window`, after NCBI's checks of
    /// the range (`what` names it); a range before the file, which those checks let through
    /// where it ends inside the file and NCBI reads before its mapping, is a format error in
    /// LOSAT's words.
    fn record_bytes<'w>(
        &self,
        at: usize,
        start: i64,
        end: i64,
        what: &str,
        window: &'w mut DataWindow,
    ) -> Result<&'w [u8], SeqDbError> {
        let file = &self.data[at];
        check_range(file, start, end, what)?;
        if start < 0 {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!("{what} (bytes {start}..{end}) starts before the file"),
            });
        }
        window.get((self.id, at), file, start as u64, (end - start) as u64)
    }
}

/// The `.nog` file of a volume (`CSeqDBGiIndex`): the GI of each OID, read as NCBI reads it
/// (`MappedFile::mapped`).
#[derive(Debug)]
struct GiIndex {
    file: MappedFile,
    size: i32,
    num_oids: i32,
}

impl GiIndex {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:52-71
    // NCBI maps the file and reads the entry size and the count at bytes 8..16 of the mapping,
    // with no check of its length: a file cut short reads as zeros to the end of its page (a
    // count of 0 sends every sequence to its definition lines; an OID's entry cut off reads as
    // GI 0), and a file of no bytes is mapped as a null pointer, which NCBI's read then faults
    // on.
    /// The `.nog` file at `path`, its header read as NCBI's first `GetSeqGI` reads it.
    fn open(path: &Path) -> Result<Self, SeqDbError> {
        let file = MappedFile::open(path, NOG_FILE)?;
        let header = match file.mapped::<8>(8)? {
            Some(header) => header,
            None => {
                return Err(SeqDbError::NcbiCrashes {
                    path: path.to_path_buf(),
                    message: format!(
                        "the {NOG_FILE} has no bytes, which NCBI maps as a null pointer"
                    ),
                })
            }
        };
        Ok(Self {
            file,
            size: i32::from_be_bytes([header[0], header[1], header[2], header[3]]),
            num_oids: i32::from_be_bytes([header[4], header[5], header[6], header[7]]),
        })
    }

    // The offset is 64-bit, where NCBI's `int` product of the OID and the entry size could
    // overflow.
    /// The GI of the sequence `vol_oid` (below the count) as a `Uint4`: big-endian at
    /// `vol_oid * size + 32` of NCBI's mapping; LOSAT's refusal where that is outside it.
    fn gi(&self, vol_oid: u32) -> Result<u32, SeqDbError> {
        let at = i64::from(vol_oid) * i64::from(self.size) + 32;
        match self.file.mapped::<4>(at)? {
            Some(bytes) => Ok(u32::from_be_bytes(bytes)),
            None => Err(SeqDbError::Unsupported {
                path: self.file.path().to_path_buf(),
                message: format!(
                    "the GI of sequence {vol_oid} of the volume is read at offset {at} of the \
                     {NOG_FILE}, which has {} bytes, outside the {} bytes of NCBI BLAST+'s \
                     mapping of it: NCBI reads memory past its mapping there, whose bytes, or \
                     whether NCBI faults, depend on its memory map",
                    self.file.len(),
                    mapped_length(self.file.len())
                ),
            }),
        }
    }
}

/// The GI masks of a database whose alias file names them (`CSeqDBGiMask`): one algorithm per
/// name of `MASKLIST`, with ids 0 to N-1 in its order.
#[derive(Debug)]
pub struct GiMasks {
    names: Vec<String>,
    search: SearchPath,
    seqtype: SeqType,
    /// The files of each algorithm, opened on first use (a failed opening is tried again, as
    /// NCBI's `x_Open` is).
    opened: Box<[OnceLock<Arc<GiMaskFiles>>]>,
    /// The `.nog` file of each volume of the database, opened on first use.
    gi_indexes: Box<[OnceLock<Option<GiIndex>>]>,
}

impl GiMasks {
    /// The masks `names` (`MASKLIST`) of a database of `volumes` volumes of `seqtype`, whose
    /// files are looked for along `search`. Nothing is read here.
    pub(super) fn new(
        names: Vec<String>,
        search: SearchPath,
        seqtype: SeqType,
        volumes: usize,
    ) -> Self {
        let opened = names.iter().map(|_| OnceLock::new()).collect();
        let gi_indexes = (0..volumes).map(|_| OnceLock::new()).collect();
        Self {
            names,
            search,
            seqtype,
            opened,
            gi_indexes,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:99-105
    /// `GetAvailableMaskAlgorithms`: 0 to N-1. No file is read.
    pub fn ids(&self) -> Vec<i32> {
        (0..self.names.len() as i32).collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:130-143
    fn available_names(&self) -> String {
        let mut out = String::from("\nAvailable filtering algorithm(s):\n\n");
        push_padded(&mut out, "Algorithm ID", 14);
        push_padded(&mut out, "Algorithm name", 40);
        out.push('\n');
        for (id, name) in self.names.iter().enumerate() {
            out.push_str("    ");
            push_padded(&mut out, &id.to_string(), 10);
            push_padded(&mut out, name, 40);
            out.push('\n');
        }
        out
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:110-120
    /// `GetAlgorithmId`: the id of the mask `name`, or the text of NCBI's exception.
    /// The name is its bytes (NCBI's `std::string`), and so is the name in the text.
    pub fn algorithm_id(&self, name: &[u8]) -> Result<i32, Vec<u8>> {
        match self.names.iter().position(|known| known.as_bytes() == name) {
            Some(id) => Ok(id as i32),
            None => Err([
                b"Filtering algorithm ".as_slice(),
                name,
                b" does not exist.\n",
                self.available_names().as_bytes(),
            ]
            .concat()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:169-178
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:115-122
    /// `x_Open` of `id`: the files of the algorithm `id`. An id outside 0 to N-1 is NCBI's text in
    /// [`SeqDbError::Range`] (a search reports it as a database error,
    /// [`SeqDbError::ncbi_message`]).
    pub(super) fn open(&self, id: i32) -> Result<Arc<GiMaskFiles>, SeqDbError> {
        let Some(at) = usize::try_from(id).ok().filter(|&at| at < self.names.len()) else {
            return Err(SeqDbError::Range(format!(
                "Filtering algorithm ID {id} does not exist.\n{}",
                self.available_names()
            )));
        };
        if let Some(files) = self.opened[at].get() {
            return Ok(Arc::clone(files));
        }
        let files = Arc::new(GiMaskFiles::open(&self.names[at], &self.search)?);
        Ok(Arc::clone(self.opened[at].get_or_init(|| files)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:56-63
    /// `GetDesc` of `id`: the stored description of the algorithm `id` and the index file that
    /// holds it. A description that is not UTF-8 is LOSAT's refusal.
    pub(super) fn description(&self, id: i32) -> Result<(String, PathBuf), SeqDbError> {
        let files = self.open(id)?;
        let path = files.index.path().to_path_buf();
        match String::from_utf8(files.desc.clone()) {
            Ok(desc) => Ok((desc, path)),
            Err(_) => Err(SeqDbError::Unsupported {
                path,
                message: "the description of the GI mask is not UTF-8 text".to_string(),
            }),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:257-268
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:57-64
    /// The `.nog` file of the volume `at` (`volume`), opened on first use; `None` where NCBI has
    /// none (no such file, or a volume without sequences).
    fn gi_index(&self, at: usize, volume: &Volume) -> Result<Option<&GiIndex>, SeqDbError> {
        if let Some(index) = self.gi_indexes[at].get() {
            return Ok(index.as_ref());
        }
        let path = with_extension(volume.base(), self.seqtype, "og");
        let index = if super::alias::exists(&path) && volume.num_oids() != 0 {
            Some(GiIndex::open(&path)?)
        } else {
            None
        };
        Ok(self.gi_indexes[at].get_or_init(|| index).as_ref())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:901-915
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:66-70
    // A `.nog` value comes from a `Uint4`, so it is never below 0 and is returned as it is (the
    // -1 that makeblastdb writes for a sequence without a GI reads as 4294967295, a cut-off
    // entry as 0). An OID at or past the `.nog`'s count (none in a file cut inside its header)
    // takes the first GI of its definition lines. The mask lookup takes the GI as a `Uint4`
    // (`GI_TO`), so `INVALID_GI` is 4294967295 too (survey D3).
    /// The key of the sequence `vol_oid` of the volume `at` (`volume`) in the GI masks:
    /// the GI of `x_GetSeqGI` of the OID, as a `Uint4`.
    pub(super) fn seq_gi_key(
        &self,
        at: usize,
        volume: &Volume,
        vol_oid: u32,
    ) -> Result<u32, SeqDbError> {
        if let Some(index) = self.gi_index(at, volume)? {
            if i64::from(vol_oid) < i64::from(index.num_oids) {
                return index.gi(vol_oid);
            }
        }
        let set = volume.defline_set(vol_oid)?;
        let gi = set.ids().find_map(|id| id.gi());
        Ok(gi.map_or(u32::MAX, |gi| gi as u32))
    }
}

#[cfg(test)]
#[path = "gi_mask_tests.rs"]
mod tests;

//! The ISAM files of a volume (`CWriteDB_Isam`, `CWriteDB_IsamIndex`, `CWriteDB_IsamData`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp`) and its OID to GI file
//! (`CWriteDB_GiIndex`, `writedb_volume.hpp`): the string index of the Seq-ids (`.psi`/`.psd`,
//! version 4), the numeric index of the GIs (`.pni`/`.pnd`), the string index of the sequence
//! hashes (`.phi`/`.phd`) and the first GI of every OID (`.pog`); `n` for nucleotide.
//!
//! NCBI keeps the string records in a `CWriteDB_PackedSemiTree` (64 KB blocks of the
//! records less their first 6 bytes, a map from those 6 bytes to a vector of pointers):
//! about the record's bytes + 8 per record. Here a record is a 16-byte entry (its first 8
//! bytes as an ordered integer, its OID, the position of the rest of its key) and the key's
//! bytes past the first 8 in an arena of 1 MiB chunks; the `\x02<oid>\n` end of a record is
//! not stored, it follows from the OID. The entries are kept in chunks of 64 Ki; when the
//! volume closes each chunk is sorted in place and the chunks are merged while the records
//! are written, so no second copy of the entries is made. The chunks (entries and tails)
//! are kept, emptied, for the next volume.

use rustc_hash::FxHashSet;

use super::files::WriteDbFile;
use super::WriteDbError;
use crate::seqid::{label, ParsedSeqId};

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_isam.hpp:156-165
// ```c++
//     enum {
//         eKeyOffset       = 9*4,  ///< Offset of the key offset table.
//         eKeyDelim        = 2,    ///< Byte indicating end of key.
//         eRecordDelim     = 10,   ///< Byte indicating end of data.
//         eMaxStringLine   = 4096, ///< Maximum line size for string.
//         eIsamNumericType = 0,    ///< Numeric ISAM file with Int4 key.
//         eIsamStringType  = 2,    ///< String ISAM file.
//         eIsamNumericLong = 5     ///< Numeric ISAM file with Int8 key.
//     };
// ```
const KEY_OFFSET: i32 = 9 * 4;
const KEY_DELIM: u8 = 2;
const RECORD_DELIM: u8 = 10;
const MAX_STRING_LINE: i32 = 4096;
const ISAM_NUMERIC_TYPE: i32 = 0;
const ISAM_STRING_TYPE: i32 = 2;

/// `EWriteDBIsamType` of the indices makeblastdb builds from FASTA input (`eTrace` is never
/// built: makeblastdb does not set `eAddTrace`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IsamType {
    /// `ePig`: protein identifier groups, never filled from FASTA input.
    Pig,
    /// `eGi`: the numeric index of the GIs.
    Gi,
    /// `eAcc`: the string index of the Seq-ids.
    Acc,
    /// `eHash`: the string index of the sequence hashes.
    Hash,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:56-96
// ```c++
// static string
// s_IsamExtension(EWriteDBIsamType itype,
//                 bool             protein,
//                 bool             is_index)
// {
//     char type_ch = '?';
//
//     switch(itype) {
//     case ePig:
//         type_ch = 'p';
//         break;
//
//     case eGi:
//         type_ch = 'n';
//         break;
//
//     case eAcc:
//         type_ch = 's';
//         break;
//
//     case eTrace:
//         type_ch = 't';
//         break;
//
//     case eHash:
//         type_ch = 'h';
//         break;
//     ...
//     string extn("???");
//     extn[0] = protein ? 'p' : 'n';
//     extn[1] = type_ch;
//     extn[2] = is_index ? 'i' : 'd';
//
//     return extn;
// }
// ```
/// `s_IsamExtension`.
fn isam_extension(itype: IsamType, protein: bool, is_index: bool) -> &'static str {
    let names: [&'static str; 8] = match itype {
        IsamType::Pig => ["ppi", "ppd", "npi", "npd", "", "", "", ""],
        IsamType::Gi => ["pni", "pnd", "nni", "nnd", "", "", "", ""],
        IsamType::Acc => ["psi", "psd", "nsi", "nsd", "", "", "", ""],
        IsamType::Hash => ["phi", "phd", "nhi", "nhd", "", "", "", ""],
    };
    names[usize::from(!protein) * 2 + usize::from(!is_index)]
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_general.hpp:53-56
// ```c++
// inline int s_DivideRoundUp(int value, int blocksize)
// {
//     return ((value + (blocksize - 1)) / blocksize);
// }
// ```
fn divide_round_up(value: i32, blocksize: i32) -> i32 {
    value.wrapping_add(blocksize - 1) / blocksize
}

/// The two files of one ISAM index (`m_IFile` and `m_DFile`), created at the volume's close
/// when the index has entries.
#[derive(Debug)]
struct IsamFiles {
    index: WriteDbFile,
    data: WriteDbFile,
}

impl IsamFiles {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:98-118
    // ```c++
    // CWriteDB_Isam::CWriteDB_Isam(EIsamType      itype,
    //                              const string & dbname,
    //                              bool           protein,
    //                              int            index,
    //                              Uint8          max_file_size,
    //                              bool           sparse)
    // {
    //     m_DFile.Reset(new CWriteDB_IsamData(itype,
    //                                         dbname,
    //                                         protein,
    //                                         index,
    //                                         max_file_size));
    //
    //     m_IFile.Reset(new CWriteDB_IsamIndex(itype,
    //                                          dbname,
    //                                          protein,
    //                                          index,
    //                                          m_DFile,
    //                                          sparse,
    //                                          max_file_size));
    // }
    // ```
    // Both files are `CWriteDB_File`s made with `always_create` false
    // (writedb_isam.cpp:181-186, 771-781).
    fn new(itype: IsamType, dbname: &str, protein: bool, index: u32, max_file_size: u64) -> Self {
        Self {
            index: WriteDbFile::lazy(
                dbname,
                isam_extension(itype, protein, true),
                index,
                max_file_size,
            ),
            data: WriteDbFile::lazy(
                dbname,
                isam_extension(itype, protein, false),
                index,
                max_file_size,
            ),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:220-266
    // ```c++
    //     WriteInt4(isam_version);
    //     WriteInt4(isam_type);
    //     WriteInt4((int)m_DataFileSize);
    //
    //     WriteInt4(num_terms);
    //     WriteInt4(samples);
    //     WriteInt4(m_PageSize);
    //
    //     WriteInt4(max_line_size);
    //     WriteInt4(m_Sparse ? 1 : 0);
    //     WriteInt4(0);
    // ```
    /// `Create()` of both files and `x_WriteHeader`.
    fn create_with_header(
        &mut self,
        isam_type: i32,
        data_file_size: u64,
        num_terms: i32,
        page_size: i32,
        max_line_size: i32,
    ) -> Result<(), WriteDbError> {
        self.index.create_now()?;
        self.data.create_now()?;
        let samples = divide_round_up(num_terms, page_size);
        for value in [
            1,
            isam_type,
            data_file_size as i32,
            num_terms,
            samples,
            page_size,
            max_line_size,
            0,
            0,
        ] {
            self.index.write_int4(value)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:144-152
    // ```c++
    // void CWriteDB_Isam::Close()
    // {
    //     // Index must be closed first, because ISAM indices are built in
    //     // memory until the volume is closed, and only then is anything
    //     // written to disk.
    //
    //     m_IFile->Close();
    //     m_DFile->Close();
    // }
    // ```
    fn close(&mut self) -> Result<(), WriteDbError> {
        self.index.close()?;
        self.data.close()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:154-164
    // ```c++
    // void CWriteDB_Isam::RenameSingle()
    // {
    //     m_IFile->RenameSingle();
    //     m_DFile->RenameSingle();
    // }
    //
    // void CWriteDB_Isam::RenameFileIndex(unsigned int num_digits)
    // {
    //     m_IFile->RenameFileIndex(num_digits);
    //     m_DFile->RenameFileIndex(num_digits);
    // }
    // ```
    fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.index.rename_single()?;
        self.data.rename_single()
    }

    fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        self.index.rename_file_index(num_digits)?;
        self.data.rename_file_index(num_digits)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:819-825
    // ```c++
    // void CWriteDB_Isam::ListFiles(vector<string> & files) const
    // {
    //     if (! m_IFile->Empty()) {
    //         files.push_back(m_IFile->GetFilename());
    //         files.push_back(m_DFile->GetFilename());
    //     }
    // }
    // ```
    // `Empty()` is false once the files were created (writedb_isam.cpp:828-835); the
    // entries not written yet are not on disk, so only created files are listed.
    fn list_files(&self, files: &mut Vec<String>) {
        if self.index.created() {
            files.push(self.index.filename().to_string());
            files.push(self.data.filename().to_string());
        }
    }

    fn abandon(&mut self) {
        self.index.abandon();
        self.data.abandon();
    }
}

/// `CWriteDB_IsamIndex::CanFit` with the file's limit.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:807-810
// ```c++
// bool CWriteDB_IsamIndex::CanFit(int num)
// {
//     return (m_DataFileSize + (num+1) * m_BytesPerElem) < m_MaxFileSize;
// }
// ```
fn can_fit(data_file_size: u64, num: usize, bytes_per_elem: u64, max_file_size: u64) -> bool {
    data_file_size + (num as u64 + 1) * bytes_per_elem < max_file_size
}

/// The bytes of the entries chunk of a string index.
const ENTRY_CHUNK: usize = 1 << 16;
/// The bytes of one arena chunk of key tails.
const TAIL_CHUNK_BITS: u32 = 20;
const TAIL_CHUNK: usize = 1 << TAIL_CHUNK_BITS;
/// The tail field of a key of at most 8 bytes, whose bytes are all in the prefix:
/// `SHORT_KEY + key length` (0 to 8). Arena positions are below it.
const SHORT_KEY: u32 = u32::MAX - 8;
/// The longest key `x_AddStringData` copies into its `char buf[256]`.
pub(crate) const MAX_STRING_KEY: usize = 256;

/// A string record as stored: `[prefix high, prefix low, tail position, OID]`, the prefix
/// being the record's first 8 bytes in an order-preserving form (`record_prefix`).
type Entry = [u32; 4];

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_general.hpp:281-296
// ```c++
//     int Cmp(const CArrayString<STR_SIZE> & other) const
//     {
//         for(int i = 0; i<STR_SIZE; i++) {
//             char ch1 = m_Data[i], ch2 = other.m_Data[i];
//
//             if (ch1 < ch2)
//                 return -1;
//
//             if (ch1 > ch2)
//                 return 1;
//
//             if ((! ch1) && (! ch2))
//                 break;
//         }
//
//         return 0;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_general.cpp:56-78
// ```c++
// void CWriteDB_PackedSemiTree::Insert(const char * x, int L)
// {
//     if (L <= PREFIX) {
//         CArrayString<PREFIX> pre(x, L);
//         CRef<TPacked> & packed = m_Packed[pre];
//         ...
//         packed->Insert("", 0);
//     } else {
//         CArrayString<PREFIX> pre(x, PREFIX);
//         CRef<TPacked> & packed = m_Packed[pre];
//         ...
//         packed->Insert(x + PREFIX, L-PREFIX);
//     }
//     m_Size++;
// }
// ```
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_general.hpp:65-81
// ```c++
// class CWriteDB_PackedStringsCompare {
// public:
//     ...
//     bool operator()(const char * a, const char * b) const
//     {
//         return strcmp(a, b) < 0;
//     }
// };
// ```
// The records are sorted by their first 6 bytes as `char`, signed on the reference platform
// (x86-64: a byte of 0x80 or more sorts before every ASCII byte), then by the rest as
// `strcmp` compares it (unsigned bytes). A record ends with its only line feed and holds no
// NUL, so no record is the start of another and the NUL padding of a short prefix never
// decides. The 8-byte prefix here takes bytes 0-5 with their sign bit flipped (signed order)
// and bytes 6-7 as they are, so integer order is NCBI's order of the first 8 bytes.
/// The order-preserving integer of a record's first 8 bytes (zero-padded).
fn record_prefix(record_head: &[u8]) -> u64 {
    let mut bytes = [0u8; 8];
    let n = record_head.len().min(8);
    bytes[..n].copy_from_slice(&record_head[..n]);
    for byte in &mut bytes[..6] {
        *byte ^= 0x80;
    }
    u64::from_be_bytes(bytes)
}

/// The bytes of a record's first 8 bytes from its prefix integer (with the zero padding).
fn prefix_bytes(prefix: u64) -> [u8; 8] {
    let mut bytes = prefix.to_be_bytes();
    for byte in &mut bytes[..6] {
        *byte ^= 0x80;
    }
    bytes
}

/// Appends `value` in decimal (`sprintf("%d")`).
fn push_decimal(out: &mut Vec<u8>, value: i64) {
    let mut digits = [0u8; 20];
    let mut n = value.unsigned_abs();
    let mut at = digits.len();
    loop {
        at -= 1;
        digits[at] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if value < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(&digits[at..]);
}

/// The decimal digits of `value` (`sprintf("%d")` of an OID) in `buf`.
fn decimal(value: u32, buf: &mut [u8; 10]) -> &[u8] {
    let mut n = value;
    let mut at = buf.len();
    loop {
        at -= 1;
        buf[at] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    &buf[at..]
}

/// The order of the record ends `<x>\n` and `<y>\n` of two OIDs: their digits as `strcmp`
/// compares them (`\n` sorts before every digit, so a shorter number that starts the other
/// sorts first, as with slices).
fn cmp_decimal(x: u32, y: u32) -> std::cmp::Ordering {
    let (mut bx, mut by) = ([0u8; 10], [0u8; 10]);
    decimal(x, &mut bx).cmp(decimal(y, &mut by))
}

/// The merge of sorted chunks of entries (a binary heap of the chunks' next entries, the
/// least on top).
struct Merge<'a> {
    chunks: &'a [Vec<Entry>],
    /// (chunk, position of its next entry).
    heap: Vec<(u32, u32)>,
}

impl<'a> Merge<'a> {
    fn new(
        chunks: &'a [Vec<Entry>],
        mut cmp: impl FnMut(&Entry, &Entry) -> std::cmp::Ordering,
    ) -> Self {
        let mut merge = Self {
            chunks,
            heap: Vec::with_capacity(chunks.len()),
        };
        for (index, chunk) in chunks.iter().enumerate() {
            if !chunk.is_empty() {
                let cursor = (index as u32, 0);
                let mut at = merge.heap.len();
                merge.heap.push(cursor);
                while at > 0 {
                    let parent = (at - 1) / 2;
                    if cmp(merge.entry(cursor), merge.entry(merge.heap[parent])).is_lt() {
                        merge.heap[at] = merge.heap[parent];
                        at = parent;
                    } else {
                        break;
                    }
                }
                merge.heap[at] = cursor;
            }
        }
        merge
    }

    fn entry(&self, (chunk, pos): (u32, u32)) -> &'a Entry {
        &self.chunks[chunk as usize][pos as usize]
    }

    /// The least entry left, or `None` at the end.
    fn next(
        &mut self,
        mut cmp: impl FnMut(&Entry, &Entry) -> std::cmp::Ordering,
    ) -> Option<&'a Entry> {
        let &(chunk, pos) = self.heap.first()?;
        let least = self.entry((chunk, pos));
        // The chunk's next entry, or the heap's last cursor, goes in from the top: the hole
        // sinks along the lesser children to a leaf, then the cursor rises to its place.
        let cursor = if (pos as usize) + 1 < self.chunks[chunk as usize].len() {
            (chunk, pos + 1)
        } else {
            let last = self.heap.pop().expect("not empty");
            if self.heap.is_empty() {
                return Some(least);
            }
            last
        };
        let len = self.heap.len();
        let mut hole = 0;
        loop {
            let left = 2 * hole + 1;
            if left >= len {
                break;
            }
            let child = if left + 1 < len
                && cmp(self.entry(self.heap[left + 1]), self.entry(self.heap[left])).is_lt()
            {
                left + 1
            } else {
                left
            };
            self.heap[hole] = self.heap[child];
            hole = child;
        }
        while hole > 0 {
            let parent = (hole - 1) / 2;
            if cmp(self.entry(cursor), self.entry(self.heap[parent])).is_lt() {
                self.heap[hole] = self.heap[parent];
                hole = parent;
            } else {
                break;
            }
        }
        self.heap[hole] = cursor;
        Some(least)
    }
}

/// The keys of the current OID (`m_OidStringData`, records of one OID differ only in their
/// key): a scan for a few keys, a hash set for many.
#[derive(Debug, Default)]
struct OidKeys {
    bytes: Vec<u8>,
    ends: Vec<usize>,
    set: Option<FxHashSet<Vec<u8>>>,
}

impl OidKeys {
    fn clear(&mut self) {
        self.bytes.clear();
        self.ends.clear();
        self.set = None;
    }

    /// `m_OidStringData.insert(record).second`.
    fn insert(&mut self, key: &[u8]) -> bool {
        if let Some(set) = &mut self.set {
            return set.insert(key.to_vec());
        }
        let mut start = 0;
        for &end in &self.ends {
            if &self.bytes[start..end] == key {
                return false;
            }
            start = end;
        }
        self.bytes.extend_from_slice(key);
        self.ends.push(self.bytes.len());
        if self.ends.len() > 32 {
            let mut set = FxHashSet::default();
            let mut start = 0;
            for &end in &self.ends {
                set.insert(self.bytes[start..end].to_vec());
                start = end;
            }
            self.set = Some(set);
        }
        true
    }
}

/// A string ISAM index (`CWriteDB_Isam` of `eAcc` or `eHash`) of one volume.
#[derive(Debug)]
pub(crate) struct StringIsam {
    files: IsamFiles,
    /// `m_DataFileSize`: the bytes of the records inserted.
    data_file_size: u64,
    /// The entries, in chunks of `ENTRY_CHUNK`.
    entries: Vec<Vec<Entry>>,
    count: usize,
    /// The key bytes past the first 8 of each key longer than 8 bytes: a LEB128 length of
    /// the whole key, then the bytes.
    tails: Vec<Vec<u8>>,
    tail_current: usize,
    /// `m_Oid`.
    oid: i64,
    /// `m_OidStringData`.
    oid_keys: OidKeys,
    /// The lower-cased key being added.
    key: Vec<u8>,
}

/// What a string index keeps between volumes: its entry and tail chunks (emptied at the
/// close), so that the next volume grows into memory the allocator already gave.
#[derive(Debug, Default)]
pub(crate) struct StringIsamScratch {
    entries: Vec<Vec<Entry>>,
    tails: Vec<Vec<u8>>,
}

impl StringIsam {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:166-213
    // ```c++
    // CWriteDB_IsamIndex::CWriteDB_IsamIndex(EWriteDBIsamType        itype,
    //                                        ...
    //     : CWriteDB_File  (dbname,
    //                       s_IsamExtension(itype, protein, true),
    //                       index,
    //                       max_file_size,
    //                       false),
    //       ...
    //       m_Oid          (-1)
    // {
    //     ...
    //     if (itype == eAcc || itype == eHash) {
    //         ...
    //         m_BytesPerElem = 1024;
    //         m_PageSize = 64;
    //     } else {
    //         m_BytesPerElem = 8;
    //         m_PageSize = 256;
    //     }
    // }
    // ```
    pub(crate) fn new(
        itype: IsamType,
        dbname: &str,
        protein: bool,
        index: u32,
        max_file_size: u64,
        scratch: StringIsamScratch,
    ) -> Self {
        let StringIsamScratch {
            mut entries,
            mut tails,
        } = scratch;
        entries.iter_mut().for_each(Vec::clear);
        tails.iter_mut().for_each(Vec::clear);
        Self {
            files: IsamFiles::new(itype, dbname, protein, index, max_file_size),
            data_file_size: 0,
            entries,
            count: 0,
            tails,
            tail_current: 0,
            oid: -1,
            oid_keys: OidKeys::default(),
            key: Vec::new(),
        }
    }

    /// `CanFit(num)`: 1024 bytes for each Seq-id and one more.
    pub(crate) fn can_fit(&self, num: usize) -> bool {
        can_fit(
            self.data_file_size,
            num,
            1024,
            self.files.index.max_file_size(),
        )
    }

    /// The entry chunks and arena chunks, for the next volume.
    pub(crate) fn take_scratch(&mut self) -> StringIsamScratch {
        StringIsamScratch {
            entries: std::mem::take(&mut self.entries),
            tails: std::mem::take(&mut self.tails),
        }
    }

    /// Stores the tail of a key longer than 8 bytes; its position.
    fn push_tail(&mut self, key: &[u8]) -> Result<u32, WriteDbError> {
        let mut head = [0u8; 2];
        let head_len = if key.len() < 0x80 {
            head[0] = key.len() as u8;
            1
        } else {
            head[0] = (key.len() & 0x7f) as u8 | 0x80;
            head[1] = (key.len() >> 7) as u8;
            2
        };
        let size = head_len + key.len() - 8;
        if self.tails.is_empty() {
            self.tails.push(Vec::with_capacity(TAIL_CHUNK));
        }
        if self.tails[self.tail_current].len() + size > TAIL_CHUNK {
            self.tail_current += 1;
            if self.tail_current == self.tails.len() {
                self.tails.push(Vec::with_capacity(TAIL_CHUNK));
            }
        }
        let pos = ((self.tail_current as u64) << TAIL_CHUNK_BITS)
            | self.tails[self.tail_current].len() as u64;
        if pos >= u64::from(SHORT_KEY) {
            return Err(WriteDbError::Creation(
                "the Seq-id keys of one volume's string index exceed 4 GiB".to_string(),
            ));
        }
        let chunk = &mut self.tails[self.tail_current];
        chunk.extend_from_slice(&head[..head_len]);
        chunk.extend_from_slice(&key[8..]);
        Ok(pos as u32)
    }

    fn push_entry(&mut self, entry: Entry) {
        if self.count / ENTRY_CHUNK == self.entries.len() {
            self.entries.push(Vec::with_capacity(ENTRY_CHUNK));
        }
        self.entries[self.count / ENTRY_CHUNK].push(entry);
        self.count += 1;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:721-769
    // ```c++
    // void CWriteDB_IsamIndex::x_AddStringData(int oid, const char * sbuf, int ssize)
    // {
    //     ...
    //     char buf[256];
    //
    //     int sz = ssize;
    //     memcpy(buf, sbuf, sz);
    //     _ASSERT(sz);
    //
    //     // lowercase the 'key' portion
    //     for(int i = 0; i < sz; i++) {
    //         buf[i] = tolower(buf[i]);
    //     }
    //
    //     buf[sz++] = (char) eKeyDelim;
    //     sz += sprintf(buf + sz, "%d", oid);
    //     buf[sz++] = (char) eRecordDelim;
    //
    //     // fix for SB-218, SB-819
    //     if (oid != m_Oid) {
    //         m_Oid = oid;
    //         m_OidStringData.clear();
    //     }
    //
    //     string tmp(buf, sz);
    //     pair< set<string>::iterator, bool> rv = m_OidStringData.insert(tmp);
    //     if (rv.second) {
    //         m_StringSort.Insert(buf, sz);
    //         m_DataFileSize += sz;
    //     }
    // }
    // ```
    // `tolower` of the C locale changes 'A'-'Z' only. A key of more than 256 bytes does not
    // fit `buf`: makeblastdb 2.17.0 stops with "*** buffer overflow detected ***" (SIGABRT,
    // the fortified `memcpy`); LOSAT refuses it with an error naming the key (makedb-v5 S3b
    // handoff, decision candidate; no NCBI crash reproduced). Keys of up to 256 bytes give
    // the records makeblastdb writes (probe of 252-256 bytes, s3b/probe).
    /// `x_AddStringData`: adds the record of `key` for `oid`, once per OID.
    pub(crate) fn add_string_data(&mut self, oid: u32, key: &[u8]) -> Result<(), WriteDbError> {
        if key.len() > MAX_STRING_KEY {
            return Err(WriteDbError::Creation(format!(
                "Seq-id key '{}' of the sequence with OID {oid} of its volume is {} bytes long; \
                 makeblastdb's string index holds keys of at most {MAX_STRING_KEY} bytes \
                 (makeblastdb 2.17.0 stops with a buffer overflow)",
                String::from_utf8_lossy(key),
                key.len()
            )));
        }
        let mut lower = std::mem::take(&mut self.key);
        lower.clear();
        lower.extend(key.iter().map(u8::to_ascii_lowercase));
        if i64::from(oid) != self.oid {
            self.oid = i64::from(oid);
            self.oid_keys.clear();
        }
        let result = if self.oid_keys.insert(&lower) {
            self.insert_record(oid, &lower)
        } else {
            Ok(())
        };
        self.key = lower;
        result
    }

    /// `m_StringSort.Insert(record)` and `m_DataFileSize += record size`.
    fn insert_record(&mut self, oid: u32, key: &[u8]) -> Result<(), WriteDbError> {
        let mut head = Vec::with_capacity(8 + 12);
        head.extend_from_slice(&key[..key.len().min(8)]);
        if key.len() < 8 {
            head.push(KEY_DELIM);
            push_decimal(&mut head, i64::from(oid));
            head.push(RECORD_DELIM);
        }
        let prefix = record_prefix(&head);
        let tail = if key.len() > 8 {
            self.push_tail(key)?
        } else {
            SHORT_KEY + key.len() as u32
        };
        self.push_entry([(prefix >> 32) as u32, prefix as u32, tail, oid]);
        let mut digits = Vec::with_capacity(10);
        push_decimal(&mut digits, i64::from(oid));
        self.data_file_size += (key.len() + 2 + digits.len()) as u64;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:503-509
    // ```c++
    // void CWriteDB_IsamIndex::AddHash(int oid, int hash)
    // {
    //     char buf[256];
    //     int sz = sprintf(buf, "%u", (unsigned)hash);
    //
    //     x_AddStringData(oid, buf, sz);
    // }
    // ```
    /// `AddHash`: the sequence hash of `oid` as a decimal key.
    pub(crate) fn add_hash(&mut self, oid: u32, hash: u32) -> Result<(), WriteDbError> {
        let mut key = Vec::with_capacity(10);
        push_decimal(&mut key, i64::from(hash));
        self.add_string_data(oid, &key)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:556-603
    // ```c++
    // void CWriteDB_IsamIndex::x_AddStringIds(int oid, const TIdList & idlist)
    // {
    //     // Build all sub-string objects and add those.
    //
    //     ITERATE(TIdList, iter, idlist) {
    //         const CSeq_id & seqid = **iter;
    //
    //         switch(seqid.Which()) {
    //
    //         case CSeq_id::e_Gi:
    //             break;
    //
    //         case CSeq_id::e_Pdb:
    //             x_AddPdb(oid, seqid);
    //             break;
    //
    //         case CSeq_id::e_Local:
    //             x_AddLocal(oid, seqid);
    //             break;
    //
    //         case CSeq_id::e_Patent:
    //             x_AddPatent(oid, seqid);
    //             break;
    //
    //         case CSeq_id::e_General:
    //             if (! m_Sparse) {
    //                 x_AddStdString(oid, seqid.AsFastaString());
    //                 const CDbtag & dbt = seqid.GetGeneral();
    //                 if (dbt.CanGetTag() && dbt.GetTag().IsStr()) {
    //                     x_AddStdString(oid, dbt.GetTag().GetStr());
    //                 }
    //             }
    //             break;
    //
    //         // default processing:
    //         default:
    //             {
    //                 const CTextseq_id * textid  = seqid.GetTextseq_Id();
    //                 if (textid) {
    //                     x_AddTextId(oid, *textid);
    //                 } else {
    //                     string acc = seqid.AsFastaString();
    //                     x_AddStringData(oid, acc.data(), acc.size());
    //                 }
    //             }
    //         }
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:605-624
    // ```c++
    // void CWriteDB_IsamIndex::x_AddLocal(int             oid,
    //                                     const CSeq_id & seqid)
    // {
    //     const CObject_id & objid = seqid.GetLocal();
    //
    //     if (! m_Sparse) {
    //         x_AddStdString(oid, seqid.AsFastaString());
    //     }
    //     if (objid.IsStr()) {
    //         x_AddStdString(oid, objid.GetStr());
    //     }
    // }
    //
    // void CWriteDB_IsamIndex::x_AddPatent(int             oid,
    //                                      const CSeq_id & seqid)
    // {
    //     if (! m_Sparse) {
    //         x_AddStdString(oid, seqid.AsFastaString());
    //     }
    // }
    // ```
    // makeblastdb never sets `eSparseIndex`, so `m_Sparse` is false.
    /// `AddIds` of the string index of the Seq-ids (`x_AddStringIds`).
    pub(crate) fn add_ids(&mut self, oid: u32, ids: &[ParsedSeqId]) -> Result<(), WriteDbError> {
        let mut fasta = Vec::new();
        for id in ids {
            fasta.clear();
            match id {
                ParsedSeqId::Gi(_) => {}
                ParsedSeqId::Pdb(_) => self.add_pdb(oid, id)?,
                ParsedSeqId::Local(tag) => {
                    label::fasta(id, &mut fasta);
                    self.add_string_data(oid, &fasta)?;
                    if let crate::seqid::model::ObjectId::Str(text) = tag {
                        self.add_string_data(oid, text)?;
                    }
                }
                ParsedSeqId::Patent(_) => {
                    label::fasta(id, &mut fasta);
                    self.add_string_data(oid, &fasta)?;
                }
                ParsedSeqId::General(dbtag) => {
                    label::fasta(id, &mut fasta);
                    self.add_string_data(oid, &fasta)?;
                    if let crate::seqid::model::ObjectId::Str(text) = &dbtag.tag {
                        self.add_string_data(oid, text)?;
                    }
                }
                _ => match id.textseq() {
                    Some(text) => self.add_text_id(oid, text)?,
                    None => {
                        label::fasta(id, &mut fasta);
                        self.add_string_data(oid, &fasta)?;
                    }
                },
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:626-673
    // ```c++
    // void CWriteDB_IsamIndex::x_AddPdb(int             oid,
    //                                   const CSeq_id & seqid)
    // {
    //     const CPDB_seq_id & pdb = seqid.GetPdb();
    //     ...
    //     CTempString mol;
    //     if (pdb.CanGetMol()) {
    //         mol = pdb.GetMol().Get();
    //     }
    //     if (! mol.size()) {
    //         NCBI_THROW(CWriteDBException,
    //                    eArgErr,
    //                    "Empty molecule string in pdb Seq-id.");
    //     }
    //     x_AddStringData(oid, mol);
    //
    //     string full_id = seqid.AsFastaString();
    //     _ASSERT(full_id.size() > 4);
    //     if (! m_Sparse) {
    //         x_AddStdString(oid, full_id);
    //     }
    //
    //     string short_id(full_id, 4);
    //     x_AddStdString(oid, short_id);
    //
    //     // mol is always 4 characters long. If chain exists, add a sting in which the
    //     // '|' delimiter is changed to a space.
    //     if (short_id[4] == '|')
    //         short_id[4] = ' ';
    //     x_AddStdString(oid, short_id);
    // }
    // ```
    // The parser takes a mol of 4 characters or more (shorter ones make "Could not construct
    // seq-id"), so `short_id` has a fifth byte: the '|' after a 4-character mol, a letter of a
    // longer one (then the same key is added twice and kept once).
    fn add_pdb(&mut self, oid: u32, id: &ParsedSeqId) -> Result<(), WriteDbError> {
        let ParsedSeqId::Pdb(pdb) = id else {
            unreachable!("a pdb id");
        };
        if pdb.mol.is_empty() {
            return Err(WriteDbError::Creation(
                "Empty molecule string in pdb Seq-id.".to_string(),
            ));
        }
        self.add_string_data(oid, &pdb.mol)?;
        let mut full_id = Vec::new();
        label::fasta(id, &mut full_id);
        self.add_string_data(oid, &full_id)?;
        let mut short_id = full_id[4.min(full_id.len())..].to_vec();
        self.add_string_data(oid, &short_id)?;
        if short_id.get(4) == Some(&b'|') {
            short_id[4] = b' ';
        }
        self.add_string_data(oid, &short_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:675-719
    // ```c++
    // bool s_NoCaseEqual(CTempString & a, CTempString & b)
    // {
    //     if (a.size() != b.size())
    //         return false;
    //
    //     return 0 == NStr::strncasecmp(a.data(), b.data(), a.size());
    // }
    //
    // void CWriteDB_IsamIndex::x_AddTextId(int                 oid,
    //                                      const CTextseq_id & id)
    // {
    //     CTempString acc, nm;
    //     ...
    //     if (id.CanGetAccession()) {
    //         acc = id.GetAccession();
    //     }
    //
    //     if (id.CanGetName()) {
    //         nm = id.GetName();
    //     }
    //
    //     if (! acc.empty()) {
    //         x_AddStringData(oid, acc);
    //     }
    //
    //     if (! m_Sparse) {
    //         // Skip name if it is empty or if it is the same as 'acc' when
    //         // case is ignored.
    //
    //         if (! (nm.empty() || s_NoCaseEqual(acc, nm))) {
    //             x_AddStringData(oid, nm);
    //         }
    //
    //         int ver = id.CanGetVersion() ? id.GetVersion() : 0;
    //
    //         if (ver && acc.size()) {
    //             x_AddString(oid, acc, ver);
    //         }
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:771-784
    // ```c++
    // void CWriteDB_IsamIndex::x_AddString(int oid, const CTempString & acc, int ver)
    // {
    //     _ASSERT(! m_Sparse);
    //
    //     if (acc.size() && ver) {
    //         char buf[256];
    //         memcpy(buf, acc.data(), acc.size());
    //
    //         int sz = acc.size();
    //         sz += sprintf(buf + sz, ".%d", ver);
    //
    //         x_AddStringData(oid, buf, sz);
    //     }
    // }
    // ```
    // The accession of a parsed id has at most 30 characters (`CFastaIdValidate`), so `buf`
    // holds `acc.ver`.
    fn add_text_id(
        &mut self,
        oid: u32,
        text: &crate::seqid::model::TextseqId,
    ) -> Result<(), WriteDbError> {
        let acc: &[u8] = text.accession.as_deref().unwrap_or(b"");
        let nm: &[u8] = text.name.as_deref().unwrap_or(b"");
        if !acc.is_empty() {
            self.add_string_data(oid, acc)?;
        }
        if !(nm.is_empty() || acc.eq_ignore_ascii_case(nm)) {
            self.add_string_data(oid, nm)?;
        }
        let ver = text.version.unwrap_or(0);
        if ver != 0 && !acc.is_empty() {
            let mut key = acc.to_vec();
            key.push(b'.');
            push_decimal(&mut key, i64::from(ver));
            self.add_string_data(oid, &key)?;
        }
        Ok(())
    }

    /// The length of the key of `entry` and its bytes past the first 8 (empty for a key of
    /// at most 8 bytes, which lies in the prefix).
    fn key_tail(&self, entry: &Entry) -> (usize, &[u8]) {
        if entry[2] >= SHORT_KEY {
            return ((entry[2] - SHORT_KEY) as usize, &[]);
        }
        let pos = entry[2];
        let chunk = &self.tails[(pos >> TAIL_CHUNK_BITS) as usize];
        let mut at = (pos as usize) & (TAIL_CHUNK - 1);
        let mut key_len = usize::from(chunk[at] & 0x7f);
        if chunk[at] & 0x80 != 0 {
            at += 1;
            key_len |= usize::from(chunk[at]) << 7;
        }
        at += 1;
        (key_len, &chunk[at..at + key_len - 8])
    }

    /// NCBI's order of two records (`key\x02<oid>\n`, the order described at
    /// `record_prefix`), without building them: the 8-byte prefixes, then for equal ones the
    /// key bytes past the first 8 as they are stored, and the OIDs' decimal digits only for
    /// equal keys. Two records with equal first 8 bytes compare from byte 8 as they compare
    /// in full, so the comparison may run over the whole keys. A key that is a proper start
    /// of the other continues with `\x02`; when the other key holds `\x02` at that place
    /// (an id with `|` may), the records are built and compared.
    fn cmp_entries(&self, x: &Entry, y: &Entry) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let px = (u64::from(x[0]) << 32) | u64::from(x[1]);
        let py = (u64::from(y[0]) << 32) | u64::from(y[1]);
        match px.cmp(&py) {
            Ordering::Equal => {}
            order => return order,
        }
        let (lx, tx) = self.key_tail(x);
        let (ly, ty) = self.key_tail(y);
        if lx == ly && lx <= 8 {
            // The same key, wholly in the equal prefixes.
            return cmp_decimal(x[3], y[3]);
        }
        if lx >= 8 && ly >= 8 {
            // The first 8 key bytes are the equal prefixes.
            let n = tx.len().min(ty.len());
            match tx[..n].cmp(&ty[..n]) {
                Ordering::Equal => {}
                order => return order,
            }
            if tx.len() == ty.len() {
                return cmp_decimal(x[3], y[3]);
            }
            // `\x02` of the shorter record against the longer key's next byte.
            let (next, shorter) = if tx.len() < ty.len() {
                (ty[n], Ordering::Less)
            } else {
                (tx[n], Ordering::Greater)
            };
            if next != KEY_DELIM {
                return if next > KEY_DELIM {
                    shorter
                } else {
                    shorter.reverse()
                };
            }
        }
        // A key that holds `\x02` where the other one ends: compare the records.
        let (mut a, mut b) = (Vec::new(), Vec::new());
        self.record(x, &mut a);
        self.record(y, &mut b);
        a.cmp(&b)
    }

    /// The full record of `entry` (`key\x02<oid>\n`) into `out`.
    fn record(&self, entry: &Entry, out: &mut Vec<u8>) {
        out.clear();
        let head = prefix_bytes((u64::from(entry[0]) << 32) | u64::from(entry[1]));
        if entry[2] >= SHORT_KEY {
            // The key is in the prefix (a key may hold any byte but the space, so its
            // length is stored, not found).
            out.extend_from_slice(&head[..(entry[2] - SHORT_KEY) as usize]);
        } else {
            let pos = entry[2];
            let chunk = &self.tails[(pos >> TAIL_CHUNK_BITS) as usize];
            let mut at = (pos as usize) & (TAIL_CHUNK - 1);
            let mut key_len = usize::from(chunk[at] & 0x7f);
            if chunk[at] & 0x80 != 0 {
                at += 1;
                key_len |= usize::from(chunk[at]) << 7;
            }
            at += 1;
            out.extend_from_slice(&head);
            out.extend_from_slice(&chunk[at..at + key_len - 8]);
        }
        out.push(KEY_DELIM);
        push_decimal(out, i64::from(entry[3]));
        out.push(RECORD_DELIM);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:268-381
    // ```c++
    // void CWriteDB_IsamIndex::x_FlushStringIndex()
    // {
    //     ...
    //     int data_pos = 0;
    //     unsigned count = m_StringSort.Size();
    //
    //     unsigned nsamples = s_DivideRoundUp(count, m_PageSize);
    //
    //     string key_buffer;
    //     vector<int> key_off;
    //     ...
    //     int output_count = 0;
    //     int index = 0;
    //
    //     m_StringSort.Sort();
    //     ...
    //     string element, prev_elem;
    //     ...
    //     element.resize(1);
    //     element[0] = char(0);
    //
    //     while(iter != end_iter) {
    //         prev_elem.swap(element);
    //         iter.Get(element);
    //
    //         if (prev_elem == element) {
    //             ++iter;
    //             continue;
    //         }
    //         ...
    //         if ((output_count & (m_PageSize-1)) == 0) {
    //             // Write the data file position to the index file.
    //
    //             WriteInt4(data_pos);
    //             ...
    //             key_off.push_back((int) key_buffer.size());
    //             ...
    //             key_buffer.append(element.data(), element.length()-1);
    //             key_buffer.append(NUL);
    //         }
    //         output_count ++;
    //
    //         data_pos = m_DataFile->Write(element);
    //         index ++;
    //
    //         ++iter;
    //     }
    //
    //     // Write the final data position.
    //
    //     WriteInt4(data_pos);
    //
    //     // Push back the final buffer offset.
    //
    //     key_off.push_back((int) key_buffer.size());
    //
    //     int key_off_start = eKeyOffset + (nsamples + 1) * 8;
    //
    //     // Write index file offsets of keys.
    //
    //     for(i = 0; i < key_off.size(); i++) {
    //         WriteInt4(key_off[i] + key_off_start);
    //     }
    //
    //     // Write buffer of keys.
    //
    //     Write(key_buffer);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:457-479
    // ```c++
    // void CWriteDB_IsamIndex::x_Flush()
    // {
    //     if (m_NumberTable.size() || m_StringSort.Size()) {
    //         Create();
    //         m_DataFile->Create();
    //
    //         // Step 1: Write header data.
    //         x_WriteHeader();
    //         ...
    //         if (m_Type == eAcc || m_Type == eHash) {
    //             x_FlushStringIndex();
    //         } else {
    //             x_FlushNumericIndex();
    //         }
    //     }
    //
    //     x_Free();
    // }
    // ```
    /// `Close`: writes both files when the index has records.
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if self.count != 0 {
            let num_terms = self.count as i32;
            self.files.create_with_header(
                ISAM_STRING_TYPE,
                self.data_file_size,
                num_terms,
                64,
                MAX_STRING_LINE,
            )?;
            // Each chunk sorted in place, then the chunks merged while the records are
            // written (`m_StringSort.Sort()` and the walk of `x_FlushStringIndex`): no
            // second copy of the entries.
            let chunks = (self.count + ENTRY_CHUNK - 1) / ENTRY_CHUNK;
            let mut entries = std::mem::take(&mut self.entries);
            for chunk in &mut entries[..chunks] {
                chunk.sort_unstable_by(|x, y| self.cmp_entries(x, y));
            }
            let nsamples = divide_round_up(num_terms, 64);
            let mut key_buffer: Vec<u8> = Vec::new();
            let mut key_off: Vec<i32> = Vec::with_capacity(nsamples.max(0) as usize + 1);
            let mut data_pos: u32 = 0;
            let mut output_count: u32 = 0;
            let mut element = Vec::with_capacity(300);
            let mut prev_elem: Vec<u8> = vec![0];
            let mut merge = Merge::new(&entries[..chunks], |x, y| self.cmp_entries(x, y));
            while let Some(entry) = merge.next(|x, y| self.cmp_entries(x, y)) {
                self.record(entry, &mut element);
                if prev_elem == element {
                    continue;
                }
                if output_count & 63 == 0 {
                    self.files.index.write_int4(data_pos as i32)?;
                    key_off.push(key_buffer.len() as i32);
                    key_buffer.extend_from_slice(&element[..element.len() - 1]);
                    key_buffer.push(0);
                }
                output_count += 1;
                data_pos = self.files.data.write(&element)?;
                std::mem::swap(&mut prev_elem, &mut element);
            }
            entries.iter_mut().for_each(Vec::clear);
            self.entries = entries;
            self.files.index.write_int4(data_pos as i32)?;
            key_off.push(key_buffer.len() as i32);
            let key_off_start = KEY_OFFSET + (nsamples + 1) * 8;
            for off in key_off {
                self.files
                    .index
                    .write_int4(off.wrapping_add(key_off_start))?;
            }
            self.files.index.write(&key_buffer)?;
        }
        // `x_Free`: the records are dropped; the chunks stay for the next volume.
        self.count = 0;
        self.tails.iter_mut().for_each(Vec::clear);
        self.tail_current = 0;
        self.files.close()
    }

    pub(crate) fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.files.rename_single()
    }

    pub(crate) fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        self.files.rename_file_index(num_digits)
    }

    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        self.files.list_files(files);
    }

    pub(crate) fn abandon(&mut self) {
        self.files.abandon();
    }

    /// The bytes of the records held (entries and tails, without the free ends of the last
    /// chunks), for the memory figures.
    #[cfg(test)]
    pub(crate) fn record_bytes(&self) -> usize {
        self.count * 16 + self.tails.iter().map(Vec::len).sum::<usize>()
    }
}

/// A numeric ISAM index (`CWriteDB_Isam` of `eGi`, or of `ePig`, which FASTA input never
/// fills) of one volume.
#[derive(Debug)]
pub(crate) struct NumericIsam {
    files: IsamFiles,
    /// `m_DataFileSize`.
    data_file_size: u64,
    /// `m_NumberTable`: (id, OID).
    table: Vec<(i64, i32)>,
}

impl NumericIsam {
    pub(crate) fn new(
        itype: IsamType,
        dbname: &str,
        protein: bool,
        index: u32,
        max_file_size: u64,
    ) -> Self {
        Self {
            files: IsamFiles::new(itype, dbname, protein, index, max_file_size),
            data_file_size: 0,
            table: Vec::new(),
        }
    }

    /// `CanFit(num)`: 8 bytes for each id and one more.
    pub(crate) fn can_fit(&self, num: usize) -> bool {
        can_fit(
            self.data_file_size,
            num,
            8,
            self.files.index.max_file_size(),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:511-522
    // ```c++
    // void CWriteDB_IsamIndex::x_AddGis(int oid, const TIdList & idlist)
    // {
    //     ITERATE(TIdList, iter, idlist) {
    //         const CSeq_id & seqid = **iter;
    //
    //         if (seqid.IsGi()) {
    //             SIdOid row(GI_TO(Int8, seqid.GetGi()), oid);
    //             m_NumberTable.push_back(row);
    //             m_DataFileSize += 8;
    //         }
    //     }
    // }
    // ```
    /// `AddIds` of the numeric index of the GIs (`x_AddGis`).
    pub(crate) fn add_gis(&mut self, oid: u32, ids: &[ParsedSeqId]) {
        for gi in ids.iter().filter_map(ParsedSeqId::gi) {
            self.table.push((gi, oid as i32));
            self.data_file_size += 8;
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_isam.cpp:383-455
    // ```c++
    // void CWriteDB_IsamIndex::x_FlushNumericIndex()
    // {
    //     _ASSERT(m_NumberTable.size());
    //
    //     int row_index = 0;
    //
    //     sort(m_NumberTable.begin(), m_NumberTable.end());
    //
    //     int count = (int) m_NumberTable.size();
    //
    //     const SIdOid * prevp = 0;
    //     ...
    //     } else {
    //         for(int i = 0; i < count; i++) {
    //             const SIdOid & elem = m_NumberTable[i];
    //
    //             if (prevp && (*prevp == elem)) {
    //                 continue;
    //             } else {
    //                 prevp = & elem;
    //             }
    //
    //             if ((row_index & (m_PageSize-1)) == 0) {
    //                 WriteInt4(elem.id());
    //                 WriteInt4(elem.oid());
    //             }
    //
    //             m_DataFile->WriteInt4(elem.id());
    //             m_DataFile->WriteInt4(elem.oid());
    //             row_index ++;
    //         }
    //
    //         // 32 bit numeric files end in (max-uint4, 0).
    //
    //         WriteInt4(-1);
    //         WriteInt4(0);
    //     }
    // }
    // ```
    // `m_UseInt8` is set by trace ids only, which makeblastdb does not index. The rows are
    // sorted as (Int8 id, OID) pairs and the id is written as an `Int4`: a GI of 2^31 or more
    // comes after the smaller ones and is written negative. The header counts the rows
    // before the duplicates are skipped (`x_WriteHeader`: `num_terms` is the table's size).
    /// `Close`: writes both files when the index has rows.
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if !self.table.is_empty() {
            self.files.create_with_header(
                ISAM_NUMERIC_TYPE,
                self.data_file_size,
                self.table.len() as i32,
                256,
                0,
            )?;
            self.table.sort_unstable();
            let mut prev: Option<(i64, i32)> = None;
            let mut row_index: u32 = 0;
            for &(id, oid) in &self.table {
                if prev == Some((id, oid)) {
                    continue;
                }
                prev = Some((id, oid));
                if row_index & 255 == 0 {
                    self.files.index.write_int4(id as i32)?;
                    self.files.index.write_int4(oid)?;
                }
                self.files.data.write_int4(id as i32)?;
                self.files.data.write_int4(oid)?;
                row_index += 1;
            }
            self.files.index.write_int4(-1)?;
            self.files.index.write_int4(0)?;
        }
        // `x_Free`: `vector<SIdOid> tmp; m_NumberTable.swap(tmp);`.
        self.table = Vec::new();
        self.files.close()
    }

    pub(crate) fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.files.rename_single()
    }

    pub(crate) fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        self.files.rename_file_index(num_digits)
    }

    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        self.files.list_files(files);
    }

    pub(crate) fn abandon(&mut self) {
        self.files.abandon();
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.hpp:53-98
// ```c++
// class CWriteDB_GiIndex : public CWriteDB_File {
// public:
//     CWriteDB_GiIndex(const string & dbname,
//                      bool           protein,
//                      int            index,
//                      Uint8          max_fsize)
//     : CWriteDB_File  (dbname, (protein ? "pog" : "nog"), index, max_fsize, false){ }
//     ...
//     void AddGi(TGi gi) {
//         m_Gi.push_back(gi);
//     }
//
// private:
//     void x_Flush() {
//
//         Int4 nGi = m_Gi.size();
//
//         if (! nGi) return;
//
//         Create();
//         WriteInt4(kVersion);
//         WriteInt4(kFileType);
//         WriteInt4(kGiSize);
//         WriteInt4(nGi);
//
//         for (Int4 i=0; i<4; i++) {
//             WriteInt4(0);
//         }
//
//         for (Int4 i=0; i<nGi; i++) {
//             WriteInt4(GI_TO(Int4, m_Gi[i]));
//         }
//
//         vector<TGi> tmp;
//         m_Gi.swap(tmp);
//     }
//
//     static const int kVersion = 1;
//     static const int kFileType = 0;
//     static const int kGiSize = 4;
//     vector<TGi> m_Gi;
// };
// ```
// The GIs are kept as the `Int4` they are written as (4 bytes per OID instead of 8).
/// `CWriteDB_GiIndex`: the `.pog` / `.nog` file of a volume.
#[derive(Debug)]
pub(crate) struct GiIndex {
    file: WriteDbFile,
    gis: Vec<i32>,
}

impl GiIndex {
    pub(crate) fn new(dbname: &str, protein: bool, index: u32, max_file_size: u64) -> Self {
        Self {
            file: WriteDbFile::lazy(
                dbname,
                if protein { "pog" } else { "nog" },
                index,
                max_file_size,
            ),
            gis: Vec::new(),
        }
    }

    /// `AddGi`, truncated to an `Int4`.
    pub(crate) fn add_gi(&mut self, gi: i64) {
        self.gis.push(gi as i32);
    }

    /// `Close` (`x_Flush`).
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if !self.gis.is_empty() {
            self.file.create_now()?;
            for value in [1, 0, 4, self.gis.len() as i32, 0, 0, 0, 0] {
                self.file.write_int4(value)?;
            }
            let mut buffer = Vec::with_capacity(1 << 14);
            for chunk in self.gis.chunks(1 << 12) {
                buffer.clear();
                for gi in chunk {
                    buffer.extend_from_slice(&gi.to_be_bytes());
                }
                self.file.write(&buffer)?;
            }
            self.gis = Vec::new();
        }
        self.file.close()
    }

    pub(crate) fn file_mut(&mut self) -> &mut WriteDbFile {
        &mut self.file
    }

    pub(crate) fn filename(&self) -> &str {
        self.file.filename()
    }
}

#[cfg(test)]
#[path = "isam_tests.rs"]
mod tests;

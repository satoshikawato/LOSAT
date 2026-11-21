//! Positioned reads of one database file.
//!
//! NCBI maps every database file into memory and hands out pointers into the
//! mapping. LOSAT copies the requested byte range with a positioned read
//! instead (`std` only, no memory mapping), so that the same reader code
//! serves native builds and WebAssembly. Only this file differs by platform.

use std::fs::File;
use std::path::{Path, PathBuf};

use super::SeqDbError;

/// The memory page of NCBI's mapping of a database file: 4 KiB, the page of the oracle platform
/// (Linux x86-64) and the smallest of the systems NCBI BLAST+ runs on. NCBI maps a file whole
/// from offset 0, and the system fills the part of the last page after the end of the file with
/// zeros (POSIX `mmap`), so NCBI reads zero bytes there (`taxdb.btd` past its end). A larger page
/// only gives NCBI more zeros to read past the file, where LOSAT then refuses or reads none.
pub(crate) const MAP_PAGE: u64 = 4096;

/// The bytes that NCBI's mapping of a file of `len` bytes holds (`MAP_PAGE`): whole pages, the
/// part of the last one after the end of the file zeros. The one form of `taxdb.btd` and
/// RPS-BLAST's `.rps` (`rps::RpsVolumeFiles::row_after_rows`).
pub(crate) fn mapped_length(len: u64) -> u64 {
    len.div_ceil(MAP_PAGE).saturating_mul(MAP_PAGE)
}

/// Where `DbFile::read_range_in` puts the bytes it reads.
enum RangeBuffer<'a> {
    /// A new buffer (`read_range`).
    New,
    /// EXPERIMENT (LOSAT_X_BNMAPPED): `x_read_range_mapped`'s buffer.
    Mapped,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): the worker's buffer (`x_read_range_into`).
    Given(&'a mut Vec<u8>),
}

/// One open database file (`.pin`, `.psq`, `.phr`, ... ) read by offset.
#[derive(Debug)]
pub struct DbFile {
    path: PathBuf,
    len: u64,
    #[cfg(any(unix, windows))]
    file: File,
    #[cfg(not(any(unix, windows)))]
    file: std::sync::Mutex<File>,
}

impl DbFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:123-129
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:93-105
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5932-5946
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:709-720
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5600-5656
    // NCBI first asks whether the file exists (`GetFileSizeL`: `stat`
    // succeeds and names a regular file, cited at `alias::exists`) and gives
    // "not found" otherwise; it then opens the file to map it, and an `open`
    // that fails is its memory-mapping error. LOSAT takes the same two steps
    // with a plain `open` in place of the mapping. NCBI maps the whole file
    // there, and a mapping that fails is an error of that place too: LOSAT
    // reads the first byte of a file that is not empty (`opened`), so a file
    // whose bytes the system does not give fails where NCBI's fails (line WQA
    // D18, 2026-11-21).
    /// Opens `path` for reading. A path that is no regular file is NCBI's
    /// "Error: File (...) not found."; a regular file that cannot be opened is
    /// [`SeqDbError::Open`], one whose first byte cannot be read
    /// [`SeqDbError::Map`].
    pub fn open(path: &Path) -> Result<Self, SeqDbError> {
        if !super::alias::exists(path) {
            return Err(SeqDbError::NotFound(
                format!("Error: File ({}) not found.", path.display()).into(),
            ));
        }
        let file = File::open(path).map_err(|source| SeqDbError::Open {
            path: path.to_path_buf(),
            source,
        })?;
        let len = file
            .metadata()
            .map_err(|source| SeqDbError::Io {
                path: path.to_path_buf(),
                source,
            })?
            .len();
        #[cfg(all(test, target_os = "linux"))]
        let file = refused_reads::handle(path, file);
        Self::opened(path, file, len)
    }

    /// The file `path` of `len` bytes, open on `file`: the one constructor (`open`, and the
    /// tests' `with_handle`). It reads the first byte of a file that is not empty, as NCBI maps
    /// the whole file where it opens it; a read that the system refuses there is
    /// [`SeqDbError::Map`], NCBI's failed mapping of the file.
    fn opened(path: &Path, file: File, len: u64) -> Result<Self, SeqDbError> {
        let opened = Self {
            path: path.to_path_buf(),
            len,
            #[cfg(any(unix, windows))]
            file,
            #[cfg(not(any(unix, windows)))]
            file: std::sync::Mutex::new(file),
        };
        if len > 0 {
            opened
                .read_exact_at_impl(&mut [0u8; 1], 0)
                .map_err(|source| SeqDbError::Map {
                    path: path.to_path_buf(),
                    len,
                    source,
                })?;
        }
        Ok(opened)
    }

    /// The file `path` of `len` bytes open on the handle `file`, through the constructor of
    /// `open` (the tests: a handle whose reads the system refuses, a directory's).
    #[cfg(all(test, target_os = "linux"))]
    pub(crate) fn with_handle(path: &Path, file: File, len: u64) -> Result<Self, SeqDbError> {
        Self::opened(path, file, len)
    }

    /// Path this file was opened with.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// File length in bytes at the time it was opened.
    pub fn len(&self) -> u64 {
        self.len
    }

    /// Whether the file held no bytes when it was opened.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1802-1817
    // NCBI returns a pointer at `offset` of the mapped file; LOSAT fills `buf`
    // from the same offset. A range that passes the end of the file, which
    // NCBI would read past its mapping, is an explicit format error. A read
    // that the system refuses here, after the file's first byte was read where
    // it was opened, is where NCBI reads mapped bytes that the system no
    // longer gives and ends by SIGBUS: LOSAT's own error,
    // [`SeqDbError::Read`] ("cannot be read"; line WQA D18, 2026-11-21).
    /// Fills `buf` with the bytes at `offset`.
    pub fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), SeqDbError> {
        let end = offset.checked_add(buf.len() as u64);
        if end.is_none_or(|end| end > self.len) {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!(
                    "byte range {}..{} passes the end of the file ({} bytes)",
                    offset,
                    offset.saturating_add(buf.len() as u64),
                    self.len
                ),
            });
        }
        self.read_exact_at_impl(buf, offset)
            .map_err(|source| SeqDbError::Read {
                path: self.path.clone(),
                source,
            })
    }

    /// The bytes of `start..end` as a new buffer.
    pub fn read_range(&self, start: u64, end: u64) -> Result<Vec<u8>, SeqDbError> {
        self.read_range_in(start, end, RangeBuffer::New)
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_range` into a buffer taken from or sized
    /// by `utils::x_mapped::X_BYTES.x_filled` (a mapped buffer from 1 MiB on, given back with
    /// `x_mapped::X_BYTES.give`; a buffer of the pool is zeroed only past the bytes it holds).
    /// No NCBI counterpart: only where the bytes are kept.
    pub fn x_read_range_mapped(&self, start: u64, end: u64) -> Result<Vec<u8>, SeqDbError> {
        self.read_range_in(start, end, RangeBuffer::Mapped)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): `read_range` into the buffer that a
    /// worker keeps from one read to the next, `kept` (taken, its capacity kept when it holds the
    /// range; a range of 1 MiB or more that it does not hold goes to a buffer of
    /// `utils::x_mapped::X_BYTES` instead and `kept` stays): `x_mapped::X_BYTES.x_filled_kept`.
    /// The buffer goes back with `X_BYTES.x_keep`. No NCBI counterpart: only where the bytes are
    /// kept.
    pub fn x_read_range_into(
        &self,
        start: u64,
        end: u64,
        kept: &mut Vec<u8>,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.read_range_in(start, end, RangeBuffer::Given(kept))
    }

    fn read_range_in(
        &self,
        start: u64,
        end: u64,
        buffer: RangeBuffer<'_>,
    ) -> Result<Vec<u8>, SeqDbError> {
        if end < start {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!("byte range {start}..{end} ends before it starts"),
            });
        }
        // Checked before the buffer is allocated: a corrupt offset must not
        // ask for gigabytes.
        if end > self.len {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!(
                    "byte range {start}..{end} passes the end of the file ({} bytes)",
                    self.len
                ),
            });
        }
        let len = usize::try_from(end - start).map_err(|_| SeqDbError::Format {
            path: self.path.clone(),
            message: format!("byte range {start}..{end} does not fit in memory"),
        })?;
        // A buffer of the pool or of the worker keeps the bytes its last read left in it, and only
        // the bytes past them are zeroed before the read writes them all.
        let mut buf = match buffer {
            RangeBuffer::New => vec![0u8; len],
            RangeBuffer::Mapped => crate::utils::x_mapped::X_BYTES.x_filled(len, 0u8),
            RangeBuffer::Given(kept) => {
                crate::utils::x_mapped::X_BYTES.x_filled_kept(kept, len, 0u8)
            }
        };
        self.read_exact_at(&mut buf, start)?;
        Ok(buf)
    }

    /// The checks of `read_range(start, end)` before its buffer
    /// is allocated, with the same errors in the same order, without the read.
    pub fn x_check_range(&self, start: u64, end: u64) -> Result<(), SeqDbError> {
        if end < start {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!("byte range {start}..{end} ends before it starts"),
            });
        }
        if end > self.len {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!(
                    "byte range {start}..{end} passes the end of the file ({} bytes)",
                    self.len
                ),
            });
        }
        usize::try_from(end - start).map_err(|_| SeqDbError::Format {
            path: self.path.clone(),
            message: format!("byte range {start}..{end} does not fit in memory"),
        })?;
        Ok(())
    }

    /// The whole file.
    pub fn read_all(&self) -> Result<Vec<u8>, SeqDbError> {
        self.read_range(0, self.len)
    }

    #[cfg(unix)]
    fn read_exact_at_impl(&self, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
        use std::os::unix::fs::FileExt;
        self.file.read_exact_at(buf, offset)
    }

    #[cfg(windows)]
    fn read_exact_at_impl(&self, mut buf: &mut [u8], mut offset: u64) -> std::io::Result<()> {
        use std::os::windows::fs::FileExt;
        while !buf.is_empty() {
            match self.file.seek_read(buf, offset) {
                Ok(0) => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        "failed to fill whole buffer",
                    ));
                }
                Ok(n) => {
                    buf = &mut buf[n..];
                    offset += n as u64;
                }
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
        Ok(())
    }

    // WASI and other targets: `std` has no stable positioned read there, so
    // the seek and the read are paired under a lock.
    #[cfg(not(any(unix, windows)))]
    fn read_exact_at_impl(&self, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
        use std::io::{Read, Seek, SeekFrom};
        let mut file = self
            .file
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        file.seek(SeekFrom::Start(offset))?;
        file.read_exact(buf)
    }
}

/// The tests' files whose bytes the system does not give: `DbFile::open` of a path listed here
/// reads, in place of the file, through a handle that the test names (a directory's, whose read
/// fails with EISDIR; a file shorter than the listed one, whose read past its end fails), with
/// the listed file's length. A test reaches the reads of `DbFile` through the readers that open
/// the file (`SeqDb::has_parsed_seqids`, `SeqDb::deflines`).
#[cfg(all(test, target_os = "linux"))]
#[path = "file_refused_reads_tests.rs"]
pub(crate) mod refused_reads;

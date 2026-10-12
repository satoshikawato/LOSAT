//! Positioned reads of one database file.
//!
//! NCBI maps every database file into memory and hands out pointers into the
//! mapping. LOSAT copies the requested byte range with a positioned read
//! instead (`std` only, no memory mapping), so that the same reader code
//! serves native builds and WebAssembly. Only this file differs by platform.

use std::fs::File;
use std::path::{Path, PathBuf};

use super::SeqDbError;

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
    // ```c++
    //     if (! m_File.Open(CSeqDB_Path(m_FileName))) {
    //         //m_Atlas.Unlock(locked);
    //
    //         string msg = string("Error: File (") + m_FileName + ") not found.";
    //
    //         NCBI_THROW(CSeqDBException, eFileErr, msg);
    //     }
    // ```
    /// Opens `path` for reading. A file that cannot be opened is NCBI's
    /// "Error: File (...) not found.".
    pub fn open(path: &Path) -> Result<Self, SeqDbError> {
        let file = File::open(path).map_err(|_| {
            SeqDbError::NotFound(format!("Error: File ({}) not found.", path.display()))
        })?;
        let len = file
            .metadata()
            .map_err(|source| SeqDbError::Io {
                path: path.to_path_buf(),
                source,
            })?
            .len();
        Ok(Self {
            path: path.to_path_buf(),
            len,
            #[cfg(any(unix, windows))]
            file,
            #[cfg(not(any(unix, windows)))]
            file: std::sync::Mutex::new(file),
        })
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
    // ```c++
    //         *buffer = m_Seq->GetFileDataPtr(start_offset-1) + 1;
    //         ...
    //         *buffer = m_Seq->GetFileDataPtr(start_offset);
    // ```
    // NCBI returns a pointer at `offset` of the mapped file; LOSAT fills `buf`
    // from the same offset. A range that passes the end of the file, which
    // NCBI would read past its mapping, is an explicit format error.
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
            .map_err(|source| SeqDbError::Io {
                path: self.path.clone(),
                source,
            })
    }

    /// The bytes of `start..end` as a new buffer.
    pub fn read_range(&self, start: u64, end: u64) -> Result<Vec<u8>, SeqDbError> {
        if end < start {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: format!("byte range {start}..{end} ends before it starts"),
            });
        }
        let len = usize::try_from(end - start).map_err(|_| SeqDbError::Format {
            path: self.path.clone(),
            message: format!("byte range {start}..{end} does not fit in memory"),
        })?;
        let mut buf = vec![0u8; len];
        self.read_exact_at(&mut buf, start)?;
        Ok(buf)
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

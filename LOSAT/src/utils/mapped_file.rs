//! Read-only views of whole files, as NCBI's `CMemoryFile` gives them to RPS-BLAST: the
//! profile database's `.loo` and `.rps` files are mapped, never read into the heap, on unix and
//! Windows; other targets (wasm32) read the file into a 4-byte-aligned buffer.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:252-262
//! ```c++
//! CRpsMmappedFile::CRpsMmappedFile(const string& filename)
//! {
//!     try { m_MmappedFile.reset(new CMemoryFile(filename)); }
//!     catch (const CFileException& e) {
//!         NCBI_RETHROW(e, CBlastException, eRpsInit,
//!              "Cannot memory map RPS-BLAST database file: " + filename);
//!     }
//! }
//! ```
//! A file changed while it is mapped can end the process (SIGBUS), as it can end NCBI's.

use std::io;
use std::path::Path;

/// A whole file, mapped read-only (or read into an aligned buffer where no mapping is used).
/// An empty file has no mapping and no bytes.
pub(crate) struct MappedFile {
    inner: Inner,
}

enum Inner {
    Empty,
    #[cfg(unix)]
    Unix {
        ptr: *mut std::ffi::c_void,
        len: usize,
    },
    #[cfg(windows)]
    Windows {
        view: windows_sys::Win32::System::Memory::MEMORY_MAPPED_VIEW_ADDRESS,
        mapping: windows_sys::Win32::Foundation::HANDLE,
        len: usize,
    },
    #[cfg(not(any(unix, windows)))]
    Heap {
        words: Vec<u32>,
        len: usize,
    },
}

// SAFETY: the mapping is read-only and owned by this value until `Drop`; nothing writes through
// it, so shared references to its bytes may cross threads.
unsafe impl Send for MappedFile {}
// SAFETY: as for `Send`: the bytes are never written while the value lives.
unsafe impl Sync for MappedFile {}

impl MappedFile {
    /// Maps the file at `path` read-only (its length at the time of the call).
    pub(crate) fn open(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let len = usize::try_from(file.metadata()?.len())
            .map_err(|_| io::Error::new(io::ErrorKind::OutOfMemory, "file too large to map"))?;
        if len == 0 {
            return Ok(Self {
                inner: Inner::Empty,
            });
        }
        Self::map(file, len)
    }

    #[cfg(unix)]
    fn map(file: std::fs::File, len: usize) -> io::Result<Self> {
        use rustix::mm::{mmap, MapFlags, ProtFlags};
        // SAFETY: a new shared read-only mapping of the whole file; no address is requested, so
        // nothing else is replaced. The mapping stays valid after `file` is closed.
        let ptr = unsafe {
            mmap(
                std::ptr::null_mut(),
                len,
                ProtFlags::READ,
                MapFlags::SHARED,
                &file,
                0,
            )
        }?;
        // The profile files are read at random places (the lookup table's cells of a query's
        // words, the rows around a hit): no read-ahead, no mapping of neighbouring pages at a
        // fault, so that the resident set holds the pages that were read. (Advice only: a failure
        // changes nothing.)
        // SAFETY: `ptr..ptr + len` is the mapping made above, and advice does not change it.
        let _ = unsafe { rustix::mm::madvise(ptr, len, rustix::mm::Advice::Random) };
        Ok(Self {
            inner: Inner::Unix { ptr, len },
        })
    }

    #[cfg(windows)]
    fn map(file: std::fs::File, len: usize) -> io::Result<Self> {
        use std::os::windows::io::AsRawHandle;
        use windows_sys::Win32::Foundation::CloseHandle;
        use windows_sys::Win32::System::Memory::{
            CreateFileMappingW, MapViewOfFile, FILE_MAP_READ, PAGE_READONLY,
        };
        // SAFETY: `file` is an open handle; a read-only mapping object of its whole length
        // (sizes 0 = the file's), unnamed, with default security.
        let mapping = unsafe {
            CreateFileMappingW(
                file.as_raw_handle() as _,
                std::ptr::null(),
                PAGE_READONLY,
                0,
                0,
                std::ptr::null(),
            )
        };
        if mapping.is_null() {
            return Err(io::Error::last_os_error());
        }
        // SAFETY: `mapping` was just created; a read view of the whole object.
        let view = unsafe { MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0) };
        if view.Value.is_null() {
            let error = io::Error::last_os_error();
            // SAFETY: `mapping` is a handle this function owns and no view of it exists.
            unsafe { CloseHandle(mapping) };
            return Err(error);
        }
        Ok(Self {
            inner: Inner::Windows { view, mapping, len },
        })
    }

    #[cfg(not(any(unix, windows)))]
    fn map(mut file: std::fs::File, len: usize) -> io::Result<Self> {
        use std::io::Read;
        let mut words = vec![0u32; len.div_ceil(4)];
        // SAFETY: the buffer holds `words.len() * 4 >= len` initialized bytes; viewed as bytes
        // (alignment 1) for the read.
        let bytes = unsafe { std::slice::from_raw_parts_mut(words.as_mut_ptr().cast::<u8>(), len) };
        file.read_exact(bytes)?;
        Ok(Self {
            inner: Inner::Heap { words, len },
        })
    }

    /// The file's bytes.
    pub(crate) fn bytes(&self) -> &[u8] {
        match &self.inner {
            Inner::Empty => &[],
            // SAFETY: `ptr` maps `len` readable bytes until `Drop`.
            #[cfg(unix)]
            Inner::Unix { ptr, len } => unsafe {
                std::slice::from_raw_parts(ptr.cast::<u8>(), *len)
            },
            // SAFETY: `view` maps the whole file, `len` readable bytes, until `Drop`.
            #[cfg(windows)]
            Inner::Windows { view, len, .. } => unsafe {
                std::slice::from_raw_parts(view.Value.cast::<u8>(), *len)
            },
            // SAFETY: the buffer holds at least `len` initialized bytes.
            #[cfg(not(any(unix, windows)))]
            Inner::Heap { words, len } => unsafe {
                std::slice::from_raw_parts(words.as_ptr().cast::<u8>(), *len)
            },
        }
    }

    /// The file's length in bytes.
    pub(crate) fn len(&self) -> usize {
        self.bytes().len()
    }

    /// `count` native-endian `i32`s from byte `offset`, or `None` when they do not lie in the
    /// file or `offset` is not a multiple of 4 (the start of the bytes is page or word aligned).
    pub(crate) fn i32s(&self, offset: usize, count: usize) -> Option<&[i32]> {
        let end = count.checked_mul(4)?.checked_add(offset)?;
        let bytes = self.bytes().get(offset..end)?;
        if count == 0 {
            return Some(&[]);
        }
        let ptr = bytes.as_ptr();
        if !(ptr as usize).is_multiple_of(std::mem::align_of::<i32>()) {
            return None;
        }
        // SAFETY: `ptr` is aligned for `i32` (checked), and `count * 4` bytes from it lie in the
        // mapping, which lives as long as `self`; every bit pattern is an `i32`.
        Some(unsafe { std::slice::from_raw_parts(ptr.cast::<i32>(), count) })
    }

    /// `count` rows of `N` native-endian `i32`s from byte `offset` (`i32s` as rows).
    pub(crate) fn rows<const N: usize>(&self, offset: usize, count: usize) -> Option<&[[i32; N]]> {
        let flat = self.i32s(offset, count.checked_mul(N)?)?;
        // SAFETY: `[i32; N]` has the size `N * 4` and the alignment of `i32`; `flat` holds
        // `count * N` of them, aligned (checked by `i32s`).
        Some(unsafe { std::slice::from_raw_parts(flat.as_ptr().cast::<[i32; N]>(), count) })
    }
}

impl Drop for MappedFile {
    fn drop(&mut self) {
        match &self.inner {
            Inner::Empty => {}
            #[cfg(unix)]
            Inner::Unix { ptr, len } => {
                // SAFETY: `ptr`/`len` is the mapping made by `map`, unmapped once here; no
                // borrow of its bytes outlives `self`.
                let _ = unsafe { rustix::mm::munmap(*ptr, *len) };
            }
            #[cfg(windows)]
            Inner::Windows { view, mapping, .. } => {
                // SAFETY: the view and the handle are those made by `map`, released once here.
                unsafe {
                    windows_sys::Win32::System::Memory::UnmapViewOfFile(*view);
                    windows_sys::Win32::Foundation::CloseHandle(*mapping);
                }
            }
            #[cfg(not(any(unix, windows)))]
            Inner::Heap { .. } => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_mapped_file_has_the_files_bytes_and_aligned_views() {
        let dir = std::env::temp_dir().join(format!(
            "losat-mapped-file-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("words.bin");
        let words: Vec<i32> = (0..70).map(|i| i * 3 - 7).collect();
        let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_ne_bytes()).collect();
        std::fs::write(&path, &bytes).unwrap();
        let mapped = MappedFile::open(&path).unwrap();
        assert_eq!(mapped.bytes(), &bytes[..]);
        assert_eq!(mapped.len(), 280);
        assert_eq!(mapped.i32s(8, 3).unwrap(), &words[2..5]);
        assert!(mapped.i32s(2, 1).is_none());
        assert!(mapped.i32s(276, 2).is_none());
        let rows = mapped.rows::<28>(0, 2).unwrap();
        assert_eq!(rows[1][0], words[28]);
        assert!(mapped.rows::<28>(0, 3).is_none());
        let empty = dir.join("empty.bin");
        std::fs::write(&empty, b"").unwrap();
        assert_eq!(MappedFile::open(&empty).unwrap().len(), 0);
        assert!(MappedFile::open(&dir.join("missing")).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

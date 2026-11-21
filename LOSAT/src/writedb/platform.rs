// Derived from LMDB 0.9.24 (mdb.c, lmdb.h) as bundled in NCBI BLAST+ 2.17.0,
// c++/src/util/lmdb/, which carries this notice:
//
// Copyright 2011-2019 Howard Chu, Symas Corp.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted only as authorized by the OpenLDAP
// Public License.
//
// A copy of this license is available in the file LICENSE in the
// top-level directory of the distribution or, alternatively, at
// <http://www.OpenLDAP.org/license.html>.
//
// This code is derived from btree.c written by Martin Hedenfalk.
//
// Copyright (c) 2009, 2010 Martin Hedenfalk <martin@bzero.se>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//
// (In LOSAT, the OpenLDAP Public License 2.8 is in THIRD_PARTY_NOTICES, section 2 (LMDB).)

//! The properties of the platform that reach the bytes makeblastdb writes.
//!
//! makeblastdb writes the same files on every platform except where NCBI's source asks the
//! platform: the page size of the LMDB files of version 5 (`GET_PAGESIZE`, mdb.c), their
//! default map sizes (`DEFAULT_LMDB_MAP_SIZE`, `NCBI_OS_MSWIN`), the length LMDB gives a data
//! file on Windows (`mdb_env_map`, `_WIN32`), `sizeof(long)` in the header size of an index
//! file (`x_Overhead`, versions 4 and 5) and the `size_t` that holds a map size. A build takes
//! these from the host once ([`Platform::host`]); the byte-for-byte tests, whose expected files
//! are makeblastdb's on Linux x86-64, pin [`Platform::LINUX_X86_64`] so that they hold on every
//! platform the tests run on (makedb-v5 line, PR1 review B1).
//!
//! makeblastdb on Windows and on macOS has not been run against LOSAT: what LOSAT writes there
//! follows NCBI's source as read here (`PD-LOSAT-DATABASE-SEARCH`).

use std::sync::OnceLock;

use super::WriteDbError;

/// The page size of the LMDB files: the operating system's, or a fixed one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PageSize {
    /// LMDB's `GET_PAGESIZE` on the host, looked up when a version 5 database creates its
    /// LMDB files (`host_page_size`).
    Host,
    /// The page size of another platform.
    Fixed(usize),
}

/// What NCBI's makeblastdb takes from the platform it runs on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Platform {
    /// LMDB's `me_os_psize` (`GET_PAGESIZE`).
    pub page_size: PageSize,
    /// `DEFAULT_LMDB_MAP_SIZE`: the map size of the `.pdb`/`.ndb` file when
    /// `BLASTDB_LMDB_MAP_SIZE` is not set.
    pub default_lmdb_map_size: u64,
    /// `DEFAULT_TAXID_MAP_SIZE`: the same for the `.ptf`/`.ntf` file.
    pub default_taxid_map_size: u64,
    /// LMDB sets the length of a data file to the map size whenever it maps it (Windows).
    pub file_length_is_map_size: bool,
    /// `sizeof(long)` in `CWriteDB_IndexFile::x_Overhead`.
    pub sizeof_long: u64,
    /// The largest map size the platform's `size_t` holds (`mdb_env_set_mapsize(size_t)`).
    pub max_map_size: u64,
}

/// The operating systems whose makeblastdb differs in NCBI's source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NcbiOs {
    /// `NCBI_OS_MSWIN` / `_WIN32`.
    Windows,
    /// Every other platform NCBI builds on (`NCBI_OS_UNIX`).
    Unix,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_lmdb.hpp:51-57
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3979-4003 (mdb_env_map)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:347-351
/// The platform values NCBI's source gives for the operating system `os` with a C `long` of
/// `sizeof_long` bytes and a `size_t` of at most `max_map_size`.
pub(crate) const fn ncbi_platform(
    os: NcbiOs,
    page_size: PageSize,
    sizeof_long: u64,
    max_map_size: u64,
) -> Platform {
    let windows = matches!(os, NcbiOs::Windows);
    Platform {
        page_size,
        default_lmdb_map_size: if windows { 500_000 } else { 700_000_000 },
        default_taxid_map_size: if windows { 500_000 } else { 300_000_000 },
        file_length_is_map_size: windows,
        sizeof_long,
        max_map_size,
    }
}

impl Platform {
    /// makeblastdb on Linux x86-64 (LP64, 4 KiB pages): the platform of every expected file
    /// of the tests, of the gates of the makedb-v5 line, and of the Web app's builds.
    pub const LINUX_X86_64: Platform =
        ncbi_platform(NcbiOs::Unix, PageSize::Fixed(4096), 8, u64::MAX);

    /// The platform LOSAT runs on, as NCBI's source sees it when built there: Windows or
    /// not, the C `long` of the target, its `size_t`, and the operating system's page size.
    /// wasm32 (no NCBI build, no operating-system page) writes what makeblastdb writes on
    /// Linux x86-64 (DECISIONS 8 of the makedb-v5 line).
    pub fn host() -> Platform {
        if cfg!(target_family = "wasm") {
            return Platform::LINUX_X86_64;
        }
        let os = if cfg!(windows) {
            NcbiOs::Windows
        } else {
            NcbiOs::Unix
        };
        ncbi_platform(
            os,
            PageSize::Host,
            std::mem::size_of::<std::ffi::c_long>() as u64,
            usize::MAX as u64,
        )
    }

    /// The values the LMDB writer needs, with the page size looked up.
    pub(crate) fn lmdb(&self) -> Result<LmdbPlatform, WriteDbError> {
        let os_page_size = match self.page_size {
            PageSize::Fixed(size) => size,
            PageSize::Host => host_page_size().map_err(|error| {
                WriteDbError::Other(format!(
                    "the page size of the operating system could not be found: {error}"
                ))
            })?,
        };
        Ok(LmdbPlatform {
            os_page_size,
            file_length_is_map_size: self.file_length_is_map_size,
            max_map_size: self.max_map_size,
        })
    }
}

/// The platform values of one LMDB environment ([`Platform::lmdb`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct LmdbPlatform {
    /// `me_os_psize`.
    pub os_page_size: usize,
    pub file_length_is_map_size: bool,
    pub max_map_size: u64,
}

impl LmdbPlatform {
    /// [`Platform::LINUX_X86_64`]'s.
    #[cfg(test)]
    pub(crate) const LINUX_X86_64: LmdbPlatform = LmdbPlatform {
        os_page_size: 4096,
        file_length_is_map_size: false,
        max_map_size: u64::MAX,
    };
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:321,395,3969
/// LMDB's `GET_PAGESIZE` on the host: `GetSystemInfo` on Windows, `sysconf(_SC_PAGE_SIZE)`
/// elsewhere. On Linux the C library's `sysconf` returns the kernel's `AT_PAGESZ`, which is
/// the kernel's page size; LOSAT asks the kernel for that size directly
/// ([`kernel_page_size`]), without `/proc` and without the C library (rustix reads `AT_PAGESZ`
/// from `/proc/self/auxv` or `PR_GET_AUXV` and panics where neither exists, as in a chroot
/// without `/proc` on Linux before 6.4; the C library's own copy is reachable only through
/// rustix's `use-libc-auxv`, which the pure-Rust runtime boundary leaves out). Found once per
/// process.
pub(crate) fn host_page_size() -> std::io::Result<usize> {
    static PAGE_SIZE: OnceLock<usize> = OnceLock::new();
    if let Some(&size) = PAGE_SIZE.get() {
        return Ok(size);
    }
    let size = os_page_size()?;
    Ok(*PAGE_SIZE.get_or_init(|| size))
}

#[cfg(any(target_os = "linux", target_os = "android"))]
fn os_page_size() -> std::io::Result<usize> {
    kernel_page_size()
}

#[cfg(all(unix, not(any(target_os = "linux", target_os = "android"))))]
fn os_page_size() -> std::io::Result<usize> {
    // rustix's libc backend: `sysconf(_SC_PAGESIZE)`, which cannot fail.
    Ok(rustix::param::page_size())
}

#[cfg(windows)]
fn os_page_size() -> std::io::Result<usize> {
    use windows_sys::Win32::System::SystemInformation::{GetSystemInfo, SYSTEM_INFO};
    // SAFETY: GetSystemInfo fills the SYSTEM_INFO it is given and cannot fail.
    let mut si: SYSTEM_INFO = unsafe { std::mem::zeroed() };
    unsafe { GetSystemInfo(&mut si) };
    Ok(si.dwPageSize as usize)
}

#[cfg(not(any(unix, windows)))]
fn os_page_size() -> std::io::Result<usize> {
    // Not reached: wasm32 has no operating-system page and takes `Platform::LINUX_X86_64`
    // (`Platform::host`); no other target without `unix` or `windows` is built.
    Err(std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "this platform has no operating-system page size",
    ))
}

/// The kernel's page size on Linux: the alignment `mprotect` asks of an address
/// (POSIX: `EINVAL` when the address is not a multiple of the page size), tried on a fresh
/// mapping without access rights, whose start is page-aligned, at offsets 1, 2, 4, ... until
/// one is accepted. Fails only when the kernel refuses `mmap` or `mprotect` altogether (no
/// memory for a mapping, or a seccomp filter); the error is returned, no size is guessed.
#[cfg(any(target_os = "linux", target_os = "android"))]
pub(crate) fn kernel_page_size() -> std::io::Result<usize> {
    use rustix::io::Errno;
    use rustix::mm::{mmap_anonymous, mprotect, munmap, MapFlags, MprotectFlags, ProtFlags};
    // Larger than any page size Linux has (at most 256 KiB, hexagon; 64 KiB on arm64,
    // powerpc64 and loongarch64).
    const SPAN: usize = 1 << 20;
    // SAFETY: a new private anonymous mapping that nothing else knows of; it is unmapped below.
    let base = unsafe {
        mmap_anonymous(
            std::ptr::null_mut(),
            SPAN,
            ProtFlags::empty(),
            MapFlags::PRIVATE,
        )
    }?;
    let mut result = Err(std::io::Error::other(
        "no offset below 1 MiB of a new mapping is page-aligned",
    ));
    let mut offset = 1usize;
    while offset < SPAN {
        // SAFETY: the byte at `offset` lies inside the mapping made above, and its protection
        // stays PROT_NONE.
        let probe = unsafe {
            mprotect(
                base.cast::<u8>().add(offset).cast(),
                1,
                MprotectFlags::empty(),
            )
        };
        match probe {
            Ok(()) => {
                result = Ok(offset);
                break;
            }
            Err(Errno::INVAL) => offset <<= 1,
            Err(error) => {
                result = Err(std::io::Error::from(error));
                break;
            }
        }
    }
    // SAFETY: the mapping made above, of SPAN bytes, used by nothing else.
    let unmapped = unsafe { munmap(base, SPAN) };
    let size = result?;
    unmapped?;
    Ok(size)
}

#[cfg(test)]
#[path = "platform_tests.rs"]
mod tests;

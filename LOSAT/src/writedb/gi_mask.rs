//! The GI-indexed mask files of `makeblastdb -gi_mask` (`CWriteDB_GiMask`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp` and `writedb_gimask.hpp`).
//!
//! Each masking algorithm of a build with `-gi_mask` writes six files under its name (the
//! `-gi_mask_name` entry, else the `-mask_data` file as typed, both relative to the current
//! directory, not to the database):
//!
//! * `<name>.gmd` and `<name>.gnd`: the masked ranges of the sequences, one record per sequence
//!   (`int32 n`, then `n` pairs of start and end, the end exclusive), big-endian and
//!   little-endian. Records go to data volumes `<name>.00.gmd`, `<name>.01.gmd`, ... when the
//!   next one would pass `-max_file_sz`; a single volume is renamed to `<name>.gmd`.
//! * `<name>.gmo` and `<name>.gno`: every GI of the masked sequences, sorted, with the data volume
//!   and offset of its record, in pages of 512 (the GIs of a page, then their offsets).
//! * `<name>.gmi` and `<name>.gni`: a header (the algorithm's description, `<program>:<options>`,
//!   and the date) and every 512th GI of the offset file with the last one.
//!
//! The headers are big-endian in both index files. Every file is created by its first write: an
//! algorithm whose masks reached no sequence with a GI writes none, and its name becomes empty
//! (the alias file's `MASKLIST` leaves it out, `super::writer::WriteDb::close`).
//!
//! Two algorithms whose names lead to the same files (one `-gi_mask_name` twice, or two names of
//! one path) are written as NCBI writes them while their data files are written one after the
//! other: an algorithm that reaches no sequence with a GI writes nothing, the offset and index
//! files are each created, written and closed in turn at the close, and a data volume that one
//! algorithm closed before the other created it is created anew. LOSAT refuses them when the
//! data volume of one would be created while the other's of the same path is open
//! ([`GiMasks::add_gi_mask`]): NCBI writes both through two streams into that file, and its bytes
//! then depend on how the C++ runtime buffers two streams on one file, which LOSAT does not
//! reproduce (survey D6).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::column::{write_pad_bytes_string, write_string_var};
use super::files::WriteDbFile;
use super::WriteDbError;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.hpp:145-152
// The bytes of a GI and of the volume and offset after it in the offset file, and the number of
// GIs in a page of the offset file.
const GI_SIZE: u32 = 4;
const OFFSET_SIZE: u32 = 8;
const PAGE_SIZE: usize = 512;

/// An element of `m_GiOffset` (`pair<TGi, pair<int, int>>`): the GI and the data volume and
/// offset of its record.
type GiOffset = (i64, i32, i32);

/// `WriteInt4` or `WriteInt4_LE` of `CBlastDbBlob`.
fn put_int4(out: &mut Vec<u8>, value: u32, le: bool) {
    out.extend_from_slice(&if le {
        value.to_le_bytes()
    } else {
        value.to_be_bytes()
    });
}

/// `CWriteDB_GiMaskData`: one data volume of an algorithm, in one byte order.
#[derive(Debug)]
struct GiMaskData {
    file: WriteDbFile,
    /// `m_DataLength`.
    data_length: u64,
    /// `m_UseLE`.
    le: bool,
    /// `m_Index`.
    index: u32,
}

impl GiMaskData {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:279-288
    /// `CWriteDB_GiMaskData`: data volume `index` of `name` with `extension`, created by its
    /// first write, empty, in the byte order `le` (little-endian for the `.gnd` files).
    fn new(name: &[u8], extension: &'static str, index: u32, max_file_size: u64, le: bool) -> Self {
        Self {
            file: WriteDbFile::lazy(name, extension, index, max_file_size),
            data_length: 0,
            le,
            index,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.hpp:96-99
    // `m_MaxFileSize` is the file's (`-max_file_sz`, 999,999,999 for 0). NCBI computes the size of
    // the record (twice the number of ranges plus one, times four) in `int`, which overflows
    // (undefined behaviour) past 268,435,455 ranges; LOSAT computes it without overflow.
    /// `CanFit`: whether a record of `num_masks` ranges keeps the volume below its size limit.
    fn can_fit(&self, num_masks: usize) -> bool {
        self.data_length + (num_masks as u64 * 2 + 1) * 4 < self.file.max_file_size()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.hpp:88-91
    /// `GetOffset`: the volume and the offset of the next record, as `int`s.
    fn offset(&self) -> (i32, i32) {
        (self.index as i32, self.data_length as i32)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:290-315
    /// `WriteMask`: one record, the number of ranges and the start and end of each, in the
    /// volume's byte order; an empty list writes nothing.
    fn write_mask(&mut self, mask: &[(u32, u32)]) {
        if mask.is_empty() {
            return;
        }
        if !self.file.created() {
            self.file.create_now();
        }
        let mut data = Vec::with_capacity(4 + 8 * mask.len());
        put_int4(&mut data, mask.len() as u32, self.le);
        for &(first, second) in mask {
            put_int4(&mut data, first, self.le);
            put_int4(&mut data, second, self.le);
        }
        self.file.write_unaddressed(&data);
        self.data_length += (1 + 2 * mask.len() as u64) * 4;
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:133-184
// NCBI converts a GI to `int` (`GI_TO`), which keeps its low 32 bits: a GI from 2^31 to 2^32 - 1 is written as the same
// four bytes (larger GIs never reach the masks, `deflines::large_gi_to_skip`).
/// `CWriteDB_GiMaskOffset::AddGIs`: the offset file, page by page.
fn write_offset_file(file: &mut WriteDbFile, gi_offset: &[GiOffset], le: bool) {
    if !file.created() {
        file.create_now();
    }
    let mut gis = Vec::with_capacity(PAGE_SIZE * GI_SIZE as usize);
    let mut offsets = Vec::with_capacity(PAGE_SIZE * OFFSET_SIZE as usize);
    for page in gi_offset.chunks(PAGE_SIZE) {
        gis.clear();
        offsets.clear();
        for &(gi, volume, offset) in page {
            put_int4(&mut gis, gi as i32 as u32, le);
            put_int4(&mut offsets, volume as u32, le);
            put_int4(&mut offsets, offset as u32, le);
        }
        file.write_unaddressed(&gis);
        file.write_unaddressed(&offsets);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:200-251
// The GIs sampled are the first of every page of the offset file and the last GI.
/// `CWriteDB_GiMaskIndex::AddGIs`: the index file, its header first.
fn write_index_file(
    file: &mut WriteDbFile,
    gi_offset: &[GiOffset],
    num_vols: u32,
    desc: &[u8],
    date: &[u8],
    le: bool,
) {
    if !file.created() {
        file.create_now();
    }
    let num_gis = gi_offset.len();
    let mut gis = Vec::new();
    let mut offsets = Vec::new();
    let mut num_index = 0u32;
    for (i, &(gi, volume, offset)) in gi_offset.iter().enumerate() {
        if i % PAGE_SIZE != 0 && i + 1 < num_gis {
            continue;
        }
        put_int4(&mut gis, gi as i32 as u32, le);
        put_int4(&mut offsets, volume as u32, le);
        put_int4(&mut offsets, offset as u32, le);
        num_index += 1;
    }
    file.write_unaddressed(&index_header(
        num_vols,
        num_index,
        num_gis as u32,
        desc,
        date,
    ));
    file.write_unaddressed(&gis);
    file.write_unaddressed(&offsets);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:253-275
// `kStringFmt` is `CBlastDbBlob::eSizeVar` (writedb_gimask.hpp:191-192). The header is
// big-endian in both index files.
/// `x_BuildHeaderFields`: the header of an index file.
fn index_header(num_vols: u32, num_index: u32, num_gis: u32, desc: &[u8], date: &[u8]) -> Vec<u8> {
    let mut header = Vec::with_capacity(48 + desc.len() + date.len());
    for value in [
        1,
        num_vols,
        GI_SIZE,
        OFFSET_SIZE,
        PAGE_SIZE as u32,
        num_index,
        num_gis,
        0,
    ] {
        put_int4(&mut header, value, false);
    }
    write_string_var(&mut header, desc);
    write_string_var(&mut header, date);
    write_pad_bytes_string(&mut header, 8);
    let size = header.len() as u32;
    header[28..32].copy_from_slice(&size.to_be_bytes());
    header
}

/// `CWriteDB_GiMask`: the files of one masking algorithm.
#[derive(Debug)]
pub(crate) struct GiMask {
    /// `m_MaskName`: empty after a close without GIs.
    name: Vec<u8>,
    max_file_size: u64,
    /// `m_DFile`, `m_DFile_LE`.
    data: [GiMaskData; 2],
    /// `m_OFile`, `m_OFile_LE`.
    offset_files: [WriteDbFile; 2],
    /// `m_IFile`, `m_IFile_LE`.
    index_files: [WriteDbFile; 2],
    /// `m_Desc` of the index files.
    desc: Vec<u8>,
    /// `m_Date` of the `.gmi` and of the `.gni` file.
    dates: [Vec<u8>; 2],
    /// `m_GiOffset`.
    gi_offset: Vec<GiOffset>,
    /// The data volumes closed for the next one, by name: those that a record created.
    closed_volumes: Vec<Vec<u8>>,
}

impl GiMask {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:46-57
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:188-198
    /// The files of the algorithm `name` with the description `desc`; `dates` are the dates of
    /// the `.gmi` and the `.gni` file (the current time as `CTime` text, which each of NCBI's
    /// constructors reads; the caller reads them in that order).
    pub(crate) fn new(name: &[u8], desc: &[u8], max_file_size: u64, dates: [Vec<u8>; 2]) -> Self {
        Self {
            name: name.to_vec(),
            max_file_size,
            data: [
                GiMaskData::new(name, "gmd", 0, max_file_size, false),
                GiMaskData::new(name, "gnd", 0, max_file_size, true),
            ],
            offset_files: [
                WriteDbFile::lazy_unindexed(name, "gmo", max_file_size),
                WriteDbFile::lazy_unindexed(name, "gno", max_file_size),
            ],
            index_files: [
                WriteDbFile::lazy_unindexed(name, "gmi", max_file_size),
                WriteDbFile::lazy_unindexed(name, "gni", max_file_size),
            ],
            desc: desc.to_vec(),
            dates,
            gi_offset: Vec::new(),
            closed_volumes: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.hpp:251-253
    /// `GetName`: the name, or empty after a close without GIs.
    pub(crate) fn name(&self) -> &[u8] {
        &self.name
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:70-88
    // A record that does not fit an empty volume opens the next one all the same: the volume
    // `.00` of a first record larger than `-max_file_sz` is never created (survey D8), and a
    // record larger than the limit is written whole into a volume of its own.
    /// The data volume that a record of `num_masks` ranges goes to (`AddGiMask` below), and
    /// whether the record creates its files.
    fn next_volume(&self, num_masks: usize) -> (u32, bool) {
        match self.data[0].can_fit(num_masks) {
            true => (self.data[0].index, !self.data[0].file.created()),
            false => (self.data[0].index + 1, true),
        }
    }

    /// The data volume whose files are open: created by a record and not yet closed.
    fn open_volume(&self) -> Option<u32> {
        self.data[0].file.created().then_some(self.data[0].index)
    }

    /// `AddGiMask`: the ranges of one sequence, under each of its GIs.
    pub(crate) fn add_gi_mask(&mut self, gis: &[i64], mask: &[(u32, u32)]) {
        if !self.data[0].can_fit(mask.len()) {
            let index = self.data[0].index + 1;
            for (data, (extension, le)) in self.data.iter_mut().zip([("gmd", false), ("gnd", true)])
            {
                data.file.close();
                if data.file.created() {
                    self.closed_volumes.push(data.file.filename().to_vec());
                }
                *data = GiMaskData::new(&self.name, extension, index, self.max_file_size, le);
            }
        }
        let (volume, offset) = self.data[0].offset();
        for data in &mut self.data {
            data.write_mask(mask);
        }
        self.gi_offset
            .extend(gis.iter().map(|&gi| (gi, volume, offset)));
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:90-120
    // `sort` uses the `operator<` of `std::pair`, which compares the GI, then the volume, then
    // the offset: elements that compare equal are equal, so the result does not depend on the
    // sort algorithm. Two sequences with one GI (a version 4 database allows that) keep both
    // entries, the record written first before the other.
    /// `Close`: the data files are closed and the offset and index files written.
    pub(crate) fn close(&mut self) {
        if self.gi_offset.is_empty() {
            self.name.clear();
            return;
        }
        for data in &mut self.data {
            data.file.close();
        }
        let num_vols = self.data[0].index + 1;
        if num_vols == 1 {
            for data in &mut self.data {
                data.file.rename_single();
            }
        }
        self.gi_offset.sort_unstable();
        for (file, (date, le)) in self
            .index_files
            .iter_mut()
            .zip([(&self.dates[0], false), (&self.dates[1], true)])
        {
            write_index_file(file, &self.gi_offset, num_vols, &self.desc, date, le);
            file.close();
        }
        for (file, le) in self.offset_files.iter_mut().zip([false, true]) {
            write_offset_file(file, &self.gi_offset, le);
            file.close();
        }
    }

    /// The files of the algorithm that a record or the close created, appended to `files`
    /// (`GiMasks::files`).
    fn files(&self, files: &mut Vec<Vec<u8>>) {
        files.extend(self.closed_volumes.iter().cloned());
        let open = self.data.iter().map(|data| &data.file);
        for file in open.chain(&self.offset_files).chain(&self.index_files) {
            if file.created() {
                files.push(file.filename().to_vec());
            }
        }
    }
}

/// `m_GiMasks` and `m_MaskAlgoMap` of `CWriteDB_Impl`: the GI masks of a build with `-gi_mask`,
/// in the order of registration.
#[derive(Debug, Default)]
pub(crate) struct GiMasks {
    masks: Vec<GiMask>,
    /// The place of each mask's files (`files_place`), in the order of `masks`.
    places: Vec<Option<(PathBuf, Vec<u8>)>>,
    /// `m_MaskAlgoMap`: the position in `masks` of each algorithm id.
    algo_map: BTreeMap<i32, usize>,
}

impl GiMasks {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1453-1456
    /// The registration of the algorithm `algorithm_id` (`RegisterMaskAlgorithm` with
    /// `m_UseGiMask`): its files under `name`.
    pub(crate) fn add(
        &mut self,
        algorithm_id: i32,
        name: &[u8],
        desc: &[u8],
        max_file_size: u64,
        dates: [Vec<u8>; 2],
    ) {
        self.algo_map.insert(algorithm_id, self.masks.len());
        self.places.push(files_place(name));
        self.masks
            .push(GiMask::new(name, desc, max_file_size, dates));
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1391-1397
    // Refused before the record is written when it would create a data volume whose file is open
    // for another algorithm (the same place and volume number): NCBI would write the two
    // algorithms' records through two streams into that file (see the module's documentation).
    /// `AddGiMask` of the algorithm `algorithm_id`, which is registered.
    pub(crate) fn add_gi_mask(
        &mut self,
        algorithm_id: i32,
        gis: &[i64],
        mask: &[(u32, u32)],
    ) -> Result<(), WriteDbError> {
        let Some(&index) = self.algo_map.get(&algorithm_id) else {
            return Ok(());
        };
        if let (Some(place), (volume, true)) = (
            &self.places[index],
            self.masks[index].next_volume(mask.len()),
        ) {
            let open = self.masks.iter().zip(&self.places).enumerate().find(
                |&(other, (mask, other_place))| {
                    other != index
                        && other_place.as_ref() == Some(place)
                        && mask.open_volume() == Some(volume)
                },
            );
            if let Some((_, (other, _))) = open {
                return Err(shared_files_refusal(other.name(), self.masks[index].name()));
            }
        }
        self.masks[index].add_gi_mask(gis, mask);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:338-342
    /// Closes every GI mask, in the order of registration.
    pub(crate) fn close(&mut self) {
        for mask in &mut self.masks {
            mask.close();
        }
    }

    /// `GetName` of every GI mask, in the order of registration.
    pub(crate) fn names(&self) -> impl Iterator<Item = &[u8]> {
        self.masks.iter().map(GiMask::name)
    }

    /// The files that the GI masks created, appended to `files`
    /// (`mask_data::MaskWriter::gi_mask_files`).
    pub(crate) fn files(&self, files: &mut Vec<Vec<u8>>) {
        for mask in &self.masks {
            mask.files(files);
        }
    }
}

/// Lower-cases every valid UTF-8 chunk of `bytes` while preserving invalid byte sequences
/// verbatim, so distinct non-UTF-8 file names (`a\xff`, `a\xfe`) do not collide under case folding.
#[cfg_attr(not(any(windows, target_os = "macos", test)), allow(dead_code))]
fn lowercase_utf8_bytes(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len());
    for chunk in bytes.utf8_chunks() {
        out.extend_from_slice(chunk.valid().to_lowercase().as_bytes());
        out.extend_from_slice(chunk.invalid());
    }
    out
}

/// The directory of the files of the GI mask `name` (resolved through symbolic links) and the
/// last part of the name, which with a dot and an extension names each file; on the platforms
/// whose file systems ignore case by default (Windows, macOS), both in lower case. `None` when
/// the directory does not exist: nothing is written there (survey D9).
fn files_place(name: &[u8]) -> Option<(PathBuf, Vec<u8>)> {
    let is_separator = |&b: &u8| b == b'/' || (cfg!(windows) && (b == b'\\' || b == b':'));
    let (dir, file) = match name.iter().rposition(is_separator) {
        Some(at) => (&name[..=at], &name[at + 1..]),
        None => (b".".as_slice(), name),
    };
    let dir = resolved_dir(&crate::utils::os_bytes::path_of_bytes(dir))?;
    if cfg!(any(windows, target_os = "macos")) {
        let dir = PathBuf::from(dir.to_string_lossy().to_lowercase());
        return Some((dir, lowercase_utf8_bytes(file)));
    }
    Some((dir, file.to_vec()))
}

/// The directory `dir` resolved through symbolic links, as `realpath(3)` resolves it
/// (`utils::same_file::canonical_path`: also where the platform has no `canonicalize`, LOSAT Web's
/// engine, so that two names of one directory are one place there too; line NP WW4b1,
/// 2026-11-27); `None` where it is not a directory.
fn resolved_dir(dir: &Path) -> Option<PathBuf> {
    crate::utils::same_file::canonical_path(dir).filter(|path| path.is_dir())
}

/// The refusal of the algorithm `second`, whose next record would open a GI-mask data file that
/// the algorithm `first` has open.
fn shared_files_refusal(first: &[u8], second: &[u8]) -> WriteDbError {
    let first = String::from_utf8_lossy(first);
    let second = String::from_utf8_lossy(second);
    let same = if first == second {
        format!("\"{first}\"")
    } else {
        format!("\"{first}\" and \"{second}\"")
    };
    WriteDbError::Unsupported(format!(
        "-gi_mask: two masking algorithms would write one GI-mask data file at the same time \
         (the files of {same}); NCBI makeblastdb's bytes then depend on how the C++ runtime \
         buffers two streams on one file, and LOSAT does not reproduce that (a -gi_mask_name \
         that names other files for one of them avoids it)"
    ))
}

#[cfg(test)]
#[path = "gi_mask_tests.rs"]
mod tests;

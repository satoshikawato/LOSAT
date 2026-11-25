//! The RPS-BLAST files of one profile-database volume (`makeprofiledb -dbtype rps`): the
//! auxiliary text (`.aux`: matrix, gap costs, scale, one K per profile), the prebuilt lookup
//! table (`.loo`) and the concatenated profile scores (`.rps`). The sequence part of the volume
//! (the consensus sequences and their deflines) is an ordinary protein volume read by
//! [`super::SeqDb`]. The `.loo` and `.rps` files are mapped (`utils::mapped_file`), never read
//! into the heap (but the rows of a 26-column `.rps`, which are widened to 28).
//!
//! The files are opened in NCBI's order, with NCBI's errors (`RpsOpenError::Ncbi`, which a
//! search reports as `BLAST engine error: <text>` with exit status 3); every inconsistency on
//! which NCBI would read outside its buffers is LOSAT's refusal (`RpsOpenError::Refused`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:549-575

use std::path::{Path, PathBuf};

use crate::utils::mapped_file::MappedFile;

/// `RPS_MAGIC_NUM`: the files of a 26-letter alphabet.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_rps.h:43-44
pub(crate) const RPS_MAGIC_NUM: i32 = 0x1e16;
/// `RPS_MAGIC_NUM_28`: the files of the 28-letter alphabet (makeprofiledb 2.17 writes these).
pub(crate) const RPS_MAGIC_NUM_28: i32 = 0x1e17;

/// `RPS_HITS_PER_CELL`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:330-347
pub(crate) const RPS_HITS_PER_CELL: i32 = 3;

/// The word size of the prebuilt table (`BLAST_WORDSIZE_PROT`) and its letters' bits
/// (`ilog2(alphabet_size) + 1`, 5 for 26 and 28 letters): `backbone_size = 1 << 15`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:148-154
pub(crate) const RPS_WORDSIZE: usize = 3;
pub(crate) const RPS_CHARSIZE: usize = 5;
pub(crate) const RPS_BACKBONE_SIZE: usize = 1 << (RPS_WORDSIZE * RPS_CHARSIZE);

/// `RPS_BUCKET_SIZE`: the concatenated rows per bucket of a scan's hits (`rps_search::scan`),
/// `2^RPS_BUCKET_BITS`: a position is its bucket's index times 2,048 plus a value below 2,048.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:349-351
pub(crate) const RPS_BUCKET_SIZE: usize = 1 << RPS_BUCKET_BITS;
pub(crate) const RPS_BUCKET_BITS: u32 = 11;

/// The buckets whose hits a pass of the scan makes and reads together (`RpsBuckets::fill_group`,
/// series PR, item PS12m): 64 buckets = 131,072 rows. No NCBI counterpart: a memory bound (NCBI
/// keeps every bucket of a scan, a volume of CDD 162 to 355 of them, 4,000,000 hits at most); the
/// order of the hits that the word finder reads is NCBI's whatever the number. A volume of at most
/// 64 buckets is one group. The scan needs the hits of every cell in groups in ascending order
/// (`list_is_grouped`; makeprofiledb writes every list ascending), so this module owns the number:
/// it regroups the cells of a table that does not (`regroup_cells`) at this granularity (mg5k_100
/// x CDD, 32 threads: groups of 128 buckets hold 472 MB, of 64 361 MB; groups of 32 cost 1 %
/// more instructions than 64).
pub(crate) const RPS_GROUP_BUCKETS: usize = 64;

/// The concatenated position of a table entry's word (`entry` less `BLAST_WORDSIZE_PROT`, plus one;
/// NCBI's `q_off` is `cell->entries[i]` minus `table_correction`): the one producer of it.
#[inline(always)]
pub(crate) fn entry_q_off(entry: i32) -> u32 {
    entry.wrapping_sub(RPS_WORDSIZE as i32 - 1) as u32
}

/// The bucket of a concatenated position (`q_off / RPS_BUCKET_SIZE`).
#[inline(always)]
pub(crate) fn q_off_bucket(q_off: u32) -> usize {
    (q_off >> RPS_BUCKET_BITS) as usize
}

/// The group of `group_buckets` buckets that a table entry's position lies in.
#[inline(always)]
fn entry_group(entry: i32, group_buckets: usize) -> u32 {
    (entry_q_off(entry) >> RPS_BUCKET_BITS) / group_buckets as u32
}

/// The most hits that one scan of a query collects (`max_hits` of `BlastRPSScanSubject`). A cell
/// with more hits than this would never fit, and NCBI's scan would not advance past its word.
pub(crate) const RPS_MAX_SCAN_HITS: i32 = 4_000_000;

/// The size of the `.loo` header (`BlastRPSLookupFileHeader`, ten `Int4`).
const LOO_HEADER_BYTES: usize = 40;

/// Why a volume's RPS files cannot be searched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum RpsOpenError {
    /// NCBI's own error (`CBlastException` with `eRpsInit`): the bytes after `BLAST engine
    /// error: ` (a path in them as the system's bytes).
    Ncbi(Vec<u8>),
    /// LOSAT's refusal of files that NCBI would read outside their bounds or that LOSAT does not
    /// search.
    Refused(String),
}

/// Why a volume's RPS files cannot be searched, as `RpsVolumeFiles::open_for` tells it: an
/// error of the files (`RpsOpenError`), or an `.aux` without a token. NCBI's stream then reads
/// no matrix name, and every query of the batch fails where its score block is made, after the
/// files are open (NCBI's error, which the search writes: "Error: Unknown error code -1 " per
/// query, as an unknown matrix of BLASTX's).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum RpsVolumeFault {
    Files(RpsOpenError),
    NoMatrix,
}

impl From<RpsOpenError> for RpsVolumeFault {
    fn from(error: RpsOpenError) -> Self {
        Self::Files(error)
    }
}

/// The contents of the `.aux` file.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:193-225
/// The loop pushes one more K after the last pair (the extraction that fails at the end of the
/// file); it is never indexed and not kept here.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RpsAux {
    pub(crate) matrix: String,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) scale: f64,
    pub(crate) karlin_k: Vec<f64>,
}

/// Reads the `.aux` text: the tokens separated by any white space, as the C++ stream's `>>`
/// reads them. `None` when a token is not of its kind or a pair is incomplete (NCBI's stream
/// then fails and leaves its values at 0).
pub(crate) fn parse_aux(text: &str) -> Option<RpsAux> {
    let mut tokens = text.split_ascii_whitespace();
    let matrix = tokens.next()?.to_string();
    let gap_open = tokens.next()?.parse::<i32>().ok()?;
    let gap_extend = tokens.next()?.parse::<i32>().ok()?;
    for _ in 0..2 {
        tokens.next()?.parse::<f64>().ok()?;
    }
    for _ in 0..2 {
        tokens.next()?.parse::<i32>().ok()?;
    }
    let scale = tokens.next()?.parse::<f64>().ok()?;
    let mut karlin_k = Vec::new();
    while let Some(length) = tokens.next() {
        length.parse::<i32>().ok()?;
        let k = tokens.next()?.parse::<f64>().ok()?;
        karlin_k.push(k);
    }
    Some(RpsAux {
        matrix,
        gap_open,
        gap_extend,
        scale,
        karlin_k,
    })
}

/// The path that NCBI's errors name: the volume's directory made absolute, with links and `..`
/// resolved (`CSeqDB::FindVolumePaths`), and the volume's name without extension, as the
/// system's bytes (NCBI's `std::string`).
pub(crate) fn canonical_volume_path(base: &Path) -> Vec<u8> {
    let name = base.file_name().map(PathBuf::from).unwrap_or_default();
    let dir = match base.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let path = std::fs::canonicalize(&dir)
        .map(|dir| dir.join(&name))
        .unwrap_or_else(|_| base.to_path_buf());
    crate::utils::os_bytes::path_bytes(&path).into_owned()
}

/// A `.freq` file read at offsets: the file itself where the platform has positioned reads, else
/// (WASI, which refuses rpsblast anyway) its bytes in the heap.
enum FreqSource {
    #[cfg(any(unix, windows))]
    File(std::fs::File),
    #[cfg(not(any(unix, windows)))]
    Heap(Vec<u8>),
}

impl FreqSource {
    fn open(path: &Path) -> std::io::Result<Self> {
        #[cfg(any(unix, windows))]
        {
            std::fs::File::open(path).map(FreqSource::File)
        }
        #[cfg(not(any(unix, windows)))]
        {
            std::fs::read(path).map(FreqSource::Heap)
        }
    }

    fn len(&self) -> std::io::Result<usize> {
        #[cfg(any(unix, windows))]
        {
            let FreqSource::File(file) = self;
            Ok(usize::try_from(file.metadata()?.len()).unwrap_or(usize::MAX))
        }
        #[cfg(not(any(unix, windows)))]
        {
            let FreqSource::Heap(bytes) = self;
            Ok(bytes.len())
        }
    }

    /// Fills `buf` from the bytes at `offset`.
    fn read_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            let FreqSource::File(file) = self;
            file.read_exact_at(buf, offset)
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::FileExt;
            let FreqSource::File(file) = self;
            let mut done = 0;
            while done < buf.len() {
                let n = file.seek_read(&mut buf[done..], offset + done as u64)?;
                if n == 0 {
                    return Err(std::io::ErrorKind::UnexpectedEof.into());
                }
                done += n;
            }
            Ok(())
        }
        #[cfg(not(any(unix, windows)))]
        {
            let FreqSource::Heap(bytes) = self;
            let start = usize::try_from(offset).unwrap_or(usize::MAX);
            let end = start.saturating_add(buf.len());
            buf.copy_from_slice(
                bytes
                    .get(start..end)
                    .ok_or(std::io::ErrorKind::UnexpectedEof)?,
            );
            Ok(())
        }
    }
}

/// The hits of a cell with `cell[0] > 0`: the first slice and the overflow run after it (the
/// cell's entries in the cell, or its first and its run), `None` when the run is not inside the
/// overflow array (`BlastRPSScanSubject` reads them in this order).
fn cell_hits<'a>(cell: &'a [i32; 4], overflow: &'a [i32]) -> Option<(&'a [i32], &'a [i32])> {
    if cell[0] <= RPS_HITS_PER_CELL {
        return Some((cell.get(1..=usize::try_from(cell[0]).ok()?)?, &[]));
    }
    let start = usize::try_from(cell[2]).ok()? / std::mem::size_of::<i32>();
    let run = overflow.get(start..start.checked_add(usize::try_from(cell[0]).ok()? - 1)?)?;
    Some((&cell[1..2], run))
}

/// A cell's hits (`cell_hits`) for `check_cells`: whether every entry lies in the rows (`entry - 2
/// < rows_limit` unsigned), and `ungrouped` set when the hits are not in groups of
/// `RPS_GROUP_BUCKETS` buckets in order (the answer of `list_is_grouped`, which a test holds this
/// to). Written without an early exit on a list's values, so that the compiler makes it a vector
/// loop: one compare per hit besides the range test.
#[inline]
fn hits_in_rows_and_groups(
    head: &[i32],
    run: &[i32],
    rows_limit: u32,
    ungrouped: &mut bool,
) -> bool {
    let group = |entry: i32| entry_group(entry, RPS_GROUP_BUCKETS);
    let (mut outside, mut decreasing) = (0u8, 0u8);
    for list in [head, run] {
        if let Some(&first) = list.first() {
            outside |= (entry_q_off(first) >= rows_limit) as u8;
        }
        for pair in list.windows(2) {
            outside |= (entry_q_off(pair[1]) >= rows_limit) as u8;
            decreasing |= (group(pair[1]) < group(pair[0])) as u8;
        }
    }
    if let (Some(&last), Some(&first)) = (head.last(), run.first()) {
        decreasing |= (group(first) < group(last)) as u8;
    }
    *ungrouped |= decreasing != 0;
    outside == 0
}

/// Whether the hits of a cell (`cell_hits`) are in groups of `group_buckets` buckets in
/// non-decreasing order: what `RpsBuckets::fill_group` needs of a table (it makes the hits of one
/// group of buckets at a time from the entries of a cell in turn, so the entries of a group must
/// be a prefix of what is left). A table made by makeprofiledb lists every cell ascending: its
/// lookup appends the hits of the profiles' rows in ascending order
/// (`s_AddPSSMNeighboringWords`, `BlastLookupAddWordHit`) and `x_RPS_DbClose` copies each chain in
/// order. (The groups, not the buckets: inside a group each bucket receives the entries in the
/// cell's order whatever their buckets.)
fn list_is_grouped(head: &[i32], run: &[i32], group_buckets: usize) -> bool {
    let mut previous = 0;
    head.iter().chain(run).all(|&entry| {
        let group = entry_group(entry, group_buckets);
        let ordered = group >= previous;
        previous = group;
        ordered
    })
}

/// The hits of a cell sorted by group of `group_buckets` buckets, each group's hits in the order
/// they had (a stable sort): the hits of one bucket stay in the order of the cell, which is all
/// that NCBI's scan makes of a cell's order (it appends the cell's hits to their buckets), and
/// the sort is idempotent on a grouped list.
fn regroup_cell(hits: &mut [i32], group_buckets: usize) {
    hits.sort_by_key(|&entry| entry_group(entry, group_buckets));
}

/// The cells and overflow array of a table whose cells are not all in groups of `group_buckets`
/// buckets in order, as copies that are, `None` when they are (`list_is_grouped`; the scan's
/// groups are `RPS_GROUP_BUCKETS`, tests read others). The regrouped cells keep their entries'
/// relative order inside each group (`regroup_cell`), so the hits of every bucket stay those
/// NCBI's scan gives it, in its order. A cell of at most `RPS_HITS_PER_CELL` hits is rewritten in
/// place, a longer one gets its first hit in the cell and its run appended to the copy of the
/// overflow array (the cell's offset rewritten; the old run stays where it was, for another cell
/// may share it).
///
/// The copy is the whole overflow array (600 MB for a table of CDD's size, per volume, and again
/// for every window of an rpstblastn run): the file is mapped read-only on every target
/// (`MappedFile`: a read-only mapping on unix and windows, a heap copy elsewhere), so copy-on-write
/// pages would need a private writable mapping on each of them; only a table that makeprofiledb
/// did not write is copied.
pub(crate) fn regroup_cells(
    cells: &[[i32; 4]],
    overflow: &[i32],
    group_buckets: usize,
) -> Result<Option<(Vec<[i32; 4]>, Vec<i32>)>, String> {
    let mut copy: Option<(Vec<[i32; 4]>, Vec<i32>)> = None;
    for (index, cell) in cells.iter().enumerate().take(RPS_BACKBONE_SIZE) {
        if cell[0] < 2 {
            continue;
        }
        let Some((head, run)) = cell_hits(cell, overflow) else {
            continue;
        };
        if list_is_grouped(head, run, group_buckets) {
            continue;
        }
        let mut hits: Vec<i32> = head.iter().chain(run).copied().collect();
        regroup_cell(&mut hits, group_buckets);
        let (cells_copy, overflow_copy) =
            copy.get_or_insert_with(|| (cells.to_vec(), overflow.to_vec()));
        let rewritten = &mut cells_copy[index];
        if cell[0] <= RPS_HITS_PER_CELL {
            rewritten[1..=hits.len()].copy_from_slice(&hits);
        } else {
            rewritten[1] = hits[0];
            rewritten[2] = i32::try_from(overflow_copy.len() * std::mem::size_of::<i32>())
                .map_err(|_| {
                    ".loo overflow array too large to regroup by group of buckets".to_string()
                })?;
            overflow_copy.extend_from_slice(&hits[1..]);
        }
    }
    Ok(copy)
}

/// One volume's `.aux`, `.loo` and `.rps`, checked against the volume's sequences.
pub(crate) struct RpsVolumeFiles {
    pub(crate) aux: RpsAux,
    /// 26 or 28: the alphabet of the `.loo` header, which NCBI's lookup table takes for the row
    /// stride of the `.rps` rows (`RPSLookupTableNew`).
    pub(crate) alphabet: usize,
    loo: MappedFile,
    rps: MappedFile,
    num_profiles: usize,
    /// The byte offset of the backbone in the `.loo` file (`start_of_backbone`).
    backbone_start: usize,
    /// The rows widened to 28 columns, for a 26-column `.rps` (columns 26 and 27 of a row are
    /// the next row's columns 0 and 1, as NCBI's stride-26 rows would read them; 0 for the last
    /// row). No search reads them: NCBI masks the query letters past a 26-letter alphabet to X
    /// (`Blast_MaskUnsupportedAA`, `rps_hooks::rps_search`).
    narrow: Option<Vec<[i32; 28]>>,
    /// The presence bits of the backbone's cells (`lookup->pv`).
    pv: Vec<u64>,
    /// The cells and overflow array with each cell's hits grouped by bucket, when the file's
    /// cells are not (`regroup_cells`; never for a table made by makeprofiledb): `cells` and
    /// `overflow` read them in place of the mapped file's.
    regrouped: Option<(Vec<[i32; 4]>, Vec<i32>)>,
    /// The `.freq` file (frequency ratios), opened by [`Self::open_freq`]: read by `pread`, not
    /// mapped, so that its pages never count in the resident set (a CDD volume's is 81 MB; the
    /// redo of a search touches almost all of it), with its length in bytes.
    freq: Option<(FreqSource, u64)>,
}

impl RpsVolumeFiles {
    /// rpsblast's `open_for` with an `.aux` without a token as a refusal (the tests' files).
    #[cfg(test)]
    pub(crate) fn open(base: &Path, lengths: &[u32]) -> Result<Self, RpsOpenError> {
        Self::open_for(base, lengths, "RPSBLAST").map_err(|fault| match fault {
            RpsVolumeFault::Files(error) => error,
            RpsVolumeFault::NoMatrix => RpsOpenError::Refused(
                "a profile database whose .aux file is not NCBI's text (matrix, gap costs, statistics, scale, then pairs of a length and a K) is not supported by LOSAT's RPSBLAST".to_string(),
            ),
        })
    }

    /// Opens the RPS files of the volume `base` (its path without extension) whose sequence part
    /// holds the profiles' consensus sequences of `lengths` letters (by local ordinal id).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:180-188
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:294-306
    /// (The `.rps` file's check, rps_aux.cpp:338-348, is the same.)
    /// `program` is the name in LOSAT's refusals; an `.aux` without a token is
    /// `RpsVolumeFault::NoMatrix`, NCBI's engine error of each query (review A of PS12a, N2).
    pub(crate) fn open_for(
        base: &Path,
        lengths: &[u32],
        program: &str,
    ) -> Result<Self, RpsVolumeFault> {
        Self::open_checking(base, lengths, program, true, false)
    }

    /// `open_for` of files that `open_for` has checked in this run (series PR, item PS12s, D6):
    /// NCBI's own errors of the files (a missing file, a magic number) and LOSAT's refusals of
    /// their layout come as there, but the table's hits are not read again (`check_cells`):
    /// only their presence bits are made. `regroup`: `open_for` found the cells not grouped by
    /// groups of buckets (`is_regrouped`), so they are regrouped again.
    pub(crate) fn open_checked(
        base: &Path,
        lengths: &[u32],
        program: &str,
        regroup: bool,
    ) -> Result<Self, RpsVolumeFault> {
        Self::open_checking(base, lengths, program, false, regroup)
    }

    /// `open_for`, the table's hits checked (`check_cells`) when `check_hits`; the cells regrouped
    /// (`regroup_cells`) when the check finds them not grouped, or when `regroup`.
    fn open_checking(
        base: &Path,
        lengths: &[u32],
        program: &str,
        check_hits: bool,
        regroup: bool,
    ) -> Result<Self, RpsVolumeFault> {
        if cfg!(target_endian = "big") {
            return Err(RpsOpenError::Refused(format!(
                "a profile database (its files are little endian) on a big-endian platform is not supported by LOSAT's {program}"
            ))
            .into());
        }
        let shown = canonical_volume_path(base);
        let with = |extension: &str| {
            let mut path = base.as_os_str().to_os_string();
            path.push(extension);
            PathBuf::from(path)
        };
        let aux_text = std::fs::read(with(".aux")).map_err(|_| {
            RpsOpenError::Ncbi(
                [
                    b"Cannot open RPS-BLAST auxiliary file: ",
                    &shown[..],
                    b".aux",
                ]
                .concat(),
            )
        })?;
        let map = |extension: &str| {
            MappedFile::open(&with(extension)).map_err(|_| {
                RpsOpenError::Ncbi(
                    [
                        b"Cannot memory map RPS-BLAST database file: ",
                        &shown[..],
                        extension.as_bytes(),
                    ]
                    .concat(),
                )
            })
        };
        let magic_error = |extension: &str| {
            RpsOpenError::Ncbi(
                [
                    b"RPS BLAST profile file (",
                    &shown[..],
                    extension.as_bytes(),
                    b") is either corrupt or constructed for an incompatible architecture",
                ]
                .concat(),
            )
        };
        let empty = |extension: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {extension} file is shorter than its header (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}"
            ))
        };
        let loo = map(".loo")?;
        let loo_magic = *loo
            .i32s(0, 1)
            .ok_or_else(|| empty(".loo"))?
            .first()
            .unwrap();
        if loo_magic != RPS_MAGIC_NUM && loo_magic != RPS_MAGIC_NUM_28 {
            return Err(magic_error(".loo").into());
        }
        let rps = map(".rps")?;
        let rps_magic = *rps
            .i32s(0, 1)
            .ok_or_else(|| empty(".rps"))?
            .first()
            .unwrap();
        if rps_magic != RPS_MAGIC_NUM && rps_magic != RPS_MAGIC_NUM_28 {
            return Err(magic_error(".rps").into());
        }
        let refused = |what: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {what} is not supported by LOSAT's {program}"
            ))
        };
        let aux_text = std::str::from_utf8(&aux_text).ok();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:205-209
        // An `.aux` of no token leaves the matrix name empty (the stream fails without throwing):
        // the files open, and the batch's score blocks fail (`RpsVolumeFault::NoMatrix`).
        if aux_text.is_some_and(|text| text.split_ascii_whitespace().next().is_none()) {
            return Err(RpsVolumeFault::NoMatrix);
        }
        let aux = aux_text
            .and_then(parse_aux)
            .ok_or_else(|| refused(".aux file is not NCBI's text (matrix, gap costs, statistics, scale, then pairs of a length and a K)"))?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:139-149
        if loo_magic != rps_magic {
            return Err(
                refused(".loo and .rps files are of different alphabets (magic numbers)").into(),
            );
        }
        let alphabet = if loo_magic == RPS_MAGIC_NUM { 26 } else { 28 };
        let num_profiles = rps
            .i32s(4, 1)
            .and_then(|n| usize::try_from(n[0]).ok())
            .ok_or_else(|| empty(".rps"))?;
        if num_profiles != lengths.len() {
            return Err(refused(&format!(
                ".rps file holds {num_profiles} profiles for {} sequences",
                lengths.len()
            ))
            .into());
        }
        let offsets = rps.i32s(8, num_profiles + 1).ok_or_else(|| empty(".rps"))?;
        // Each profile's rows: its columns and one sentinel row. A profile whose `.psq` letters
        // (the `.pin` length) are more or fewer than its columns opens: NCBI reads its letters
        // only in the traceback, where LOSAT asks each read of them
        // (`algorithm::rps_search::profile`; series PR, item PS10e).
        if offsets[0] != 0 {
            return Err(refused(".rps profiles do not start at row 0").into());
        }
        if let Some(i) = (0..num_profiles).find(|&i| offsets[i + 1] <= offsets[i]) {
            return Err(refused(&format!(
                ".rps profile {i} has no rows, not even its sentinel row (its offsets {} and {} do not ascend)",
                offsets[i],
                offsets[i + 1]
            ))
            .into());
        }
        let num_rows =
            usize::try_from(offsets[num_profiles]).map_err(|_| refused(".rps row count"))?;
        // NCBI reads the rows that the offsets name and nothing after them: bytes after the last
        // row are searched as NCBI searches them (not read).
        let rows_start = 8 + 4 * (num_profiles + 1);
        if rps.len()
            < num_rows
                .checked_mul(alphabet * 4)
                .and_then(|bytes| bytes.checked_add(rows_start))
                .ok_or_else(|| refused(".rps row count"))?
        {
            return Err(refused(
                ".rps file is shorter than its rows (header, offsets, one row per position)",
            )
            .into());
        }
        if aux.karlin_k.len() < num_profiles {
            return Err(refused(&format!(
                ".aux file has {} Karlin-Altschul K values for {num_profiles} profiles",
                aux.karlin_k.len()
            ))
            .into());
        }
        if aux.karlin_k[..num_profiles]
            .iter()
            .any(|k| !(k.is_finite() && *k > 0.0))
        {
            return Err(refused(".aux file has a K value that is not positive").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:154-160
        let header = loo
            .i32s(0, LOO_HEADER_BYTES / 4)
            .ok_or_else(|| empty(".loo"))?;
        let backbone_start = usize::try_from(header[8])
            .ok()
            .filter(|&start| start >= LOO_HEADER_BYTES && start % 4 == 0)
            .ok_or_else(|| refused(".loo backbone offset"))?;
        let mut files = Self {
            aux,
            alphabet,
            loo,
            rps,
            num_profiles,
            backbone_start,
            narrow: None,
            pv: vec![0u64; RPS_BACKBONE_SIZE / 64],
            regrouped: None,
            freq: None,
        };
        if files.cells().is_none() || files.overflow().is_none() {
            return Err(refused(".loo file is shorter than its backbone").into());
        }
        let ungrouped = files
            .check_cells(num_rows, check_hits)
            .map_err(|what| refused(&what))?;
        if ungrouped || regroup {
            let (cells, overflow) = (
                files.cells().expect("checked above"),
                files.overflow().expect("checked above"),
            );
            files.regrouped =
                regroup_cells(cells, overflow, RPS_GROUP_BUCKETS).map_err(|what| refused(&what))?;
        }
        if alphabet == 26 {
            files.narrow = Some(files.widen_rows(num_rows));
        }
        Ok(files)
    }

    /// Opens the volume's `.freq` file (the frequency ratios of the composition-based
    /// statistics), after the other files as NCBI does. A missing file, an empty one or a wrong
    /// magic number is NCBI's error (the file's own text is replaced by the message below, which
    /// names `rps_dbname`: `database_name`, the `-db` value of a database of one volume, else the
    /// volume's path; `RpsVolumeSetup::with_database`);
    /// a file shorter than the rows that the `.rps` offsets name (NCBI reads past it) is
    /// LOSAT's refusal.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:594-601
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:485-499
    /// The rows start after the header and `num_profiles + 1` offsets, with the profile count
    /// of the `.rps` file (`freq_header->start_offsets` plus `num_profiles` plus one,
    /// blast_traceback.c:956-957), and their stride is the `.rps` alphabet.
    pub(crate) fn open_freq(
        &mut self,
        base: &Path,
        program: &str,
        database_name: Option<&[u8]>,
    ) -> Result<(), RpsOpenError> {
        let shown = database_name.map_or_else(|| canonical_volume_path(base), <[u8]>::to_vec);
        let missing = || {
            RpsOpenError::Ncbi(
                [
                    &shown[..],
                    b" contains no frequency ratios needed for composition-based statistics.\nPlease disable composition-based statistics when searching against ",
                    &shown[..],
                    b".",
                ]
                .concat(),
            )
        };
        let mut path = base.as_os_str().to_os_string();
        path.push(".freq");
        let freq = FreqSource::open(&PathBuf::from(path)).map_err(|_| missing())?;
        let mut magic = [0u8; 4];
        freq.read_at(&mut magic, 0).map_err(|_| missing())?;
        let magic = i32::from_le_bytes(magic);
        if magic != RPS_MAGIC_NUM && magic != RPS_MAGIC_NUM_28 {
            return Err(missing());
        }
        let num_rows = self.offsets()[self.num_profiles] as usize;
        let needed = num_rows
            .checked_mul(self.alphabet * 4)
            .and_then(|bytes| bytes.checked_add(self.rows_start()));
        let length = freq.len().ok();
        if needed.is_none_or(|needed| length.is_none_or(|len| len < needed)) {
            return Err(RpsOpenError::Refused(
                format!("a profile database whose .freq file is shorter than its rows (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}"),
            ));
        }
        self.freq = Some((freq, length.unwrap_or(0) as u64));
        Ok(())
    }

    /// The sum of the positive maxima of the rows of each profile, read from the `.rps` file
    /// through a small buffer (the mapping would bring the whole file into the resident set).
    /// An alignment uses each row of its profile at most once, so no alignment of the profile
    /// scores more than this; `None` where the file cannot be read.
    ///
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h
    /// (`BlastGapAlignStruct::score` is an `Int4`, as the extensions' scores are: a sum above
    /// `INT4_MAX` overflows in C, which is undefined behavior).
    pub(crate) fn profile_score_bounds(&self, base: &Path) -> Option<Vec<i64>> {
        use std::io::{Read, Seek, SeekFrom};
        let mut path = base.as_os_str().to_os_string();
        path.push(".rps");
        let mut file = std::fs::File::open(PathBuf::from(path)).ok()?;
        file.seek(SeekFrom::Start(self.rows_start() as u64)).ok()?;
        let offsets = self.offsets();
        let row_bytes = self.alphabet * 4;
        const ROWS_PER_READ: usize = 8192;
        let mut buffer = vec![0u8; row_bytes * ROWS_PER_READ];
        let num_rows = offsets[self.num_profiles] as usize;
        let mut bounds = vec![0i64; self.num_profiles];
        let (mut profile, mut row) = (0usize, 0usize);
        while row < num_rows {
            let rows = ROWS_PER_READ.min(num_rows - row);
            file.read_exact(&mut buffer[..rows * row_bytes]).ok()?;
            for chunk in buffer[..rows * row_bytes].chunks_exact(row_bytes) {
                while row >= offsets[profile + 1] as usize {
                    profile += 1;
                }
                let best = chunk
                    .chunks_exact(4)
                    .map(|cell| i32::from_le_bytes([cell[0], cell[1], cell[2], cell[3]]))
                    .max()
                    .unwrap_or(0);
                bounds[profile] += i64::from(best.max(0));
                row += 1;
            }
        }
        Some(bounds)
    }

    /// The frequency ratios of the rows `start..start + count` into `out` (28 per row): each
    /// row takes `alphabet` integers from the `.freq` file at the `.rps` alphabet's stride,
    /// divided by `FREQ_RATIO_SCALE` and padded with 0.0 (a 26-column `.rps` with a 28-column
    /// `.freq` reads rows two integers early, as NCBI does). NCBI reads the rows in the file's
    /// mapping: bytes after the file's end read as the zeros of its last page
    /// (`file::mapped_length`; the row after the rows of a profile whose `.pin` length is above
    /// its columns, `rps_search::profile`), and `false` is returned where a row's bytes lie past
    /// that page (NCBI's read faults).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:985-1006
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_rps.h:83
    pub(crate) fn freq_ratios(
        &self,
        start: usize,
        count: usize,
        bytes: &mut Vec<u8>,
        out: &mut Vec<[f64; 28]>,
    ) -> std::io::Result<bool> {
        out.clear();
        let Some((freq, file_length)) = self.freq.as_ref() else {
            return Ok(true);
        };
        let file_length = *file_length;
        // (`bytes` only grows: a list's rows are not zeroed again for the next one.)
        let wanted = count * self.alphabet * 4;
        if bytes.len() < wanted {
            bytes.resize(wanted, 0);
        }
        let bytes = &mut bytes[..wanted];
        let offset = (self.rows_start() + start * self.alphabet * 4) as u64;
        let end = offset + wanted as u64;
        if end > super::file::mapped_length(file_length) {
            return Ok(false);
        }
        let in_file = usize::try_from(file_length.saturating_sub(offset))
            .unwrap_or(usize::MAX)
            .min(wanted);
        let (read, zeros) = bytes.split_at_mut(in_file);
        freq.read_at(read, offset)?;
        zeros.fill(0);
        // Each row is written where it stays (its zeros, then its ratios over them), not built
        // aside and copied in.
        out.reserve(count);
        for row in bytes.chunks_exact(self.alphabet * 4) {
            out.push([0.0; 28]);
            let Some(ratios) = out.last_mut() else {
                break;
            };
            for (ratio, value) in ratios.iter_mut().zip(row.as_chunks::<4>().0) {
                *ratio = f64::from(i32::from_le_bytes(*value)) / 1_000_000_000.0;
            }
        }
        Ok(true)
    }

    /// Every hit of the table names a position whose word lies in the concatenated rows
    /// (`entry - 2` from 0, `entry` before the last row), every overflow run lies in the
    /// overflow array, no cell holds more hits than a scan can take; the presence bits. Returns
    /// whether some cell lists its hits not in groups of buckets in order (`list_is_grouped`; false when
    /// `check_hits` is not set).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:162-173
    fn check_cells(&mut self, num_rows: usize, check_hits: bool) -> Result<bool, String> {
        let cells = self.cells().expect("checked by the caller");
        let overflow = self.overflow().expect("checked by the caller");
        // `entry >= RPS_WORDSIZE - 1 && entry < num_rows` as one unsigned compare of `entry - 2`
        // (a negative entry wraps above any bound under 2^31 - 2).
        let rows_limit = u32::try_from(num_rows.saturating_sub(RPS_WORDSIZE - 1))
            .unwrap_or(u32::MAX)
            .min(i32::MAX as u32 - 1);
        let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
        let mut ungrouped = false;
        for (index, cell) in cells.iter().enumerate().take(RPS_BACKBONE_SIZE) {
            let num_used = cell[0];
            if num_used <= 0 {
                continue;
            }
            pv[index / 64] |= 1u64 << (index % 64);
            if !check_hits {
                continue;
            }
            if num_used > RPS_MAX_SCAN_HITS {
                return Err(format!(
                    ".loo word {index} has {num_used} hits, more than one scan takes (NCBI BLAST+ does not end)"
                ));
            }
            // The cell's hits (`cell_hits`), each inside the rows, and whether their groups of
            // buckets never decrease (`list_is_grouped`; one compare per hit in this walk).
            let Some((head, run)) = cell_hits(cell, overflow) else {
                return Err(format!(
                    ".loo word {index} names a position outside the profiles' rows"
                ));
            };
            let all_in_rows = hits_in_rows_and_groups(head, run, rows_limit, &mut ungrouped);
            if !all_in_rows {
                return Err(format!(
                    ".loo word {index} names a position outside the profiles' rows"
                ));
            }
        }
        self.pv = pv;
        Ok(ungrouped)
    }

    /// Whether the cells are regrouped (`regroup_cells`): the table lists some cell's hits not in
    /// groups of buckets in order.
    pub(crate) fn is_regrouped(&self) -> bool {
        self.regrouped.is_some()
    }

    /// The rows of a 26-column file widened to 28 columns.
    fn widen_rows(&self, num_rows: usize) -> Vec<[i32; 28]> {
        let flat = self
            .rps
            .i32s(self.rows_start(), num_rows * 26)
            .expect("checked by `open`");
        (0..num_rows)
            .map(|i| {
                let mut row = [0i32; 28];
                row[..26].copy_from_slice(&flat[i * 26..i * 26 + 26]);
                if i + 1 < num_rows {
                    row[26..].copy_from_slice(&flat[(i + 1) * 26..(i + 1) * 26 + 2]);
                }
                row
            })
            .collect()
    }

    fn rows_start(&self) -> usize {
        8 + 4 * (self.num_profiles + 1)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:184-190 (cited at
    // `rps_search::RpsRows`)
    /// NCBI's row after the concatenated rows (`rps_pssm[num_pssm_rows]`, a pointer to the bytes
    /// that follow the last row in the mapped `.rps` file, at the row stride of the alphabet) as
    /// each letter reads it: the file's own bytes there, then the zeros of its last memory page
    /// (`file::mapped_length`), `None` for a letter whose four bytes lie past that page (NCBI's
    /// read of them faults, SIGBUS: the file ends on a page boundary). The ungapped extensions
    /// read it where the unit's concatenated length `Lc` is above the rows
    /// (`rps_search::RpsPastRows`).
    pub(crate) fn row_after_rows(&self) -> [Option<i32>; 28] {
        let bytes = self.rps.bytes();
        let mapped = super::file::mapped_length(bytes.len() as u64);
        let num_rows = self.offsets()[self.num_profiles] as usize;
        let start = (self.rows_start() + num_rows * self.alphabet * 4) as u64;
        std::array::from_fn(|letter| {
            let at = start + 4 * letter as u64;
            (at + 4 <= mapped).then(|| {
                let mut word = [0u8; 4];
                for (k, byte) in word.iter_mut().enumerate() {
                    *byte = bytes.get((at as usize) + k).copied().unwrap_or(0);
                }
                i32::from_ne_bytes(word)
            })
        })
    }

    /// The `num_profiles + 1` row offsets of the profiles (`start_offsets`).
    pub(crate) fn offsets(&self) -> &[i32] {
        self.rps
            .i32s(8, self.num_profiles + 1)
            .expect("checked by `open`")
    }

    /// The concatenated rows, 28 columns each (one per consensus letter and one sentinel row
    /// after each profile): the mapped rows of a 28-column file, the widened copy otherwise.
    pub(crate) fn rows(&self) -> &[[i32; 28]] {
        match &self.narrow {
            Some(rows) => rows,
            None => {
                let num_rows = self.offsets()[self.num_profiles] as usize;
                self.rps
                    .rows::<28>(self.rows_start(), num_rows)
                    .expect("checked by `open`")
            }
        }
    }

    /// The `backbone_size + 1` cells of the table.
    pub(crate) fn cells(&self) -> Option<&[[i32; 4]]> {
        if let Some((cells, _)) = &self.regrouped {
            return Some(cells);
        }
        self.loo
            .rows::<4>(self.backbone_start, RPS_BACKBONE_SIZE + 1)
    }

    /// The overflow array: the file after the cells.
    pub(crate) fn overflow(&self) -> Option<&[i32]> {
        if let Some((_, overflow)) = &self.regrouped {
            return Some(overflow);
        }
        let start = self.backbone_start + (RPS_BACKBONE_SIZE + 1) * 16;
        let count = self.loo.len().checked_sub(start)? / 4;
        self.loo.i32s(start, count)
    }

    /// The presence bits of the cells with hits.
    pub(crate) fn pv(&self) -> &[u64] {
        &self.pv
    }
}

#[cfg(test)]
#[path = "rps_tests.rs"]
mod tests;

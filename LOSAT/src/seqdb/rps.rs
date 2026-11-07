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
//! ```c++
//!     if (flags & fAuxInfoFile) {
//!         m_AuxFile.Reset(new CRpsAuxFile(path));
//!         ...
//!     if (flags & fLookupTableFile) {
//!         m_LutFile.Reset(new CRpsLookupTblFile(path));
//!         ...
//!     if (flags & fPssmFile) {
//!         m_PssmFile.Reset(new CRpsPssmFile(path));
//! ```

use std::path::{Path, PathBuf};

use crate::utils::mapped_file::MappedFile;

/// `RPS_MAGIC_NUM`: the files of a 26-letter alphabet.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_rps.h:43-44
/// ```c
/// #define RPS_MAGIC_NUM 0x1e16    /**< RPS blast version number */
/// #define RPS_MAGIC_NUM_28 0x1e17    /**< Version number for 28-letter alphabet */
/// ```
pub(crate) const RPS_MAGIC_NUM: i32 = 0x1e16;
/// `RPS_MAGIC_NUM_28`: the files of the 28-letter alphabet (makeprofiledb 2.17 writes these).
pub(crate) const RPS_MAGIC_NUM_28: i32 = 0x1e17;

/// `RPS_HITS_PER_CELL`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:330-344
/// ```c
/// #define RPS_HITS_PER_CELL 3 /**< maximum number of hits in an RPS backbone
///                                  cell; ...
/// typedef struct RPSBackboneCell {
///     Int4 num_used;                   /**< number of hits in this cell */
///     Int4 entries[RPS_HITS_PER_CELL]; /**< if the number of hits in this cell
///                                           is RPS_HITS_PER_CELL or less, all
///                                           hits go into this array. Otherwise,
///                                           the first hit in the list goes into
///                                           element 0 of the array, and element 1
///                                           contains the byte offset into the
///                                           overflow array where the list of the
///                                           remaining hits begins */
/// } RPSBackboneCell;
/// ```
pub(crate) const RPS_HITS_PER_CELL: i32 = 3;

/// The word size of the prebuilt table (`BLAST_WORDSIZE_PROT`) and its letters' bits
/// (`ilog2(alphabet_size) + 1`, 5 for 26 and 28 letters): `backbone_size = 1 << 15`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:148-153
/// ```c
///     lookup->wordsize = BLAST_WORDSIZE_PROT;
///     lookup->charsize = ilog2(lookup->alphabet_size) + 1;
///     lookup->backbone_size = 1 << (lookup->wordsize * lookup->charsize);
///     lookup->mask = lookup->backbone_size - 1;
/// ```
pub(crate) const RPS_WORDSIZE: usize = 3;
pub(crate) const RPS_CHARSIZE: usize = 5;
pub(crate) const RPS_BACKBONE_SIZE: usize = 1 << (RPS_WORDSIZE * RPS_CHARSIZE);

/// The most hits that one scan of a query collects (`max_hits` of `BlastRPSScanSubject`). A cell
/// with more hits than this would never fit, and NCBI's scan would not advance past its word.
pub(crate) const RPS_MAX_SCAN_HITS: i32 = 4_000_000;

/// The size of the `.loo` header (`BlastRPSLookupFileHeader`, ten `Int4`).
const LOO_HEADER_BYTES: usize = 40;

/// Why a volume's RPS files cannot be searched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum RpsOpenError {
    /// NCBI's own error (`CBlastException` with `eRpsInit`): the text after `BLAST engine
    /// error: `.
    Ncbi(String),
    /// LOSAT's refusal of files that NCBI would read outside their bounds or that LOSAT does not
    /// search.
    Refused(String),
}

/// The contents of the `.aux` file.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:193-225
/// ```c++
///     // Read the data
///     input >> matrix;
///     input >> gap_open;
///     input >> gap_extend;
///     input >> ignore_me_d; // corresponds to BlastRPSAuxInfo::ungapped_k
///     input >> ignore_me_d; // corresponds to BlastRPSAuxInfo::ungapped_h
///     input >> ignore_me_i; // corresponds to BlastRPSAuxInfo::max_db_seq_length
///     input >> ignore_me_i; // corresponds to BlastRPSAuxInfo::db_length
///     input >> scale_factor;
///     while (input) {
///         int seq_size;   // unused value
///         double k;
///         input >> seq_size;
///         input >> k;
///         karlinK.push_back(k);
///     }
/// ```
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
/// resolved (`CSeqDB::FindVolumePaths`), and the volume's name without extension.
pub(crate) fn canonical_volume_path(base: &Path) -> String {
    let name = base.file_name().map(PathBuf::from).unwrap_or_default();
    let dir = match base.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir.to_path_buf(),
        _ => PathBuf::from("."),
    };
    std::fs::canonicalize(&dir)
        .map(|dir| dir.join(&name))
        .unwrap_or_else(|_| base.to_path_buf())
        .to_string_lossy()
        .into_owned()
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
    /// the next row's columns 0 and 1, as NCBI's stride-26 rows read them; 0 for the last row).
    narrow: Option<Vec<[i32; 28]>>,
    /// The presence bits of the backbone's cells (`lookup->pv`).
    pv: Vec<u64>,
}

impl RpsVolumeFiles {
    /// Opens the RPS files of the volume `base` (its path without extension) whose sequence part
    /// holds the profiles' consensus sequences of `lengths` letters (by local ordinal id).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:180-188
    /// ```c++
    ///     const string file2open(filename_no_extn + kExtension);
    ///     CNcbiIfstream auxfile(file2open.c_str());
    ///     if (auxfile.bad() || auxfile.fail()) {
    ///         NCBI_THROW(CBlastException, eRpsInit,
    ///                    "Cannot open RPS-BLAST auxiliary file: " + file2open);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:296-306
    /// ```c++
    ///     m_Data = (BlastRPSLookupFileHeader*) m_MmappedFile->GetPtr();
    ///     if (m_Data->magic_number != RPS_MAGIC_NUM &&
    ///         m_Data->magic_number != RPS_MAGIC_NUM_28) {
    ///         m_Data = NULL;
    ///         NCBI_THROW(CBlastException, eRpsInit,
    ///                "RPS BLAST profile file (" + filename_no_extn + kExtension +
    ///                ") is either corrupt or constructed for an incompatible "
    ///                "architecture");
    ///     }
    /// ```
    /// (The `.rps` file's check, rps_aux.cpp:338-348, is the same.)
    pub(crate) fn open(base: &Path, lengths: &[u32]) -> Result<Self, RpsOpenError> {
        if cfg!(target_endian = "big") {
            return Err(RpsOpenError::Refused(
                "a profile database (its files are little endian) on a big-endian platform is not supported by LOSAT's RPSBLAST".to_string(),
            ));
        }
        let shown = canonical_volume_path(base);
        let with = |extension: &str| {
            let mut path = base.as_os_str().to_os_string();
            path.push(extension);
            PathBuf::from(path)
        };
        let aux_text = std::fs::read(with(".aux")).map_err(|_| {
            RpsOpenError::Ncbi(format!("Cannot open RPS-BLAST auxiliary file: {shown}.aux"))
        })?;
        let map = |extension: &str| {
            MappedFile::open(&with(extension)).map_err(|_| {
                RpsOpenError::Ncbi(format!(
                    "Cannot memory map RPS-BLAST database file: {shown}{extension}"
                ))
            })
        };
        let magic_error = |extension: &str| {
            RpsOpenError::Ncbi(format!(
                "RPS BLAST profile file ({shown}{extension}) is either corrupt or constructed for an incompatible architecture"
            ))
        };
        let empty = |extension: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {extension} file is shorter than its header (NCBI BLAST+ reads past it) is not supported by LOSAT's RPSBLAST"
            ))
        };
        let loo = map(".loo")?;
        let loo_magic = *loo
            .i32s(0, 1)
            .ok_or_else(|| empty(".loo"))?
            .first()
            .unwrap();
        if loo_magic != RPS_MAGIC_NUM && loo_magic != RPS_MAGIC_NUM_28 {
            return Err(magic_error(".loo"));
        }
        let rps = map(".rps")?;
        let rps_magic = *rps
            .i32s(0, 1)
            .ok_or_else(|| empty(".rps"))?
            .first()
            .unwrap();
        if rps_magic != RPS_MAGIC_NUM && rps_magic != RPS_MAGIC_NUM_28 {
            return Err(magic_error(".rps"));
        }
        let refused = |what: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {what} is not supported by LOSAT's RPSBLAST"
            ))
        };
        let aux = std::str::from_utf8(&aux_text)
            .ok()
            .and_then(parse_aux)
            .ok_or_else(|| refused(".aux file is not NCBI's text (matrix, gap costs, statistics, scale, then pairs of a length and a K)"))?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:139-146
        // ```c
        //     /* set the alphabet size. Use hardwired numbers, since we cannot rely on
        //        #define'd constants matching up to the sizes implicit in disk files */
        //     if (lookup_header->magic_number == RPS_MAGIC_NUM)
        //         lookup->alphabet_size = 26;
        //     else
        //         lookup->alphabet_size = 28;
        // ```
        if loo_magic != rps_magic {
            return Err(refused(
                ".loo and .rps files are of different alphabets (magic numbers)",
            ));
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
            )));
        }
        let offsets = rps.i32s(8, num_profiles + 1).ok_or_else(|| empty(".rps"))?;
        // Each profile's rows: its consensus letters and one sentinel row.
        if offsets[0] != 0 {
            return Err(refused(".rps profiles do not start at row 0"));
        }
        for (i, &length) in lengths.iter().enumerate() {
            let rows = i64::from(offsets[i + 1]) - i64::from(offsets[i]);
            if rows != i64::from(length) + 1 {
                return Err(refused(&format!(
                    ".rps profile {i} has {} rows for a consensus sequence of {length} letters",
                    rows - 1
                )));
            }
        }
        let num_rows =
            usize::try_from(offsets[num_profiles]).map_err(|_| refused(".rps row count"))?;
        let rows_start = 8 + 4 * (num_profiles + 1);
        if rps.len()
            != num_rows
                .checked_mul(alphabet * 4)
                .and_then(|bytes| bytes.checked_add(rows_start))
                .ok_or_else(|| refused(".rps row count"))?
        {
            return Err(refused(
                ".rps file length is not that of its rows (header, offsets, one row per position)",
            ));
        }
        if aux.karlin_k.len() < num_profiles {
            return Err(refused(&format!(
                ".aux file has {} Karlin-Altschul K values for {num_profiles} profiles",
                aux.karlin_k.len()
            )));
        }
        if aux.karlin_k[..num_profiles]
            .iter()
            .any(|k| !(k.is_finite() && *k > 0.0))
        {
            return Err(refused(".aux file has a K value that is not positive"));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:154-160
        // ```c
        //     lookup->rps_backbone = (RPSBackboneCell *) ((Uint1 *) lookup_header +
        //                                                 lookup_header->
        //                                                 start_of_backbone);
        //     lookup->overflow =
        //         (Int4 *) ((Uint1 *) lookup_header + lookup_header->start_of_backbone +
        //                   (lookup->backbone_size + 1) * sizeof(RPSBackboneCell));
        // ```
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
        };
        if files.cells().is_none() || files.overflow().is_none() {
            return Err(refused(".loo file is shorter than its backbone"));
        }
        files.check_cells(num_rows).map_err(|what| refused(&what))?;
        if alphabet == 26 {
            files.narrow = Some(files.widen_rows(num_rows));
        }
        Ok(files)
    }

    /// Every hit of the table names a position whose word lies in the concatenated rows
    /// (`entry - 2` from 0, `entry` before the last row), every overflow run lies in the
    /// overflow array, no cell holds more hits than a scan can take; the presence bits.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:162-173
    /// ```c
    ///     pv = lookup->pv = (PV_ARRAY_TYPE *) calloc(
    ///                             (lookup->backbone_size >> PV_ARRAY_BTS),
    ///                             sizeof(PV_ARRAY_TYPE));
    ///
    ///     for (i = 0; i < lookup->backbone_size; i++) {
    ///         if (lookup->rps_backbone[i].num_used > 0) {
    ///             PV_SET(pv, i, PV_ARRAY_BTS);
    ///         }
    ///     }
    /// ```
    fn check_cells(&mut self, num_rows: usize) -> Result<(), String> {
        let cells = self.cells().expect("checked by the caller");
        let overflow = self.overflow().expect("checked by the caller");
        let in_rows = |entry: i32| {
            usize::try_from(entry)
                .ok()
                .is_some_and(|entry| entry >= RPS_WORDSIZE - 1 && entry < num_rows)
        };
        let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
        for (index, cell) in cells.iter().enumerate().take(RPS_BACKBONE_SIZE) {
            let num_used = cell[0];
            if num_used <= 0 {
                continue;
            }
            pv[index / 64] |= 1u64 << (index % 64);
            if num_used > RPS_MAX_SCAN_HITS {
                return Err(format!(
                    ".loo word {index} has {num_used} hits, more than one scan takes (NCBI BLAST+ does not end)"
                ));
            }
            let hits_ok = if num_used <= RPS_HITS_PER_CELL {
                cell[1..=num_used as usize]
                    .iter()
                    .all(|&entry| in_rows(entry))
            } else {
                let start = usize::try_from(cell[2]).ok().map(|bytes| bytes / 4);
                let run = start.and_then(|start| {
                    overflow.get(start..start.checked_add(num_used as usize - 1)?)
                });
                in_rows(cell[1]) && run.is_some_and(|run| run.iter().all(|&entry| in_rows(entry)))
            };
            if !hits_ok {
                return Err(format!(
                    ".loo word {index} names a position outside the profiles' rows"
                ));
            }
        }
        self.pv = pv;
        Ok(())
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
        self.loo
            .rows::<4>(self.backbone_start, RPS_BACKBONE_SIZE + 1)
    }

    /// The overflow array: the file after the cells.
    pub(crate) fn overflow(&self) -> Option<&[i32]> {
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
mod tests {
    use super::*;

    /// A temporary directory of this test (tests of one binary run in parallel).
    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("losat-rps-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn words(values: &[i32]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_le_bytes()).collect()
    }

    /// A volume of two profiles (lengths 3 and 2), with `columns` columns, whose table holds
    /// the word ending at row 2 in a cell and three words in one overflow run.
    pub(crate) fn write_volume(dir: &Path, columns: usize, aux: &str) -> PathBuf {
        let base = dir.join("v");
        let magic = if columns == 26 {
            RPS_MAGIC_NUM
        } else {
            RPS_MAGIC_NUM_28
        };
        let offsets = [0, 4, 7];
        let mut rps = vec![magic, 2];
        rps.extend_from_slice(&offsets);
        for row in 0..7 {
            for column in 0..columns {
                rps.push((row * 100 + column) as i32);
            }
        }
        std::fs::write(base.with_extension("rps"), words(&rps)).unwrap();
        let mut loo = vec![magic, 1, 0, 0, 0, 0, 0, 0, 40, 0];
        let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
        cells[5] = [1, 2, 0, 0];
        cells[9] = [4, 2, 4, 0];
        for cell in &cells {
            loo.extend_from_slice(cell);
        }
        loo.extend_from_slice(&[99, 3, 5, 6]);
        std::fs::write(base.with_extension("loo"), words(&loo)).unwrap();
        std::fs::write(base.with_extension("aux"), aux).unwrap();
        base
    }

    const AUX: &str =
        "BLOSUM62\n11\n1\n0.000000e+00\n0.000000e+00\n0\n0\n100.000000\n3\n4.1e-02\n2\n3.9e-02\n";

    #[test]
    fn the_aux_text_is_read_as_the_streams_do() {
        let aux = parse_aux(AUX).unwrap();
        assert_eq!(aux.matrix, "BLOSUM62");
        assert_eq!((aux.gap_open, aux.gap_extend, aux.scale), (11, 1, 100.0));
        assert_eq!(aux.karlin_k, vec![0.041, 0.039]);
        // Any white space separates the tokens.
        assert_eq!(parse_aux(&AUX.replace('\n', " \t")).unwrap(), aux);
        // A pair without its K, and a token of the wrong kind.
        assert!(parse_aux(&format!("{AUX}7\n")).is_none());
        assert!(parse_aux(&AUX.replace("\n11\n", "\neleven\n")).is_none());
    }

    #[test]
    fn a_volume_of_28_columns_is_viewed_in_place() {
        let dir = scratch("v28");
        let base = write_volume(&dir, 28, AUX);
        let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        assert_eq!(files.alphabet, 28);
        assert_eq!(files.offsets(), &[0, 4, 7]);
        assert_eq!(files.rows().len(), 7);
        assert_eq!(files.rows()[5][27], 527);
        assert_eq!(files.cells().unwrap()[9], [4, 2, 4, 0]);
        assert_eq!(files.overflow().unwrap(), &[99, 3, 5, 6]);
        assert_eq!(files.pv()[0], (1 << 5) | (1 << 9));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_volume_of_26_columns_reads_the_next_rows_first_columns() {
        let dir = scratch("v26");
        let base = write_volume(&dir, 26, AUX);
        let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        assert_eq!(files.alphabet, 26);
        let rows = files.rows();
        assert_eq!(rows[2][25], 225);
        // Letters 26 and 27 read the next row's columns 0 and 1 (stride 26).
        assert_eq!((rows[2][26], rows[2][27]), (300, 301));
        assert_eq!((rows[6][26], rows[6][27]), (0, 0));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn ncbis_errors_come_in_ncbis_file_order() {
        let dir = scratch("errors");
        let base = write_volume(&dir, 28, AUX);
        let shown = canonical_volume_path(&base);
        let open = || RpsVolumeFiles::open(&base, &[3, 2]).err();
        std::fs::write(base.with_extension("rps"), words(&[7, 2])).unwrap();
        assert_eq!(
            open(),
            Some(RpsOpenError::Ncbi(format!(
                "RPS BLAST profile file ({shown}.rps) is either corrupt or constructed for an incompatible architecture"
            )))
        );
        std::fs::remove_file(base.with_extension("loo")).unwrap();
        assert_eq!(
            open(),
            Some(RpsOpenError::Ncbi(format!(
                "Cannot memory map RPS-BLAST database file: {shown}.loo"
            )))
        );
        std::fs::remove_file(base.with_extension("aux")).unwrap();
        assert_eq!(
            open(),
            Some(RpsOpenError::Ncbi(format!(
                "Cannot open RPS-BLAST auxiliary file: {shown}.aux"
            )))
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn inconsistent_files_are_refused() {
        let dir = scratch("refused");
        let base = write_volume(&dir, 28, AUX);
        // The consensus lengths of the sequence part differ.
        assert!(matches!(
            RpsVolumeFiles::open(&base, &[3, 3]),
            Err(RpsOpenError::Refused(_))
        ));
        // Too few K values.
        std::fs::write(base.with_extension("aux"), AUX.replace("2\n3.9e-02\n", "")).unwrap();
        assert!(matches!(
            RpsVolumeFiles::open(&base, &[3, 2]),
            Err(RpsOpenError::Refused(_))
        ));
        // An empty .rps (NCBI crashes).
        let base = write_volume(&dir, 28, AUX);
        std::fs::write(base.with_extension("rps"), b"").unwrap();
        assert!(matches!(
            RpsVolumeFiles::open(&base, &[3, 2]),
            Err(RpsOpenError::Refused(_))
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }
}

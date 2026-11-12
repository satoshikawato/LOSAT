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
    /// The `.freq` file (frequency ratios), opened by [`Self::open_freq`]: read by `pread`, not
    /// mapped, so that its pages never count in the resident set (a CDD volume's is 81 MB; the
    /// redo of a search touches almost all of it).
    freq: Option<FreqSource>,
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
        Self::open_for(base, lengths, "RPSBLAST").map_err(|fault| match fault {
            RpsVolumeFault::Files(error) => error,
            RpsVolumeFault::NoMatrix => RpsOpenError::Refused(
                "a profile database whose .aux file is not NCBI's text (matrix, gap costs, statistics, scale, then pairs of a length and a K) is not supported by LOSAT's RPSBLAST".to_string(),
            ),
        })
    }

    /// `open` for the program `program` (the name in LOSAT's refusals), which tells an `.aux`
    /// without a token apart (`RpsVolumeFault::NoMatrix`) from the other faults of the files.
    pub(crate) fn open_for(
        base: &Path,
        lengths: &[u32],
        program: &str,
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:205-208
        // ```c++
        //     // Read the data
        //     input >> matrix;
        //     input >> gap_open;
        //     input >> gap_extend;
        // ```
        // An `.aux` of no token leaves the matrix name empty (the stream fails without throwing):
        // the files open, and the batch's score blocks fail (`RpsVolumeFault::NoMatrix`).
        if aux_text.is_some_and(|text| text.split_ascii_whitespace().next().is_none()) {
            return Err(RpsVolumeFault::NoMatrix);
        }
        let aux = aux_text
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
        // Each profile's rows: its consensus letters and one sentinel row.
        if offsets[0] != 0 {
            return Err(refused(".rps profiles do not start at row 0").into());
        }
        for (i, &length) in lengths.iter().enumerate() {
            let rows = i64::from(offsets[i + 1]) - i64::from(offsets[i]);
            if rows != i64::from(length) + 1 {
                return Err(refused(&format!(
                    ".rps profile {i} has {} rows for a consensus sequence of {length} letters",
                    rows - 1
                ))
                .into());
            }
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
            freq: None,
        };
        if files.cells().is_none() || files.overflow().is_none() {
            return Err(refused(".loo file is shorter than its backbone").into());
        }
        files.check_cells(num_rows).map_err(|what| refused(&what))?;
        if alphabet == 26 {
            files.narrow = Some(files.widen_rows(num_rows));
        }
        Ok(files)
    }

    /// Opens the volume's `.freq` file (the frequency ratios of the composition-based
    /// statistics), after the other files as NCBI does. A missing file, an empty one or a wrong
    /// magic number is NCBI's error (the file's own text is replaced by the message below);
    /// a file shorter than the rows that the `.rps` offsets name (NCBI reads past it) is
    /// LOSAT's refusal.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:568-584
    /// ```c++
    ///     if (flags & fFreqRatiosFile) {
    ///         try {
    ///             // read frequency ratios data
    ///             m_FreqRatiosFile.Reset(new CRpsFreqRatiosFile(path));
    ///         } catch (const CBlastException& e) {
    ///         	string msg = rps_dbname + " contains no frequency ratios needed for composition-based statistics.\n" \
    ///         			     "Please disable composition-based statistics when searching against " + rps_dbname + ".";
    ///             NCBI_RETHROW(e, CBlastException, eRpsInit, msg);
    ///         }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:485-499
    /// ```c++
    ///     m_Data = (BlastRPSFreqRatiosHeader *) m_MmappedFile->GetPtr();
    ///     if (m_Data->magic_number != RPS_MAGIC_NUM &&
    ///             m_Data->magic_number != RPS_MAGIC_NUM_28) {
    ///             m_Data = NULL;
    ///             NCBI_THROW(CBlastException, eRpsInit,
    ///                    "RPS BLAST freq ratios file (" + filename_no_extn + kExtension +
    ///                    ") is either corrupt or constructed for an incompatible "
    ///                    "architecture");
    /// ```
    /// The rows start after the header and `num_profiles + 1` offsets, with the profile count
    /// of the `.rps` file (`freq_header->start_offsets + num_profiles + 1`,
    /// blast_traceback.c:956-957), and their stride is the `.rps` alphabet.
    pub(crate) fn open_freq(&mut self, base: &Path) -> Result<(), RpsOpenError> {
        let shown = canonical_volume_path(base);
        let missing = || {
            RpsOpenError::Ncbi(format!(
                "{shown} contains no frequency ratios needed for composition-based statistics.\nPlease disable composition-based statistics when searching against {shown}."
            ))
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
        if needed.is_none_or(|needed| freq.len().map_or(true, |len| len < needed)) {
            return Err(RpsOpenError::Refused(
                "a profile database whose .freq file is shorter than its rows (NCBI BLAST+ reads past it) is not supported by LOSAT's RPSBLAST".to_string(),
            ));
        }
        self.freq = Some(freq);
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
    /// `.freq` reads rows two integers early, as NCBI does).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:985-1006
    /// ```c
    ///    for (ic=0; ic<ncol; ic++) {
    ///        for (ir=0; ir<psi_matrix->pssm->nrows; ir++) {
    ///            psi_matrix->freq_ratios[ic][ir] = (double) freq[ic][ir]/ (double) FREQ_RATIO_SCALE;
    ///        }
    ///        for (;  ir<BLASTAA_SIZE; ir++) {
    ///            psi_matrix->freq_ratios[ic][ir] = 0.0;
    ///        }
    ///    }
    /// ```
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_rps.h:91
    /// ```c
    /// #define FREQ_RATIO_SCALE 1000000000
    /// ```
    pub(crate) fn freq_ratios(
        &self,
        start: usize,
        count: usize,
        bytes: &mut Vec<u8>,
        out: &mut Vec<[f64; 28]>,
    ) -> std::io::Result<()> {
        out.clear();
        let Some(freq) = self.freq.as_ref() else {
            return Ok(());
        };
        // (`bytes` only grows: a list's rows are not zeroed again for the next one.)
        let wanted = count * self.alphabet * 4;
        if bytes.len() < wanted {
            bytes.resize(wanted, 0);
        }
        let bytes = &mut bytes[..wanted];
        freq.read_at(
            bytes,
            (self.rows_start() + start * self.alphabet * 4) as u64,
        )?;
        out.extend(bytes.chunks_exact(self.alphabet * 4).map(|row| {
            let mut ratios = [0.0f64; 28];
            for (ratio, value) in ratios.iter_mut().zip(row.chunks_exact(4)) {
                *ratio = f64::from(i32::from_le_bytes([value[0], value[1], value[2], value[3]]))
                    / 1_000_000_000.0;
            }
            ratios
        }));
        Ok(())
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

    /// A `.freq` file of `columns` columns for the two profiles of `write_volume` (rows 0 to
    /// 6): the ratio of row `r`, column `c` is `r * 1000 + c` over 10^9.
    fn write_freq(base: &Path, columns: usize) {
        let magic = if columns == 26 {
            RPS_MAGIC_NUM
        } else {
            RPS_MAGIC_NUM_28
        };
        let mut freq = vec![magic, 2, 0, 4, 7];
        for row in 0..7 {
            for column in 0..columns {
                freq.push((row * 1000 + column) as i32);
            }
        }
        std::fs::write(base.with_extension("freq"), words(&freq)).unwrap();
    }

    #[test]
    fn the_freq_rows_have_the_rps_alphabets_stride() {
        let dir = scratch("freq");
        let mut ratios = Vec::new();
        let mut bytes = Vec::new();
        // 28 columns on both sides: rows as written, ÷ 10^9.
        let base = write_volume(&dir, 28, AUX);
        write_freq(&base, 28);
        let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        files.open_freq(&base).unwrap();
        files.freq_ratios(4, 2, &mut bytes, &mut ratios).unwrap();
        assert_eq!(ratios.len(), 2);
        assert_eq!(ratios[0][0], 4000.0 / 1e9);
        assert_eq!(ratios[1][27], 5027.0 / 1e9);
        // A 26-column `.rps` with a 28-column `.freq` reads rows two integers early: row r
        // starts at integer 26 r of the rows (NCBI's drift), and columns 26 and 27 are 0.0.
        let base = write_volume(&dir, 26, AUX);
        write_freq(&base, 28);
        let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        files.open_freq(&base).unwrap();
        files.freq_ratios(2, 1, &mut bytes, &mut ratios).unwrap();
        let flat: Vec<i32> = (0..7)
            .flat_map(|row| (0..28).map(move |column| row * 1000 + column))
            .collect();
        assert_eq!(ratios[0][0], f64::from(flat[52]) / 1e9);
        assert_eq!(ratios[0][25], f64::from(flat[77]) / 1e9);
        assert_eq!((ratios[0][26], ratios[0][27]), (0.0, 0.0));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_score_bound_of_a_profile_sums_its_rows_positive_maxima() {
        let dir = scratch("bounds");
        let base = write_volume(&dir, 28, AUX);
        let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        // Row r holds r * 100 + c: the maximum of a row is r * 100 + 27. Profile 0 has rows
        // 0..=3 (the last is its sentinel row), profile 1 rows 4..=6.
        let bounds = files.profile_score_bounds(&base).unwrap();
        assert_eq!(bounds, vec![27 + 127 + 227 + 327, 427 + 527 + 627]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_missing_or_bad_freq_file_is_ncbis_error_and_a_short_one_is_refused() {
        let dir = scratch("freqerr");
        let base = write_volume(&dir, 28, AUX);
        let shown = canonical_volume_path(&base);
        let text = format!(
            "{shown} contains no frequency ratios needed for composition-based statistics.\nPlease disable composition-based statistics when searching against {shown}."
        );
        let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        assert_eq!(
            files.open_freq(&base).err(),
            Some(RpsOpenError::Ncbi(text.clone()))
        );
        std::fs::write(base.with_extension("freq"), words(&[7, 2, 0, 4, 7])).unwrap();
        assert_eq!(files.open_freq(&base).err(), Some(RpsOpenError::Ncbi(text)));
        write_freq(&base, 28);
        let freq = std::fs::read(base.with_extension("freq")).unwrap();
        std::fs::write(base.with_extension("freq"), &freq[..freq.len() - 4]).unwrap();
        assert!(matches!(
            files.open_freq(&base),
            Err(RpsOpenError::Refused(_))
        ));
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
        // A .rps cut inside its rows (NCBI crashes).
        let base = write_volume(&dir, 28, AUX);
        let rps = std::fs::read(base.with_extension("rps")).unwrap();
        std::fs::write(base.with_extension("rps"), &rps[..rps.len() - 4]).unwrap();
        assert!(matches!(
            RpsVolumeFiles::open(&base, &[3, 2]),
            Err(RpsOpenError::Refused(_))
        ));
        // Bytes after the last row are not read (NCBI searches the file).
        let mut longer = rps.clone();
        longer.extend_from_slice(&[0x7f; 112]);
        std::fs::write(base.with_extension("rps"), &longer).unwrap();
        let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        assert_eq!(files.rows().len(), 7);
        // (Windows does not write a mapped file.)
        drop(files);
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

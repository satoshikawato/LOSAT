//! Unit counts files (`CSeqMaskerIstat` and its readers): the format check of
//! `CSeqMaskerIstatFactory::DiscoverStatType`, the optimized binary reader
//! (`CSeqMaskerIstatOBinary` + `CSeqMaskerUsetHash`), the binary reader
//! (`CSeqMaskerIstatBin` + `CSeqMaskerUsetArray`), the optimized ascii reader
//! (`CSeqMaskerIstatOAscii` + `CSeqMaskerUsetHash`) and the ascii reader
//! (`CSeqMaskerIstatAscii` + `CSeqMaskerUsetSimple`).
//!
//! blastn reads the two ascii formats only through `-window_masker_taxid`: `-window_masker_db`
//! refuses them (`blast_args.cpp:446-458`, "Only optimized binary windowmasker stat files are
//! supported"), but the file of a taxid's directory is read with no check of its format.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use super::ostat::StatAlgoVersion;
use super::pct::{apply_pct, PctSettings};
use super::{MaskerParams, WinMaskError};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str;
use crate::utils::lower_bound::lower_bound_by;

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_istat_factory.hpp:56-62
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatType {
    Ascii,
    Binary,
    OAscii,
    OBinary,
    Unknown,
}

/// What `DiscoverStatType` finds: the format, the metadata lines (bytes, as NCBI's
/// `std::string`s), and the bytes before the counts (the binary metadata block).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StatFile {
    pub stat_type: StatType,
    pub metadata: Vec<Vec<u8>>,
    pub skip: u64,
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:165-197
pub fn discover_stat_type(path: &Path) -> Result<StatFile, WinMaskError> {
    let name = crate::utils::os_bytes::path_bytes(path).into_owned();
    let mut file = open_named(path).map_err(|_| WinMaskError::Open { name })?;
    discover(&mut file)
}

/// The counts file `path` opened by name.
///
/// No NCBI counterpart: `CNcbiIfstream` and POSIX `open` of an empty name fail (ENOENT) on every
/// platform NCBI runs on, which is what the empty name of a `-window_masker_taxid` run without a
/// counts file must give; WASI's `path_open` takes the empty path as the preopened directory
/// and succeeds, so the empty name is "not found" here on every target.
fn open_named(path: &Path) -> std::io::Result<File> {
    if path.as_os_str().is_empty() {
        return Err(std::io::ErrorKind::NotFound.into());
    }
    File::open(path)
}

/// `discover_stat_type` on an opened file or an in-memory copy of one (LOSAT Web's counts files,
/// `web_adapter`).
pub fn discover_stat_type_in<R: Read + Seek>(r: &mut R) -> Result<StatFile, WinMaskError> {
    discover(r)
}

/// Where NCBI's counts readers throw the format errors that LOSAT's readers report: the
/// message, and the place that the `T0` line of NCBI's exception report names
/// (`WinMaskError::create_report`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormatSite {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:57-63
    FileTooShort,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:113-115
    UnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:124-126
    HashKeySize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:128-130
    HashOffset,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:132-134
    HashShift,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:196-198
    HashTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:207-209
    ValuesTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:78-79
    BinFileSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:98-99
    BinUnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:107-108
    BinLength,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:100-102
    OAsciiUnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:113-116
    OAsciiHashKeySize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:118-121
    OAsciiHashOffset,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:123-127
    OAsciiHashShift,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:144-145
    OAsciiFileTooShort,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:175-177
    OAsciiHashTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:191-193
    OAsciiValuesTable,
}

impl FormatSite {
    /// NCBI's message.
    pub fn message(self) -> &'static str {
        match self {
            Self::FileTooShort => "file too short",
            Self::UnitSize => "unit size must be in [1,16]",
            Self::HashKeySize => "hash key size must be in [1,2*unit_size - 1]",
            Self::HashOffset => "offset must by in [0,32 - hash_key_size]",
            Self::HashShift => "shift must be in [1, 32 - hash_key_size]",
            Self::HashTable => "not enough data to fill the hash table",
            Self::ValuesTable => "not enough data to fill the values table",
            Self::BinFileSize => "wrong file size",
            Self::BinUnitSize => "illegal unit size",
            Self::BinLength => "wrong length",
            Self::OAsciiUnitSize => "unit size must be in [1,16]",
            Self::OAsciiHashKeySize => "hash key size must be in [1,2*unit_size - 1]",
            Self::OAsciiHashOffset => "offset must by in [0,32 - hash_key_size]",
            Self::OAsciiHashShift => "shift must be in [1, 32 - hash_key_size]",
            Self::OAsciiFileTooShort => "file too short",
            Self::OAsciiHashTable => "not enough lines to fill the hash table",
            Self::OAsciiValuesTable => "not enough lines to fill the values table",
        }
    }

    /// The file, line, class and error code text (`GetType()::GetErrCodeString()`) and function
    /// that NCBI's report prints for the throw. The function of `readWord` is not printed (NCBI 2.17.0's Linux
    /// build, observed).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:48-52
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:49-50
    /// The optimized ascii reader's codes are those of the optimized binary one
    /// (`seq_masker_istat_oascii.cpp:50-54`). The lines are those of `NCBI_THROW` (NCBI BLAST+
    /// 2.17.0, line NP M3d).
    pub(super) fn throw_site(self) -> ThrowSite {
        const OBINARY: &str = "src/algo/winmask/seq_masker_istat_obinary.cpp";
        const OBINARY_NEW: Option<&str> =
            Some("ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary");
        const BIN: &str = "src/algo/winmask/seq_masker_istat_bin.cpp";
        const BIN_NEW: Option<&str> = Some("ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin");
        const OASCII: &str = "src/algo/winmask/seq_masker_istat_oascii.cpp";
        const OASCII_NEW: Option<&str> = Some("ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii");
        let (file, line, kind, function) = match self {
            Self::FileTooShort => (OBINARY, 63, "Exception::format error", None),
            Self::UnitSize => (OBINARY, 114, "Exception::bad parameter", OBINARY_NEW),
            Self::HashKeySize => (OBINARY, 125, "Exception::bad hash parameter", OBINARY_NEW),
            Self::HashOffset => (OBINARY, 129, "Exception::bad hash parameter", OBINARY_NEW),
            Self::HashShift => (OBINARY, 133, "Exception::bad hash parameter", OBINARY_NEW),
            Self::HashTable => (OBINARY, 197, "Exception::format error", OBINARY_NEW),
            Self::ValuesTable => (OBINARY, 208, "Exception::format error", OBINARY_NEW),
            Self::BinFileSize => (BIN, 79, "Exception::file format error", BIN_NEW),
            Self::BinUnitSize => (BIN, 99, "Exception::file format error", BIN_NEW),
            Self::BinLength => (BIN, 108, "Exception::file format error", BIN_NEW),
            Self::OAsciiUnitSize => (OASCII, 101, "Exception::bad parameter", OASCII_NEW),
            Self::OAsciiHashKeySize => (OASCII, 114, "Exception::bad hash parameter", OASCII_NEW),
            Self::OAsciiHashOffset => (OASCII, 119, "Exception::bad hash parameter", OASCII_NEW),
            Self::OAsciiHashShift => (OASCII, 124, "Exception::bad hash parameter", OASCII_NEW),
            Self::OAsciiFileTooShort => (OASCII, 145, "Exception::format error", OASCII_NEW),
            Self::OAsciiHashTable => (OASCII, 176, "Exception::format error", OASCII_NEW),
            Self::OAsciiValuesTable => (OASCII, 192, "Exception::format error", OASCII_NEW),
        };
        ThrowSite {
            file,
            line,
            kind,
            function,
        }
    }
}

/// The body of `DiscoverStatType` on an opened file.
pub(super) fn discover<R: Read + Seek>(r: &mut R) -> Result<StatFile, WinMaskError> {
    let (metadata, skip) = metadata_lines(r)?;
    seek(r, 0)?;

    // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:174-195
    let mut reader = BufReader::new(r);
    let mut readable = true;
    if skip > 0 {
        readable = skip_bytes(&mut reader, skip);
    } else {
        let mut line = Vec::new();
        for _ in 0..metadata.len() {
            line.clear();
            if reader.read_until(b'\n', &mut line).unwrap_or(0) == 0 {
                readable = false;
            }
        }
    }
    let mut stat_type = StatType::Unknown;
    if readable {
        let mut word = [0u8; 4];
        if reader.read_exact(&mut word).is_ok() {
            stat_type = match u32::from_ne_bytes(word) {
                0 => StatType::Binary,
                0x4141_4141 => StatType::OAscii,
                1 | 2 => StatType::OBinary,
                _ => StatType::Ascii,
            };
        }
    }
    Ok(StatFile {
        stat_type,
        metadata,
        skip,
    })
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:61-107
/// A short first word keeps the zero bytes of `word( 0 )`; a short length word keeps the bytes
/// of the 3 it replaces. A file shorter than its metadata block fails the skip in
/// `DiscoverStatType`, which then returns eUnknown, so the block is read only as far as the
/// file goes (NCBI's bytes past the end are uninitialized and never used).
fn metadata_lines<R: Read + Seek>(r: &mut R) -> Result<(Vec<Vec<u8>>, u64), WinMaskError> {
    seek(r, 0)?;
    let mut word = [0u8; 4];
    read_up_to(r, &mut word);
    let word = u32::from_ne_bytes(word);
    if word == 3 {
        let mut len = 3u32.to_ne_bytes();
        let got = read_up_to(r, &mut len);
        let len = u32::from_ne_bytes(len);
        let mut data = Vec::new();
        if got == 4 {
            // A short read of the block leaves the stream failed; nothing follows it here.
            let _ = r.by_ref().take(len as u64).read_to_end(&mut data);
        }
        let metadata = data
            .split(|&b| b == 0)
            .take_while(|s| !s.is_empty())
            .map(<[u8]>::to_vec)
            .collect();
        return Ok((metadata, 8 + len as u64));
    }
    if word < 3 {
        return Ok((Vec::new(), 0));
    }
    seek(r, 0)?;
    let mut reader = BufReader::new(r);
    let mut metadata = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).unwrap_or(0) == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.len() >= 2 && line[0] == b'#' && line[1] == b'#' {
            metadata.push(line.clone());
        } else {
            break;
        }
    }
    Ok((metadata, 0))
}

/// The unit counts with the thresholds of the masker (`CSeqMaskerIstat`).
pub struct UnitCounts {
    unit_size: u8,
    threshold: u32,
    textend: u32,
    max_count: u32,
    use_max_count: u32,
    min_count: u32,
    use_min_count: u32,
    /// The cache bit array of the optimized binary format (`optimization_data::cba_`).
    cba: Option<Vec<u32>>,
    uset: Uset,
    metadata: Vec<Vec<u8>>,
    warnings: Vec<String>,
    /// The format of the file.
    stat_type: StatType,
    /// The counts algorithm's version that the metadata names (`GetStatAlgoVersion`).
    algo_version: StatAlgoVersion,
}

enum Uset {
    Hash(UsetHash),
    Array(UsetArray),
}

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_uset_hash.hpp:122-128
struct UsetHash {
    unit_size: u8,
    k: u8,
    roff: u8,
    bc: u8,
    m: u32,
    cmask: u32,
    ht: Vec<u32>,
    vt: Vec<u16>,
}

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_uset_array.hpp:124-130
struct UsetArray {
    unit_size: u8,
    entries: Vec<(u32, u32)>,
}

impl UnitCounts {
    /// `CSeqMaskerIstatFactory::create` for the two formats blastn accepts.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:317-340
    /// blastn passes no `-t_*_pct` setting (-1 each, `seq_masker.hpp:178-181`); the
    /// windowmasker application's are [`UnitCounts::open_with_pct`].
    pub fn open(path: &Path, params: &MaskerParams) -> Result<Self, WinMaskError> {
        Self::open_with_pct(path, params, &PctSettings::default(), &mut |_| {})
    }

    /// [`UnitCounts::open`] on an opened file or an in-memory copy of one.
    pub fn read<R: Read + Seek>(r: &mut R, params: &MaskerParams) -> Result<Self, WinMaskError> {
        Self::read_with_pct(r, params, &PctSettings::default(), &mut |_| {})
    }

    /// [`UnitCounts::open`] with the windowmasker application's `-t_*_pct` settings `pct`
    /// (`pct::apply_pct`); `warn` gets the warning that NCBI posts, at once, when the file has
    /// no percentages.
    pub fn open_with_pct(
        path: &Path,
        params: &MaskerParams,
        pct: &PctSettings,
        warn: &mut dyn FnMut(String),
    ) -> Result<Self, WinMaskError> {
        let name = crate::utils::os_bytes::path_bytes(path).into_owned();
        let mut file = open_named(path).map_err(|_| WinMaskError::Open { name })?;
        Self::read_with_pct(&mut file, params, pct, warn)
    }

    /// [`UnitCounts::open_with_pct`] on an opened file or an in-memory copy of one.
    pub fn read_with_pct<R: Read + Seek>(
        r: &mut R,
        params: &MaskerParams,
        pct: &PctSettings,
        warn: &mut dyn FnMut(String),
    ) -> Result<Self, WinMaskError> {
        let found = discover(r)?;
        check_pct_lines(&found.metadata, found.stat_type)?;
        let params = &apply_pct(params, &found.metadata, found.stat_type, pct, warn)?;
        let mut counts = Self {
            unit_size: 0,
            threshold: params.cutoff_score,
            textend: params.textend,
            max_count: params.max_score,
            use_max_count: params.set_max_score,
            min_count: params.min_score,
            use_min_count: params.set_min_score,
            cba: None,
            uset: Uset::Array(UsetArray {
                unit_size: 0,
                entries: Vec::new(),
            }),
            metadata: found.metadata,
            warnings: Vec::new(),
            stat_type: found.stat_type,
            algo_version: StatAlgoVersion::default(),
        };
        let start_line = counts.metadata.len();
        match found.stat_type {
            StatType::OBinary => counts.read_obinary(r, found.skip)?,
            StatType::Binary => counts.read_binary(r, found.skip)?,
            StatType::Ascii => counts.read_ascii(r, start_line)?,
            StatType::OAscii => counts.read_oascii(r, start_line)?,
            StatType::Unknown => {
                return Err(WinMaskError::CreateFail {
                    message: "unrecognized unit counts format".to_string(),
                })
            }
        }
        counts.algo_version =
            check_stat_algo_version(&counts.metadata).map_err(|e| counts.exception(e))?;
        Ok(counts)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_ascii.cpp:65-176
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_simple.cpp:55-67
    /// The units are kept in increasing order, so the search of `CSeqMaskerUsetSimple::get_info`
    /// (`lower_bound`) is that of the binary format's array ([`UsetArray::get_info`]). A file
    /// without a unit size line keeps the default of 15 (`seq_masker_uset_simple.hpp:79`). The
    /// ambiguity unit is used only by the merge pass, which blastn does not run. A line is a
    /// `getline` of the file opened in text mode, as on NCBI's Linux builds: a carriage return
    /// stays in it.
    fn read_ascii<R: Read + Seek>(
        &mut self,
        r: &mut R,
        start_line: usize,
    ) -> Result<(), WinMaskError> {
        seek(r, 0)?;
        let mut lines = TextLines::new(r);
        for _ in 0..start_line {
            lines.getline();
        }
        let mut linenum = start_line as u32;
        let mut unit_size: u8 = 15;
        let mut start = true;
        let mut entries: Vec<(u32, u32)> = Vec::new();
        while let Some(line) = lines.getline() {
            linenum = linenum.wrapping_add(1);
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if line[0] == b'>' {
                let name_end = line.iter().position(|&c| c == b' ' || c == b'\t');
                let val_start = name_end.and_then(|end| {
                    line[end..]
                        .iter()
                        .position(|&c| c != b' ' && c != b'\t')
                        .map(|at| end + at)
                });
                let (Some(name_end), Some(val_start)) = (name_end, val_start) else {
                    return Err(self.ascii_error(ascii_site(96), format!("at line {linenum}")));
                };
                let name = &line[1..name_end];
                let value = &line[val_start..];
                let parse = |value: &[u8]| ncbi_str::string_to_uint_base(value, 0);
                match name {
                    b"t_threshold" if self.threshold == 0 => {
                        self.threshold = parse(value).map_err(|e| self.exception(e))?;
                    }
                    b"t_extend" if self.textend == 0 => {
                        self.textend = parse(value).map_err(|e| self.exception(e))?;
                    }
                    b"t_low" => {
                        let count = parse(value).map_err(|e| self.exception(e))?;
                        self.set_min_count(count);
                    }
                    b"t_high" if self.max_count == 0 => {
                        self.max_count = parse(value).map_err(|e| self.exception(e))?;
                    }
                    _ => {}
                }
                continue;
            }
            if start {
                start = false;
                unit_size = ncbi_str::string_to_uint(line).map_err(|e| self.exception(e))? as u8;
                continue;
            }
            let blank = |c: &u8| *c == b' ' || *c == b'\t';
            let unit_start = line.iter().position(|c| !blank(c));
            let unit_end =
                unit_start.and_then(|from| line[from..].iter().position(blank).map(|at| from + at));
            let cnt_start = unit_end.and_then(|from| {
                line[from..]
                    .iter()
                    .position(|c| !blank(c))
                    .map(|at| from + at)
            });
            let (Some(unit_start), Some(unit_end), Some(cnt_start)) =
                (unit_start, unit_end, cnt_start)
            else {
                return Err(self.ascii_error(ascii_site(137), format!("at line {linenum}")));
            };
            let unit = ncbi_str::string_to_uint_base(&line[unit_start..unit_end], 16)
                .map_err(|e| self.exception(e))?;
            let count =
                ncbi_str::string_to_uint(&line[cnt_start..]).map_err(|e| self.exception(e))?;
            if count >= self.min_count {
                if let Some(&(last, _)) = entries.last() {
                    if unit <= last {
                        return Err(self.exception(NcbiException::new(
                            USET_SIMPLE_BAD_ORDER,
                            format!("last unit: {last:x} ; adding {unit:x}"),
                        )));
                    }
                }
                entries.push((unit, count));
            }
        }
        let mut bad_param = String::new();
        for (missing, name) in [
            (self.threshold == 0, "t_threhold "),
            (self.textend == 0, "t_extend "),
            (self.max_count == 0, "t_high "),
            (self.min_count == 0, "t_low "),
        ] {
            if missing {
                bad_param.push_str(name);
            }
        }
        if !bad_param.is_empty() {
            return Err(self.ascii_error(ASCII_BAD_PARAM, bad_param));
        }
        if self.use_min_count == 0 {
            self.use_min_count = self.min_count.wrapping_add(1) / 2;
        }
        if self.use_max_count == 0 {
            self.use_max_count = self.max_count;
        }
        entries.shrink_to_fit();
        self.unit_size = unit_size;
        self.uset = Uset::Array(UsetArray { unit_size, entries });
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:72-196
    /// The numbers are glibc's `atoi` ([`ncbi_str::c_atoi`]) truncated to the field, and those of
    /// the hash table parameters are read as C++'s `>>` reads an `unsigned int`
    /// ([`extract_uints`]). A value that `>>` does not reach (the line ends, or an earlier value
    /// failed) is left uninitialized by NCBI: a check that reads one is refused
    /// ([`WinMaskError::UnsetHashParameter`]). NCBI 2.17.0's errors for such lines change from
    /// one load of the masker to the next (line NP M3d, review 2: `3694` gives the hash key
    /// size's error for the first query, then the offset's). A value that `>>` failed to read is
    /// set (0, or `UINT_MAX` past it), so `3694 x` and `3694 17 0 x` are NCBI's errors. The tables
    /// grow with the lines read, so a file that names a larger table than it holds fails without
    /// NCBI's allocation of the named size. The optimized ascii format has no cache bit array.
    fn read_oascii<R: Read + Seek>(
        &mut self,
        r: &mut R,
        start_line: usize,
    ) -> Result<(), WinMaskError> {
        const HEADER_LINES: u32 = 7;
        let file_len = r
            .seek(SeekFrom::End(0))
            .map_err(|_| WinMaskError::format(FormatSite::FileTooShort))?;
        seek(r, 0)?;
        let mut lines = TextLines::new(r);
        for _ in 0..start_line {
            lines.getline();
        }
        let mut linenum = 0u32;
        let mut unit_size: u8 = 0;
        let [mut m, mut k, mut roff, mut bc] = [0u32; 4];
        let mut thresholds = [0u32; 4];
        while let Some(line) = lines.getline() {
            linenum += 1;
            match linenum {
                2 => {
                    unit_size = ncbi_str::c_atoi(line) as u8;
                    if unit_size == 0 || unit_size > 16 {
                        return Err(self.format_error(FormatSite::OAsciiUnitSize));
                    }
                }
                3 => {
                    // The checks read k, roff and bc in this order (M is set whenever k is).
                    let set;
                    ([m, k, roff, bc], set) = extract_uints(line);
                    let unset = |field: &'static str| WinMaskError::UnsetHashParameter { field };
                    if set < 2 {
                        return Err(unset("k"));
                    }
                    if k == 0 || k > 2 * u32::from(unit_size) - 1 {
                        return Err(self.format_error(FormatSite::OAsciiHashKeySize));
                    }
                    if set < 3 {
                        return Err(unset("roff"));
                    }
                    if roff > 32 - k {
                        return Err(self.format_error(FormatSite::OAsciiHashOffset));
                    }
                    if set < 4 {
                        return Err(unset("bc"));
                    }
                    if bc == 0 || bc > 32 - k {
                        return Err(self.format_error(FormatSite::OAsciiHashShift));
                    }
                }
                4..=7 => thresholds[linenum as usize - 4] = ncbi_str::c_atoi(line) as u32,
                _ => {}
            }
            if linenum == HEADER_LINES {
                break;
            }
        }
        if linenum < HEADER_LINES {
            return Err(self.format_error(FormatSite::OAsciiFileTooShort));
        }
        let [t_low, t_extend, t_thres, t_high] = thresholds;
        self.unit_size = unit_size;
        self.set_thresholds(t_low, t_extend, t_thres, t_high);

        // Every line read has at least one byte: no more lines than the file has bytes.
        let room = usize::try_from(file_len).unwrap_or(usize::MAX);
        let ht_size = 1usize << k;
        let mut ht = Vec::with_capacity(ht_size.min(room));
        while ht.len() < ht_size {
            let Some(line) = lines.getline() else {
                return Err(self.format_error(FormatSite::OAsciiHashTable));
            };
            ht.push(ncbi_str::c_atoi(line) as u32);
        }
        let mut vt = Vec::with_capacity((m as usize).min(room));
        while vt.len() < m as usize {
            let Some(line) = lines.getline() else {
                return Err(self.format_error(FormatSite::OAsciiValuesTable));
            };
            vt.push(ncbi_str::c_atoi(line) as u16);
        }
        self.uset = Uset::Hash(UsetHash {
            unit_size,
            k: k as u8,
            roff: roff as u8,
            bc: bc as u8,
            m,
            cmask: ((1u64 << bc) - 1) as u32,
            ht,
            vt,
        });
        Ok(())
    }

    /// A format error of a reader, with the warnings posted before it.
    fn format_error(&self, site: FormatSite) -> WinMaskError {
        WinMaskError::Format {
            site,
            warnings: self.warnings.clone(),
        }
    }

    /// An exception thrown inside the reader (the factory rethrows it), with the warnings
    /// posted before it.
    fn exception(&self, exception: NcbiException) -> WinMaskError {
        WinMaskError::Exception {
            exception: Box::new(exception),
            warnings: self.warnings.clone(),
        }
    }

    /// The ascii reader's own exception at `site`.
    fn ascii_error(&self, site: ThrowSite, message: String) -> WinMaskError {
        self.exception(NcbiException::new(site, message))
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:70-216
    /// The words are read in the machine's byte order, as NCBI reads them, through [`Stream`]:
    /// `readWord` throws "file too short" when an earlier read came up short, and the word of a
    /// short read keeps 0 in the bytes it did not get (NCBI's are uninitialized). So a file that
    /// ends inside the last header word (`t_high`, or the flag word of format 2) gives NCBI's
    /// "not enough data to fill the hash table". A short cache bit array posts NCBI's warning
    /// and then fails the hash table read, as NCBI does; the warning goes with the error.
    fn read_obinary<R: Read + Seek>(&mut self, r: &mut R, skip: u64) -> Result<(), WinMaskError> {
        let mut s = Stream::new(r)?;
        s.skip(skip);
        let read_word = |s: &mut Stream<'_, R>| -> Result<u32, WinMaskError> {
            if !s.good() {
                return Err(WinMaskError::format(FormatSite::FileTooShort));
            }
            let mut word = 0u32;
            s.read_into(&mut word);
            Ok(word)
        };
        let use_opt = read_word(&mut s)? != 1;
        let unit_size = read_word(&mut s)? as u8;
        if unit_size == 0 || unit_size > 16 {
            return Err(WinMaskError::format(FormatSite::UnitSize));
        }
        let m = read_word(&mut s)?;
        let k = read_word(&mut s)?;
        let roff = read_word(&mut s)?;
        let bc = read_word(&mut s)?;
        if k == 0 || k > (2 * unit_size as u32 - 1) {
            return Err(WinMaskError::format(FormatSite::HashKeySize));
        }
        if roff > 32 - k {
            return Err(WinMaskError::format(FormatSite::HashOffset));
        }
        if bc == 0 || bc > 32 - k {
            return Err(WinMaskError::format(FormatSite::HashShift));
        }
        let t_low = read_word(&mut s)?;
        let t_extend = read_word(&mut s)?;
        let t_thres = read_word(&mut s)?;
        let t_high = read_word(&mut s)?;
        self.unit_size = unit_size;
        self.set_thresholds(t_low, t_extend, t_thres, t_high);

        let mut warnings = Vec::new();
        if use_opt {
            let divisor = read_word(&mut s)?;
            if divisor > 0 {
                // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:164-186
                // `arg_use_ba` is true for blastn and the application (`winmask_filter.cpp:121`,
                // `win_mask_config.cpp:361`), so the array is kept. Units of 1 or 2 letters give
                // an array of 0 words that is still in use (a `new` of 0 `Uint4` words); see `SeqMasker::mask`.
                let cba_size = (1u64 << (2 * unit_size as u32)) / 32;
                match s.read_u32s(cba_size) {
                    Some(cba) => self.cba = Some(cba),
                    None => warnings.push(
                        "file read failed: bit array optimizations are not used.".to_string(),
                    ),
                }
            }
        }

        let ht = s.read_u32s(1u64 << k).ok_or_else(|| WinMaskError::Format {
            site: FormatSite::HashTable,
            warnings: warnings.clone(),
        })?;
        let vt = s.read_u16s(m as u64).ok_or_else(|| WinMaskError::Format {
            site: FormatSite::ValuesTable,
            warnings: warnings.clone(),
        })?;
        self.warnings.extend(warnings);

        // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:53-62
        self.uset = Uset::Hash(UsetHash {
            unit_size,
            k: k as u8,
            roff: roff as u8,
            bc: bc as u8,
            m,
            cmask: ((1u64 << bc) - 1) as u32,
            ht,
            vt,
        });
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:76-122
    /// One `data` word takes every read, so a short read keeps the earlier value in the bytes it
    /// did not get. `datalen % 8` is unsigned (`sizeof`), so a data part of -8 or -16 bytes
    /// passes and the trailer comes from the bytes after the header. DiscoverStatType has read
    /// the first word, so `skip` and the format word are in the file; the data read is never
    /// short (the file holds the data and the trailer behind the header).
    fn read_binary<R: Read + Seek>(&mut self, r: &mut R, skip: u64) -> Result<(), WinMaskError> {
        const HEADER_LEN: i64 = 8;
        const TRAILER_LEN: i64 = 16;
        let mut s = Stream::new(r)?;
        let iflen = s.left as i64;
        if iflen < HEADER_LEN + TRAILER_LEN {
            return Err(WinMaskError::format(FormatSite::BinFileSize));
        }
        s.skip(skip);
        let adj = if s.good() { skip as i64 } else { -1 };
        let mut data = 0u32;
        s.read_into(&mut data);
        s.read_into(&mut data);
        let unit_size = data as u8;
        if unit_size == 0 || unit_size > 16 {
            return Err(WinMaskError::format(FormatSite::BinUnitSize));
        }
        let datalen = iflen - HEADER_LEN - TRAILER_LEN - adj;
        if !(datalen as u64).is_multiple_of(8) {
            return Err(WinMaskError::format(FormatSite::BinLength));
        }
        let mut entries = Vec::new();
        if datalen > 0 {
            entries = s
                .read_u32_pairs(datalen as u64 / 8)
                .ok_or_else(|| WinMaskError::format(FormatSite::BinLength))?;
        }
        s.read_into(&mut data);
        let t_low = data;
        s.read_into(&mut data);
        let t_extend = data;
        s.read_into(&mut data);
        let t_thres = data;
        s.read_into(&mut data);
        let t_high = data;
        self.unit_size = unit_size;
        self.set_thresholds(t_low, t_extend, t_thres, t_high);
        self.uset = Uset::Array(UsetArray { unit_size, entries });
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:141-157
    /// The binary reader does the same (`seq_masker_istat_bin.cpp:120-142`).
    fn set_thresholds(&mut self, t_low: u32, t_extend: u32, t_thres: u32, t_high: u32) {
        self.set_min_count(t_low);
        if self.textend == 0 {
            self.textend = t_extend;
        }
        if self.threshold == 0 {
            self.threshold = t_thres;
        }
        if self.max_count == 0 {
            self.max_count = t_high;
        }
        if self.use_min_count == 0 {
            self.use_min_count = self.min_count.wrapping_add(1) / 2;
        }
        if self.use_max_count == 0 {
            self.use_max_count = self.max_count;
        }
    }

    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_istat.hpp:260-273
    /// The warning is kept for the caller to print ([`UnitCounts::warnings`]).
    fn set_min_count(&mut self, arg_min_count: u32) {
        if self.min_count != 0 && self.min_count < arg_min_count {
            self.warnings.push(format!(
                "Requested value of t_low ({}) is less than the one stored with the N-mer counts ({}).The value {} will be used.",
                self.min_count, arg_min_count, arg_min_count
            ));
            self.min_count = arg_min_count;
        } else if self.min_count == 0 {
            self.min_count = arg_min_count;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:219-227
    /// `CSeqMaskerIstatBin::at` is the same (`seq_masker_istat_bin.cpp:151-159`).
    #[inline]
    pub(super) fn at(&self, unit: u32) -> Result<u32, WinMaskError> {
        let res = match &self.uset {
            Uset::Hash(h) => h.get_info(unit)?,
            Uset::Array(a) => a.get_info(unit),
        };
        Ok(if res == 0 || res < self.min_count {
            self.use_min_count
        } else if res > self.max_count {
            self.use_max_count
        } else {
            res
        })
    }

    /// The count stored for a unit, without the T_low and T_high replacement (`trueat`).
    pub fn true_count(&self, unit: u32) -> Result<u32, WinMaskError> {
        match &self.uset {
            Uset::Hash(h) => h.get_info(unit),
            Uset::Array(a) => Ok(a.get_info(unit)),
        }
    }

    pub fn unit_size(&self) -> u8 {
        self.unit_size
    }

    /// T_threshold (`get_threshold`).
    pub fn threshold(&self) -> u32 {
        self.threshold
    }

    /// T_extend (`get_textend`).
    pub fn textend(&self) -> u32 {
        self.textend
    }

    /// T_high (`get_max_count`).
    pub fn max_count(&self) -> u32 {
        self.max_count
    }

    /// T_low (`get_min_count`).
    pub fn min_count(&self) -> u32 {
        self.min_count
    }

    /// The count used above T_high (`get_use_max_count`).
    pub fn use_max_count(&self) -> u32 {
        self.use_max_count
    }

    /// The count used below T_low (`get_use_min_count`).
    pub fn use_min_count(&self) -> u32 {
        self.use_min_count
    }

    /// The metadata lines of the file (`##...`), as bytes.
    pub fn metadata(&self) -> &[Vec<u8>] {
        &self.metadata
    }

    /// The warnings NCBI posts while reading the file (`ERR_POST` with `Warning`).
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// The counts algorithm's version that the file names (`GetStatAlgoVersion`).
    pub fn stat_algo_version(&self) -> StatAlgoVersion {
        self.algo_version
    }

    /// The text of the file's first `##note:` line (`GetMetaData`, `ExtractMetaDataStr`): what
    /// follows the first `:` after the line's `##` of the first line whose name is `note`; empty
    /// for none.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:149-162
    pub fn note(&self) -> &[u8] {
        for line in &self.metadata {
            let Some(at) = line.iter().skip(2).position(|&b| b == b':') else {
                continue;
            };
            let pos = at + 2;
            if &line[2..pos] == b"note" {
                return &line[pos + 1..];
            }
        }
        &[]
    }

    /// The percentage of each count from 0 that the file's `##pct:` lines hold (`GetCountMap`),
    /// read as the factory reads them; empty for none.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:237-276
    pub fn count_map(&self) -> Result<Vec<f64>, WinMaskError> {
        super::pct::count_map(&self.metadata, self.stat_type)
    }

    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_istat.hpp:155-156
    pub(super) fn cache_bit_array(&self) -> Option<&[u32]> {
        self.cba.as_deref()
    }
}

impl UsetHash {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:73-111
    // Always inlined: the masker reaches it through `at` for every unit of every window, and
    // with `true_count` (`-convert`) as a second caller LLVM kept it out of line, which cost
    // blastn `-window_masker_db` 0.2-0.5 % more instructions.
    #[inline(always)]
    fn get_info(&self, unit: u32) -> Result<u32, WinMaskError> {
        let runit = reverse_complement(unit, self.unit_size);
        let unit = if runit < unit { runit } else { unit };
        let (key, check) = hash_code(unit, self.k, self.roff);
        let hval = self.ht[key as usize];
        let coll = hval & self.cmask;
        if coll == 0 {
            Ok(0)
        } else if coll == 1 {
            if check as u32 != hval >> 24 {
                Ok(0)
            } else {
                Ok((hval >> self.bc) & 0xFFF)
            }
        } else {
            // In 64 bits: a sum that wraps `Uint4` would read outside NCBI's table.
            let first = (hval >> self.bc) as u64;
            if first + coll as u64 > self.m as u64 {
                return Err(WinMaskError::BadIndex { key, value: hval });
            }
            let values = &self.vt[first as usize..(first + coll as u64) as usize];
            Ok(values
                .iter()
                .find(|&&v| (v >> 9) as u32 == check as u32)
                .map_or(0, |&v| (v & 0x1FF) as u32))
        }
    }
}

impl UsetArray {
    /// `CSeqMaskerUsetArray::get_info` of a unit: the count of the unit or of its reverse
    /// complement, whichever is smaller, 0 when the array has no entry for it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_array.cpp:64-82
    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_uset_array.hpp:151-153
    /// The units of a binary file are not checked for order (`add_info`), so in a file whose
    /// units are not increasing the entry found is the one NCBI's `lower_bound` finds
    /// (`utils::lower_bound`). No data (`unit_data` null) is an empty array here: 0 either way.
    // Always inlined into `at`, as `UsetHash::get_info` is.
    #[inline(always)]
    fn get_info(&self, unit: u32) -> u32 {
        let runit = reverse_complement(unit, self.unit_size);
        let unit = if runit < unit { runit } else { unit };
        let entries = &self.entries;
        let at = lower_bound_by(0, entries.len(), |i| entries[i].0 < unit);
        match entries.get(at) {
            Some(&(u, count)) if u == unit => count,
            _ => 0,
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_util.cpp:54-65
#[inline]
pub(super) fn reverse_complement(seq: u32, size: u8) -> u32 {
    let mut result = 0u32;
    for i in 0..size as u32 {
        let letter = !(((seq >> (2 * i)) & 0x3) | !0x3);
        result = (result << 2) | letter;
    }
    result
}

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_util.hpp:81-87
/// With units of 16 letters `roff + k` can be 32, a shift that C++ leaves undefined. NCBI's
/// x86-64 and arm64 builds shift by the amount modulo 32, and the counts files they write
/// (`CSeqMaskerOstatOpt::doFinalize`, `seq_masker_ostat_opt.cpp:289-313`) hold the check bits
/// computed that way, so the reader computes them the same way to find the stored counts.
#[inline]
pub(super) fn hash_code(unit: u32, k: u8, roff: u8) -> (u32, u8) {
    let key = (unit >> roff) & ((1u32 << k).wrapping_sub(1));
    let low = unit & ((1u32 << roff).wrapping_sub(1));
    let check = (unit.wrapping_shr(roff as u32 + k as u32) << roff).wrapping_add(low) as u8;
    (key, check)
}

fn seek<R: Seek>(r: &mut R, pos: u64) -> Result<(), WinMaskError> {
    r.seek(SeekFrom::Start(pos))
        .map(|_| ())
        .map_err(|_| WinMaskError::format(FormatSite::FileTooShort))
}

/// Reads as many bytes as the stream has, up to the length of `buf` (`istream::read` that
/// stops at the end of the file).
fn read_up_to<R: Read>(r: &mut R, buf: &mut [u8]) -> usize {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(got) => n += got,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(_) => break,
        }
    }
    n
}

fn skip_bytes<R: Read>(r: &mut R, n: u64) -> bool {
    std::io::copy(&mut r.take(n), &mut std::io::sink()).is_ok_and(|got| got == n)
}

/// The reads of NCBI's readers on one `CNcbiIfstream`: a read that comes up short sets
/// `failbit`, and every later read then reads nothing.
struct Stream<'a, R: Read> {
    r: BufReader<&'a mut R>,
    failed: bool,
    /// The bytes from the read position to the end of the file.
    left: u64,
}

impl<'a, R: Read + Seek> Stream<'a, R> {
    fn new(r: &'a mut R) -> Result<Self, WinMaskError> {
        let len = r
            .seek(SeekFrom::End(0))
            .map_err(|_| WinMaskError::format(FormatSite::FileTooShort))?;
        seek(r, 0)?;
        Ok(Self {
            r: BufReader::with_capacity(1 << 16, r),
            failed: false,
            left: len,
        })
    }
}

impl<R: Read> Stream<'_, R> {
    /// `operator bool()` of the stream.
    fn good(&self) -> bool {
        !self.failed
    }

    /// `is.read( buf, n )`: the bytes read replace the first bytes of `buf`.
    fn read(&mut self, buf: &mut [u8]) -> usize {
        if self.failed {
            return 0;
        }
        let got = read_up_to(&mut self.r, buf);
        self.left -= got as u64;
        if got < buf.len() {
            self.failed = true;
        }
        got
    }

    fn skip(&mut self, n: u64) {
        if self.failed {
            return;
        }
        if skip_bytes(&mut self.r, n) {
            self.left -= n;
        } else {
            self.failed = true;
            self.left = 0;
        }
    }

    /// `is.read` of one `Uint4` into `word`, in the machine's byte order.
    fn read_into(&mut self, word: &mut u32) {
        let mut bytes = word.to_ne_bytes();
        self.read(&mut bytes);
        *word = u32::from_ne_bytes(bytes);
    }

    /// `is.read` of `n` `Uint4` words into a new array; `None` when the read
    /// fails (checked against the file's length before the array is allocated).
    fn read_u32s(&mut self, n: u64) -> Option<Vec<u32>> {
        self.read_array(n, u32::from_ne_bytes)
    }

    /// `is.read` of `n` `Uint2` words into a new array; `None` when it fails.
    fn read_u16s(&mut self, n: u64) -> Option<Vec<u16>> {
        self.read_array(n, u16::from_ne_bytes)
    }

    /// `n` pairs of words (unit, count), each pair read into its entry at once.
    fn read_u32_pairs(&mut self, n: u64) -> Option<Vec<(u32, u32)>> {
        self.read_array(n, |pair: [u8; 8]| {
            let [a, b, c, d, e, f, g, h] = pair;
            (
                u32::from_ne_bytes([a, b, c, d]),
                u32::from_ne_bytes([e, f, g, h]),
            )
        })
    }

    /// `n` values of `W` bytes each, the vector of exactly `n` filled through a block
    /// buffer; `None` (the stream failed) when the file has fewer bytes.
    fn read_array<T, const W: usize>(
        &mut self,
        n: u64,
        value: impl Fn([u8; W]) -> T,
    ) -> Option<Vec<T>> {
        let bytes = n.checked_mul(W as u64)?;
        if self.failed || bytes > self.left {
            self.failed = true;
            self.left = 0;
            return None;
        }
        let mut out = Vec::with_capacity(n as usize);
        let mut buf = vec![0u8; 1 << 16];
        let mut todo = bytes as usize;
        while todo > 0 {
            let take = todo.min(buf.len() / W * W);
            if self.read(&mut buf[..take]) < take {
                return None;
            }
            out.extend(
                buf[..take]
                    .chunks_exact(W)
                    .map(|w| value(w.try_into().expect("a chunk of W bytes"))),
            );
            todo -= take;
        }
        Some(out)
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:239-276
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:361-366
/// The binary formats read the first `##pct:` line only (the loop breaks after it), the ascii
/// ones every such line. A line of exactly those six letters throws `std::out_of_range` from
/// `substr` at position 7; the text is libstdc++'s, as in NCBI's Linux builds. A binary line whose
/// largest count reads as `u32::MAX` (also a larger number, which the stream stores as that) is
/// refused whether or not the percentages are used: NCBI's table has room for the count plus
/// one, 0 in 32 bits, and its loop reads past it (undefined behaviour in C++). The percentages
/// themselves are used only by the application's `-t_*_pct` settings (`pct::apply_pct`) and
/// `-convert` (`pct::count_map`).
fn check_pct_lines(metadata: &[Vec<u8>], stat_type: StatType) -> Result<(), WinMaskError> {
    let mut pct = metadata.iter().filter(|line| line.starts_with(b"##pct:"));
    let (short, first) = match stat_type {
        StatType::Binary | StatType::OBinary => {
            let first = pct.next();
            (first.is_some_and(|line| line.len() == 6), first)
        }
        StatType::Ascii | StatType::OAscii => (pct.any(|line| line.len() == 6), None),
        StatType::Unknown => (false, None),
    };
    if short {
        return Err(WinMaskError::CreateFail {
            message: "could not create a unit counts container: basic_string::substr: __pos (which is 7) > this->size() (which is 6)".to_string(),
        });
    }
    if let Some(line) = first {
        let mut at = 0;
        if extract_uint(&line[7..], &mut at).is_some_and(|(max, _)| max == u32::MAX) {
            return Err(WinMaskError::PctLine {
                line: String::from_utf8_lossy(line).into_owned(),
                why: super::pct::PctWhy::MaxCountWraps,
            });
        }
    }
    Ok(())
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:111-146
/// `CSeqMaskerIstatFactory::create` calls it after the reader (`:343-348`), inside its `try`:
/// a part of digits that `StringToInt` cannot convert (empty, or beyond `int`) is its
/// exception, which the factory rethrows. The masker does not use the version; `-convert`
/// writes it into the converted file (1.0.0 when no line names one).
fn check_stat_algo_version(metadata: &[Vec<u8>]) -> Result<StatAlgoVersion, NcbiException> {
    const NAME: &[u8] = b"windowmasker-statistics-algorithm-version";
    let digits = |t: &[u8]| t.iter().all(u8::is_ascii_digit);
    let find = |line: &[u8], c: u8, from: usize| {
        line.get(from..)
            .and_then(|rest| rest.iter().position(|&b| b == c))
            .map(|at| from + at)
    };
    for line in metadata {
        let line = line.as_slice();
        let Some(p) = find(line, b':', 2) else {
            continue;
        };
        if &line[2..p] != NAME {
            continue;
        }
        let pp = p + 1;
        let Some(p) = find(line, b'.', pp) else {
            continue;
        };
        let t = &line[pp..p];
        if !digits(t) {
            continue;
        }
        let major = ncbi_str::string_to_int(t, 0)?;
        let pp = p + 1;
        let Some(p) = find(line, b'.', pp) else {
            continue;
        };
        let t = &line[pp..p];
        if !digits(t) {
            continue;
        }
        let minor = ncbi_str::string_to_int(t, 0)?;
        let t = &line[p + 1..];
        if !digits(t) {
            continue;
        }
        let patch = ncbi_str::string_to_int(t, 0)?;
        return Ok(StatAlgoVersion {
            major,
            minor,
            patch,
        });
    }
    Ok(StatAlgoVersion::default())
}

/// The function of the ascii reader's throws.
const ASCII_NEW: Option<&str> = Some("ncbi::CSeqMaskerIstatAscii::CSeqMaskerIstatAscii");

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_ascii.cpp:40-50
/// The ascii reader's syntax errors: line 96 (a `>` line) and 137 (a unit line).
const fn ascii_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/algo/winmask/seq_masker_istat_ascii.cpp",
        line,
        kind: "Exception::syntax error",
        function: ASCII_NEW,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_ascii.cpp:168-169
const ASCII_BAD_PARAM: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker_istat_ascii.cpp",
    line: 169,
    kind: "Exception::bad parameter value",
    function: ASCII_NEW,
};

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_simple.cpp:44-52
const USET_SIMPLE_BAD_ORDER: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker_uset_simple.cpp",
    line: 62,
    kind: "Exception::bad unit order",
    function: Some("ncbi::CSeqMaskerUsetSimple::add_info"),
};

/// The lines of a text counts file as NCBI's readers take them with `getline` on an
/// `ifstream` opened in text mode (Linux: no translation of line ends): each line without its
/// line feed, the last one without one too; `None` once the stream has failed (the end of the
/// file was reached with nothing read, or a read error), and after that.
struct TextLines<R: Read> {
    r: BufReader<R>,
    line: Vec<u8>,
    failed: bool,
}

impl<R: Read> TextLines<R> {
    fn new(r: R) -> Self {
        Self {
            r: BufReader::with_capacity(1 << 16, r),
            line: Vec::new(),
            failed: false,
        }
    }

    fn getline(&mut self) -> Option<&[u8]> {
        if self.failed {
            return None;
        }
        self.line.clear();
        match self.r.read_until(b'\n', &mut self.line) {
            Ok(0) | Err(_) => {
                self.failed = true;
                None
            }
            Ok(_) => {
                if self.line.last() == Some(&b'\n') {
                    self.line.pop();
                }
                Some(&self.line)
            }
        }
    }
}

/// Reads up to `N` unsigned integers from `line` as C++'s `std::istringstream` extracts them with
/// `>>` into `unsigned int` (default flags: decimal, white space skipped), and returns them with
/// the number of values that the stream stored.
///
/// A value that was not stored (the line ended before it, or an earlier value failed) is 0 in the
/// returned array and is not counted: NCBI leaves such a variable as it was. A value that failed is
/// stored and counted (0 when no number was found, `u32::MAX` when it was too large).
///
/// The rules are those of the C++ standard's `num_get` stages 1 to 3 for an unsigned integer
/// (specifier `%u`) and `strtoull`, as the table `tests/fixtures/winmask_uints/oascii_hash_line.tsv`
/// observed them (NCBI 2.17.0's C++ library, 3,038 lines, 4 fields each):
///
/// * Each `>>` first skips white space (space, tab, line feed, vertical tab, form feed, carriage
///   return); if the text ends there nothing is stored and the stream has failed.
/// * Stage 2 accumulates an optional `+` or `-` and then decimal digits, and stops at any other
///   character (`0x10` is `0`; `1.5` is `1`; a sign followed by a space has no digits).
/// * Stage 3 converts the field. No digit: 0, and the stream has failed. Otherwise the digits give
///   a magnitude; one above `u32::MAX` stores `u32::MAX` and the stream has failed, whatever the
///   sign (`-4294967296`, and `-18446744073709551615`, where `strtoull` alone would wrap in 64
///   bits, fail alike). A magnitude that fits is stored as it is, or negated in 32 bits when the
///   field has a `-` (`-1` is `u32::MAX`, `-4294967295` is 1, `-0` is 0), and the stream has not
///   failed.
/// * A stream that has failed, or that reached the end of the text while accumulating a field,
///   stores nothing more.
pub(super) fn extract_uints<const N: usize>(line: &[u8]) -> ([u32; N], usize) {
    let mut values = [0u32; N];
    let mut stored = 0;
    let mut at = 0;
    while stored < N {
        let Some((value, failed)) = extract_uint(line, &mut at) else {
            break;
        };
        values[stored] = value;
        stored += 1;
        if failed || at == line.len() {
            break;
        }
    }
    (values, stored)
}

/// The white space that `>>` skips (the "C" locale's `isspace`).
pub(super) fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// One `>>` of an `unsigned int` from `line` at `*at` (by the rules of [`extract_uints`]), which
/// moves past what it read: `None` when the text ends at white space (nothing stored, the
/// stream has failed); else the value stored and whether the stream has failed. A stream that
/// has failed, or that has reached the end of the text, stores nothing more.
pub(super) fn extract_uint(line: &[u8], at: &mut usize) -> Option<(u32, bool)> {
    while *at < line.len() && is_c_space(line[*at]) {
        *at += 1;
    }
    if *at == line.len() {
        return None;
    }
    let negative = line[*at] == b'-';
    if matches!(line[*at], b'-' | b'+') {
        *at += 1;
    }
    let first_digit = *at;
    let mut magnitude = 0u64;
    let mut too_large = false;
    while *at < line.len() && line[*at].is_ascii_digit() {
        magnitude = magnitude * 10 + u64::from(line[*at] - b'0');
        if magnitude > u64::from(u32::MAX) {
            // Past this the magnitude only has to stay above the limit.
            too_large = true;
            magnitude = u64::from(u32::MAX) + 1;
        }
        *at += 1;
    }
    Some(if *at == first_digit {
        (0, true)
    } else if too_large {
        (u32::MAX, true)
    } else if negative {
        ((magnitude as u32).wrapping_neg(), false)
    } else {
        (magnitude as u32, false)
    })
}

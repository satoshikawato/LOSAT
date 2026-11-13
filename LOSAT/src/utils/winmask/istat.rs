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

use super::{MaskerParams, WinMaskError};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str;

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_istat_factory.hpp:56-62
/// ```c
///     enum EStatType {
///         eAscii,
///         eBinary,
///         eOAscii,
///         eOBinary,
///         eUnknown
///     };
/// ```
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
/// ```c
/// CSeqMaskerIstatFactory::EStatType CSeqMaskerIstatFactory::DiscoverStatType(
///         string const & name, vector< string > & md, size_t & skip ) {
///     GetMetaDataLines( name, md, skip );
///     CNcbiIfstream check( name.c_str() );
///
///     if( !check ) {
///         NCBI_THROW( Exception, eOpen, string( "could not open " ) + name );
///     }
/// ...
/// ```
pub fn discover_stat_type(path: &Path) -> Result<StatFile, WinMaskError> {
    let name = path.to_string_lossy().into_owned();
    let mut file = File::open(path).map_err(|_| WinMaskError::Open { name })?;
    discover(&mut file)
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
    /// ```c
    /// Uint4 CSeqMaskerIstatOBinary::readWord( CNcbiIstream & is ) const
    /// {
    ///     Uint4 result;
    ///
    ///     if( !is )
    ///         NCBI_THROW( Exception, eFormat, "file too short" );
    /// ```
    FileTooShort,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:113-115
    /// ```c
    ///     if( unit_size == 0 || unit_size > 16 )
    ///         NCBI_THROW( Exception, eBadParam,
    ///             "unit size must be in [1,16]" );
    /// ```
    UnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:124-126
    /// ```c
    ///     if( k == 0U || k > (Uint4)(2*unit_size - 1) )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "hash key size must be in [1,2*unit_size - 1]" );
    /// ```
    HashKeySize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:128-130
    /// ```c
    ///     if( roff > 32 - k )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "offset must by in [0,32 - hash_key_size]" );
    /// ```
    HashOffset,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:132-134
    /// ```c
    ///     if( bc == 0 || bc > 32 - k )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "shift must be in [1, 32 - hash_key_size]" );
    /// ```
    HashShift,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:196-198
    /// ```c
    ///     if( !input_stream.read( (char *)ht, ht_size*sizeof( Uint4 ) ) )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough data to fill the hash table" );
    /// ```
    HashTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:207-209
    /// ```c
    ///     if( !input_stream.read( (char *)vt, M*sizeof( Uint2 ) ) )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough data to fill the values table" );
    /// ```
    ValuesTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:78-79
    /// ```c
    ///         if( iflen < HEADER_LEN + TRAILER_LEN )
    ///             NCBI_THROW( Exception, eFormat, "wrong file size" );
    /// ```
    BinFileSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:98-99
    /// ```c
    ///         if( unit_size == 0 || unit_size > 16 )
    ///             NCBI_THROW( Exception, eFormat, "illegal unit size" );
    /// ```
    BinUnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:107-108
    /// ```c
    ///         if( datalen%(2*sizeof( Uint4 )) != 0 )
    ///             NCBI_THROW( Exception, eFormat, "wrong length" );
    /// ```
    BinLength,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:100-102
    /// ```c
    ///                 if( unit_size == 0 || unit_size > 16 )
    ///                     NCBI_THROW( Exception, eBadParam,
    ///                                 "unit size must be in [1,16]" );
    /// ```
    OAsciiUnitSize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:113-116
    /// ```c
    ///                     if( k == 0U || k > (Uint4)(2*unit_size - 1) )
    ///                         NCBI_THROW(
    ///                             Exception, eBadHashParam,
    ///                             "hash key size must be in [1,2*unit_size - 1]" );
    /// ```
    OAsciiHashKeySize,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:118-121
    /// ```c
    ///                     if( roff > 32 - k )
    ///                         NCBI_THROW(
    ///                             Exception, eBadHashParam,
    ///                             "offset must by in [0,32 - hash_key_size]" );
    /// ```
    OAsciiHashOffset,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:123-127
    /// ```c
    ///                     if( bc == 0 || bc > 32 - k )
    ///                         NCBI_THROW(
    ///                             Exception, eBadHashParam,
    ///                             "shift must be in "
    ///                             "[1, 32 - hash_key_size]" );
    /// ```
    OAsciiHashShift,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:144-145
    /// ```c
    ///     if( linenum < HEADER_LINES )
    ///         NCBI_THROW( Exception, eFormat, "file too short" );
    /// ```
    OAsciiFileTooShort,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:175-177
    /// ```c
    ///     if( linenum < ht_size )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough lines to fill the hash table" );
    /// ```
    OAsciiHashTable,
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_oascii.cpp:191-193
    /// ```c
    ///     if( linenum < M )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough lines to fill the values table" );
    /// ```
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
    /// ```c
    ///         case eStreamOpenFail:   return "open failed";
    ///         case eBadHashParam:     return "bad hash parameter";
    ///         case eBadParam:         return "bad parameter";
    ///         case eFormat:           return "format error";
    ///         case eAlloc:            return "allocation failure";
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_bin.cpp:49-50
    /// ```c
    ///         case eStreamOpenFail:   return "open failed";
    ///         case eFormat:           return "file format error";
    /// ```
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

    // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:174-194
    // ```c
    //     if( skip > 0 ) {
    //         char * buf( new char[skip] );
    //         check.read( buf, skip );
    //         delete[] buf;
    //     }
    //     else if( !md.empty() ) {
    //         string line;
    //         for( size_t i( 0 ); i < md.size(); ++i ) getline( check, line );
    //     }
    //
    //     if( check ) {
    //         Uint4 data = 1;
    //
    //         if( check.read( (char *)&data, sizeof( Uint4 ) ) ) {
    //             if( data == 0 ) return eBinary;
    //             if( data == 0x41414141 ) return eOAscii;
    //             if( data == 1 || data == 2 ) return eOBinary;
    //             return eAscii;
    //         }
    //     }
    //
    //     return eUnknown;
    // ```
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
/// ```c
/// static bool CheckBin(
///         string const & name, vector< string > & md, size_t & skip ) {
///     skip = 0;
///     CNcbiIfstream input_stream( name.c_str(), IOS_BASE::binary );
///
///     if( input_stream ) {
///         Uint4 word( 0 );
///         input_stream.read( (char *)&word, sizeof( word ) );
///
///         if( word == 3 ) {
///             skip = 8;
///             input_stream.read( (char *)&word, sizeof( word ) );
///             skip += word;
///             char * data( new char[word + 2] );
///             data[word] = data[word + 1] = 0;
///             input_stream.read( data, word );
///
///             for( char * d( data ); *d != 0; d += strlen( d ) + 1 ) {
///                 md.push_back( string( d ) );
///             }
/// ...
///             return true;
///         }
///         else if( word < 3 ) return true;
///         else return false;
///     }
/// ...
/// static void GetMetaDataLines(
///         string const & name, vector< string > & md, size_t & skip ) {
///     if( !CFile( name ).Exists() ) return;
///     if( CheckBin( name, md, skip ) ) return;
///     CNcbiIfstream check( name.c_str() );
///     string line;
///
///     while( check ) {
///         getline( check, line );
///
///         if( line.size() >= 2 && line[0] == '#' && line[1] == '#' ) {
///             md.push_back( line );
///         }
///         else break;
///     }
/// }
/// ```
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
}

enum Uset {
    Hash(UsetHash),
    Array(UsetArray),
}

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_uset_hash.hpp:122-128
/// ```c
///         Uint1 unit_size;    /**<\internal Unit size in bases. */
///         Uint1 k;            /**<\internal Hash key size in bits. */
///         Uint1 roff;         /**<\internal Right offset in bits. */
///         Uint1 bc;           /**<\internal Bit size of "number of collisions" field. */
///         Uint4 M;            /**<\internal Size of the secondary table (in 2-byte words). */
///         Uint4 cmask;        /**<\internal Mask used to extract the "number of collisions" field. */
/// ```
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
/// ```c
///         struct entry
///         {
///             Uint4 u;    /**\internal The unit value. */
///             Uint4 c;    /**\internal The unit count. */
///         };
///
///         Uint1 unit_size;            /**<\internal The unit size. */
/// ```
struct UsetArray {
    unit_size: u8,
    entries: Vec<(u32, u32)>,
}

impl UnitCounts {
    /// `CSeqMaskerIstatFactory::create` for the two formats blastn accepts.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:317-340
    /// ```c
    ///             case eBinary: res = new CSeqMaskerIstatBin(
    ///                                   name, threshold, textend,
    ///                                   max_count, use_max_count,
    ///                                   min_count, use_min_count,
    ///                                   skip );
    ///                           break;
    /// ...
    ///             case eOBinary: res = new CSeqMaskerIstatOBinary(
    ///                                    name,
    ///                                    threshold, textend,
    ///                                    max_count, use_max_count,
    ///                                    min_count, use_min_count,
    ///                                    use_ba, skip );
    ///                            break;
    ///
    /// 			default: NCBI_THROW(Exception, eCreateFail, "unrecognized unit counts format");
    /// ```
    /// The `-t_*_pct` settings of the windowmasker application are not ported (blastn passes
    /// -1, `seq_masker.hpp:178-181`).
    pub fn open(path: &Path, params: &MaskerParams) -> Result<Self, WinMaskError> {
        let name = path.to_string_lossy().into_owned();
        let mut file = File::open(path).map_err(|_| WinMaskError::Open { name })?;
        Self::read(&mut file, params)
    }

    /// [`UnitCounts::open`] on an opened file or an in-memory copy of one.
    pub fn read<R: Read + Seek>(r: &mut R, params: &MaskerParams) -> Result<Self, WinMaskError> {
        let found = discover(r)?;
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
        };
        check_pct_lines(&counts.metadata, found.stat_type)?;
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
        check_stat_algo_version(&counts.metadata).map_err(|e| counts.exception(e))?;
        Ok(counts)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_ascii.cpp:65-176
    /// ```c
    ///     bool start = true;
    ///     Uint4 linenum = start_line;
    ///     Uint4 ambig_len = kMax_UI4;
    ///     string line;
    ///     for( Uint4 i( 0 ); i < start_line; ++i ) getline( input_stream, line );
    ///
    ///     while( input_stream )
    ///     {
    ///         line.erase();
    ///         getline( input_stream, line );
    ///         ++linenum;
    ///         if( line.length() == 0 ) continue;
    ///         if( line[0] == '#' ) continue;
    ///
    ///         // Check if we have a precomputed parameter.
    ///         if( line[0] == '>' )
    ///         {
    ///             SIZE_TYPE name_end = line.find_first_of( " \t", 0 );
    ///             SIZE_TYPE val_start = line.find_first_not_of( " \t", name_end );
    ///
    ///             if( name_end == NPOS || val_start == NPOS )
    ///             {
    ///                 CNcbiOstrstream str;
    ///                 str << "at line " << linenum;
    ///                 string msg = CNcbiOstrstreamToString(str);
    ///                 NCBI_THROW( Exception, eSyntax, msg);
    ///             }
    ///
    ///             string name = line.substr( 1, name_end - 1 );
    ///
    ///             if( name == "t_threshold" && get_threshold() == 0 )
    ///                 set_threshold(
    ///                     NStr::StringToUInt(line.substr(val_start, NPOS), 0, 0));
    ///
    ///             if( name == "t_extend" && get_textend() == 0 )
    ///                 set_textend(
    ///                     NStr::StringToUInt(line.substr(val_start, NPOS), 0, 0));
    ///
    ///             if( name == "t_low" )
    ///                 set_min_count(
    ///                     NStr::StringToUInt(line.substr(val_start, NPOS), 0, 0));
    ///
    ///             if( name == "t_high" && get_max_count() == 0 )
    ///                 set_max_count(
    ///                     NStr::StringToUInt(line.substr(val_start, NPOS), 0, 0));
    ///
    ///             continue;
    ///         }
    ///
    ///         if( start )
    ///         {
    ///             start = false;
    ///             uset.set_unit_size(
    ///                 static_cast< Uint1 >( NStr::StringToUInt( line ) ) );
    ///             continue;
    ///         }
    ///
    ///         SIZE_TYPE unit_start = line.find_first_not_of( " \t", 0 );
    ///         SIZE_TYPE unit_end   = line.find_first_of( " \t", unit_start );
    ///         SIZE_TYPE cnt_start  = line.find_first_not_of( " \t", unit_end );
    ///
    ///         if( unit_start == NPOS || unit_end == NPOS || cnt_start == NPOS )
    ///         {
    /// ...
    ///             NCBI_THROW( Exception, eSyntax, msg );
    ///         }
    ///
    ///         Uint4 unit = NStr::StringToUInt(line.substr(unit_start,
    ///                                                     unit_end - unit_start),
    ///                                         0, 16);
    ///         Uint4 cnt = NStr::StringToUInt(line.substr(cnt_start));
    ///
    ///         if( cnt < ambig_len ) {
    ///             ambig_len = cnt;
    ///             set_ambig_unit( unit );
    ///         }
    ///
    ///         if( cnt >= get_min_count() )
    ///             uset.add_info( unit, cnt );
    ///     }
    ///
    ///     string bad_param;
    ///
    ///     if( get_threshold() == 0 )
    ///         bad_param += "t_threhold ";
    /// ...
    ///     if( !bad_param.empty() )
    ///         NCBI_THROW( Exception, eParam, bad_param );
    ///
    ///     if( get_use_min_count() == 0 )
    ///       set_use_min_count( (get_min_count() + 1)/2 );
    ///
    ///     if( get_use_max_count() == 0 )
    ///       set_use_max_count( get_max_count() );
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_simple.cpp:55-67
    /// ```c
    /// void CSeqMaskerUsetSimple::add_info( Uint4 unit, Uint4 count )
    /// {
    ///     if( !units.empty() && unit <= units[units.size() - 1] )
    ///     {
    ///         ostringstream s;
    ///         s << "last unit: " << hex << units[units.size() - 1]
    ///           << " ; adding " << hex << unit;
    ///         NCBI_THROW( Exception, eBadOrder, s.str() );
    ///     }
    /// ```
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
    /// ```c
    ///     Uint4 linenum = 0;
    ///     string line;
    ///     for( Uint4 i( 0 ); i < start_line; ++i ) getline( input_stream, line );
    ///     Uint1 unit_size( 0 );
    ///     Uint4 k, roff, bc;
    ///     Uint4 t_low    = 0,
    ///           t_extend = 0,
    ///           t_thres  = 0,
    ///           t_high   = 0,
    ///           M;
    ///
    ///     while( getline( input_stream, line ) )
    ///     {
    ///         ++linenum;
    ///
    ///         switch( linenum )
    ///         {
    ///             case 1: break; //skip the file format identifier
    ///             case 2: // unit size
    ///
    ///                 unit_size = (Uint1)atoi( line.c_str() );
    /// ...
    ///             case 3: // hash table parameters
    ///
    ///                 {
    ///                     istringstream i( line );
    ///                     i >> M >> k >> roff >> bc;
    /// ...
    ///             case 4: t_low    = atoi( line.c_str() ); break;
    ///             case 5: t_extend = atoi( line.c_str() ); break;
    ///             case 6: t_thres  = atoi( line.c_str() ); break;
    ///             case 7: t_high   = atoi( line.c_str() ); break;
    ///
    ///             default: break;
    ///         }
    ///
    ///         if( linenum == HEADER_LINES )
    ///             break;
    ///     }
    ///
    ///     if( linenum < HEADER_LINES )
    ///         NCBI_THROW( Exception, eFormat, "file too short" );
    ///
    ///     set_min_count( t_low );
    /// ...
    ///     Uint4 ht_size = (1ULL<<k);
    ///     Uint4 * ht = new Uint4[ht_size];
    /// ...
    ///     for( linenum = 0;
    ///          linenum < ht_size && getline( input_stream, line );
    ///          ++linenum )
    ///         ht[linenum] = atoi( line.c_str() );
    /// ...
    ///     Uint2 * vt = new Uint2[M];
    /// ...
    ///     for( linenum = 0;
    ///          linenum < M && getline( input_stream, line );
    ///          ++linenum )
    ///         vt[linenum] = (Uint2)atoi( line.c_str() );
    /// ```
    /// The numbers are glibc's `atoi` ([`ncbi_str::c_atoi`]) truncated to the field, and those of
    /// the hash table parameters are read as libstdc++'s `>>` reads an `unsigned int`
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
    /// ```c
    ///     word = readWord( input_stream );
    ///
    ///     if( word == 1 )
    ///         use_opt = false;
    ///
    ///     word = readWord( input_stream );
    ///     unit_size = (Uint1)word;
    ///
    ///     if( unit_size == 0 || unit_size > 16 )
    ///         NCBI_THROW( Exception, eBadParam,
    ///             "unit size must be in [1,16]" );
    ///
    ///     uset.setUnitSize( unit_size );
    ///
    ///     M    = readWord( input_stream );
    ///     k    = readWord( input_stream );
    ///     roff = readWord( input_stream );
    ///     bc   = readWord( input_stream );
    ///
    ///     if( k == 0U || k > (Uint4)(2*unit_size - 1) )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "hash key size must be in [1,2*unit_size - 1]" );
    ///
    ///     if( roff > 32 - k )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "offset must by in [0,32 - hash_key_size]" );
    ///
    ///     if( bc == 0 || bc > 32 - k )
    ///         NCBI_THROW( Exception, eBadHashParam,
    ///                     "shift must be in [1, 32 - hash_key_size]" );
    ///
    ///     t_low    = readWord( input_stream );
    ///     t_extend = readWord( input_stream );
    ///     t_thres  = readWord( input_stream );
    ///     t_high   = readWord( input_stream );
    ///
    ///     set_min_count( t_low );
    /// ...
    ///     if( use_opt )
    ///     {
    ///         Uint4 divisor = readWord( input_stream );
    ///
    ///         if( divisor > 0 )
    ///         {
    ///             Uint8 total = (1ULL<<(2*unit_size));
    ///             Uint4 cba_size = (Uint4)(total/(8*sizeof( Uint4 )));
    /// ...
    ///             else if( !input_stream.read( (char *)cba, cba_size*sizeof( Uint4 ) ) )
    ///             {
    ///                 ERR_POST( Warning << "file read failed: "
    ///                                   << "bit array optimizations are not used." );
    /// ...
    ///     Uint4 ht_size = (1ULL<<k);
    /// ...
    ///     if( !input_stream.read( (char *)ht, ht_size*sizeof( Uint4 ) ) )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough data to fill the hash table" );
    /// ...
    ///     if( !input_stream.read( (char *)vt, M*sizeof( Uint2 ) ) )
    ///         NCBI_THROW( Exception, eFormat,
    ///                     "not enough data to fill the values table" );
    /// ```
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
                // ```c
                //             Uint8 total = (1ULL<<(2*unit_size));
                //             Uint4 cba_size = (Uint4)(total/(8*sizeof( Uint4 )));
                // ...
                //             else if( !input_stream.read( (char *)cba, cba_size*sizeof( Uint4 ) ) )
                //             {
                //                 ERR_POST( Warning << "file read failed: "
                //                                   << "bit array optimizations are not used." );
                //                 delete[] cba;
                //                 cba = 0;
                //             }
                // ...
                //             optimization_data opt_data( 8*sizeof( Uint4 ), cba );
                //             set_optimization_data( opt_data );
                // ```
                // `arg_use_ba` is true for blastn and the application (`winmask_filter.cpp:121`,
                // `win_mask_config.cpp:361`), so the array is kept. Units of 1 or 2 letters give
                // an array of 0 words that is still in use (`new Uint4[0]`); see `SeqMasker::mask`.
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
        // ```c
        //     k = arg_k;
        //     roff = arg_roff;
        //     bc = arg_bc;
        //     cmask = (1ULL<<bc) - 1;
        // ```
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
    /// ```c
    ///         iflen = (streamsize)input_file.GetLength();
    ///
    ///         if( iflen < HEADER_LEN + TRAILER_LEN )
    ///             NCBI_THROW( Exception, eFormat, "wrong file size" );
    /// ...
    ///     streamsize adj( in_stream.tellg() );
    ///     Uint4 data;
    ///     in_stream.read( (char *)&data, sizeof( Uint4 ) );
    ///
    ///     {
    ///         in_stream.read( (char *)&data, sizeof( Uint4 ) );
    ///         Uint1 us = (Uint1)data;
    ///
    ///         if( us == 0 || us > 16 )
    ///             NCBI_THROW( Exception, eFormat, "illegal unit size" );
    /// ...
    ///         streamsize datalen = iflen - HEADER_LEN - TRAILER_LEN - adj;
    ///
    ///         if( datalen%(2*sizeof( Uint4 )) != 0 )
    ///             NCBI_THROW( Exception, eFormat, "wrong length" );
    /// ...
    ///         if( datalen > 0 )
    ///         {
    ///             cdata = new Uint4[datalen/sizeof( Uint4 )];
    ///             in_stream.read( (char *)cdata, datalen );
    ///             uset.add_info( cdata, datalen/sizeof( Uint4 ) );
    ///         }
    ///     }
    ///
    ///     in_stream.read( (char *)&data, sizeof( Uint4 ) );
    ///
    ///     set_min_count( data );
    /// ```
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
            let words = s
                .read_u32s(datalen as u64 / 4)
                .ok_or_else(|| WinMaskError::format(FormatSite::BinLength))?;
            entries = words.chunks_exact(2).map(|e| (e[0], e[1])).collect();
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
    /// ```c
    ///     set_min_count( t_low );
    ///
    ///     if( get_textend() == 0 )
    ///         set_textend( t_extend );
    ///
    ///     if( get_threshold() == 0 )
    ///         set_threshold( t_thres );
    ///
    ///     if( get_max_count() == 0 )
    ///         set_max_count( t_high );
    ///
    ///     if( get_use_min_count() == 0 )
    ///       set_use_min_count( (get_min_count() + 1)/2 );
    ///
    ///     if( get_use_max_count() == 0 )
    ///       set_use_max_count( get_max_count() );
    /// ```
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
    /// ```c
    ///     void set_min_count( Uint4 arg_min_count )
    ///     {
    ///         if( min_count != 0 && min_count < arg_min_count ) {
    ///             ERR_POST( Warning << "Requested value of t_low ("
    ///                               << min_count
    ///                               << ") is less than the one stored with the "
    ///                               << "N-mer counts (" << arg_min_count << ")."
    ///                               << "The value " << arg_min_count
    ///                               << " will be used." );
    ///             min_count = arg_min_count;
    ///         }
    ///         else if( min_count == 0 ) {
    ///             min_count = arg_min_count;
    ///         }
    ///     }
    /// ```
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
    /// ```c
    /// Uint4 CSeqMaskerIstatOBinary::at( Uint4 unit ) const
    /// {
    ///     Uint4 res = uset.get_info( unit );
    ///
    ///     if( res == 0 || res < get_min_count() )
    ///         return get_use_min_count();
    ///
    ///     return (res > get_max_count()) ? get_use_max_count() : res;
    /// }
    /// ```
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

    /// The warnings NCBI posts while reading the file (`ERR_POST( Warning << ...)`).
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_istat.hpp:155-156
    /// ```c
    ///     const optimization_data * get_optimization_data() const
    ///     { return opt_data_.cba_ == 0 ? 0 : &opt_data_; }
    /// ```
    pub(super) fn cache_bit_array(&self) -> Option<&[u32]> {
        self.cba.as_deref()
    }
}

impl UsetHash {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:73-111
    /// ```c
    ///     Uint4 runit = CSeqMaskerUtil::reverse_complement( unit, unit_size );
    ///
    ///     if( runit < unit )
    ///         unit = runit;
    ///
    ///     pair< Uint4, Uint1 > hash = CSeqMaskerUtil::hash_code( unit, k, roff );
    ///     Uint4 hval = htp[hash.first];
    ///     Uint4 coll = hval&cmask;
    ///
    ///     if( coll == 0 )
    ///         return 0;
    ///     else if( coll == 1 )
    ///     {
    ///         if( hash.second != (hval>>24) )
    ///             return 0;
    ///         else return (hval>>bc)&0xFFF;
    ///     }
    ///     else
    ///     {
    ///         if( (hval>>bc) + coll > M )
    ///         {
    ///             ostringstream r;
    ///             r << "bad index at key " << hash.first
    ///               << " : " << htp[hash.first];
    ///             NCBI_THROW( Exception, eBadIndex, r.str() );
    ///         }
    ///
    ///         const Uint2 * start = vtp + (hval>>bc);
    ///         const Uint2 * end = start + coll;
    ///
    ///         for( ; start < end; ++start )
    ///             if( ((*start)>>9) == hash.second )
    ///                 return (*start)&0x1FF;
    ///
    ///         return 0;
    ///     }
    /// ```
    #[inline]
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
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_array.cpp:64-82
    /// ```c
    ///     const entry * ud = unit_data.get();
    ///
    ///     if( ud == 0 )
    ///         return 0;
    ///
    ///     Uint4 runit = CSeqMaskerUtil::reverse_complement( unit, unit_size );
    ///
    ///     if( runit < unit )
    ///         unit = runit;
    ///
    ///     entry target = { unit, 0 };
    ///     const entry * r = lower_bound( ud, ud + asize, target, less< entry >() );
    ///
    ///     if( r == ud + asize || r->u != unit )
    ///         return 0;
    ///     else return r->c;
    /// ```
    /// `lower_bound` is libstdc++'s halving loop (`stl_algobase.h`, `__lower_bound`), so a file
    /// whose units are out of order is probed at the same entries as by NCBI's Linux builds.
    #[inline]
    fn get_info(&self, unit: u32) -> u32 {
        if self.entries.is_empty() {
            return 0;
        }
        let runit = reverse_complement(unit, self.unit_size);
        let unit = if runit < unit { runit } else { unit };
        let mut first = 0usize;
        let mut len = self.entries.len();
        while len > 0 {
            let half = len >> 1;
            let middle = first + half;
            if self.entries[middle].0 < unit {
                first = middle + 1;
                len = len - half - 1;
            } else {
                len = half;
            }
        }
        match self.entries.get(first) {
            Some(&(u, c)) if u == unit => c,
            _ => 0,
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_util.cpp:54-65
/// ```c
/// Uint4 CSeqMaskerUtil::reverse_complement( Uint4 seq, Uint1 size )
/// {
///     Uint4 result( 0 );
///
///     for( Uint1 i( 0 ); i < size; ++i )
///     {
///         Uint4 letter( ~(((seq>>(2*i))&0x3)|(~0x3)) );
///         result = (result<<2)|letter;
///     }
///
///     return result;
/// }
/// ```
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
/// ```c
///     static pair< Uint4, Uint1 > hash_code( Uint4 unit,
///                                            Uint1 k, Uint1 roff )
///     {
///         return make_pair( (unit>>roff)&((((Uint4)1)<<k) - 1),
///                           (Uint1)(((unit>>(roff + k))<<roff)
///                             + (unit&((1<<roff) - 1))) );
///     }
/// ```
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

    /// `is.read( (char *)&word, sizeof( Uint4 ) )`, in the machine's byte order.
    fn read_into(&mut self, word: &mut u32) {
        let mut bytes = word.to_ne_bytes();
        self.read(&mut bytes);
        *word = u32::from_ne_bytes(bytes);
    }

    /// `is.read( (char *)array, n*sizeof( Uint4 ) )` into a new array; `None` when the read
    /// fails (checked against the file's length before the array is allocated).
    fn read_u32s(&mut self, n: u64) -> Option<Vec<u32>> {
        let bytes = n.checked_mul(4)?;
        if self.failed || bytes > self.left {
            self.failed = true;
            self.left = 0;
            return None;
        }
        let mut out = Vec::with_capacity(n as usize);
        let mut buf = vec![0u8; 1 << 16];
        let mut todo = bytes as usize;
        while todo > 0 {
            let take = todo.min(buf.len());
            if self.read(&mut buf[..take]) < take {
                return None;
            }
            out.extend(
                buf[..take]
                    .chunks_exact(4)
                    .map(|w| u32::from_ne_bytes([w[0], w[1], w[2], w[3]])),
            );
            todo -= take;
        }
        Some(out)
    }

    /// `is.read( (char *)array, n*sizeof( Uint2 ) )` into a new array; `None` when it fails.
    fn read_u16s(&mut self, n: u64) -> Option<Vec<u16>> {
        let bytes = n.checked_mul(2)?;
        if self.failed || bytes > self.left {
            self.failed = true;
            self.left = 0;
            return None;
        }
        let mut out = Vec::with_capacity(n as usize);
        let mut buf = vec![0u8; 1 << 16];
        let mut todo = bytes as usize;
        while todo > 0 {
            let take = todo.min(buf.len());
            if self.read(&mut buf[..take]) < take {
                return None;
            }
            out.extend(
                buf[..take]
                    .chunks_exact(2)
                    .map(|w| u16::from_ne_bytes([w[0], w[1]])),
            );
            todo -= take;
        }
        Some(out)
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:239-276
/// ```c
///         if( stat_type == eAscii || stat_type == eOAscii )
///         {
///             for( auto const & line : md )
///             {
///                 if( !line.empty() && line.substr( 0, 6 ) == "##pct:" )
///                 {
///                     double pct( 0.0 );
///                     Uint4 count;
///                     std::istringstream iss(
///                         line.substr( 7, line.size() - 7 ) );
///                     iss >> count >> pct;
///                     count_map.push_back( pct );
///                 }
///             }
/// ...
///         if( stat_type == eBinary || stat_type == eOBinary )
///         {
///             for( auto const & line : md )
///             {
///                 if( !line.empty() && line.substr( 0, 6 ) == "##pct:" )
///                 {
///                     std::istringstream iss(
///                         line.substr( 7, line.size() - 7 ) );
/// ...
///                     break;
///                 }
///             }
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:361-366
/// ```c
///     catch( std::exception & e )
///     {
///         NCBI_THROW( Exception, eCreateFail,
///                     std::string( "could not create a unit counts container: " ) +
///                     e.what() );
///     }
/// ```
/// The binary formats read the first `##pct:` line only (the loop breaks after it), the ascii
/// ones every such line. A line of exactly those six letters throws `std::out_of_range` from
/// `substr( 7 )`; the text is libstdc++'s, as in NCBI's Linux builds. The percentages themselves
/// are used only by the application's `-t_*_pct` settings, which are not ported.
fn check_pct_lines(metadata: &[Vec<u8>], stat_type: StatType) -> Result<(), WinMaskError> {
    let mut pct = metadata.iter().filter(|line| line.starts_with(b"##pct:"));
    let short = match stat_type {
        StatType::Binary | StatType::OBinary => pct.next().is_some_and(|line| line.len() == 6),
        StatType::Ascii | StatType::OAscii => pct.any(|line| line.len() == 6),
        StatType::Unknown => false,
    };
    if short {
        return Err(WinMaskError::CreateFail {
            message: "could not create a unit counts container: basic_string::substr: __pos (which is 7) > this->size() (which is 6)".to_string(),
        });
    }
    Ok(())
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_factory.cpp:111-146
/// ```c
/// void ExtractStatAlgoVersion(
///         vector< string > const & md, CSeqMaskerVersion & ver ) {
///     static char const * DIGITS = "0123456789";
///     string t;
///
///     for( vector< string >::const_iterator i( md.begin() );
///             i != md.end(); ++i ) {
///         string::size_type p( i->find( ':', 2 ) ), pp;
///
///         if( p != string::npos ) {
///             if( i->substr( 2, p -2 ) ==
///                     CSeqMaskerOstat::STAT_ALGO_COMPONENT_NAME ) {
///                 int major = 0, minor = 0, patch = 0;
///                 pp = p + 1;
///                 p = i->find( '.', pp );
///                 if( p == string::npos ) continue;
///                 t = i->substr( pp, p - pp );
///                 if( t.find_first_not_of( DIGITS ) != string::npos ) continue;
///                 major = NStr::StringToInt( t );
///                 pp = p + 1;
///                 p = i->find( '.', pp );
///                 if( p == string::npos ) continue;
///                 t = i->substr( pp, p - pp );
///                 if( t.find_first_not_of( DIGITS ) != string::npos ) continue;
///                 minor = NStr::StringToInt( t );
///                 t = i->substr( p + 1 );
///                 if( t.find_first_not_of( DIGITS ) != string::npos ) continue;
///                 patch = NStr::StringToInt( t );
///                 ver = CSeqMaskerVersion(
///                         CSeqMaskerOstat::STAT_ALGO_COMPONENT_NAME,
///                         major, minor, patch );
///                 return;
/// ```
/// `CSeqMaskerIstatFactory::create` calls it after the reader (`:343-348`), inside its `try`:
/// a part of digits that `StringToInt` cannot convert (empty, or beyond `int`) is its
/// exception, which the factory rethrows. The version itself is not used by the masker.
fn check_stat_algo_version(metadata: &[Vec<u8>]) -> Result<(), NcbiException> {
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
        ncbi_str::string_to_int(t, 0)?;
        let pp = p + 1;
        let Some(p) = find(line, b'.', pp) else {
            continue;
        };
        let t = &line[pp..p];
        if !digits(t) {
            continue;
        }
        ncbi_str::string_to_int(t, 0)?;
        let t = &line[p + 1..];
        if !digits(t) {
            continue;
        }
        ncbi_str::string_to_int(t, 0)?;
        return Ok(());
    }
    Ok(())
}

/// The function of the ascii reader's throws.
const ASCII_NEW: Option<&str> = Some("ncbi::CSeqMaskerIstatAscii::CSeqMaskerIstatAscii");

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_ascii.cpp:40-50
/// ```c
/// CSeqMaskerIstatAscii::Exception::GetErrCodeString() const
/// {
///     switch( GetErrCode() )
///     {
///         case eStreamOpenFail:   return "open failed";
///         case eSyntax:           return "syntax error";
///         case eParam:            return "bad parameter value";
/// ```
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
/// ```c
///     if( !bad_param.empty() )
///         NCBI_THROW( Exception, eParam, bad_param );
/// ```
const ASCII_BAD_PARAM: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker_istat_ascii.cpp",
    line: 169,
    kind: "Exception::bad parameter value",
    function: ASCII_NEW,
};

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_simple.cpp:44-52
/// ```c
/// const char * CSeqMaskerUsetSimple::Exception::GetErrCodeString() const
/// {
///     switch( GetErrCode() )
///     {
///         case eBadOrder:     return "bad unit order";
/// ```
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

/// `istringstream( line ) >> v0 >> v1 >> ...` into `unsigned int`s, as libstdc++ reads them
/// (`num_get::_M_extract_int` with the stream's default flags): white space, then a sign and the
/// decimal digits up to the first other byte; no digit is 0 and a failure, a value beyond
/// `UINT_MAX` is `UINT_MAX` and a failure, and a `-` gives the value's negation modulo 2^32.
/// After a failure, and once the line has ended, the stream reads nothing more: the values and
/// how many of them `>>` set (the first ones; a failed one is set), the others 0 here (NCBI
/// leaves them uninitialized; see [`UnitCounts::read_oascii`]).
pub(super) fn extract_uints<const N: usize>(line: &[u8]) -> ([u32; N], usize) {
    let mut values = [0u32; N];
    let mut set = 0;
    let mut pos = 0;
    for value in values.iter_mut() {
        // The sentry: skip white space; at the end of the line it fails.
        while pos < line.len() && ncbi_str::is_space(line[pos]) {
            pos += 1;
        }
        if pos == line.len() {
            break;
        }
        set += 1;
        let negative = line[pos] == b'-';
        if negative || line[pos] == b'+' {
            pos += 1;
        }
        let digits_from = pos;
        let mut result: u32 = 0;
        let mut overflow = false;
        while pos < line.len() && line[pos].is_ascii_digit() {
            let digit = u32::from(line[pos] - b'0');
            if result > u32::MAX / 10 {
                overflow = true;
            } else {
                result *= 10;
                overflow |= result > u32::MAX - digit;
                result = result.wrapping_add(digit);
            }
            pos += 1;
        }
        if pos == digits_from {
            break;
        }
        if overflow {
            *value = u32::MAX;
            break;
        }
        *value = if negative {
            result.wrapping_neg()
        } else {
            result
        };
        if pos == line.len() {
            break;
        }
    }
    (values, set)
}

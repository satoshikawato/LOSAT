//! `x_ProcessMaskData` of makeblastdb (`c++/src/app/blastdb/makeblastdb.cpp`): the options
//! `-mask_data`, `-mask_id`, `-mask_desc` and `-gi_mask_name` are split into lists, every mask
//! file is read with [`MaskFileReader`], and the masks are registered with the `BuildDatabase`
//! and collected for the sequences it will write.
//!
//! # The order in which `LOSAT makedb` must call things when `-mask_data` is wired
//!
//! Line NP, M4: nothing calls [`process_mask_data`] yet (`LOSAT makedb` still refuses
//! `-mask_data`). `makeblastdb`'s `x_BuildDatabase` does, in this order
//! (makeblastdb.cpp:1231-1296; `Run`, makeblastdb.cpp:1322-1331):
//!
//! ```text
//! 1. BuildDatabase::new(..)          CBuildDatabase constructor: the head of the log
//!                                    ("Building a new DB, current time", "New DB name:", ...)
//!                                    to the log stream (standard output, or -logfile)
//! 2. set_max_file_size(bytes)        after "Maximum file size: " (log, standard output)
//! 3. process_mask_data(..)           x_ProcessMaskData: one "Mask file: <name>" log line per file
//!                                    read; "Ignoring mask file ..." to the standard error stream
//! 4. add_fasta(..)                   x_ProcessInputData: the checks of the input file ("File X
//!                                    does not exist", "is empty") come after the mask files
//! 5. end_build(..)                   m_DB->EndBuild(): the closing lines of the log
//! 6. (CATCH_ALL)                     a failure of 1-5: the message and the status of
//!                                    `MaskDataError::into_makedb_error` /
//!                                    `MakedbError::exit_code`, one line on the standard error
//! 7. post_destruction_diagnostics    ~CBuildDatabase, when the application object goes: after
//!                                    the exception message and after the metadata file
//! ```
//!
//! A failure of step 3 ends the run at once: steps 4 and 5 are not run (no database is written
//! and the input file is not looked at), step 6 prints the message, and step 7 still posts the
//! "No sequences matched any of the masks provided." message when a mask list with masks was
//! read before the failure (`BuildDatabase::mask_data_source_mut` is set). A file that ends where
//! `more TRUE` announced a continuation ([`MaskDataError::EofInContinuation`]) is such a failure
//! with no message and exit status 0: makeblastdb leaves without a database and without an error
//! (`CATCH_ALL` ignores a `CIOException` that is not a flush failure).
//!
//! Streams: the log (steps 1, 2, 3 "Mask file:", 5) is makeblastdb's standard output; the
//! diagnostics (`ERR_POST`: "Ignoring mask file", the warnings of the FASTA reader, the message of
//! step 7) and the message of step 6 (`LOG_POST`) are on the standard error stream, in the order
//! of the events.
//!
//! The options `-mask_id`, `-mask_desc`, `-gi_mask_name` depend on others in the argument parser
//! (makeblastdb.cpp:276-288): `-mask_id` needs `-mask_data`, `-mask_desc` needs `-mask_id`,
//! `-gi_mask_name` needs `-gi_mask`, which needs `-parse_seqids` and excludes `-mask_id`. LOSAT
//! has no `-parse_seqids`, so the CLI refuses `-gi_mask` and `-gi_mask_name` (NCBI's parser
//! prints the usage text for `-gi_mask_name` without `-gi_mask`, status 1, before any of this
//! runs), and `gi_mask_name` of [`MaskDataOptions`] is always `None` there.

use std::fmt;
use std::io::Write;
use std::path::{Path, PathBuf};

use super::MaskedRangeSet;
use crate::writedb::build_db::BuildDatabase;
use crate::writedb::fasta_source::Warn;
use crate::writedb::makeblastdb::MakedbError;
use crate::writedb::mask_info::{MaskFileError, MaskFileReader};
use crate::writedb::WriteDbError;

/// The values of the options that `x_ProcessMaskData` reads: `None` for an option that was not
/// given (`CArgValue::HasValue()` false), `Some("")` for one given with an empty value.
#[derive(Debug, Clone, Copy, Default)]
pub struct MaskDataOptions<'a> {
    /// `-mask_data`.
    pub mask_data: Option<&'a str>,
    /// `-mask_id`.
    pub mask_id: Option<&'a str>,
    /// `-mask_desc`.
    pub mask_desc: Option<&'a str>,
    /// `-gi_mask_name`.
    pub gi_mask_name: Option<&'a str>,
    /// The directory in which a relative mask file name is looked up (`None`: the current
    /// directory, as makeblastdb does). The names in the log, the messages and the registry stay
    /// as given. A hook for tests, which cannot change the directory of the process.
    pub root: Option<&'a Path>,
}

/// Why `x_ProcessMaskData` stopped. Every one of them ends the run of makeblastdb (the exception
/// leaves `x_BuildDatabase`); [`MaskDataError::into_makedb_error`] gives the message and the exit
/// status.
#[derive(Debug)]
pub enum MaskDataError {
    /// An exception whose text is NCBI's: `CInvalidDataException` ([`WriteDbError::Options`]) of
    /// `x_ProcessMaskData` and of the format guess, `CWriteDBException`
    /// ([`WriteDbError::Creation`]) of the registry and of the masked ranges.
    Ncbi(WriteDbError),
    /// A failure of the mask file reader that NCBI reports in its own words, or does not
    /// report: [`MaskFileError::Malformed`] (NCBI: a `CSerialException`, status 255, or
    /// `std::bad_alloc`, status 4), [`MaskFileError::Unsupported`],
    /// [`MaskFileError::EofInContinuation`].
    File(MaskFileError),
}

impl From<MaskFileError> for MaskDataError {
    fn from(error: MaskFileError) -> Self {
        match error {
            MaskFileError::Ncbi(error) => MaskDataError::Ncbi(error),
            other => MaskDataError::File(other),
        }
    }
}

impl From<WriteDbError> for MaskDataError {
    fn from(error: WriteDbError) -> Self {
        MaskDataError::Ncbi(error)
    }
}

impl fmt::Display for MaskDataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaskDataError::Ncbi(error) => fmt::Display::fmt(error, f),
            MaskDataError::File(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for MaskDataError {}

impl MaskDataError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
    // ```c++
    // #define CATCH_ALL(exit_code)                                                \
    //     catch (const CInvalidDataException& e) {                                \
    //         LOG_POST(Error << "BLAST options error: " << e.GetMsg());           \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:212-215
    // ```c++
    //     catch (const CWriteDBException& e) {                                    \
    //         LOG_POST(Error << "BLAST Database creation error: " << e.GetMsg()); \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:236-247
    // ```c++
    //     catch (const CIOException& e) {                                         \
    //         if (e.GetErrCode() == CIOException::eFlush) {                       \
    //             LOG_POST(Error << "BLAST failed to write output: " << e.GetMsg());\
    //             exit_code = BLAST_OUTPUT_ERROR;                                 \
    //         }                                                                   \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    /// Where each outcome of `x_ProcessMaskData` becomes makeblastdb's (exit status, message);
    /// the only place that decides it. `None` is a run that ends with status 0 and no message.
    ///
    /// * [`MaskDataError::Ncbi`]: `CATCH_ALL` as NCBI (`BLAST options error: ...`, status 1;
    ///   `BLAST Database creation error: ...`, status 1; the message of a failure of the file
    ///   system, "Error: ...", status 255).
    /// * [`MaskFileError::Malformed`]: a file that NCBI's serial reader refuses. NCBI ends with
    ///   status 255 and "Error: NCBI C++ Exception: ..." (a `CSerialException`, whose text names
    ///   the directory of NCBI's build and is not reproduced); LOSAT prints "Error: mask data in
    ///   '<file>': <its own text>" with status 255. NCBI's `std::bad_alloc` for a string whose
    ///   length cannot be allocated (status 4, "BLAST ran out of memory") depends on the memory
    ///   of the host and is not told apart.
    /// * [`MaskFileError::Unsupported`]: NCBI reads the file and LOSAT does not
    ///   (`MakedbError::Unsupported`, status 1).
    /// * [`MaskFileError::EofInContinuation`]: `None`, exactly as NCBI (status 0, nothing is
    ///   printed here; the destructor's message follows).
    pub fn into_makedb_error(self) -> Option<MakedbError> {
        match self {
            MaskDataError::Ncbi(error) => Some(MakedbError::Build(error)),
            MaskDataError::File(MaskFileError::EofInContinuation { .. }) => None,
            MaskDataError::File(MaskFileError::Unsupported(what)) => Some(
                MakedbError::Unsupported(format!("mask data that LOSAT does not read: {what}")),
            ),
            MaskDataError::File(error @ MaskFileError::Malformed { .. }) => {
                Some(MakedbError::Build(WriteDbError::Other(error.to_string())))
            }
            // The reader turns this into `MaskDataError::Ncbi`.
            MaskDataError::File(MaskFileError::Ncbi(error)) => Some(MakedbError::Build(error)),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3146-3185
// ```c++
// template <typename TStr>
// TStr s_TruncateSpaces(const TStr& str, NStr::ETrunc where,
//                       const TStr& empty_str)
// {
//     ...
//         while ( isspace((unsigned char) str[beg]) ) {
//     ...
//         while (isspace((unsigned char) str[--end])) {
// ```
/// `NStr::TruncateSpaces(str)`: white space (`isspace` of the C locale) at both ends goes.
fn truncate_spaces(text: &str) -> &str {
    text.trim_matches(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r'))
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-275,282-301
// ```c++
//     void Do(TContainer&     target,
//             TPosContainer&  token_pos,
//             const TString&  empty_str = TString())
//     {
//         auto target_initial_size = target.size();
//
//         // Special cases
//         if (m_Str.empty()) {
//             return;
//         } else if (m_Delim.empty()) {
//             target.push_back(m_Str);
//             token_pos.push_back(0);
// ...
//
//         CTempStringList part_collector(m_Storage);
//         SIZE_TYPE prev_pos;
//         SIZE_TYPE delim_pos = NPOS;
//         m_Pos = 0;
//         do {
//             Advance(&part_collector, &prev_pos, &delim_pos);
//             target.push_back(empty_str); // reserve space for value added next line
//             part_collector.Join(&target.back());
//             part_collector.Clear();
//             token_pos.push_back(prev_pos);
//         } while ( !AtEnd() );
//
//         if ( (m_Flags & NStr::fSplit_Truncate_End) == 0 ) {
//             // account training delimiter
//             if (delim_pos != NPOS) {
//                 // add empty token after last delimiter
//                 target.push_back(empty_str);
//                 token_pos.push_back(delim_pos + 1);
//             }
// ```
/// `NStr::Split(text, ",", list)` (no flags): the empty string gives no token, any other string
/// one token for each comma plus one, empty tokens included ("," gives two, "a," gives "a" and
/// ""). The tokens are not trimmed.
fn split_commas(text: &str) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(',').collect()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:897-909
// ```c++
//         // make sure this is not a numeric id
//         for (unsigned int i = 0; i < id_list.size(); ++i) {
//             Int4 nid(-1);
//             if (NStr::StringToNumeric(id_list[i], &nid, NStr::fConvErr_NoThrow, 10)) {
//                 NCBI_THROW(CInvalidDataException, eInvalidInput,
//                    "mask_id can not be numeric.");
//             }
//         }
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:3611-3616
// ```c++
//     template <typename TNumeric>
//     static typename enable_if< is_integral<TNumeric>::value && is_signed<TNumeric>::value && (sizeof(TNumeric) == sizeof(int) && !is_same<TNumeric, long>::value), bool>::type
//     x_StringToNumeric(const CTempString str, TNumeric* value, TStringToNumFlags flags, int base)
//     {
//         *value = StringToInt(str, flags, base);
//         return (*value || !errno);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-873
// ```c++
// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
// {
//     ...
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     ...
//     while (char ch = str[pos]) {
//         ...
//         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
//             break;
//         }
//         // Overflow check
//         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
//             S2N_CONVERT_ERROR_OVERFLOW(Int8);
//         }
//     ...
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     ...
//     CHECK_ENDPTR(Int8);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
// ```c++
// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Int8 value = StringToInt8(str, flags, base);
//     if ( value < kMin_Int  ||  value > kMax_Int ) {
//         S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
//     }
//     return (int) value;
// }
// ```
// With `fConvErr_NoThrow` and no other flag the text converts when it is an optional sign and one
// or more decimal digits and nothing else (no spaces, no commas, no "0x"), and its value is in
// the range of Int4 (a value outside it, or beyond Int8, sets errno and gives 0). A zero that
// leaves errno clear ("0", "-0", "007") is numeric.
/// `NStr::StringToNumeric<Int4>(text, &nid, fConvErr_NoThrow, 10)`: whether `text` converts.
fn is_numeric_int4(text: &str) -> bool {
    let digits = text
        .strip_prefix('-')
        .or_else(|| text.strip_prefix('+'))
        .unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return false;
    }
    let negative = text.starts_with('-');
    // Int4 holds -2^31 ..= 2^31 - 1; every other value, and an overflow of Int8 on the way,
    // is an error as well.
    let mut value: i64 = 0;
    for b in digits.bytes() {
        value = value * 10 + i64::from(b - b'0');
        if value > i64::from(i32::MAX) + 1 {
            return false;
        }
    }
    if negative {
        value = -value;
    }
    (i64::from(i32::MIN)..=i64::from(i32::MAX)).contains(&value)
}

fn options_error(message: &str) -> MaskDataError {
    MaskDataError::Ncbi(WriteDbError::Options(message.to_string()))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:876-1000
// ```c++
// void CMakeBlastDBApp::x_ProcessMaskData()
// {
//     const CArgs & args = GetArgs();
//
//     const CArgValue & files    = args["mask_data"];
//     const CArgValue & ids      = args["mask_id"];
//     const CArgValue & descs    = args["mask_desc"];
//     const CArgValue & gi_names = args["gi_mask_name"];
//
//     vector<string> file_list;
//     vector<string> id_list;
//     vector<string> desc_list;
//     vector<string> gi_mask_names;
//
//     if (! files.HasValue()) return;
//     NStr::Split(NStr::TruncateSpaces(files.AsString()), ",", file_list);
//     if (! file_list.size()) {
//         NCBI_THROW(CInvalidDataException, eInvalidInput,
//                 "mask_data option found, but no files were specified.");
//     }
//
//     if (ids.HasValue()) {
//         NStr::Split(NStr::TruncateSpaces(ids.AsString()), ",", id_list);
//         if (file_list.size() != id_list.size()) {
//             NCBI_THROW(CInvalidDataException, eInvalidInput,
//                 "the size of mask_id does not match that of mask_data.");
//         }
//         // make sure this is not a numeric id
//         ...
//     }
//
//     if (descs.HasValue()) {
//         NStr::Split(NStr::TruncateSpaces(descs.AsString()), ",", desc_list);
//         if (file_list.size() != desc_list.size()) {
//             NCBI_THROW(CInvalidDataException, eInvalidInput,
//                 "the size of mask_desc does not match that of mask_data.");
//         }
//     } else {
//         // description is optional
//         vector<string> default_desc(id_list.size(), "");
//         desc_list.swap(default_desc);
//     }
//
//     if (gi_names.HasValue()) {
//         NStr::Split(NStr::TruncateSpaces(gi_names.AsString()), ",", gi_mask_names);
//         if (file_list.size() != gi_mask_names.size()) {
//             NCBI_THROW(CInvalidDataException, eInvalidInput,
//                 "the size of gi_mask_name does not match that of mask_data.");
//         }
//     }
//
//     for (unsigned int i = 0; i < file_list.size(); ++i) {
//         if ( !CFile(file_list[i]).Exists() ) {
//             ERR_POST(Error << "Ignoring mask file '" << file_list[i]
//                            << "' as it does not exist.");
//             continue;
//         }
//
//         CNcbiIfstream mask_file(file_list[i].c_str(), ios::binary);
//         CFormatGuess::EFormat mask_file_format = CFormatGuess::eUnknown;
//         {{
//              CFormatGuess fg(mask_file);
//              ...
//              mask_file_format = fg.GuessFormat();
//          }}
//
//         int algo_id = -1;
//         while (true) {
//             CRef<CBlast_db_mask_info> first_obj;
//
//             try {
//                 s_ReadObject(mask_file, mask_file_format, first_obj, "mask data in '" + file_list[i] + "'");
//             }
//             catch (CEofException&) {
//                 // must be end of file
//                 break;
//             }
//
//             if (algo_id < 0) {
//                 *m_LogFile << "Mask file: " << file_list[i] << endl;
//                 string opts = first_obj->GetAlgo_options();
//                 if (id_list.size())  {
//                     algo_id = m_DB->RegisterMaskingAlgorithm(id_list[i], desc_list[i], opts);
//                 } else {
//                     EBlast_filter_program prog_id =
//                         static_cast<EBlast_filter_program>(first_obj->GetAlgo_program());
//                     string name = gi_mask_names.size() ? gi_mask_names[i] : file_list[i];
//                     algo_id = m_DB->RegisterMaskingAlgorithm(prog_id, opts, name);
//                 }
//             }
//
//             CRef<CBlast_mask_list> masks(& first_obj->SetMasks());
//             first_obj.Reset();
//
//             while(1) {
//                 if (m_Ranges.Empty() && ! masks->GetMasks().empty()) {
//                     m_Ranges.Reset(new CMaskedRangeSet);
//                     m_DB->SetMaskDataSource(*m_Ranges);
//                 }
//
//                 ITERATE(CBlast_mask_list::TMasks, iter, masks->GetMasks()) {
//                     CConstRef<CSeq_id> seqid((**iter).GetId());
//
//                     if (seqid.Empty()) {
//                         NCBI_THROW(CInvalidDataException, eInvalidInput,
//                                      "Cannot get masked range Seq-id");
//                     }
//
//                     m_Ranges->Insert(algo_id, *seqid, **iter);
//                 }
//
//                 if (! masks->GetMore())
//                     break;
//
//                 s_ReadObject(mask_file, mask_file_format, masks, "mask data (continuation)");
//             }
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp (include/corelib/ncbifile.hpp:4039-4042)
// ```c++
// inline
// bool CFile::Exists(void) const
// {
//     return IsFile();
// }
// ```
// `CFile::Exists` is true for a regular file (a link is followed); a directory, a device and a
// FIFO "do not exist" (the oracle prints "Ignoring mask file '/dev/null' as it does not
// exist."). A file that exists and cannot be opened is a failed stream, which the format guess
// finds neither binary nor text: "Unknown encoding for mask data in '<file>'" (the oracle, a
// file with mode 000).
/// `x_ProcessMaskData`. `db` is the `CBuildDatabase` (`m_DB`; the log `m_LogFile` is its log) and
/// `warn` takes the lines that NCBI posts with `ERR_POST` (without the final newline).
///
/// Returns at once, without a trace, when `-mask_data` was not given. On an error nothing is
/// undone (NCBI does not either): what was registered and read stays in `db`.
pub fn process_mask_data(
    db: &mut BuildDatabase<'_>,
    options: &MaskDataOptions<'_>,
    warn: Warn<'_>,
) -> Result<(), MaskDataError> {
    let Some(files) = options.mask_data else {
        return Ok(());
    };
    let file_list = split_commas(truncate_spaces(files));
    if file_list.is_empty() {
        return Err(options_error(
            "mask_data option found, but no files were specified.",
        ));
    }

    let mut id_list: Vec<&str> = Vec::new();
    let mut desc_list: Vec<&str> = Vec::new();
    let mut gi_mask_names: Vec<&str> = Vec::new();

    if let Some(ids) = options.mask_id {
        id_list = split_commas(truncate_spaces(ids));
        if file_list.len() != id_list.len() {
            return Err(options_error(
                "the size of mask_id does not match that of mask_data.",
            ));
        }
        // make sure this is not a numeric id
        if id_list.iter().any(|id| is_numeric_int4(id)) {
            return Err(options_error("mask_id can not be numeric."));
        }
    }

    if let Some(descs) = options.mask_desc {
        desc_list = split_commas(truncate_spaces(descs));
        if file_list.len() != desc_list.len() {
            return Err(options_error(
                "the size of mask_desc does not match that of mask_data.",
            ));
        }
    } else {
        // description is optional
        desc_list = vec![""; id_list.len()];
    }

    if let Some(names) = options.gi_mask_name {
        gi_mask_names = split_commas(truncate_spaces(names));
        if file_list.len() != gi_mask_names.len() {
            return Err(options_error(
                "the size of gi_mask_name does not match that of mask_data.",
            ));
        }
    }

    for (i, &file) in file_list.iter().enumerate() {
        let path: PathBuf = match options.root {
            Some(root) => root.join(file),
            None => PathBuf::from(file),
        };
        if !std::fs::metadata(&path).is_ok_and(|meta| meta.is_file()) {
            warn(
                format!("Error: [makeblastdb] Ignoring mask file '{file}' as it does not exist.")
                    .as_bytes(),
            );
            continue;
        }

        // The name in the messages is the one given, not the path it was looked up under.
        let mut reader = open_mask_file(&path, Path::new(file))?;

        let mut algo_id: Option<i32> = None;
        while let Some(first_obj) = reader.next_info()? {
            let algo_id = match algo_id {
                Some(id) => id,
                None => {
                    let _ = writeln!(db.log(), "Mask file: {file}");
                    let opts = &first_obj.algo_options;
                    let id = if !id_list.is_empty() {
                        db.register_masking_algorithm_id(
                            id_list[i].as_bytes(),
                            desc_list[i].as_bytes(),
                            opts,
                        )?
                    } else {
                        let name = if gi_mask_names.is_empty() {
                            file
                        } else {
                            gi_mask_names[i]
                        };
                        db.register_masking_algorithm(
                            first_obj.algo_program,
                            opts,
                            name.as_bytes(),
                        )?
                    };
                    algo_id = Some(id);
                    id
                }
            };

            let mut masks = first_obj.masks;
            loop {
                if db.mask_data_source_mut().is_none() && !masks.masks.is_empty() {
                    db.set_mask_data_source(MaskedRangeSet::new());
                }

                for mask in &masks.masks {
                    let Some(seq_id) = mask.id.as_ref() else {
                        return Err(options_error("Cannot get masked range Seq-id"));
                    };
                    if let Some(set) = db.mask_data_source_mut() {
                        set.insert(algo_id, seq_id, mask)?;
                    }
                }

                if !masks.more {
                    break;
                }

                masks = reader.next_list()?;
            }
        }
    }
    Ok(())
}

/// The reader of the mask file at `path`, whose name in messages is `name`. A file that cannot
/// be opened or read is "Unknown encoding" as in NCBI (`MaskFileReader::open_named`).
fn open_mask_file(path: &Path, name: &Path) -> Result<MaskFileReader, MaskDataError> {
    Ok(MaskFileReader::open_named(path, name)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spaces_and_commas_are_split_as_ncbi_does() {
        assert_eq!(truncate_spaces(" \t a b \n"), "a b");
        assert_eq!(truncate_spaces(" \u{a0}"), "\u{a0}");
        assert_eq!(truncate_spaces("  "), "");
        assert!(split_commas("").is_empty());
        assert_eq!(split_commas(","), ["", ""]);
        assert_eq!(split_commas("a,"), ["a", ""]);
        assert_eq!(split_commas("a,,b"), ["a", "", "b"]);
        assert_eq!(split_commas("a, b"), ["a", " b"]);
    }

    /// The answers of makeblastdb 2.17.0 (`-mask_id` with one or two mask files).
    #[test]
    fn numeric_ids_are_those_that_stringtonumeric_converts() {
        for numeric in [
            "0",
            "5",
            "007",
            "+5",
            "-0",
            "-5",
            "2147483647",
            "-2147483648",
        ] {
            assert!(is_numeric_int4(numeric), "{numeric}");
        }
        for not_numeric in [
            "",
            "+",
            "-",
            " 5",
            "5 ",
            "5a",
            "a5",
            "1e3",
            "0x10",
            "٣",
            "+-5",
            "2147483648",
            "-2147483649",
            "99999999999999999999",
            "9223372036854775808",
        ] {
            assert!(!is_numeric_int4(not_numeric), "{not_numeric}");
        }
    }
}

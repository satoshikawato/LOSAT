//! The DELTA-BLAST files and checks of makeprofiledb (`-dbtype delta`): `x_CheckDelta`,
//! `x_ValidateCd`, `x_UpdateDelta` and `x_WrapUpDelta`, writing the weighted residue frequencies
//! (`.wcounts`) and run-length compressed effective observations (`.obsr`) for each volume.
//!
//! Each profile's `.wcounts` and `.obsr` records are streamed to unnamed temporary scratch files
//! in the temporary directory during the volume pass (matching NCBI's `CTmpFile`s while keeping
//! memory bounded on full CDD volumes, falling back to memory on WASI where no temporary
//! directory exists). The final `<name>.wcounts` and `<name>.obsr` files are created only at
//! `finish` (`x_WrapUpDelta`), after Pass 2 completes and Pass-3 ASN.1 posts / errors are
//! replayed, so an input profile sharing an output name is not truncated early and a `.wcounts`
//! open error preserves the protein database files as in NCBI.

use std::borrow::Cow;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};

use crate::blastinput::input_files::{temporary_directory, unnamed_temporary_file};
use crate::core::blast_util::BLASTAA_SIZE;
use crate::ncbi_app::Diag;
use crate::report::serial::asn_types::Post;
use crate::report::serial::pssm::Pssm;
use crate::seqdb::rps::RPS_MAGIC_NUM_28;

use super::rps_files::RpsError;

/// `RPS_MAGIC_NUM_28`, written as the header magic of `.wcounts` and `.obsr`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1841-1844
const DELTA_MAGIC: i32 = RPS_MAGIC_NUM_28;

/// `kFixedPointScaleFactor`: multiplier applied before converting delta floats to `Uint4`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:112-113
const FIXED_POINT_SCALE_FACTOR: f64 = 1000.0;

/// `kEpsylon`: allowed distance of a column's weighted-frequency sum from 1.0 in `x_ValidateCd`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:114
const EPSILON: f64 = 0.0001;

/// Run-length merge tolerance for consecutive `numIndeptObsr` values in `x_UpdateDelta`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1665
const OBSR_MERGE_EPS: f64 = 1e-4;

/// Buffer capacity for the temporary scratch files and final `.wcounts`/`.obsr` outputs.
const DELTA_IO_BUFFER_BYTES: usize = 1 << 16;

/// Result of `x_CheckDelta` for one profile in `x_CreateDeltaList`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DeltaCheck {
    /// Profile passes the observation threshold and (unless `-exclude_invalid false`) validation.
    Included,
    /// Maximum independent observations plus one fell below `-obsr_threshold`.
    ExcludedObsr,
    /// Profile failed `x_ValidateCd` and `-exclude_invalid` is true.
    ExcludedInvalidCd,
}

/// Errors from delta checking or writing, mapped to makeprofiledb's `RunError`.
#[derive(Debug)]
pub(super) enum DeltaError {
    /// `CInputException`: "INPUT ERROR: ", exit status 1.
    Input(Vec<u8>),
    /// `CException`: "ERROR: ", exit status 255.
    Exception(Vec<u8>),
    /// `CSeqDBException` / I/O / refusal shared with `RpsError`.
    Rps(RpsError),
}

impl From<RpsError> for DeltaError {
    fn from(error: RpsError) -> Self {
        DeltaError::Rps(error)
    }
}

impl From<std::io::Error> for DeltaError {
    fn from(error: std::io::Error) -> Self {
        DeltaError::Rps(RpsError::Io(error))
    }
}

/// Builds a rule-3 refusal error for `-dbtype delta`.
fn refused(reason: impl Into<String>) -> DeltaError {
    DeltaError::Rps(RpsError::Refused(reason.into()))
}

/// A `double` converted to `Uint4` (`uint32_t`): `None` where C++'s floating-to-unsigned
/// conversion is undefined (NaN, infinity, or a truncated value outside `0..=u32::MAX`).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1674,1688
pub(super) fn uint4_of(value: f64) -> Option<u32> {
    (value > -1.0 && value < 4_294_967_296.0).then_some(value as u32)
}

/// Converts a scaled delta float to `Uint4` via [`uint4_of`], returning a rule-3 refusal when out
/// of C++'s defined range.
fn encode_uint4(scaled: f64, field: &str, expr: impl std::fmt::Display) -> Result<u32, DeltaError> {
    uint4_of(scaled).ok_or_else(|| {
        refused(format!(
            "a {field} for -dbtype delta (NCBI's conversion of {expr} to a 32-bit unsigned \
             integer has no defined result)"
        ))
    })
}

/// Advances a 32-bit delta offset counter (`CurrFreqOffset` or `CurrObsrOffset`) by `added`
/// entries, adding the buffer length as unsigned and converting back to two's-complement `Int4`
/// (wrapping modulo 2^32).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1702-1703
pub(super) fn advance_offset(current: i32, added: usize) -> i32 {
    current.wrapping_add(added as i32)
}

/// Converts `CurrFreqOffset` (an `Int4` entry offset) to the column offset stored in `.wcounts`
/// by dividing by `BLASTAA_SIZE` with signed truncation toward zero.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1682,1825
pub(super) fn freq_record_offset(curr_freq_offset: i32) -> i32 {
    curr_freq_offset / BLASTAA_SIZE as i32
}

/// Extracts and validates the `weightedResFreqsPerPos` and `numIndeptObsr` layout for
/// `x_CheckDelta` and `x_UpdateDelta`, refusing the shapes where NCBI copies into an empty
/// `std::list::begin()` or dereferences an empty `numIndeptObsr` list.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1549-1560,1602-1649,1755-1803
fn delta_fields<'a>(
    pssm: &'a Pssm,
    seq_size: i32,
    filename: &[u8],
) -> Result<(Cow<'a, [f64]>, &'a [f64]), DeltaError> {
    let input = |what: &[u8]| DeltaError::Input([filename, what].concat());
    let inter = pssm.intermediate_data.as_ref();
    let Some(orig_freqs) = inter.and_then(|d| d.weighted_res_freqs_per_pos.as_deref()) else {
        return Err(input(
            b" contains no weighted residue frequencies for building delta database",
        ));
    };
    let Some(obsr) = inter.and_then(|d| d.num_indept_obsr.as_deref()) else {
        return Err(input(
            b" contains no observations information for building delta database",
        ));
    };
    let alphabet_size = pssm.num_rows as usize;
    let by_row = pssm.by_row.unwrap_or(false);
    // NCBI's delta check and update transpose weighted frequencies into an empty list's begin
    // iterator (and index over `seq_size`), overwriting the list header node or reading past
    // `orig_freqs`.
    if by_row && (seq_size > 0 || !orig_freqs.is_empty()) {
        return Err(refused(
            "a profile stored by row (byRow TRUE) for -dbtype delta (NCBI's makeprofiledb copies \
             the transposed weighted residue frequencies into an empty list's begin iterator)",
        ));
    }
    // NCBI's padding branch also copies `orig_freqs` into an empty list's begin iterator whenever
    // `orig_freqs` is non-empty.
    let freqs = if alphabet_size < BLASTAA_SIZE {
        if !orig_freqs.is_empty() {
            return Err(refused(format!(
                "a profile with numRows {alphabet_size} (fewer than {BLASTAA_SIZE}) and weighted \
                 residue frequencies for -dbtype delta (NCBI's makeprofiledb copies them into an \
                 empty list's begin iterator)"
            )));
        }
        let pad = (BLASTAA_SIZE - alphabet_size).saturating_mul(seq_size.max(0) as usize);
        Cow::Owned(vec![0.0; pad])
    } else {
        Cow::Borrowed(orig_freqs)
    };
    if obsr.is_empty() {
        return Err(refused(
            "a profile with an empty numIndeptObsr list for -dbtype delta (NCBI's makeprofiledb \
             dereferences the end iterator of the empty list)",
        ));
    }
    Ok((freqs, obsr))
}

/// `x_ValidateCd`: checks column count agreement (`freqs.len() / 28 == obsr.len()`), each full
/// 28-residue column's sum within `kEpsylon` of `1.0`, and every observation `>= 1.0`. Refuses a
/// partial trailing column (`len % 28 != 0`) only when every full column sums within `kEpsylon`
/// of `1.0`, which is the exact condition under which NCBI's loop reaches the partial tail and
/// increments past `freqs.end()`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1715-1749
fn validate_cd(freqs: &[f64], obsr: &[f64]) -> Result<bool, DeltaError> {
    if freqs.len() / BLASTAA_SIZE != obsr.len() {
        return Err(DeltaError::Exception(
            b"Number of frequency and observations columns do not match".to_vec(),
        ));
    }
    for col in freqs.chunks_exact(BLASTAA_SIZE) {
        let mut sum = 0.0f64;
        for &val in col {
            sum += val;
        }
        if (sum - 1.0).abs() > EPSILON {
            return Ok(false);
        }
    }
    if !freqs.len().is_multiple_of(BLASTAA_SIZE) {
        return Err(refused(format!(
            "a profile whose weightedResFreqsPerPos has {} entries (not a multiple of \
             {BLASTAA_SIZE}) and every full column sums within {EPSILON} of 1.0 for -dbtype \
             delta (NCBI's makeprofiledb increments past the end of the list)",
            freqs.len()
        )));
    }
    if obsr.iter().any(|&v| v < 1.0) {
        return Ok(false);
    }
    Ok(true)
}

/// `x_CheckDelta` (plus `x_CreateDeltaList`'s field presence checks): returns whether the profile
/// is included, excluded for too few independent observations, or excluded as an invalid CD.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1549-1565,1752-1818
pub(super) fn check_delta(
    pssm: &Pssm,
    seq_size: i32,
    filename: &[u8],
    obsr_threshold: f64,
    exclude_invalid: bool,
) -> Result<DeltaCheck, DeltaError> {
    let (freqs, obsr) = delta_fields(pssm, seq_size, filename)?;
    // Matches C++ max-element scan: keeps the first element unless a later element is strictly
    // greater.
    let mut max_val = obsr[0];
    for &v in &obsr[1..] {
        if max_val < v {
            max_val = v;
        }
    }
    let max_obsr = max_val + 1.0;
    if max_obsr < obsr_threshold {
        return Ok(DeltaCheck::ExcludedObsr);
    }
    // `x_ValidateCd` runs even when `exclude_invalid` is false because NCBI evaluates it on the
    // left-hand side of the exclusion check, so column-count mismatch still throws `CException`.
    if !validate_cd(&freqs, obsr)? && exclude_invalid {
        return Ok(DeltaCheck::ExcludedInvalidCd);
    }
    Ok(DeltaCheck::Included)
}

/// Scratch storage for `.wcounts` or `.obsr` payload bytes until `finish` (`x_WrapUpDelta`): an
/// unnamed temporary file in `TMPDIR` on native targets (as NCBI's `CTmpFile`), or an in-memory
/// buffer on WASI where no temporary directory exists.
enum DeltaScratch {
    File(BufWriter<File>),
    Mem(Vec<u8>),
}

impl DeltaScratch {
    fn file(file: File) -> Self {
        DeltaScratch::File(BufWriter::with_capacity(DELTA_IO_BUFFER_BYTES, file))
    }

    fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        match self {
            DeltaScratch::File(writer) => writer.write_all(bytes),
            DeltaScratch::Mem(buf) => {
                buf.extend_from_slice(bytes);
                Ok(())
            }
        }
    }

    fn copy_to(&mut self, out: &mut BufWriter<File>) -> std::io::Result<()> {
        match self {
            DeltaScratch::File(writer) => {
                writer.flush()?;
                let file = writer.get_mut();
                file.seek(SeekFrom::Start(0))?;
                std::io::copy(file, out)?;
                Ok(())
            }
            DeltaScratch::Mem(buf) => out.write_all(buf),
        }
    }
}

/// The buffered `.wcounts` and `.obsr` state of one delta volume (`x_UpdateDelta` /
/// `x_WrapUpDelta`), held in temporary scratch files until `finish` creates the output files.
pub(super) struct DeltaVolume {
    db_name: Vec<u8>,
    tmp_freq: Option<DeltaScratch>,
    tmp_obsr: Option<DeltaScratch>,
    /// Failure creating the temporary scratch files in `TMPDIR`, held until the start of `finish`
    /// so any Pass-2 error takes precedence, matching NCBI's `CTmpFile` creation at the start of
    /// `x_UpdateDelta`.
    tmp_open_error: Option<DeltaError>,
    freq_offsets: Vec<i32>,
    obsr_offsets: Vec<i32>,
    curr_freq_offset: i32,
    curr_obsr_offset: i32,
    scratch: Vec<u8>,
    /// ASN.1 reader posts captured during the volume pass, replayed at `finish` to match NCBI's
    /// third read pass inside `x_UpdateDelta`.
    deferred_posts: Vec<(usize, Vec<Post>)>,
    /// First Pass-3 error (e.g. float-to-`Uint4` UB refusal or scratch write I/O error), held
    /// until `finish` so any earlier Pass-2 error is reported first as in NCBI.
    deferred_error: Option<(usize, DeltaError)>,
}

/// Prepares the scratch buffers for one delta volume (`x_UpdateDelta`). Does not create
/// `<name>.wcounts` or `<name>.obsr` yet; those are opened in `finish` (`x_WrapUpDelta`).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1570-1580
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3660-3705,5389-5396,6561-6637
pub(super) fn new_volume(name: &[u8], num_seqs: usize) -> DeltaVolume {
    let (tmp_obsr, tmp_freq, tmp_open_error) = match temporary_directory() {
        Some(dir) => match unnamed_temporary_file(&dir, "delta-obsr")
            .and_then(|obsr| Ok((obsr, unnamed_temporary_file(&dir, "delta-wcounts")?)))
        {
            Ok((obsr, freq)) => (
                Some(DeltaScratch::file(obsr)),
                Some(DeltaScratch::file(freq)),
                None,
            ),
            Err(error) => (
                None,
                None,
                Some(DeltaError::Rps(RpsError::Refused(format!(
                    "a temporary directory ('{}') where temporary files for -dbtype delta cannot \
                     be created ({error}; NCBI's makeprofiledb stops with a message containing \
                     its process ID)",
                    dir.display()
                )))),
            ),
        },
        None => {
            let mem = || Some(DeltaScratch::Mem(Vec::new()));
            (mem(), mem(), None)
        }
    };
    DeltaVolume {
        db_name: name.to_vec(),
        tmp_freq,
        tmp_obsr,
        tmp_open_error,
        freq_offsets: Vec::with_capacity(num_seqs + 1),
        obsr_offsets: Vec::with_capacity(num_seqs + 1),
        curr_freq_offset: 0,
        curr_obsr_offset: 0,
        scratch: Vec::new(),
        deferred_posts: Vec::new(),
        deferred_error: None,
    }
}

/// Run-length encodes `obsr` into `(value, count)` `u32` pairs as `x_UpdateDelta` writes them to
/// `.obsr`, returning `Err` if any `(current + 1.0) * 1000` falls outside C++'s defined `Uint4`
/// range. Run counts wrap modulo 2^32 as NCBI's unsigned 32-bit run counter does.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1654-1677
fn encode_obsr(obsr: &[f64]) -> Result<Vec<u32>, DeltaError> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < obsr.len() {
        let current = obsr[i];
        let mut num = 1u32;
        i += 1;
        while i < obsr.len() && (obsr[i] - current).abs() < OBSR_MERGE_EPS {
            i += 1;
            num = num.wrapping_add(1);
        }
        let scaled = (current + 1.0) * FIXED_POINT_SCALE_FACTOR;
        let encoded = encode_uint4(
            scaled,
            &format!("numIndeptObsr value of {current}"),
            format_args!("({current} + 1) * 1000"),
        )?;
        out.push(encoded);
        out.push(num);
    }
    Ok(out)
}

/// Writes the header (magic, record count, offset table) and copies `scratch` into `out`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1841-1863
fn write_delta_file(
    out: &mut BufWriter<File>,
    offsets: &[i32],
    scratch: &mut DeltaScratch,
) -> std::io::Result<()> {
    let num_records = (offsets.len() - 1) as i32;
    out.write_all(&DELTA_MAGIC.to_le_bytes())?;
    out.write_all(&num_records.to_le_bytes())?;
    for offset in offsets {
        out.write_all(&offset.to_le_bytes())?;
    }
    out.flush()?;
    scratch.copy_to(out)?;
    out.flush()
}

impl DeltaVolume {
    /// Records the current `.wcounts` column offset and `.obsr` word offset in the offset tables.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1652,1682,1825-1826
    fn push_current_offsets(&mut self) {
        self.freq_offsets
            .push(freq_record_offset(self.curr_freq_offset));
        self.obsr_offsets.push(self.curr_obsr_offset);
    }

    /// Encodes one profile's `.wcounts` and `.obsr` records into the scratch buffers.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1598-1707
    fn encode_and_write_profile(
        &mut self,
        pssm: &Pssm,
        seq_size: i32,
        filename: &[u8],
    ) -> Result<(), DeltaError> {
        let (weighted_freqs, obsr) = delta_fields(pssm, seq_size, filename)?;
        let obsr_words = encode_obsr(obsr)?;
        if obsr.len() != weighted_freqs.len() / BLASTAA_SIZE {
            return Err(DeltaError::Exception(
                [
                    b"Number of frequencies and observations columns do not match in ".as_slice(),
                    filename,
                ]
                .concat(),
            ));
        }
        self.scratch.clear();
        self.scratch
            .reserve((weighted_freqs.len() + BLASTAA_SIZE) * 4);
        for &val in weighted_freqs.as_ref() {
            let scaled = val * FIXED_POINT_SCALE_FACTOR;
            let encoded = encode_uint4(
                scaled,
                &format!("weightedResFreqsPerPos value of {val}"),
                format_args!("{val} * 1000"),
            )?;
            self.scratch.extend_from_slice(&encoded.to_le_bytes());
        }
        // Additional column of zeros for compatibility with RPS database.
        self.scratch.extend_from_slice(&[0u8; BLASTAA_SIZE * 4]);

        self.tmp_freq
            .as_mut()
            .expect("tmp_freq open")
            .write_all(&self.scratch)?;
        self.scratch.clear();
        self.scratch.reserve(obsr_words.len() * 4);
        for word in &obsr_words {
            self.scratch.extend_from_slice(&word.to_le_bytes());
        }
        self.tmp_obsr
            .as_mut()
            .expect("tmp_obsr open")
            .write_all(&self.scratch)?;

        self.push_current_offsets();
        self.curr_freq_offset =
            advance_offset(self.curr_freq_offset, weighted_freqs.len() + BLASTAA_SIZE);
        self.curr_obsr_offset = advance_offset(self.curr_obsr_offset, obsr_words.len());
        Ok(())
    }

    /// Appends one profile's `.wcounts` and `.obsr` records to the scratch buffers during the
    /// volume pass, deferring any Pass-3 ASN.1 posts or errors to `finish`.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1582-1707
    pub(super) fn add_profile(
        &mut self,
        profile_index: usize,
        pssm: &Pssm,
        seq_size: i32,
        filename: &[u8],
        posts: Vec<Post>,
    ) {
        if !posts.is_empty() {
            self.deferred_posts.push((profile_index, posts));
        }
        if self.tmp_open_error.is_some() || self.deferred_error.is_some() {
            return;
        }
        if let Err(error) = self.encode_and_write_profile(pssm, seq_size, filename) {
            self.deferred_error = Some((profile_index, error));
        }
    }

    /// Completes `x_UpdateDelta` / `x_WrapUpDelta`: checks temporary-file creation, replays Pass-3
    /// ASN.1 posts, raises any deferred Pass-3 error, then creates `<name>.wcounts` and
    /// `<name>.obsr` (in NCBI's order) and writes their headers, offset tables, and records.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1572-1575,1709-1712,1822-1863
    pub(super) fn finish(
        mut self,
        diag: &mut Diag,
        post_fn: impl Fn(&mut Diag, &[Post]),
    ) -> Result<(), DeltaError> {
        if let Some(error) = self.tmp_open_error.take() {
            return Err(error);
        }
        if let Some((fail_index, error)) = self.deferred_error.take() {
            for (idx, posts) in &self.deferred_posts {
                if *idx <= fail_index {
                    post_fn(diag, posts);
                }
            }
            return Err(error);
        }
        for (_, posts) in &self.deferred_posts {
            post_fn(diag, posts);
        }
        self.push_current_offsets();

        let open = |extension: &str, message: &'static str| {
            let path = [&self.db_name[..], b".", extension.as_bytes()].concat();
            File::create(crate::utils::os_bytes::path_of_bytes(&path))
                .map_err(|_| DeltaError::Rps(RpsError::SeqDb(message)))
        };
        let mut wcounts = BufWriter::with_capacity(
            DELTA_IO_BUFFER_BYTES,
            open("wcounts", "Failed to open output .wcounts file")?,
        );
        let mut obsr = BufWriter::with_capacity(
            DELTA_IO_BUFFER_BYTES,
            open("obsr", "Failed to open output .obsr file")?,
        );
        write_delta_file(
            &mut wcounts,
            &self.freq_offsets,
            self.tmp_freq.as_mut().expect("tmp_freq open"),
        )?;
        write_delta_file(
            &mut obsr,
            &self.obsr_offsets,
            self.tmp_obsr.as_mut().expect("tmp_obsr open"),
        )?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "delta_tests.rs"]
mod tests;

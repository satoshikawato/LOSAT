//! The DELTA-BLAST files and checks of makeprofiledb (`-dbtype delta`): `x_CheckDelta`,
//! `x_ValidateCd`, `x_UpdateDelta` and `x_WrapUpDelta`, writing the weighted residue frequencies
//! (`.wcounts`) and run-length compressed effective observations (`.obsr`) for each volume.
//!
//! Each profile's `.wcounts` and `.obsr` records are streamed to disk during the volume pass so
//! memory stays bounded on full CDD volumes, and the offset tables at byte 8 of both files are
//! written at `finish` (on any error makeprofiledb's destructor removes both files). Any ASN.1
//! reader posts from the volume pass are replayed at `finish` to match NCBI's third read pass in
//! `x_UpdateDelta`.

use std::borrow::Cow;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};

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

/// A `double` converted to `Uint4` (`uint32_t`): `None` where C++'s floating-to-unsigned
/// conversion is undefined (NaN, infinity, or a truncated value outside `0..=u32::MAX`).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1674,1688
pub(super) fn uint4_of(value: f64) -> Option<u32> {
    (value > -1.0 && value < 4_294_967_296.0).then_some(value as u32)
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
    let Some(inter) = pssm.intermediate_data.as_ref() else {
        return Err(input(
            b" contains no weighted residue frequencies for building delta database",
        ));
    };
    let Some(orig_freqs) = inter.weighted_res_freqs_per_pos.as_deref() else {
        return Err(input(
            b" contains no weighted residue frequencies for building delta database",
        ));
    };
    let Some(obsr) = inter.num_indept_obsr.as_deref() else {
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
        return Err(DeltaError::Rps(RpsError::Refused(
            "a profile stored by row (byRow TRUE) for -dbtype delta (NCBI's makeprofiledb copies \
             the transposed weighted residue frequencies into an empty list's begin iterator)"
                .to_string(),
        )));
    }
    // NCBI's padding branch also copies `orig_freqs` into an empty list's begin iterator whenever
    // `orig_freqs` is non-empty.
    let freqs = if alphabet_size < BLASTAA_SIZE {
        if !orig_freqs.is_empty() {
            return Err(DeltaError::Rps(RpsError::Refused(format!(
                "a profile with numRows {alphabet_size} (fewer than {BLASTAA_SIZE}) and weighted \
                 residue frequencies for -dbtype delta (NCBI's makeprofiledb copies them into an \
                 empty list's begin iterator)"
            ))));
        }
        let pad = (BLASTAA_SIZE - alphabet_size).saturating_mul(seq_size.max(0) as usize);
        Cow::Owned(vec![0.0; pad])
    } else {
        Cow::Borrowed(orig_freqs)
    };
    if obsr.is_empty() {
        return Err(DeltaError::Rps(RpsError::Refused(
            "a profile with an empty numIndeptObsr list for -dbtype delta (NCBI's makeprofiledb \
             dereferences the end iterator of the empty list)"
                .to_string(),
        )));
    }
    Ok((freqs, obsr))
}

/// `x_ValidateCd`: checks column count agreement (`freqs.len() / 28 == obsr.len()`), each
/// 28-residue column's sum within `kEpsylon` of `1.0`, and every observation `>= 1.0`. Refuses a
/// partial trailing column whose length is not a multiple of 28, where NCBI increments past
/// `freqs.end()`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1715-1749
fn validate_cd(freqs: &[f64], obsr: &[f64]) -> Result<bool, DeltaError> {
    if freqs.len() / BLASTAA_SIZE != obsr.len() {
        return Err(DeltaError::Exception(
            b"Number of frequency and observations columns do not match".to_vec(),
        ));
    }
    if !freqs.len().is_multiple_of(BLASTAA_SIZE) {
        return Err(DeltaError::Rps(RpsError::Refused(format!(
            "a profile whose weightedResFreqsPerPos has {} entries (not a multiple of \
             {BLASTAA_SIZE}) for -dbtype delta (NCBI's makeprofiledb increments past the end of \
             the list)",
            freqs.len()
        ))));
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

/// The open `.wcounts` and `.obsr` files of one delta volume (`x_UpdateDelta` / `x_WrapUpDelta`).
pub(super) struct DeltaVolume {
    wcounts: BufWriter<File>,
    obsr: BufWriter<File>,
    freq_offsets: Vec<i32>,
    obsr_offsets: Vec<i32>,
    curr_freq_offset: i32,
    curr_obsr_offset: i32,
    scratch: Vec<u8>,
    /// ASN.1 reader posts captured during the volume pass, replayed at `finish` to match NCBI's
    /// third read pass inside `x_UpdateDelta`.
    deferred_posts: Vec<(usize, Vec<Post>)>,
    /// First Pass-3 error (e.g. float-to-`Uint4` UB refusal or 32-bit offset overflow), held until
    /// `finish` so any earlier Pass-2 error is reported first as in NCBI.
    deferred_error: Option<(usize, DeltaError)>,
}

/// Opens `<name>.wcounts` and `<name>.obsr` and writes their magic, record count, and placeholder
/// offset tables (replaced at `finish`).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1828-1858
pub(super) fn open_volume(name: &[u8], num_seqs: usize) -> Result<DeltaVolume, DeltaError> {
    let open = |extension: &str, message: &'static str| {
        let path = [name, b".", extension.as_bytes()].concat();
        File::create(crate::utils::os_bytes::path_of_bytes(&path))
            .map_err(|_| DeltaError::Rps(RpsError::SeqDb(message)))
    };
    let mut wcounts = BufWriter::new(open("wcounts", "Failed to open output .wcounts file")?);
    let mut obsr = BufWriter::new(open("obsr", "Failed to open output .obsr file")?);
    let count = i32::try_from(num_seqs)
        .ok()
        .filter(|&count| count < i32::MAX)
        .ok_or_else(|| {
            DeltaError::Rps(RpsError::Refused(format!(
                "a delta volume of {num_seqs} profiles (NCBI's makeprofiledb counts its offsets \
                 past its 32-bit integers)"
            )))
        })?;
    let mut header = Vec::with_capacity(8 + 4 * (num_seqs + 1));
    header.extend_from_slice(&DELTA_MAGIC.to_le_bytes());
    header.extend_from_slice(&count.to_le_bytes());
    for i in 0..=count {
        header.extend_from_slice(&i.to_le_bytes());
    }
    wcounts.write_all(&header)?;
    obsr.write_all(&header)?;
    Ok(DeltaVolume {
        wcounts,
        obsr,
        freq_offsets: Vec::with_capacity(num_seqs + 1),
        obsr_offsets: Vec::with_capacity(num_seqs + 1),
        curr_freq_offset: 0,
        curr_obsr_offset: 0,
        scratch: Vec::new(),
        deferred_posts: Vec::new(),
        deferred_error: None,
    })
}

/// Run-length encodes `obsr` into `(value, count)` `u32` pairs as `x_UpdateDelta` writes them to
/// `.obsr`, returning `Err` if any `(current + 1.0) * 1000` falls outside C++'s defined `Uint4`
/// range.
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
            num = num.checked_add(1).ok_or_else(|| {
                DeltaError::Rps(RpsError::Refused(
                    "a profile with 2^32 or more consecutive equal observations for -dbtype delta"
                        .to_string(),
                ))
            })?;
        }
        let scaled = (current + 1.0) * FIXED_POINT_SCALE_FACTOR;
        let Some(encoded) = uint4_of(scaled) else {
            return Err(DeltaError::Rps(RpsError::Refused(format!(
                "a numIndeptObsr value of {current} for -dbtype delta (NCBI's conversion of \
                 ({current} + 1) * 1000 to a 32-bit unsigned integer has no defined result)"
            ))));
        };
        out.push(encoded);
        out.push(num);
    }
    Ok(out)
}

impl DeltaVolume {
    /// Appends one profile's `.wcounts` and `.obsr` records during the volume pass, deferring any
    /// Pass-3 ASN.1 posts or errors to `finish`.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1582-1707
    pub(super) fn add_profile(
        &mut self,
        profile_index: usize,
        pssm: &Pssm,
        seq_size: i32,
        filename: &[u8],
        posts: Vec<Post>,
    ) -> Result<(), DeltaError> {
        if !posts.is_empty() {
            self.deferred_posts.push((profile_index, posts));
        }
        if self.deferred_error.is_some() {
            return Ok(());
        }
        let (freqs, obsr) = match delta_fields(pssm, seq_size, filename) {
            Ok(fields) => fields,
            Err(error) => {
                self.deferred_error = Some((profile_index, error));
                return Ok(());
            }
        };
        let obsr_words = match encode_obsr(obsr) {
            Ok(words) => words,
            Err(error) => {
                self.deferred_error = Some((profile_index, error));
                return Ok(());
            }
        };
        if obsr.len() != freqs.len() / BLASTAA_SIZE {
            self.deferred_error = Some((
                profile_index,
                DeltaError::Exception(
                    [
                        b"Number of frequencies and observations columns do not match in "
                            .as_slice(),
                        filename,
                    ]
                    .concat(),
                ),
            ));
            return Ok(());
        }
        self.scratch.clear();
        self.scratch.reserve((freqs.len() + BLASTAA_SIZE) * 4);
        for &val in freqs.as_ref() {
            let scaled = val * FIXED_POINT_SCALE_FACTOR;
            let Some(encoded) = uint4_of(scaled) else {
                self.deferred_error = Some((
                    profile_index,
                    DeltaError::Rps(RpsError::Refused(format!(
                        "a weightedResFreqsPerPos value of {val} for -dbtype delta (NCBI's \
                         conversion of {val} * 1000 to a 32-bit unsigned integer has no defined \
                         result)"
                    ))),
                ));
                return Ok(());
            };
            self.scratch.extend_from_slice(&encoded.to_le_bytes());
        }
        // Additional column of zeros for compatibility with RPS database.
        self.scratch.extend_from_slice(&[0u8; BLASTAA_SIZE * 4]);

        let freq_entries = freqs.len() + BLASTAA_SIZE;
        let Some(next_freq_offset) = i32::try_from(freq_entries)
            .ok()
            .and_then(|n| self.curr_freq_offset.checked_add(n))
        else {
            self.deferred_error = Some((
                profile_index,
                DeltaError::Rps(RpsError::Refused(
                    "a delta volume whose .wcounts entries pass 2^31 (NCBI's makeprofiledb counts \
                     them past its 32-bit integers)"
                        .to_string(),
                )),
            ));
            return Ok(());
        };
        let Some(next_obsr_offset) = i32::try_from(obsr_words.len())
            .ok()
            .and_then(|n| self.curr_obsr_offset.checked_add(n))
        else {
            self.deferred_error = Some((
                profile_index,
                DeltaError::Rps(RpsError::Refused(
                    "a delta volume whose .obsr entries pass 2^31 (NCBI's makeprofiledb counts \
                     them past its 32-bit integers)"
                        .to_string(),
                )),
            ));
            return Ok(());
        };

        self.freq_offsets
            .push(self.curr_freq_offset / BLASTAA_SIZE as i32);
        self.obsr_offsets.push(self.curr_obsr_offset);
        self.wcounts.write_all(&self.scratch)?;
        self.scratch.clear();
        self.scratch.reserve(obsr_words.len() * 4);
        for word in obsr_words {
            self.scratch.extend_from_slice(&word.to_le_bytes());
        }
        self.obsr.write_all(&self.scratch)?;
        self.curr_freq_offset = next_freq_offset;
        self.curr_obsr_offset = next_obsr_offset;
        Ok(())
    }

    /// Completes `x_UpdateDelta` / `x_WrapUpDelta`: replays Pass-3 ASN.1 posts, raises any
    /// deferred Pass-3 error, and writes the final offset tables into `.wcounts` and `.obsr`.
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1709-1712,1822-1863
    pub(super) fn finish(
        mut self,
        diag: &mut Diag,
        post_fn: impl Fn(&mut Diag, &[Post]),
    ) -> Result<(), DeltaError> {
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
        self.freq_offsets
            .push(self.curr_freq_offset / BLASTAA_SIZE as i32);
        self.obsr_offsets.push(self.curr_obsr_offset);

        let write_offsets = |file: &mut BufWriter<File>, offsets: &[i32]| -> std::io::Result<()> {
            file.flush()?;
            let mut table = Vec::with_capacity(offsets.len() * 4);
            for offset in offsets {
                table.extend_from_slice(&offset.to_le_bytes());
            }
            file.get_mut().seek(SeekFrom::Start(8))?;
            file.get_mut().write_all(&table)?;
            Ok(())
        };
        write_offsets(&mut self.wcounts, &self.freq_offsets)?;
        write_offsets(&mut self.obsr, &self.obsr_offsets)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "delta_tests.rs"]
mod tests;

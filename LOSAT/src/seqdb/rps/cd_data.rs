//! The `.wcounts` and `.obsr` files of a delta database (`makeprofiledb -dbtype delta`): the
//! weighted residue counts and the numbers of independent observations of each conserved domain,
//! which DELTA-BLAST's PSSM of the domain hits reads (`cdd_pssm_input.cpp`). Like the `.freq`
//! file they are read at offsets (`FreqSource`), the rows of the hit domains only, never mapped
//! (NCBI maps them whole).
//!
//! Both files start as the `.rps` file does: the magic number, the number of profiles and
//! `num_profiles + 1` offsets, then the data of every profile. The offsets of `.wcounts` count
//! rows of 28 counts (`Uint4`), those of `.obsr` count `Uint4` values.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:361-459
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:804-889

use std::path::{Path, PathBuf};

use super::{canonical_volume_path, FreqSource, RpsOpenError, RPS_MAGIC_NUM, RPS_MAGIC_NUM_28};

/// `kAlphabetSize` of the conserved-domain input: the columns of a `.wcounts` row.
pub(crate) const CD_ALPHABET: usize = 28;

/// `kRpsScaleFactor`: the observations are stored times 1000.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/cdd_pssm_input.hpp:104-106
const OBSERVATION_SCALE: f64 = 1000.0;

/// One of the two files: its positioned reader and its length.
struct CdFile {
    source: FreqSource,
    /// The file's length in bytes when it was opened.
    len: u64,
    /// The file's extension (`.wcounts` or `.obsr`).
    extension: &'static str,
    /// The file's name in messages (`<volume path><extension>`).
    shown: String,
}

impl CdFile {
    /// Opens `<base><extension>`: NCBI's errors for a file that cannot be mapped and for a wrong
    /// magic number; an empty (0-byte) file (NCBI dereferences a null mapping and crashes) is
    /// LOSAT's refusal. For a 1-to-3-byte file, POSIX mmap zero-fills the remainder of the last
    /// page so NCBI's 4-byte `magic_number` read is defined; LOSAT zero-pads the available bytes
    /// to 4 and applies the same magic check. The profile count and offset table are read only
    /// when a hit domain asks for them (`profile_offsets`), matching NCBI's open.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:255-264
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:386-400
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:427-449
    fn open(base: &Path, extension: &'static str, program: &str) -> Result<Self, RpsOpenError> {
        let volume = canonical_volume_path(base);
        let shown_bytes = [&volume[..], extension.as_bytes()].concat();
        let shown = String::from_utf8_lossy(&shown_bytes).into_owned();
        let mut path = base.as_os_str().to_os_string();
        path.push(extension);
        let map_error = || {
            RpsOpenError::Ncbi(
                [
                    b"Cannot memory map RPS-BLAST database file: ",
                    &shown_bytes[..],
                ]
                .concat(),
            )
        };
        let source = FreqSource::open(&PathBuf::from(path)).map_err(|_| map_error())?;
        let len = source.len().map_err(|_| map_error())? as u64;
        if len == 0 {
            return Err(RpsOpenError::Refused(format!(
                "a profile database whose {extension} file is empty (NCBI BLAST+ crashes on a null mapping) is not supported by LOSAT's {program}"
            )));
        }
        let mut magic_bytes = [0u8; 4];
        let to_read = usize::try_from(len.min(4)).unwrap_or(4);
        source
            .read_at(&mut magic_bytes[..to_read], 0)
            .map_err(|_| map_error())?;
        let magic = i32::from_le_bytes(magic_bytes);
        if magic != RPS_MAGIC_NUM && magic != RPS_MAGIC_NUM_28 {
            return Err(RpsOpenError::Ncbi(
                [
                    b"RPS BLAST profile file (",
                    &shown_bytes[..],
                    b") is either corrupt or constructed for an incompatible architecture",
                ]
                .concat(),
            ));
        }
        Ok(Self {
            source,
            len,
            extension,
            shown,
        })
    }

    /// Reads the profile count and checks the offset table of `<base><extension>`, returning the
    /// profile count and byte offset where profile data starts.
    fn header_layout(&self, program: &str) -> Result<(usize, u64), RpsOpenError> {
        let too_short = || {
            RpsOpenError::Refused(format!(
                "a profile database whose {} file is shorter than its header (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}",
                self.extension,
            ))
        };
        let mut count_bytes = [0u8; 4];
        self.source
            .read_at(&mut count_bytes, 4)
            .map_err(|_| too_short())?;
        let num_profiles = i32::from_le_bytes(count_bytes);
        let num_profiles = usize::try_from(num_profiles).map_err(|_| too_short())?;
        let data_start = 8 + 4 * (num_profiles as u64 + 1);
        if self.len < data_start {
            return Err(too_short());
        }
        Ok((num_profiles, data_start))
    }

    /// The data start and offsets of the profile `oid` and of the next: `(data_start, first, next)`.
    fn profile_offsets(&self, oid: usize, program: &str) -> Result<(u64, u64, u64), RpsOpenError> {
        let (num_profiles, data_start) = self.header_layout(program)?;
        let inconsistent = || {
            RpsOpenError::Refused(format!(
                "a profile database whose {} file does not hold the profile of a hit domain (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}",
                self.extension,
            ))
        };
        if oid >= num_profiles {
            return Err(inconsistent());
        }
        let mut pair = [0u8; 8];
        self.source
            .read_at(&mut pair, 8 + 4 * oid as u64)
            .map_err(|_| inconsistent())?;
        let first = i32::from_le_bytes([pair[0], pair[1], pair[2], pair[3]]);
        let next = i32::from_le_bytes([pair[4], pair[5], pair[6], pair[7]]);
        match (u64::try_from(first), u64::try_from(next)) {
            (Ok(first), Ok(next)) if first <= next => Ok((data_start, first, next)),
            _ => Err(inconsistent()),
        }
    }

    /// `count` little-endian `Uint4` values from value index `from` of the data.
    fn read_u32s(
        &self,
        data_start: u64,
        from: u64,
        count: usize,
        program: &str,
    ) -> Result<Vec<u32>, RpsOpenError> {
        let mut bytes = vec![0u8; count * 4];
        self.source
            .read_at(&mut bytes, data_start + from * 4)
            .map_err(|_| self.shorter(program))?;
        Ok(bytes
            .chunks_exact(4)
            .map(|value| u32::from_le_bytes([value[0], value[1], value[2], value[3]]))
            .collect())
    }

    /// The number of values of the data that the file holds from value index `from` on.
    fn values_from(&self, data_start: u64, from: u64) -> u64 {
        ((self.len.saturating_sub(data_start)) / 4).saturating_sub(from)
    }

    /// LOSAT's refusal of a file shorter than the profiles that a read needs.
    fn shorter(&self, program: &str) -> RpsOpenError {
        RpsOpenError::Refused(format!(
            "a profile database whose {} is shorter than its profiles (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}",
            self.shown
        ))
    }
}

/// The two files of a delta database's first volume.
pub(crate) struct CdProfileFiles {
    wcounts: CdFile,
    obsr: CdFile,
    program: String,
}

/// The data of columns `from..to` of a conserved domain.
pub(crate) struct CdColumns {
    /// Each column's counts over its sum (`x_FillResidueCounts`).
    pub(crate) wfreqs: Vec<[f64; CD_ALPHABET]>,
    /// Each column's independent observations (`x_FillObservations`).
    pub(crate) iobsr: Vec<f64>,
}

impl CdProfileFiles {
    /// Opens `.wcounts`, then `.obsr`, of the volume `base` as `CBlastRPSInfo` does for
    /// DELTA-BLAST (the frequencies file and the observations file).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:576-592
    pub(crate) fn open(base: &Path, program: &str) -> Result<Self, RpsOpenError> {
        Ok(Self {
            wcounts: CdFile::open(base, ".wcounts", program)?,
            obsr: CdFile::open(base, ".obsr", program)?,
            program: program.to_string(),
        })
    }

    /// The reads of one hit of the profile `oid` whose segments end at column `end` at most
    /// (`CdHitData`).
    pub(crate) fn hit(&self, oid: usize, end: usize) -> CdHitData<'_> {
        CdHitData {
            files: self,
            oid,
            end,
            offsets: None,
            counts: Vec::new(),
            counts_from: 0,
            observations: None,
        }
    }
}

/// The data of one hit's profile as its segments read it (`CdHitData::columns`, segment after
/// segment): the profile's offsets in both files and its observations read once, and its counts
/// from the first segment's column to the hit's end (`end`) in one read, so that a hit costs a
/// few positioned reads, not four per segment. Each segment's checks and refusals come as they
/// would from reads of its own columns alone, in the same order.
pub(crate) struct CdHitData<'a> {
    files: &'a CdProfileFiles,
    oid: usize,
    end: usize,
    /// The data start, first and next row of the profile in `.wcounts`.
    offsets: Option<(u64, u64, u64)>,
    /// The counts of the columns `counts_from..` (`CD_ALPHABET` per column).
    counts: Vec<u32>,
    counts_from: usize,
    /// The profile's run-length pairs of `.obsr`.
    observations: Option<Vec<u32>>,
}

impl CdHitData<'_> {
    /// The columns `from..to` of the hit's profile.
    ///
    /// The counts of a column sum in `Uint4` (wrapping, as C's does) and each count divides by the
    /// sum (yielding NaN when the sum is 0, as NCBI's division does); the observations come from
    /// the run-length pairs (a value and a number of columns) that the profile's data holds.
    /// NCBI's range checks are assertions that a release build drops: a range beyond the profile's
    /// data (reading outside it in NCBI) is LOSAT's refusal.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:804-846
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:852-889
    pub(crate) fn columns(&mut self, from: usize, to: usize) -> Result<CdColumns, RpsOpenError> {
        let files = self.files;
        let program = files.program.as_str();
        let beyond = |file: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {file} file holds fewer columns than a domain hit covers (NCBI BLAST+ reads past them) is not supported by LOSAT's {program}"
            ))
        };
        let (data_start, first, next) = match self.offsets {
            Some(offsets) => offsets,
            None => {
                let offsets = files.wcounts.profile_offsets(self.oid, program)?;
                self.offsets = Some(offsets);
                offsets
            }
        };
        // The last row of a profile is a column of zero counts that the sequence does not have.
        let rows = next - first;
        if rows == 0 || to as u64 > rows - 1 {
            return Err(beyond(".wcounts"));
        }
        let held = self.counts.len() / CD_ALPHABET;
        if from < self.counts_from || to > self.counts_from + held {
            // The columns `from..to` as one read of them reads them: refused where the file ends
            // before `to`; with them, the hit's later columns up to its end that the profile and
            // the file hold.
            let in_file = files
                .wcounts
                .values_from(data_start, (first + from as u64) * CD_ALPHABET as u64)
                / CD_ALPHABET as u64;
            if (to - from) as u64 > in_file {
                return Err(files.wcounts.shorter(program));
            }
            let last = (self.end.min((rows - 1) as usize).max(to) - from).min(in_file as usize);
            self.counts = files.wcounts.read_u32s(
                data_start,
                (first + from as u64) * CD_ALPHABET as u64,
                last * CD_ALPHABET,
                program,
            )?;
            self.counts_from = from;
        }
        let start = (from - self.counts_from) * CD_ALPHABET;
        let counts = &self.counts[start..start + (to - from) * CD_ALPHABET];
        let mut wfreqs = Vec::with_capacity(to - from);
        for row in counts.chunks_exact(CD_ALPHABET) {
            let sum = row.iter().fold(0u32, |sum, &count| sum.wrapping_add(count));
            let mut freqs = [0.0; CD_ALPHABET];
            for (freq, &count) in freqs.iter_mut().zip(row) {
                *freq = f64::from(count) / f64::from(sum);
            }
            wfreqs.push(freqs);
        }
        if self.observations.is_none() {
            let (obsr_start, first, next) = files.obsr.profile_offsets(self.oid, program)?;
            self.observations = Some(files.obsr.read_u32s(
                obsr_start,
                first,
                usize::try_from(next - first).unwrap_or(usize::MAX),
                program,
            )?);
        }
        let pairs = self.observations.as_deref().unwrap_or_default();
        let mut iobsr = Vec::with_capacity(to - from);
        let mut start = 0usize;
        for pair in pairs.chunks_exact(2) {
            // The run's second word, read as a signed 32-bit number, is its number of columns: a
            // negative one has none.
            let columns = usize::try_from(pair[1] as i32).unwrap_or(0);
            let end = start.saturating_add(columns);
            for column in start.max(from)..end.min(to) {
                debug_assert_eq!(column - from, iobsr.len());
                iobsr.push(f64::from(pair[0]) / OBSERVATION_SCALE);
            }
            start = end;
            if start >= to {
                break;
            }
        }
        if iobsr.len() != to - from {
            return Err(beyond(".obsr"));
        }
        Ok(CdColumns { wfreqs, iobsr })
    }
}

#[cfg(test)]
#[path = "cd_data_tests.rs"]
mod tests;

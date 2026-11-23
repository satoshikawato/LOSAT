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

/// One of the two files: its positioned reader and its profile count.
struct CdFile {
    source: FreqSource,
    num_profiles: usize,
    /// The file's name in messages (`<volume path><extension>`).
    shown: String,
}

impl CdFile {
    /// Opens `<base><extension>`: NCBI's errors for a file that cannot be mapped and for a wrong
    /// magic number; a file shorter than its header (NCBI reads past it) is LOSAT's refusal.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:255-264
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:386-400
    fn open(base: &Path, extension: &str, program: &str) -> Result<Self, RpsOpenError> {
        let volume = canonical_volume_path(base);
        let shown = format!("{volume}{extension}");
        let mut path = base.as_os_str().to_os_string();
        path.push(extension);
        let source = FreqSource::open(&PathBuf::from(path)).map_err(|_| {
            RpsOpenError::Ncbi(format!("Cannot memory map RPS-BLAST database file: {shown}"))
        })?;
        let too_short = || {
            RpsOpenError::Refused(format!(
                "a profile database whose {extension} file is shorter than its header (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}"
            ))
        };
        let mut header = [0u8; 8];
        source.read_at(&mut header, 0).map_err(|_| too_short())?;
        let magic = i32::from_le_bytes([header[0], header[1], header[2], header[3]]);
        if magic != RPS_MAGIC_NUM && magic != RPS_MAGIC_NUM_28 {
            return Err(RpsOpenError::Ncbi(format!(
                "RPS BLAST profile file ({shown}) is either corrupt or constructed for an incompatible architecture"
            )));
        }
        let num_profiles = i32::from_le_bytes([header[4], header[5], header[6], header[7]]);
        let num_profiles = usize::try_from(num_profiles).map_err(|_| too_short())?;
        let file = Self {
            source,
            num_profiles,
            shown,
        };
        let offsets_end = 8 + 4 * (num_profiles as u64 + 1);
        if file.source.len().map_or(true, |len| (len as u64) < offsets_end) {
            return Err(too_short());
        }
        Ok(file)
    }

    /// The offsets of the profile `oid` and of the next: its data's first and one past its last.
    fn profile_offsets(&self, oid: usize, program: &str) -> Result<(u64, u64), RpsOpenError> {
        let inconsistent = || {
            RpsOpenError::Refused(format!(
                "a profile database whose {} file does not hold the profile of a hit domain (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}",
                self.shown.rsplit('.').next().map_or(String::new(), |e| format!(".{e}")),
            ))
        };
        if oid >= self.num_profiles {
            return Err(inconsistent());
        }
        let mut pair = [0u8; 8];
        self.source
            .read_at(&mut pair, 8 + 4 * oid as u64)
            .map_err(|_| inconsistent())?;
        let first = i32::from_le_bytes([pair[0], pair[1], pair[2], pair[3]]);
        let next = i32::from_le_bytes([pair[4], pair[5], pair[6], pair[7]]);
        match (u64::try_from(first), u64::try_from(next)) {
            (Ok(first), Ok(next)) if first <= next => Ok((first, next)),
            _ => Err(inconsistent()),
        }
    }

    /// The byte offset where the data of the profiles starts.
    fn data_start(&self) -> u64 {
        8 + 4 * (self.num_profiles as u64 + 1)
    }

    /// `count` little-endian `Uint4` values from value index `from` of the data.
    fn read_u32s(&self, from: u64, count: usize, program: &str) -> Result<Vec<u32>, RpsOpenError> {
        let mut bytes = vec![0u8; count * 4];
        self.source
            .read_at(&mut bytes, self.data_start() + from * 4)
            .map_err(|_| {
                RpsOpenError::Refused(format!(
                    "a profile database whose {} is shorter than its profiles (NCBI BLAST+ reads past it) is not supported by LOSAT's {program}",
                    self.shown
                ))
            })?;
        Ok(bytes
            .chunks_exact(4)
            .map(|value| u32::from_le_bytes([value[0], value[1], value[2], value[3]]))
            .collect())
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
    /// DELTA-BLAST (`fFrequenciesFile | fObservationsFile`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:576-592
    pub(crate) fn open(base: &Path, program: &str) -> Result<Self, RpsOpenError> {
        Ok(Self {
            wcounts: CdFile::open(base, ".wcounts", program)?,
            obsr: CdFile::open(base, ".obsr", program)?,
            program: program.to_string(),
        })
    }

    /// The columns `from..to` of the profile `oid`.
    ///
    /// The counts of a column sum in `Uint4` (wrapping, as C's does) and each count divides by the
    /// sum; the observations come from the run-length pairs (a value and a number of columns)
    /// that the profile's data holds. NCBI's range checks are assertions that a release build
    /// drops: a range beyond the profile's data (reading outside it in NCBI) is LOSAT's refusal.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:804-846
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:852-889
    pub(crate) fn columns(
        &self,
        oid: usize,
        from: usize,
        to: usize,
    ) -> Result<CdColumns, RpsOpenError> {
        let program = self.program.as_str();
        let beyond = |file: &str| {
            RpsOpenError::Refused(format!(
                "a profile database whose {file} file holds fewer columns than a domain hit covers (NCBI BLAST+ reads past them) is not supported by LOSAT's {program}"
            ))
        };
        let (first, next) = self.wcounts.profile_offsets(oid, program)?;
        // The last row of a profile is a column of zero counts that the sequence does not have.
        let rows = next - first;
        if rows == 0 || to as u64 > rows - 1 {
            return Err(beyond(".wcounts"));
        }
        let counts = self.wcounts.read_u32s(
            (first + from as u64) * CD_ALPHABET as u64,
            (to - from) * CD_ALPHABET,
            program,
        )?;
        let mut wfreqs = Vec::with_capacity(to - from);
        for row in counts.chunks_exact(CD_ALPHABET) {
            let sum = row.iter().fold(0u32, |sum, &count| sum.wrapping_add(count));
            if sum == 0 {
                // NCBI divides by zero and goes on with NaN frequencies.
                return Err(RpsOpenError::Refused(format!(
                    "a profile database whose .wcounts file has a column of zero counts inside a profile is not supported by LOSAT's {program}"
                )));
            }
            let mut freqs = [0.0; CD_ALPHABET];
            for (freq, &count) in freqs.iter_mut().zip(row) {
                *freq = f64::from(count) / f64::from(sum);
            }
            wfreqs.push(freqs);
        }
        let (first, next) = self.obsr.profile_offsets(oid, program)?;
        let pairs = self
            .obsr
            .read_u32s(first, usize::try_from(next - first).unwrap_or(usize::MAX), program)?;
        let mut iobsr = Vec::with_capacity(to - from);
        let mut start = 0usize;
        for pair in pairs.chunks_exact(2) {
            // `(Int4)data[i + 1]` as the number of columns of the run: a negative one has none.
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

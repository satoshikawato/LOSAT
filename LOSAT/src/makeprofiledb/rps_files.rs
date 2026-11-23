//! The RPS files of one volume of makeprofiledb (`CRPS_DbInfo`): the profiles' scores (`.rps`),
//! their frequency ratios (`.freq`), the table of their words (`.loo`) and the text of their
//! parameters (`.aux`), in the 28-letter form that rpsblast reads (`crate::seqdb::rps`).
//!
//! The scores and the frequency ratios of each profile are written as it is added (NCBI writes
//! them at the end of the files too); the offsets at the start of both files are written when
//! the volume closes, where NCBI writes each as its profile is added: on any error NCBI's
//! application removes these files, so only the bytes of a complete volume are seen. The table
//! of words is the thin backbone of every profile of the volume, finalized and repacked at the
//! close (`x_RPS_DbClose`).

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};

use crate::algorithm::pssm_search::lookup::ThinBackbone;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::seqdb::rps::{RPS_BACKBONE_SIZE, RPS_HITS_PER_CELL, RPS_MAGIC_NUM_28, RPS_WORDSIZE};

/// `RPS_DATABASE_VERSION` of makeprofiledb 2.17.0, which writes the 28-letter alphabet.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:117-122
const RPS_DATABASE_VERSION: i32 = RPS_MAGIC_NUM_28;

/// `FREQ_RATIO_SCALE`.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_rps.h:83
const FREQ_RATIO_SCALE: f64 = 1_000_000_000.0;

/// `BLAST_SCORE_MAX` (`INT2_MAX`).
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

/// `INT2_MIN`: the score of the rows past a profile's alphabet.
const INT2_MIN: i32 = i16::MIN as i32;

/// The byte size of the `.loo` header (`BlastRPSLookupFileHeader`).
const LOO_HEADER_BYTES: i32 = 40;

/// The byte size of a backbone cell (`RPSBackboneCell`: the count and three hits).
const CELL_BYTES: i32 = 16;

/// Why the RPS files of a volume cannot go on.
#[derive(Debug)]
pub(super) enum RpsError {
    /// NCBI's `CSeqDBException` (a file that does not open): exit status 2.
    SeqDb(&'static str),
    /// NCBI's `CBlastException` (`eCoreBlastError`): exit status 1.
    Blast(&'static str),
    /// A `std::exception` that leaves NCBI's `Run` ("Unknown exception", exit status 255): the
    /// `runtime_error` of `GetFreqRatios` for an empty list of frequency ratios, the
    /// `bad_array_new_length` of a negative length's matrix (`negative_length_error`).
    Unknown,
    /// LOSAT's refusal of what NCBI does not do in a defined way.
    Refused(String),
    /// A write that failed: NCBI's streams are not tested and the files end short; LOSAT says so.
    Io(std::io::Error),
}

impl From<std::io::Error> for RpsError {
    fn from(error: std::io::Error) -> Self {
        RpsError::Io(error)
    }
}

/// The parameters of the volume that its first profile fixes (`CRPS_DbInfo`'s `gap_open`,
/// `gap_extend`, `matrix` and `scale_factor`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct VolumeParams {
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) matrix: Vec<u8>,
    pub(super) scale_factor: i32,
}

/// The open RPS files of one volume (`CRPS_DbInfo`).
pub(super) struct RpsVolume {
    /// `num_seqs`: the profiles the volume is made for.
    num_seqs: usize,
    pssm: BufWriter<File>,
    lookup_file: File,
    aux: BufWriter<File>,
    freq: BufWriter<File>,
    /// The offset of each profile added (`curr_seq_offset` when it was added).
    offsets: Vec<i32>,
    /// `curr_seq_offset`: the rows of the profiles added, one more for each.
    curr_seq_offset: i32,
    /// The table of words (`lookup`), made by the first profile.
    lookup: Option<ThinBackbone>,
    /// The threshold of the profile being added (`lookup->threshold`, `word_threshold`).
    threshold: i64,
    pub(super) params: VolumeParams,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:725-783
/// `x_InitRPSDbInfo`: opens the four files of the volume `name` in NCBI's order (a file that
/// does not open stops the volume, the files opened before it stay) and writes the magic
/// number, the number of profiles and the place of their offsets (0 to `num_seqs`, which the
/// close replaces); the parameters start as the command line's.
pub(super) fn open_volume(
    name: &str,
    num_seqs: usize,
    params: VolumeParams,
) -> Result<RpsVolume, RpsError> {
    let open = |extension: &str, message: &'static str| {
        File::create(format!("{name}.{extension}")).map_err(|_| RpsError::SeqDb(message))
    };
    let pssm = open("rps", "Failed to open output .rps file ")?;
    let lookup_file = open("loo", "Failed to open output .loo file")?;
    let aux = open("aux", "Failed to open output .aux file")?;
    let freq = open("freq", "Failed to open output .freq file")?;
    let count = i32::try_from(num_seqs)
        .ok()
        .filter(|&count| count < i32::MAX)
        .ok_or_else(|| {
            RpsError::Refused(format!(
                "a volume of {num_seqs} profiles (NCBI's makeprofiledb counts its offsets past \
                 its 32-bit integers)"
            ))
        })?;
    let mut volume = RpsVolume {
        num_seqs,
        pssm: BufWriter::new(pssm),
        lookup_file,
        aux: BufWriter::new(aux),
        freq: BufWriter::new(freq),
        offsets: Vec::with_capacity(num_seqs),
        curr_seq_offset: 0,
        lookup: None,
        threshold: 0,
        params,
    };
    let mut header = Vec::with_capacity(8 + 4 * (num_seqs + 1));
    header.extend_from_slice(&RPS_DATABASE_VERSION.to_le_bytes());
    header.extend_from_slice(&count.to_le_bytes());
    for i in 0..=count {
        header.extend_from_slice(&i.to_le_bytes());
    }
    volume.pssm.write_all(&header)?;
    volume.freq.write_all(&header)?;
    Ok(volume)
}

/// A `double` `value` cast to `Int4`: `None` where C's conversion is undefined (NaN, or a value
/// whose integer part an `Int4` does not hold).
pub(super) fn int4_of(value: f64) -> Option<i32> {
    (value > f64::from(i32::MIN) - 1.0 && value < f64::from(i32::MAX) + 1.0).then_some(value as i32)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
/// `BLAST_Nint` of `value` cast to `Int4`: the nearest `long` (halves away from zero), then its low
/// 32 bits (GCC's conversion to `Int4`); `None` where the conversion to `long` is undefined.
fn nint_int4(value: f64) -> Option<i32> {
    const LONG_LIMIT: f64 = 9_223_372_036_854_775_808.0;
    let x = value + if value >= 0.0 { 0.5 } else { -0.5 };
    (-LONG_LIMIT..LONG_LIMIT)
        .contains(&x)
        .then_some(x as i64 as i32)
}

/// C's `%e` (a stream in `scientific` notation): six digits after the point, the exponent with its
/// sign and at least two digits; `inf`, `-inf`, `nan` and `-nan` as glibc writes them. Rust's
/// precision formatting rounds the exact binary value half to even, as glibc does.
pub(super) fn c_scientific(value: f64) -> String {
    if !value.is_finite() {
        let text = if value.is_nan() { "nan" } else { "inf" };
        return if value.is_sign_negative() {
            format!("-{text}")
        } else {
            text.to_string()
        };
    }
    let e_style = format!("{value:.6e}");
    let (mantissa, exponent) = e_style
        .split_once('e')
        .expect("Rust exponent formatting has an 'e'");
    let exponent: i32 = exponent
        .parse()
        .expect("Rust exponent formatting has an integer exponent");
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exponent.unsigned_abs())
}

impl RpsVolume {
    /// Whether the first profile has made the table (`lookup` is not null).
    pub(super) fn has_lookup(&self) -> bool {
        self.lookup.is_some()
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:909-966
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1129-1172
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:236-247
    /// The end of `x_RPSAddFirstSequence`: the table of words with the threshold `threshold`
    /// (the scale factor times the word score threshold: a negative value is 0, a positive one
    /// is kept, 0 and NaN leave blastp's default 11, then `(Int4)`), and the header of the
    /// `.aux` file.
    pub(super) fn add_first(&mut self, threshold: f64) -> Result<(), RpsError> {
        const BLAST_WORD_THRESHOLD_BLASTP: f64 = 11.0;
        let mut option = BLAST_WORD_THRESHOLD_BLASTP;
        if threshold < 0.0 {
            option = 0.0;
        }
        if threshold > 0.0 {
            option = threshold;
        }
        self.threshold = word_threshold(option)?;
        self.lookup = Some(ThinBackbone::new(RPS_WORDSIZE));
        let params = &self.params;
        self.aux.write_all(&params.matrix)?;
        writeln!(self.aux)?;
        writeln!(self.aux, "{}", params.gap_open)?;
        writeln!(self.aux, "{}", params.gap_extend)?;
        writeln!(self.aux, "{}", c_scientific(0.0))?;
        writeln!(self.aux, "{}", c_scientific(0.0))?;
        writeln!(self.aux, "0")?;
        writeln!(self.aux, "0")?;
        writeln!(self.aux, "{:.6}", f64::from(params.scale_factor))?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1087-1092
    /// `lookup->threshold` of a later profile: `scale_factor * threshold` converted to `Int4`
    /// (`word_threshold`).
    pub(super) fn set_threshold(&mut self, threshold: f64) -> Result<(), RpsError> {
        self.threshold = word_threshold(threshold)?;
        Ok(())
    }

    /// `curr_seq_offset` after a profile of `seq_size` letters, or LOSAT's refusal where NCBI's
    /// `Int4` overflows.
    fn next_offset(&self, seq_size: usize) -> Result<i32, RpsError> {
        i32::try_from(seq_size)
            .ok()
            .and_then(|size| size.checked_add(1))
            .and_then(|rows| self.curr_seq_offset.checked_add(rows))
            .ok_or_else(|| {
                RpsError::Refused(
                    "a volume whose profiles hold 2^31 positions or more (NCBI's makeprofiledb \
                     counts them past its 32-bit integers)"
                        .to_string(),
                )
            })
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1112-1184
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1088-1112
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1009-1057
    /// `x_RPSUpdatePSSM`, `x_RPSUpdateLookup`, `x_UpdateFreqRatios` and the profile's lines of
    /// the `.aux` file: `scores` and `freq` in the order of the file (`num_rows` rows of
    /// `seq_size` positions, by row or by position), `freq_matrix` the frequency ratios as
    /// `GetFreqRatios` holds them for a profile by row (28 rows of `seq_size` values, row by
    /// row), `kappa` the profile's.
    pub(super) fn add_profile(&mut self, profile: &ProfileData<'_>) -> Result<(), RpsError> {
        let seq_size = profile.seq_size;
        let next_offset = self.next_offset(seq_size)?;
        let rows = pssm_rows(profile.scores, profile.num_rows, profile.by_row, seq_size)?;
        let lookup = self
            .lookup
            .as_mut()
            .expect("the first profile made the table");
        lookup.add_pssm_words(
            &rows,
            self.curr_seq_offset,
            self.threshold,
            &[(0, seq_size as i64 - 1)],
        );
        let mut bytes = Vec::with_capacity(rows.rows.len() * BLASTAA_SIZE * 4);
        for row in &rows.rows {
            for score in row {
                bytes.extend_from_slice(&score.to_le_bytes());
            }
        }
        self.pssm.write_all(&bytes)?;
        bytes.clear();
        freq_rows(profile, &mut bytes)?;
        self.freq.write_all(&bytes)?;
        writeln!(self.aux, "{seq_size}")?;
        writeln!(self.aux, "{}", c_scientific(profile.kappa))?;
        self.offsets.push(self.curr_seq_offset);
        self.curr_seq_offset = next_offset;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1190-1301
    /// `x_RPS_DbClose`: the offsets of both files (each profile's, then the total), then the
    /// table of words: "Empty database" without a profile; else finalized (`eBackbone`), each
    /// hit moved to the end of its word, the overflow of the cells of more than three hits
    /// packed, and written after the header with the backbone filled to 32,769 cells.
    pub(super) fn close(mut self) -> Result<(), RpsError> {
        let mut table = Vec::with_capacity(4 * (self.offsets.len() + 1));
        for offset in self.offsets.iter().chain([&self.curr_seq_offset]) {
            table.extend_from_slice(&offset.to_le_bytes());
        }
        debug_assert!(self.offsets.len() <= self.num_seqs);
        for file in [&mut self.pssm, &mut self.freq] {
            file.flush()?;
            file.get_mut().seek(SeekFrom::Start(8))?;
            file.get_mut().write_all(&table)?;
        }
        self.aux.flush()?;
        let Some(lookup) = self.lookup.take() else {
            return Err(RpsError::Blast("Empty database"));
        };
        // The table's own threshold is not written (the header's fields are 0).
        let table = lookup.finalize(self.curr_seq_offset as usize, 0);
        let mut overflow = table.overflow;
        let mut backbone = table.backbone;
        let word_end = RPS_WORDSIZE as i32 - 1;
        let mut cursor = 0usize;
        for cell in backbone.iter_mut() {
            let num_used = cell.num_used;
            if num_used == 0 {
                continue;
            }
            if num_used <= RPS_HITS_PER_CELL {
                for entry in &mut cell.entries[..num_used as usize] {
                    *entry += word_end;
                }
            } else {
                let old_cursor = cell.cursor();
                cell.entries = [overflow[old_cursor] + word_end, (cursor * 4) as i32, 0];
                for i in 1..num_used as usize {
                    overflow[cursor] = overflow[old_cursor + i] + word_end;
                    cursor += 1;
                }
            }
        }
        let cells = RPS_BACKBONE_SIZE as i32 + 1;
        let end_of_overflow = i32::try_from(cursor)
            .ok()
            .and_then(|cursor| cursor.checked_mul(4))
            .and_then(|bytes| (LOO_HEADER_BYTES + cells * CELL_BYTES).checked_add(bytes))
            .ok_or_else(|| {
                RpsError::Refused(
                    "a table of words of 2^31 bytes or more (NCBI's makeprofiledb counts its \
                     size past its 32-bit integers)"
                        .to_string(),
                )
            })?;
        let mut out = BufWriter::new(&mut self.lookup_file);
        let mut header = [0i32; 10];
        header[0] = RPS_DATABASE_VERSION;
        header[8] = LOO_HEADER_BYTES;
        header[9] = end_of_overflow;
        for value in header {
            out.write_all(&value.to_le_bytes())?;
        }
        for cell in &backbone {
            out.write_all(&cell.num_used.to_le_bytes())?;
            for entry in cell.entries {
                out.write_all(&entry.to_le_bytes())?;
            }
        }
        let empty = [0u8; CELL_BYTES as usize];
        for _ in backbone.len()..cells as usize {
            out.write_all(&empty)?;
        }
        drop(backbone);
        for value in &overflow[..cursor] {
            out.write_all(&value.to_le_bytes())?;
        }
        out.flush()?;
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:245
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1087-1092
// NCBI converts the scaled threshold to `Int4` (truncation toward zero), which C leaves
// undefined outside `Int4` (x86: `INT_MIN`, every word a hit; cases V25, S_cwstbig). LOSAT keeps
// the truncated value in 64 bits and compares the words' scores with it (`ThinBackbone::
// add_pssm_words`): the same threshold inside `Int4`, the mathematical one outside (beyond the
// 64 bits, where no word's score reaches or every one does, the saturated value). A NaN has no
// value to compare with and is refused.
/// `lookup->threshold`: `threshold` cast to `Int4`, kept in 64 bits.
fn word_threshold(threshold: f64) -> Result<i64, RpsError> {
    if threshold.is_nan() {
        return Err(RpsError::Refused(
            "a word score threshold that is not a number after the scaling (NCBI's conversion \
             of it to a 32-bit integer has no defined result)"
                .to_string(),
        ));
    }
    Ok(threshold as i64)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1106-1167
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:144-156
// A query length past `INT_MAX` is a negative `int`, and passes `x_CheckInputScoremat` when
// `numColumns` is that number. `Create` of one more than `seq_size` allocates an array of that
// many `Int4` row pointers: from -2 down a negative array size, whose new-expression throws
// `std::bad_array_new_length` (the handler's unknown exception); at -1 no row, so the loops over
// the positions do not run: by position the scores are found short ("PSSM was truncated early"), by
// row they are found too many ("PSSM too large for this sequence"), and by row without scores NCBI
// writes its extra row at index -1, before the array (undefined behaviour), which LOSAT refuses.
/// The error of `x_RPSUpdatePSSM` for a profile of the negative length `seq_size` whose scores
/// `scores` are stored by row when `by_row`.
pub(super) fn negative_length_error(seq_size: i32, by_row: bool, scores: &[i32]) -> RpsError {
    debug_assert!(seq_size < 0);
    if seq_size < -1 {
        RpsError::Unknown
    } else if !by_row {
        RpsError::Blast("PSSM was truncated early")
    } else if !scores.is_empty() {
        RpsError::Blast("PSSM too large for this sequence")
    } else {
        RpsError::Refused(
            "a profile by row without scores whose query length is 4294967295 (NCBI's \
             makeprofiledb writes the row past its last position before its matrix)"
                .to_string(),
        )
    }
}

/// One profile as the RPS files read it.
pub(super) struct ProfileData<'a> {
    pub(super) seq_size: usize,
    /// `numRows`.
    pub(super) num_rows: usize,
    /// `byRow`.
    pub(super) by_row: bool,
    /// `finalData.scores`, or the engine's (28 per position) for a profile of frequency ratios.
    pub(super) scores: &'a [i32],
    /// `intermediateData.freqRatios` as read.
    pub(super) freq_ratios: &'a [f64],
    /// `finalData.kappa`.
    pub(super) kappa: f64,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1112-1170
/// `x_RPSUpdatePSSM`'s matrix: `seq_size + 1` rows of 28 scores read from `scores` by position
/// (`num_rows` per position, the rest `INT2_MIN`) or by row, with the extra row of
/// `-BLAST_SCORE_MAX`; "PSSM was truncated early" when the scores end first, "PSSM too large for
/// this sequence" when some are left.
fn pssm_rows(
    scores: &[i32],
    num_rows: usize,
    by_row: bool,
    seq_size: usize,
) -> Result<QueryPssm, RpsError> {
    let mut rows = vec![[INT2_MIN; BLASTAA_SIZE]; seq_size + 1];
    let mut values = scores.iter().copied();
    let complete = if !by_row {
        let mut complete = true;
        'positions: for row in rows.iter_mut().take(seq_size) {
            for cell in row.iter_mut().take(num_rows) {
                match values.next() {
                    Some(score) => *cell = score,
                    None => {
                        complete = false;
                        break 'positions;
                    }
                }
            }
        }
        complete
    } else {
        let mut complete = true;
        'letters: for j in 0..num_rows {
            for row in rows.iter_mut().take(seq_size) {
                match values.next() {
                    Some(score) => row[j] = score,
                    None => {
                        complete = false;
                        break 'letters;
                    }
                }
            }
        }
        complete
    };
    if !complete {
        return Err(RpsError::Blast("PSSM was truncated early"));
    }
    if values.next().is_some() {
        return Err(RpsError::Blast("PSSM too large for this sequence"));
    }
    rows[seq_size] = [-BLAST_SCORE_MAX; BLASTAA_SIZE];
    Ok(QueryPssm { rows })
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1009-1057
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:215-237,263-282
// NCBI reference (598d8ae6): c++/include/util/math/matrix.hpp:370-376
// By position, each position takes the next `numRows` ratios while they last (the rest 0). By
// row, NCBI reads `GetFreqRatios`'s matrix (28 rows of `seq_size` columns, stored row by row) at
// `(i, j)` for position `i` and letter `j`: the value at `i * seq_size + j` of its storage, which
// lies inside it only for profiles of up to 28 positions; NCBI reads past the matrix for longer
// ones, and past the list of ratios for a list shorter than the matrix's rows, which LOSAT refuses.
/// `x_UpdateFreqRatios`: the profile's rows of the `.freq` file, then a row of zeros, appended
/// to `out`.
fn freq_rows(profile: &ProfileData<'_>, out: &mut Vec<u8>) -> Result<(), RpsError> {
    let seq_size = profile.seq_size;
    let push = |value: f64, out: &mut Vec<u8>| -> Result<(), RpsError> {
        let scaled = nint_int4(value * FREQ_RATIO_SCALE).ok_or_else(|| {
            RpsError::Refused(format!(
                "a frequency ratio of {value} (NCBI's conversion of it to an integer has no \
                 defined result)"
            ))
        })?;
        out.extend_from_slice(&scaled.to_le_bytes());
        Ok(())
    };
    if !profile.by_row {
        let mut ratios = profile.freq_ratios.iter().copied();
        for _ in 0..seq_size {
            let mut j = 0;
            while j < profile.num_rows {
                let Some(ratio) = ratios.next() else {
                    break;
                };
                push(ratio, out)?;
                j += 1;
            }
            for _ in j..BLASTAA_SIZE {
                out.extend_from_slice(&0i32.to_le_bytes());
            }
        }
    } else {
        let matrix = freq_ratio_matrix(profile.freq_ratios, profile.num_rows, seq_size, true)?;
        if seq_size > BLASTAA_SIZE {
            return Err(RpsError::Refused(format!(
                "a profile by row of {seq_size} positions (more than {BLASTAA_SIZE}; NCBI's \
                 makeprofiledb reads its frequency ratios past their matrix)"
            )));
        }
        for i in 0..seq_size {
            for j in 0..BLASTAA_SIZE {
                push(matrix[i * seq_size + j], out)?;
            }
        }
    }
    out.extend_from_slice(&[0u8; BLASTAA_SIZE * 4]);
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:215-237,263-282
/// `CScorematPssmConverter::GetFreqRatios`: `RpsError::Unknown` for an empty list (its
/// `runtime_error`), else the matrix of 28 rows of `num_columns` ratios,
/// stored row by row, the rows past `num_rows` 0, filled from `ratios` by row or by column.
/// NCBI reads past the list where it holds fewer than `num_rows * num_columns` ratios, which
/// LOSAT refuses; ratios past those are not read.
pub(super) fn freq_ratio_matrix(
    ratios: &[f64],
    num_rows: usize,
    num_columns: usize,
    by_row: bool,
) -> Result<Vec<f64>, RpsError> {
    if ratios.is_empty() {
        return Err(RpsError::Unknown);
    }
    if ratios.len() < num_rows * num_columns {
        return Err(RpsError::Refused(format!(
            "a profile of {num_columns} positions and {num_rows} rows with {} frequency ratios \
             (NCBI's makeprofiledb reads past the list)",
            ratios.len()
        )));
    }
    let mut matrix = vec![0.0f64; BLASTAA_SIZE * num_columns];
    let mut values = ratios.iter().copied();
    if by_row {
        for r in 0..num_rows {
            for c in 0..num_columns {
                matrix[r * num_columns + c] = values.next().expect("counted");
            }
        }
    } else {
        for c in 0..num_columns {
            for r in 0..num_rows {
                matrix[r * num_columns + c] = values.next().expect("counted");
            }
        }
    }
    Ok(matrix)
}

#[cfg(test)]
#[path = "rps_files_tests.rs"]
mod tests;

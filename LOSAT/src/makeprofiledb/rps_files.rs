//! The RPS files of one volume of makeprofiledb (`CRPS_DbInfo`): the profiles' scores (`.rps`),
//! their frequency ratios (`.freq`), the table of their words (`.loo`) and the text of their
//! parameters (`.aux`), in the 28-letter form that rpsblast reads (`crate::seqdb::rps`).
//!
//! Each file is NCBI's `CNcbiOfstream`, the file stream of the database's writer
//! (`crate::writedb::ofstream::OfStream`), written with NCBI's calls: the header of `.rps` and
//! `.freq` as one write per integer; for each profile its rows at the end of the file (a
//! `seekp` to the end), then its offset at its place in the header (a `seekp` there), the `.aux`
//! text one value and one newline at a time; the last offset and the table of words at the close.
//! Where the bytes reach the files is the stream's rule, so a run that ends at a closed pipe
//! (`crate::writedb::ofstream::end_at_closed_pipe`) leaves the bytes NCBI's process, ended by
//! SIGPIPE there, leaves: the rows of the profiles added, the offsets but the last one's, the
//! `.aux` text only past a full buffer, no table of words (line RG item T8c). The table of words is
//! the thin backbone of every profile of the volume, finalized and repacked at the close
//! (`x_RPS_DbClose`).

use crate::algorithm::pssm_search::lookup::ThinBackbone;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::seqdb::rps::{RPS_BACKBONE_SIZE, RPS_HITS_PER_CELL, RPS_MAGIC_NUM_28, RPS_WORDSIZE};
use crate::utils::os_bytes::path_of_bytes;
use crate::writedb::ofstream::OfStream;

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
    /// A write that failed: NCBI's streams are not tested and the files end short; LOSAT says so
    /// (the RPS files' streams: the file whose write failed; the delta files: the error).
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
    /// `db_name`: the volume's name, which the files' names extend.
    name: Vec<u8>,
    /// `num_seqs`: the profiles the volume is made for.
    num_seqs: usize,
    /// `pssm_file`: the scores (`.rps`).
    pssm: OfStream,
    /// `lookup_file`: the table of words (`.loo`), written at the close.
    lookup_file: OfStream,
    /// `aux_file`: the text of the parameters and of each profile's length and kappa (`.aux`).
    aux: OfStream,
    /// `freq_file`: the frequency ratios (`.freq`).
    freq: OfStream,
    /// The end of `.rps` and of `.freq`, where NCBI's `seekp` to the end puts the next profile's
    /// rows: their header, then the rows written.
    pssm_end: u64,
    freq_end: u64,
    /// The profiles added: the `seq_index` of the next one.
    added: usize,
    /// `curr_seq_offset`: the rows of the profiles added, one more for each.
    curr_seq_offset: i32,
    /// The table of words (`lookup`), made by the first profile.
    lookup: Option<ThinBackbone>,
    /// The threshold of the profile being added (`lookup->threshold`, `word_threshold`).
    threshold: i64,
    pub(super) params: VolumeParams,
}

/// The place of the offset of the profile `seq_index` in `.rps` and `.freq`: after the magic
/// number and the count, four bytes for each profile before it.
fn offset_place(seq_index: usize) -> u64 {
    8 + 4 * seq_index as u64
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:725-783
/// `x_InitRPSDbInfo`: opens the four files of the volume `name` in NCBI's order (a file that
/// does not open stops the volume, the files opened before it stay) and writes the magic
/// number, the number of profiles and the place of their offsets (0 to `num_seqs`, which each
/// profile and the close replace), each integer to both files in turn; the parameters start as
/// the command line's.
pub(super) fn open_volume(
    name: &[u8],
    num_seqs: usize,
    params: VolumeParams,
) -> Result<RpsVolume, RpsError> {
    let open = |extension: &str, message: &'static str| {
        let stream = OfStream::create(&path_of_bytes(&[name, b".", extension.as_bytes()].concat()));
        if stream.fail() {
            Err(RpsError::SeqDb(message))
        } else {
            Ok(stream)
        }
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
    let header_end = offset_place(num_seqs + 1);
    let mut volume = RpsVolume {
        name: name.to_vec(),
        num_seqs,
        pssm,
        lookup_file,
        aux,
        freq,
        pssm_end: header_end,
        freq_end: header_end,
        added: 0,
        curr_seq_offset: 0,
        lookup: None,
        threshold: 0,
        params,
    };
    for value in [RPS_DATABASE_VERSION, count].into_iter().chain(0..=count) {
        volume.pssm.write(&value.to_le_bytes());
        volume.freq.write(&value.to_le_bytes());
    }
    volume.check()?;
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
    /// `.aux` file, each value and each newline a write of its own (NCBI's `<<`).
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
        let lines = [
            self.params.matrix.clone(),
            self.params.gap_open.to_string().into_bytes(),
            self.params.gap_extend.to_string().into_bytes(),
            c_scientific(0.0).into_bytes(),
            c_scientific(0.0).into_bytes(),
            b"0".to_vec(),
            b"0".to_vec(),
            format!("{:.6}", f64::from(self.params.scale_factor)).into_bytes(),
        ];
        for line in &lines {
            self.aux_line(line);
        }
        self.check()
    }

    /// A line of `.aux` as NCBI inserts it: the value, then a newline, each a write of its own.
    fn aux_line(&mut self, value: &[u8]) {
        self.aux.write(value);
        self.aux.write(b"\n");
    }

    /// `Err(RpsError::Io)` naming the first file whose stream failed (a write, a seek or the
    /// close that its system call refused), which NCBI does not test.
    fn check(&self) -> Result<(), RpsError> {
        let files = [
            (&self.pssm, "rps"),
            (&self.lookup_file, "loo"),
            (&self.aux, "aux"),
            (&self.freq, "freq"),
        ];
        match files.iter().find(|(stream, _)| stream.fail()) {
            Some((_, extension)) => Err(RpsError::Io(std::io::Error::other(format!(
                "{}.{extension}",
                String::from_utf8_lossy(&self.name)
            )))),
            None => Ok(()),
        }
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1478-1479
    // NCBI writes the rows one by one; written here as one write, they reach the file at the same
    // point (the `seekp` after them writes what the stream holds), and nothing can end the run
    // between.
    /// `x_RPSUpdatePSSM`, `x_RPSUpdateLookup`, `x_UpdateFreqRatios` and the profile's lines of
    /// the `.aux` file: `scores` and `freq` in the order of the file (`num_rows` rows of
    /// `seq_size` positions, by row or by position), `freq_matrix` the frequency ratios as
    /// `GetFreqRatios` holds them for a profile by row (28 rows of `seq_size` values, row by
    /// row), `kappa` the profile's. In `.rps` and `.freq` the rows go to the end of the file and
    /// the profile's offset (`curr_seq_offset`) to its place in the header.
    pub(super) fn add_profile(&mut self, profile: &ProfileData<'_>) -> Result<(), RpsError> {
        let seq_size = profile.seq_size;
        let next_offset = self.next_offset(seq_size)?;
        let rows = pssm_rows(profile.scores, profile.num_rows, profile.by_row, seq_size)?;
        let place = offset_place(self.added);
        let offset = self.curr_seq_offset.to_le_bytes();
        let mut bytes = Vec::with_capacity(rows.rows.len() * BLASTAA_SIZE * 4);
        for row in &rows.rows {
            for score in row {
                bytes.extend_from_slice(&score.to_le_bytes());
            }
        }
        self.pssm.seekp(self.pssm_end);
        self.pssm.write(&bytes);
        self.pssm_end += bytes.len() as u64;
        self.pssm.seekp(place);
        self.pssm.write(&offset);
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
        bytes.clear();
        freq_rows(profile, &mut bytes)?;
        self.freq.seekp(self.freq_end);
        self.freq.write(&bytes);
        self.freq_end += bytes.len() as u64;
        self.freq.seekp(place);
        self.freq.write(&offset);
        self.aux_line(seq_size.to_string().as_bytes());
        self.aux_line(c_scientific(profile.kappa).as_bytes());
        self.added += 1;
        self.curr_seq_offset = next_offset;
        self.check()
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1190-1301
    // NCBI writes the backbone and the overflow each as one write; written here in pieces, they
    // reach the file whole by its close, and nothing can end the run between.
    /// `x_RPS_DbClose`: the last offset of both files (the total, at the place after the
    /// profiles'), then the table of words: "Empty database" without a profile; else finalized
    /// (`eBackbone`), each hit moved to the end of its word, the overflow of the cells of more
    /// than three hits packed, and written after the header with the backbone filled to 32,769
    /// cells; then the files are flushed and closed in NCBI's order.
    pub(super) fn close(mut self) -> Result<(), RpsError> {
        let place = offset_place(self.num_seqs);
        let total = self.curr_seq_offset.to_le_bytes();
        self.pssm.seekp(place);
        self.pssm.write(&total);
        self.freq.seekp(place);
        self.freq.write(&total);
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
        let mut header = [0i32; 10];
        header[0] = RPS_DATABASE_VERSION;
        header[8] = LOO_HEADER_BYTES;
        header[9] = end_of_overflow;
        write_int4s(&mut self.lookup_file, header.into_iter());
        write_int4s(
            &mut self.lookup_file,
            backbone
                .iter()
                .flat_map(|cell| std::iter::once(cell.num_used).chain(cell.entries)),
        );
        let empty = [0u8; CELL_BYTES as usize];
        for _ in backbone.len()..cells as usize {
            self.lookup_file.write(&empty);
        }
        drop(backbone);
        write_int4s(&mut self.lookup_file, overflow[..cursor].iter().copied());
        for stream in [
            &mut self.lookup_file,
            &mut self.pssm,
            &mut self.aux,
            &mut self.freq,
        ] {
            stream.flush();
            stream.close();
        }
        self.check()
    }
}

/// The bytes of a piece of `write_int4s`.
const INT4_PIECE: usize = 64 * 1024;

/// Writes the `Int4`s `values` (little-endian) to `stream` in pieces of `INT4_PIECE` bytes,
/// each past the stream's buffer, so written at once.
fn write_int4s(stream: &mut OfStream, values: impl Iterator<Item = i32>) {
    let mut piece = Vec::with_capacity(INT4_PIECE);
    for value in values {
        piece.extend_from_slice(&value.to_le_bytes());
        if piece.len() == INT4_PIECE {
            stream.write(&piece);
            piece.clear();
        }
    }
    if !piece.is_empty() {
        stream.write(&piece);
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

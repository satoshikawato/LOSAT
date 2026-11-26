//! The TBLASTX reports from the final HSP list of a run: the displayed rows of every HSP
//! (outfmt 0, the hit records, and the identity counts of outfmt 6 and 7), the per-query
//! blocks of outfmt 7, and the data of the pairwise report.

use std::io::{self, Write};
use std::sync::Arc;

use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::{record_frame, Placements};
use crate::config::ScoringMatrix;
use crate::core::blast_hits::trim_by_max_hsps;
use crate::report::outfmt6::{write_hit_fields_bytes, OutputConfig};
use crate::report::pairwise::PairwiseHit;
use crate::sequence::packed_subject::PackedNucleotide;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::protein_display_score;

use super::blast_engine::{TblastxHsp, TblastxQueryStats};
#[cfg(test)]
use super::hit_collector::QueryLists;
use super::hit_collector::TblastxHspList;
use super::scoring::{Blosum62Score, TbxScore};

/// The 4-bit IUPAC code of a nucleotide as the display's translation table reads it
/// (`A` 1, `C` 2, `G` 4, `T` 8, ambiguity letters their union, `U` as `T`, `X` as `N`);
/// another byte is a gap (0), which translates to `X`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:94,97-112
#[cfg(test)]
fn display_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'U' => 8,
        b'X' => 15,
        upper => b"-ACMGRSVTWYHKDBN"
            .iter()
            .position(|&code| code == upper)
            .map_or(0, |code| code as u8),
    }
}

// No NCBI counterpart (the arrangement): `display_base` of every byte, on the plus strand and
// complemented for the minus strand (the bits reversed, as in `display_residue_scan`), for the
// display table (`display_table`). NCBI's `sm_BaseToIdx` (above) is the plus-strand table.
const fn display_base_table(minus: bool) -> [u8; 256] {
    let codes = b"-ACMGRSVTWYHKDBN";
    let mut table = [0u8; 256];
    let mut byte = 0;
    while byte < 256 {
        let upper = (byte as u8).to_ascii_uppercase();
        let mut mask = 0;
        if upper == b'U' {
            mask = 8;
        } else if upper == b'X' {
            mask = 15;
        } else {
            let mut code = 0;
            while code < codes.len() {
                if codes[code] == upper {
                    mask = code as u8;
                    break;
                }
                code += 1;
            }
        }
        if minus {
            mask = ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3);
        }
        table[byte] = mask;
        byte += 1;
    }
    table
}

const DISPLAY_BASE_PLUS: [u8; 256] = display_base_table(false);
const DISPLAY_BASE_MINUS: [u8; 256] = display_base_table(true);

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:94-112
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:145-148
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:170-209

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:155-160
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:228-231
pub(crate) fn display_codon(masks: [u8; 3], code: &crate::utils::genetic_code::GeneticCode) -> u8 {
    let expansions = [(1, 2), (2, 1), (4, 3), (8, 0)];
    let mut aa = 0;
    for (x, i) in expansions {
        if masks[0] & x == 0 {
            continue;
        }
        for (y, j) in expansions {
            if masks[1] & y == 0 {
                continue;
            }
            for (z, k) in expansions {
                if masks[2] & z == 0 {
                    continue;
                }
                let ch = code.table[16 * i + 4 * j + k];
                if aa == 0 {
                    aa = ch;
                } else if aa != ch {
                    aa = match (aa, ch) {
                        (b'B' | b'D' | b'N', b'D' | b'N') => b'B',
                        (b'Z' | b'E' | b'Q', b'E' | b'Q') => b'Z',
                        (b'J' | b'I' | b'L', b'I' | b'L') => b'J',
                        _ => b'X',
                    };
                }
            }
        }
    }
    if aa == 0 {
        b'X'
    } else {
        aa
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:155-160,228-231 (cited at
// `display_codon` above): NCBI's display reads the residue of a codon's state from a
// table (`m_AminoAcid`) that it fills once per code; `alnvec.cpp:911-918` (at
// `display_residue`) steps the state by the three bases.
/// The display's residue of each codon of 4-bit codes, at `(first << 8) | (second << 4) |
/// third`: `display_codon` of the 4,096 codons of `code`, made at the first use of the code and
/// kept for the process (`GeneticCode::display_slot`).
fn display_table(code: &GeneticCode) -> &'static [u8; 4096] {
    code.display_slot().get_or_init(|| {
        Box::new(std::array::from_fn(|codon| {
            let masks = [
                (codon >> 8) as u8,
                ((codon >> 4) & 15) as u8,
                (codon & 15) as u8,
            ];
            display_codon(masks, code)
        }))
    })
}

/// A genetic code as the display translates with it (`display_residue`): the code's table
/// (`display_table`; graduated 2026-10-27, the switch was `LOSAT_X_DISPTABLE`), looked up once for
/// many residues.
#[derive(Clone, Copy)]
pub(crate) struct DisplayCode {
    table: &'static [u8; 4096],
}

impl DisplayCode {
    pub(crate) fn new(code: &GeneticCode) -> Self {
        Self {
            table: display_table(code),
        }
    }

    /// `display_residue` with this code.
    pub(crate) fn residue(&self, sequence: &[u8], frame: i8, offset: usize) -> u8 {
        let table = self.table;
        let first = 3 * offset + frame.unsigned_abs() as usize - 1;
        let base = |i: usize| -> usize {
            if frame > 0 {
                DISPLAY_BASE_PLUS[sequence[first + i] as usize] as usize
            } else {
                DISPLAY_BASE_MINUS[sequence[sequence.len() - first - i - 1] as usize] as usize
            }
        };
        table[(base(0) << 8) | (base(1) << 4) | base(2)]
    }

    /// `CAlnVec::TranslateNAToAA` of the codons of `nucleotides` (IUPAC letters, whole codons), or
    /// of their reverse complement's when `minus` (the letters of a minus-strand row,
    /// `CAlnVec::GetSeqString`): `residue` of each codon (the query-anchored reports,
    /// `report::query_anchored`).
    pub(crate) fn translate(&self, nucleotides: &[u8], minus: bool, out: &mut Vec<u8>) {
        let frame = if minus { -1 } else { 1 };
        out.extend(
            (0..nucleotides.len() / 3).map(|offset| self.residue(nucleotides, frame, offset)),
        );
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:292-302
// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:292-302
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-121
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:911-918
/// The displayed residue of codon `offset` of `frame` (+1..+3, -1..-3) of `sequence`: the
/// codon's bases read on the frame's strand (the IUPAC complement on the minus strand) and
/// translated with the display's table (`display_codon`, which merges D/N, E/Q and I/L).
#[cfg(test)]
pub(crate) fn display_residue(sequence: &[u8], frame: i8, offset: usize, code: &GeneticCode) -> u8 {
    DisplayCode::new(code).residue(sequence, frame, offset)
}

/// `display_residue` without the display table (what the display did before `LOSAT_X_DISPTABLE`
/// graduated): each base's code found in the IUPAC letters and the codon translated by
/// `display_codon`; the reference of the unit test of the table.
#[cfg(test)]
fn display_residue_scan(sequence: &[u8], frame: i8, offset: usize, code: &GeneticCode) -> u8 {
    let first = 3 * offset + frame.unsigned_abs() as usize - 1;
    let masks = std::array::from_fn(|i| {
        if frame > 0 {
            display_base(sequence[first + i])
        } else {
            // The complement reverses the bits of the 4-bit code (A <-> T, C <-> G).
            let mask = display_base(sequence[sequence.len() - first - i - 1]);
            ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
        }
    });
    display_codon(masks, code)
}

/// The displayed query and subject rows of an HSP, in uppercase: the HSP's codons of each
/// sequence translated again from the nucleotides, the query with `-query_gencode` and the
/// subject with `-db_gencode` (`query_code`, `db_code`). The HSP's residue offsets within
/// its frames are those of the search (`Hit::sort_*`); `subject` holds the subject's letters
/// from the codon of residue `subject_first` on (`SubjectLetters`; 0 for a whole record).
fn displayed_rows(
    hsp: &TblastxHsp,
    query: &[u8],
    subject: &[u8],
    subject_first: usize,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) -> (Vec<u8>, Vec<u8>) {
    let (query_row, subject_row) =
        displayed_residues(hsp, query, subject, subject_first, query_code, db_code);
    (query_row.collect(), subject_row.collect())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1948-1953
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
/// The residues of the displayed rows of `hsp` (`displayed_rows`), query and subject, in row
/// order, each translated as it is read.
fn displayed_residues<'r>(
    hsp: &'r TblastxHsp,
    query: &'r [u8],
    subject: &'r [u8],
    subject_first: usize,
    query_code: &'r GeneticCode,
    db_code: &'r GeneticCode,
) -> (impl Iterator<Item = u8> + 'r, impl Iterator<Item = u8> + 'r) {
    let hit = &hsp.hit;
    let query_frame = hit.query_frame as i8;
    let (query_display, subject_display) =
        (DisplayCode::new(query_code), DisplayCode::new(db_code));
    let query_row = (hit.sort_query_offset..hit.sort_query_end)
        .map(move |offset| query_display.residue(query, query_frame, offset));
    let subject_row = (hit.sort_subject_offset..hit.sort_subject_end).map(move |offset| {
        subject_display.residue(subject, hsp.subject_frame, offset - subject_first)
    });
    (query_row, subject_row)
}

// No NCBI counterpart: NCBI counts on the rows it formats (tabular.cpp:1003-1021, at
// `pair_counts`). The identities pass needs only the counts, so it reads the residues pair by
// pair and makes no rows. On LOSAT Web's threaded build (wasi-libc's dlmalloc, one lock for
// the whole heap) the rows of this pass and of `PairwiseRows::hit` were the allocator's time
// of the report stage on the pool (TBXRPTPAR): d04 at 8 threads 138-186 ms per search, of
// which 60 ms under this pass (v4-W s8).
/// `row_counts` of the displayed rows of `hsp` (`displayed_rows`), without making them. The
/// HSPs are ungapped: both rows have the HSP's length.
fn displayed_counts<S: TbxScore>(
    hsp: &TblastxHsp,
    query: &[u8],
    subject: &[u8],
    subject_first: usize,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) -> (usize, usize) {
    let hit = &hsp.hit;
    debug_assert_eq!(
        hit.sort_query_end - hit.sort_query_offset,
        hit.sort_subject_end - hit.sort_subject_offset,
        "ungapped HSP"
    );
    let (query_row, subject_row) =
        displayed_residues(hsp, query, subject, subject_first, query_code, db_code);
    pair_counts(sc, query_row.zip(subject_row))
}

/// `displayed_rows` into two buffers (cleared first), which the BLAST XML reuses from HSP to HSP.
#[allow(clippy::too_many_arguments)]
fn displayed_rows_into(
    hsp: &TblastxHsp,
    query: &[u8],
    subject: &[u8],
    subject_first: usize,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) {
    let (query_residues, subject_residues) =
        displayed_residues(hsp, query, subject, subject_first, query_code, db_code);
    query_row.clear();
    query_row.extend(query_residues);
    subject_row.clear();
    subject_row.extend(subject_residues);
}

/// The letters of the subjects that the displayed rows read, by subject index.
#[derive(Clone, Copy)]
pub(crate) enum SubjectLetters<'a> {
    /// The records of `-subject` (`run_local`'s, or the command line's with LOSAT_X_TBXSTORE=0):
    /// their letters.
    Records(&'a [FastaRecord]),
    /// (LOSAT_X_TBXSTORE) The packed records of `-subject` (`ordinals`, the record of each
    /// subject index, `0..n`, into `stores`): a row reads the letters of its HSP's codons only,
    /// made from the stored codes (`codon_window`).
    Stored {
        stores: &'a [PackedNucleotide],
        ordinals: &'a [u32],
    },
    /// The reported sequences of a `-db` search by rank (`database::ReportedSubjects`): their
    /// lengths, and their stores, read where their batches were reported (empty when no report
    /// reads letters); a row reads the letters of its HSP's codons as `Stored`, the engine's
    /// counts those of the re-evaluation (`database::ReportedStore`).
    Database {
        lengths: &'a [u32],
        stores: &'a [super::database::ReportedStore],
    },
}

impl<'a> SubjectLetters<'a> {
    /// The number of bases of subject `s_idx`.
    pub(crate) fn length(&self, s_idx: usize) -> usize {
        match self {
            Self::Records(records) => records[s_idx].seq().len(),
            Self::Stored { stores, ordinals } => stores[ordinals[s_idx] as usize].len(),
            Self::Database { lengths, .. } => lengths[s_idx] as usize,
        }
    }

    /// Appends the searched letters `[from, to)` of subject `s_idx` in upper case (the
    /// query-anchored sections, `anchored::AnchoredSequences`).
    pub(crate) fn letters_of(&self, s_idx: usize, from: usize, to: usize, out: &mut Vec<u8>) {
        let store = match self {
            Self::Records(records) => {
                out.extend(
                    records[s_idx].seq()[from..to]
                        .iter()
                        .map(u8::to_ascii_uppercase),
                );
                return;
            }
            Self::Stored { stores, ordinals } => &stores[ordinals[s_idx] as usize],
            Self::Database { stores, .. } => &*stores[s_idx].shown,
        };
        let mut letters = Vec::new();
        store.letters_into(from, to, &mut letters);
        out.extend_from_slice(&letters);
    }

    /// The letters that the rows of `hsp` read and the residue offset they start at: its
    /// record's letters from 0, or a stored subject's codons of the HSP made into `window`.
    fn letters<'w>(&self, hsp: &TblastxHsp, window: &'w mut Vec<u8>) -> (&'w [u8], usize)
    where
        'a: 'w,
    {
        let s_idx = hsp.hit.s_idx as usize;
        match self {
            Self::Records(records) => (records[s_idx].seq(), 0),
            Self::Stored { stores, ordinals } => {
                stored_letters(&stores[ordinals[s_idx] as usize], hsp, window)
            }
            Self::Database { stores, .. } => stored_letters(&stores[s_idx].shown, hsp, window),
        }
    }

    /// `letters` as the search's re-evaluation read them, for the engine's counts: a database
    /// sequence's re-evaluated store (`database::ReportedStore::engine`, `B` over the masked
    /// bases of `-db_hard_mask`), the rows' letters otherwise.
    fn engine_letters<'w>(&self, hsp: &TblastxHsp, window: &'w mut Vec<u8>) -> (&'w [u8], usize)
    where
        'a: 'w,
    {
        match self {
            Self::Database { stores, .. } => {
                stored_letters(&stores[hsp.hit.s_idx as usize].engine, hsp, window)
            }
            Self::Records(_) | Self::Stored { .. } => self.letters(hsp, window),
        }
    }

    /// The engine counts of `hsp` (`hsp_engine_counts`) against its subject's letters as the
    /// re-evaluation read them (`engine_letters`); `window` is the buffer of a stored subject's
    /// letters.
    fn engine_counts<S: TbxScore>(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
        sc: S,
    ) -> (usize, usize) {
        let (subject, subject_first) = self.engine_letters(hsp, window);
        hsp_engine_counts(hsp, query, subject, subject_first, query_code, db_code, sc)
    }

    /// The displayed rows of `hsp` (`displayed_rows`) against its subject; `window` is the
    /// buffer of a stored subject's letters.
    fn rows(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
    ) -> (Vec<u8>, Vec<u8>) {
        let (subject, subject_first) = self.letters(hsp, window);
        displayed_rows(hsp, query, subject, subject_first, query_code, db_code)
    }

    // No NCBI counterpart for the attribute: on wasm32 the counts of an HSP stay a function of
    // their own, called once per HSP, so that V8 runs the later HSPs in its optimized code (as
    // `PairwiseRows::hit`); inlined, the serial pass over a whole list is one call in V8's
    // baseline code. Other targets keep the compiler's choice.
    /// The counts of the displayed rows of `hsp` (`displayed_counts`) against its subject;
    /// `window` is the buffer of a stored subject's letters.
    #[cfg_attr(target_arch = "wasm32", inline(never))]
    fn counts<S: TbxScore>(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
        sc: S,
    ) -> (usize, usize) {
        let (subject, subject_first) = self.letters(hsp, window);
        displayed_counts(hsp, query, subject, subject_first, query_code, db_code, sc)
    }
}

impl SubjectLetters<'_> {
    /// `rows` into two buffers (cleared first): the BLAST XML builds the rows of one HSP at a
    /// time into buffers it reuses.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn rows_into(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
        query_row: &mut Vec<u8>,
        subject_row: &mut Vec<u8>,
    ) {
        let (subject, subject_first) = self.letters(hsp, window);
        displayed_rows_into(
            hsp,
            query,
            subject,
            subject_first,
            query_code,
            db_code,
            query_row,
            subject_row,
        );
    }

    /// The engine's identities of `hsp` (`hsp_engine_counts`; the Seq-align's `num_ident`,
    /// which the custom tabular reports read where they show `nident` without the rows,
    /// `tabular`). The identities do not depend on the matrix.
    pub(crate) fn engine_identities(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
    ) -> usize {
        self.engine_counts(hsp, query, query_code, db_code, window, Blosum62Score)
            .0
    }
}

/// The codons of the residues of `hsp` in its subject's packed letters `store`, made into
/// `window`, and the residue offset they start at (`codon_window`).
fn stored_letters<'w>(
    store: &PackedNucleotide,
    hsp: &TblastxHsp,
    window: &'w mut Vec<u8>,
) -> (&'w [u8], usize) {
    let (from, to) = (hsp.hit.sort_subject_offset, hsp.hit.sort_subject_end);
    let (start, length) = codon_window(store.len(), hsp.subject_frame, from, to);
    store.letters_into(start, start + length, window);
    (window, from)
}

/// The engine's identities of `hsp` (`SubjectLetters::engine_identities`) against `store`, its
/// subject's packed letters as the re-evaluation read them (`database::ReportedStore::engine`):
/// the batch loop's, which names the subjects by ordinal id before the reports name them by rank.
pub(crate) fn stored_engine_identities(
    hsp: &TblastxHsp,
    store: &PackedNucleotide,
    query: &[u8],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    window: &mut Vec<u8>,
) -> usize {
    let (subject, subject_first) = stored_letters(store, hsp, window);
    hsp_engine_counts(
        hsp,
        query,
        subject,
        subject_first,
        query_code,
        db_code,
        Blosum62Score,
    )
    .0
}

/// The nucleotides `(start, number)` of the codons of residues `[from, to)` of `frame` in a
/// subject of `length` bases. Within them, the codon of residue `o` is the one
/// `display_residue` reads for `o - from`: from `3 * (o - from) + |frame| - 1` on the plus
/// strand, and as far from the window's end as from the subject's on the minus strand (the
/// window of TBLASTN's report rows, tblastn/stage_e_report.rs `aligned_sequences`).
///
/// No NCBI counterpart: NCBI's formatter reads the subject through its object manager
/// (`CAlnVec`); LOSAT makes only these letters of a stored subject.
fn codon_window(length: usize, frame: i8, from: usize, to: usize) -> (usize, usize) {
    let shift = usize::from(frame.unsigned_abs());
    let number = 3 * (to - from) + shift - 1;
    let start = if frame > 0 {
        3 * from
    } else {
        length + 1 - 3 * to - shift
    };
    (start, number)
}

/// The identities and the positives (identities included) of two displayed rows.
fn row_counts<S: TbxScore>(sc: S, query_row: &[u8], subject_row: &[u8]) -> (usize, usize) {
    pair_counts(
        sc,
        query_row.iter().copied().zip(subject_row.iter().copied()),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2131-2150
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1003-1021
/// The identities and the positives (identities included) of the letter pairs of two
/// displayed rows, in row order (`row_counts`, `displayed_counts`), by the display matrix of
/// the search's matrix (BLOSUM62's instance reads the constant matrix).
fn pair_counts<S: TbxScore>(sc: S, pairs: impl Iterator<Item = (u8, u8)>) -> (usize, usize) {
    let mut identities = 0;
    let mut positives = 0;
    for (q, s) in pairs {
        if q == s {
            identities += 1;
            positives += 1;
        } else if protein_display_score(sc.matrix(), q, s) > 0 {
            positives += 1;
        }
    }
    (identities, positives)
}

/// The identity counts of outfmt 6 and 7, which NCBI takes from the displayed rows
/// (`tblastx` does not set the formatter's no-fetch mode, blast_format.cpp:791-792).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1460-1465
pub(crate) fn set_displayed_identities<S: TbxScore>(
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) {
    let mut window = Vec::new();
    for hsp in hits {
        let (identities, positives) = subjects.counts(
            hsp,
            queries[hsp.hit.q_idx as usize].seq(),
            query_code,
            db_code,
            &mut window,
            sc,
        );
        let hit = &mut hsp.hit;
        hit.num_ident = identities;
        hit.num_positives = positives;
        hit.mismatch = (hit.length - identities) as i64;
        hit.identity = if hit.length > 0 {
            identities as f64 / hit.length as f64 * 100.0
        } else {
            0.0
        };
    }
}

/// NCBI's engine identities and positives (identities included) of an HSP, which its
/// Seq-aligns show (`num_ident`, `num_positives`): counted on the residues the search
/// translated (`translation::frame_residue`: the unmasked query frame with `-query_gencode`,
/// the subject frame with `-db_gencode`, as the re-evaluation's buffers, not the display's
/// table), a pair that is not an identity being positive by the search matrix. `subject`
/// holds the subject's letters from the codon of residue `subject_first` on (as
/// `displayed_rows`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2717
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:770-792
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:831-833
fn hsp_engine_counts<S: TbxScore>(
    hsp: &TblastxHsp,
    query: &[u8],
    subject: &[u8],
    subject_first: usize,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) -> (usize, usize) {
    use super::translation::frame_residue;
    let hit = &hsp.hit;
    let query_frame = hit.query_frame as i8;
    let mut identities = 0;
    let mut positives = 0;
    for (q_offset, s_offset) in (hit.sort_query_offset..hit.sort_query_end)
        .zip(hit.sort_subject_offset..hit.sort_subject_end)
    {
        let q = frame_residue(query, query_frame, q_offset, query_code);
        let s = frame_residue(
            subject,
            hsp.subject_frame,
            s_offset - subject_first,
            db_code,
        );
        if q == s {
            identities += 1;
        } else if sc.score(q, s) > 0 {
            positives += 1;
        }
    }
    (identities, identities + positives)
}

/// The engine counts of the Seq-aligns of `hits` (outfmt 8, 12), by HSP index: the engine's
/// `num_ident` (the search's count, which `set_displayed_identities` replaces by the display
/// count; read before it) and its `num_positives`, which the search does not count
/// (`hsp_engine_counts` on the letters the re-evaluation read, `SubjectLetters::engine_letters`,
/// with `B` over the masked bases of `-db_hard_mask`: the same letters give the same identities,
/// checked in debug builds). Made only when a Seq-align format is requested.
pub(crate) fn engine_counts<S: TbxScore>(
    hits: &[TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) -> Vec<(u32, u32)> {
    let mut window = Vec::new();
    hits.iter()
        .map(|hsp| {
            let (identities, positives) = subjects.engine_counts(
                hsp,
                queries[hsp.hit.q_idx as usize].seq(),
                query_code,
                db_code,
                &mut window,
                sc,
            );
            debug_assert_eq!(identities, hsp.hit.num_ident, "engine identities");
            (hsp.hit.num_ident as u32, positives as u32)
        })
        .collect()
}

/// The SEG masks of a query in nucleotide coordinates, by frame.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1256
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:908-931
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:143-160
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:923-929
/// The minus-frame conversion is NCBI's: it starts one codon after the masked run, so the
/// last residue of a run is not lowercased (and a run of one residue is empty).
pub(super) fn query_dna_masks(
    stats: &TblastxQueryStats,
    query_length: usize,
) -> Vec<(i8, Vec<(i32, i32)>)> {
    let protein_length = |context: usize| {
        if query_length == 0 || query_length <= context % 3 {
            0
        } else {
            (query_length - context % 3) / 3
        }
    };
    let forward: [usize; 3] = std::array::from_fn(protein_length);
    let reverse: [usize; 3] = std::array::from_fn(|i| protein_length(i + 3));
    let dna_length = 2 + if forward[0] == 0 { reverse } else { forward }
        .iter()
        .sum::<usize>();
    stats
        .seg_masks
        .iter()
        .map(|(frame, masks)| {
            let ranges = masks
                .iter()
                .map(|&(start, end)| (start as i32, end as i32 - 1))
                .collect();
            let combined = crate::algorithm::common::masks::combine_masks(ranges);
            (
                *frame,
                crate::algorithm::common::masks::protein_to_dna_masks(
                    &combined, dna_length, *frame,
                ),
            )
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:936-951
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4305-4326
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:551-558,593-595
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2500-2521
/// The masks that the report shows for a query: NCBI's `Map` puts each mask of the searched
/// letters (nucleotide coordinates of its frame) into the query's record by the start of
/// its interval (`offset`, 0 without `-query_loc`), cut at the interval's end, and drops a
/// mask that covers the whole interval (`kTarget`). Each keeps the frame of the search
/// (`BLAST_ContextToFrame`), which counts from the interval's ends.
fn shown_query_masks(
    masks: Vec<(i8, Vec<(i32, i32)>)>,
    query_length: usize,
    offset: usize,
) -> Vec<(i8, Vec<(i32, i32)>)> {
    masks
        .into_iter()
        .map(|(frame, ranges)| {
            let shown = ranges
                .into_iter()
                .filter_map(|(from, to)| {
                    crate::blastinput::seq_range::record_mask(from, to, query_length, offset)
                })
                .collect();
            (frame, shown)
        })
        .collect()
}

/// The masks of a query as the query-anchored sections read them (`query_anchored::Mask`): the
/// shown masks of every frame (`shown_query_masks`, those of outfmt 0's rows), each with the frame
/// of its context.
pub(crate) fn query_frame_masks(
    stats: &TblastxQueryStats,
    query_length: usize,
    offset: usize,
) -> Vec<crate::report::query_anchored::Mask> {
    shown_query_masks(query_dna_masks(stats, query_length), query_length, offset)
        .into_iter()
        .flat_map(|(frame, ranges)| {
            ranges
                .into_iter()
                .map(move |(from, to)| (from as u32, to as u32, frame))
        })
        .collect()
}

/// The masked regions of a query as NCBI hands them to the BLAST XML writer
/// (`Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959): the SEG intervals of each frame in
/// nucleotides, joined as `BlastSeqLocCombine` of `filter_out` with 0 joins them (only overlapping ones,
/// `query_dna_masks`), moved into the query's record and without the whole interval
/// (`shown_query_masks`), in context order, each with the frame of its context. `out` is cleared
/// first.
pub(crate) fn xml_query_masks(
    stats: &TblastxQueryStats,
    query_length: usize,
    offset: usize,
    out: &mut Vec<crate::report::blast_xml::QueryMask>,
) {
    out.clear();
    let shown = shown_query_masks(query_dna_masks(stats, query_length), query_length, offset);
    for (frame, ranges) in shown {
        for (from, to) in ranges {
            out.push(crate::report::blast_xml::QueryMask {
                from,
                to,
                frame: i32::from(frame),
            });
        }
    }
}

/// The masks of a displayed query row `row`, set as row 0 of `row_masks` (`HspMasks`, the
/// rows' builder's, which keeps its buffers from row to row): the shown masks (record
/// coordinates, `shown_query_masks`, as the query-anchored sections show them) labelled with the
/// row's frame `frame` (`s_GetStdsegMasterFrame`, the frame `x_OutputSeq` matches) that meet the
/// row, at its columns (`HspMasks::set_row_met`, `x_FillLocList`: width 3, `[low, high]` the
/// 0-based nucleotide range of the row on the plus strand of the record, a residue for each
/// codon), in the masks' order; and the row with them in lower case (`HspMasks::lower_row`), as
/// the hit records of the library API and the Web adapter show it (and outfmt 0, which gets no
/// lists, `query_anchored::output_seq`). `false` when no mask meets the row. With
/// `-query_loc` a mask's label is the frame of the searched letters, which can differ from the
/// frame the record gives the same letters: NCBI then shows a mask of another frame on the row
/// (`locFrame` equal to `frame`), and so does LOSAT.
///
/// The scan of a frame's ranges (`MaskSeek::for_each_meeting`; graduated 2026-10-27, the switch
/// was `LOSAT_X_TBXMASKSEEK`) starts at the first range whose end, or the end of a range before
/// it, reaches `low`, and stops at the first range that starts after `high`; the ranges it leaves
/// out are the ones a scan of all the ranges skips (that scan is the reference of the unit test
/// `mask_seek_lowercases_the_letters_of_the_full_scan`). The ranges it keeps go back into the
/// order of the shown masks (each range's place in the list, which `MaskSeek` holds,
/// `HspMasks::add_meeting`).
fn query_row_masks(
    row: &mut [u8],
    frame: i8,
    low: i32,
    high: i32,
    masks: &[FrameMasks],
    row_masks: &mut crate::report::query_anchored::HspMasks,
) -> bool {
    for frame_masks in masks.iter().filter(|masks| masks.frame == frame) {
        row_masks.add_meeting(&frame_masks.seek, low, high, |from, to| {
            (from.max(0) as u32, to.max(0) as u32, frame)
        });
    }
    row_masks.set_row_met(
        0,
        frame,
        row,
        |_| true,
        (low.max(0) as u32, high.max(0) as u32),
        frame > 0,
        3,
    );
    if row_masks.row(0).is_none() {
        return false;
    }
    // A row of TBLASTX's ungapped HSP holds one residue for each codon of its range.
    debug_assert_eq!(i64::from(high) - i64::from(low) + 1, 3 * row.len() as i64);
    row_masks.lower_row(0, row);
    true
}

/// Every HSP of the final hit list (in the order of all formats) with its displayed rows,
/// frames, counts, linked-set size and subject: the rows made up front (`run_impl::rows_up_front`,
/// the Web adapter's hit records); otherwise outfmt 0 makes the rows of one window of queries at a
/// time with the same `PairwiseRows::hit` (`run_impl::write_tblastx_outputs`).
// Kept out of line: this runs only for outfmt 0 or hit records.
#[inline(never)]
///
/// With `-query_loc` and `-subject_loc` the rows are read from the searched letters, and
/// the hits move into the records (`shift_to_records`) with the records' frames and
/// lengths. The subjects' titles come from `subjects`, their letters and lengths from
/// `subject_letters`. `seqalign_identities` holds the Seq-align's `num_ident` of each HSP of
/// `hits` when the report reads them (`seqalign_identities`). Under `-html` (`html_masks`) each
/// row's masks are added to `html_masks` in the order of `hits` (`report::pairwise::HtmlRowMasks`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn pairwise_hits(
    hits: &[TblastxHsp],
    seqalign_identities: Option<&[u32]>,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    query_stats: &[TblastxQueryStats],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    query_placements: &Placements,
    subject_placements: &Placements,
    mut html_masks: Option<&mut Vec<crate::report::query_anchored::HspMasks>>,
) -> Vec<PairwiseHit> {
    let rows = PairwiseRows::new(
        queries,
        subjects,
        subject_letters,
        query_stats,
        query_code,
        db_code,
        query_placements,
        subject_placements,
    );
    debug_assert!(seqalign_identities.is_none_or(|identities| identities.len() == hits.len()));
    let mut window = Vec::new();
    let mut row_masks = crate::report::query_anchored::HspMasks::default();
    hits.iter()
        .enumerate()
        .map(|(index, hsp)| {
            let hit = rows.hit(
                hsp,
                seqalign_identities.map(|identities| identities[index]),
                &mut window,
                &mut row_masks,
            );
            if let Some(table) = html_masks.as_deref_mut() {
                table.push(row_masks.copy_rows());
            }
            hit
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203
/// The Seq-align's `num_ident` of each HSP of `hits`, by HSP index (`PairwiseHit::seqalign_num_ident`,
/// which the outfmt 0 report reads with `-sorthits` and `-sorthsps`): the engine's identities,
/// which `hits` hold until `set_displayed_identities` replaces them by the display counts (the
/// first value of `engine_counts`, without counting again). Made before that, when outfmt 0's
/// rows or the hit records are made and the report reads them (`-sorthits`, or a `-sorthsps`
/// above 0: `ReportSort::reads_hsp_scores`); a default report reads none and copies none.
pub(crate) fn seqalign_identities(hits: &[TblastxHsp]) -> Vec<u32> {
    hits.iter().map(|hsp| hsp.hit.num_ident as u32).collect()
}

/// The windows of the queries `0..num_queries` whose outfmt 0 rows are made and written
/// together (`x_pairwise_par::query_windows` of the estimates of their HSPs' lengths in a report
/// of `line_length` columns a line, the report's `PairwiseConfig::line_length`, before any row is
/// made, with the budget of `x_pool`, `x_pairwise_par::window_bytes`), as BLASTN's
/// (`blastn::pairwise::pairwise_windows`). `hits` is the final list in report order (queries in
/// input order).
#[inline(never)]
pub(crate) fn pairwise_windows(
    hits: &[TblastxHsp],
    num_queries: usize,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    line_length: usize,
) -> Vec<std::ops::Range<usize>> {
    use crate::report::x_pairwise_par::{query_estimate_of, query_windows, window_bytes};
    let mut estimates = vec![query_estimate_of(std::iter::empty(), line_length); num_queries];
    for query_hits in hits.chunk_by(|a, b| a.hit.q_idx == b.hit.q_idx) {
        if let Some(estimate) = estimates.get_mut(query_hits[0].hit.q_idx as usize) {
            *estimate = query_estimate_of(query_hits.iter().map(|hsp| hsp.hit.length), line_length);
        }
    }
    query_windows(&estimates, window_bytes(x_pool))
}

/// The positives of the outfmt 0 rows `hits` by `matrix`, a matrix other than BLOSUM62, counted
/// again from the rows (`row_counts`) after `PairwiseRows::hit` counted BLOSUM62's: the query
/// row's masked letters are lower case there, and the letters of the rows are upper case. One
/// `hit` for every matrix keeps the helpers that it calls once (`SubjectLetters::rows`,
/// `query_row_masks`) in its body.
pub(crate) fn set_row_positives(hits: &mut [PairwiseHit], matrix: ScoringMatrix) {
    for hit in hits {
        let query_row = hit
            .query_seq
            .as_deref()
            .unwrap_or_default()
            .as_bytes()
            .to_ascii_uppercase();
        let subject_row = hit.subject_seq.as_deref().unwrap_or_default().as_bytes();
        let (_, positives) = row_counts(matrix, &query_row, subject_row);
        hit.hit.num_positives = positives;
        hit.positives = Some(positives);
    }
}

/// The shown SEG masks of one frame of a query (`shown_query_masks`), arranged for the scans of
/// the rows (`query_anchored::MaskSeek`, which holds each mask's place in the frame's shown list).
struct FrameMasks {
    frame: i8,
    seek: crate::report::query_anchored::MaskSeek<i32>,
}

impl FrameMasks {
    fn new(frame: i8, ranges: Vec<(i32, i32)>) -> Self {
        Self {
            frame,
            seek: crate::report::query_anchored::MaskSeek::new(ranges),
        }
    }
}

/// What the outfmt 0 row of an HSP (`PairwiseRows::hit`) reads besides the HSP: the inputs of
/// `pairwise_hits`, with the queries' shown masks made once.
pub(crate) struct PairwiseRows<'a> {
    queries: &'a [FastaRecord],
    subjects: &'a [FastaRecord],
    subject_letters: SubjectLetters<'a>,
    query_masks: Vec<Vec<FrameMasks>>,
    query_code: &'a GeneticCode,
    db_code: &'a GeneticCode,
    query_placements: &'a Placements,
    subject_placements: &'a Placements,
}

impl<'a> PairwiseRows<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        queries: &'a [FastaRecord],
        subjects: &'a [FastaRecord],
        subject_letters: SubjectLetters<'a>,
        query_stats: &[TblastxQueryStats],
        query_code: &'a GeneticCode,
        db_code: &'a GeneticCode,
        query_placements: &'a Placements,
        subject_placements: &'a Placements,
    ) -> Self {
        let query_masks: Vec<Vec<FrameMasks>> = query_stats
            .iter()
            .zip(queries)
            .enumerate()
            .map(|(q_idx, (stats, query))| {
                shown_query_masks(
                    query_dna_masks(stats, query.seq().len()),
                    query.seq().len(),
                    query_placements.offset(q_idx),
                )
                .into_iter()
                .map(|(frame, ranges)| FrameMasks::new(frame, ranges))
                .collect()
            })
            .collect();
        Self {
            queries,
            subjects,
            subject_letters,
            query_masks,
            query_code,
            db_code,
            query_placements,
            subject_placements,
        }
    }

    // No NCBI counterpart for the attribute: on wasm32 (LOSAT Web) the row of an HSP stays a
    // function of its own, called once per HSP, so that V8 runs the later HSPs in its optimized
    // code (it switches a function to that code only at a later call, never inside a running
    // one); inlined into the map, the rows of a whole list (or of a pool task's part of it) were
    // one call in V8's baseline code. Other targets keep the compiler's choice.
    /// The outfmt 0 row of `hsp` (its displayed rows, frames, counts, linked-set size and
    /// subject), whose Seq-align's `num_ident` is `seqalign_num_ident` (`seqalign_identities`)
    /// when the report reads it; without it the row keeps the HSP's identities there, which
    /// nothing reads. `window` is the buffer of a stored subject's letters, `row_masks` the
    /// masks of the rows (`query_row_masks`), both kept from HSP to HSP: after the call it holds
    /// this HSP's (the query row's, set at every HSP; the subject row has none), which `-html`
    /// copies (`HspMasks::copy_rows`).
    #[cfg_attr(target_arch = "wasm32", inline(never))]
    pub(crate) fn hit(
        &self,
        hsp: &TblastxHsp,
        seqalign_num_ident: Option<u32>,
        window: &mut Vec<u8>,
        row_masks: &mut crate::report::query_anchored::HspMasks,
    ) -> PairwiseHit {
        let (queries, subjects, subject_letters) =
            (self.queries, self.subjects, self.subject_letters);
        let (query_code, db_code) = (self.query_code, self.db_code);
        let (query_placements, subject_placements) =
            (self.query_placements, self.subject_placements);
        let query = &queries[hsp.hit.q_idx as usize];
        let subject = &subjects[hsp.hit.s_idx as usize];
        let (mut query_row, subject_row) =
            subject_letters.rows(hsp, query.seq(), query_code, db_code, window);
        let (identities, positives) = row_counts(Blosum62Score, &query_row, &subject_row);
        let (q_idx, s_idx) = (hsp.hit.q_idx as usize, hsp.hit.s_idx as usize);
        let mut hit = hsp.hit.clone();
        shift_hit(&mut hit, query_placements, subject_placements);
        let query_length = query_placements.length(q_idx, query.seq().len());
        let subject_length = subject_placements.length(s_idx, subject_letters.length(s_idx));
        let (low, high) = (
            hit.q_start.min(hit.q_end) as i32 - 1,
            hit.q_start.max(hit.q_end) as i32 - 1,
        );
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:430-440
        // The frames that the report prints come from the records' coordinates and
        // lengths (`seq_range::record_frame`); without ranges they are the search's.
        let query_frame = record_frame(hit.query_frame > 0, hit.q_start - 1, query_length) as i8;
        let subject_frame =
            record_frame(hsp.subject_frame > 0, hit.s_start - 1, subject_length) as i8;
        // The row holds its masks in lower case; the lists are left to the next row
        // (`query_anchored::output_seq`).
        query_row_masks(
            &mut query_row,
            query_frame,
            low,
            high,
            &self.query_masks[q_idx],
            row_masks,
        );
        hit.num_ident = identities;
        hit.num_positives = positives;
        PairwiseHit {
            hit,
            query_seq: Some(row_string(query_row)),
            subject_seq: Some(row_string(subject_row)),
            query_frame: Some(query_frame),
            subject_frame: Some(subject_frame),
            positives: Some(positives),
            gaps: Some(0),
            subject_length: Some(subject_length),
            subject_title: crate::blastinput::fasta_reader::hit_subject_title(&subject.title),
            comp_adjust_method: None,
            sum_n: Some(hsp.num),
            subject_stored_gaps: Vec::new(),
            query_stored_gaps: Vec::new(),
            query_case: Default::default(),
            seqalign_num_ident: seqalign_num_ident.map_or(hsp.hit.num_ident, |n| n as usize),
            edge: crate::report::pairwise::PairwiseEdge::Normal,
        }
    }
}

// No NCBI counterpart: a row becomes its `String` in its own buffer, with no copy and no
// second allocation (see `displayed_counts` for why that counts on LOSAT Web).
/// A displayed row as a `String`: valid UTF-8 (the rows are ASCII letters) keeps the row's
/// buffer; a row with an invalid byte is the text of `String::from_utf8_lossy`.
fn row_string(row: Vec<u8>) -> String {
    String::from_utf8(row)
        .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned())
}

/// Moves a hit's reported coordinates from the searched letters into the records: by the
/// start of the query's interval (`RemapToQueryLoc`) and of the subject's interval
/// (`s_RemapToSubjectLoc`), in nucleotides, on either strand.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
/// A subject interval has the strand both, so `RemapAlignToLoc` (seq_align_util.cpp:71-78)
/// moves its row by the interval's start.
fn shift_hit(hit: &mut crate::common::Hit, queries: &Placements, subjects: &Placements) {
    let q_shift = queries.offset(hit.q_idx as usize);
    let s_shift = subjects.offset(hit.s_idx as usize);
    hit.q_start += q_shift;
    hit.q_end += q_shift;
    hit.s_start += s_shift;
    hit.s_end += s_shift;
}

/// `shift_hit` for every HSP of the final hit list (the tabular rows). Without ranges
/// nothing moves.
pub(crate) fn shift_to_records(
    hits: &mut [TblastxHsp],
    queries: &Placements,
    subjects: &Placements,
) {
    if !queries.ranged() && !subjects.ranged() {
        return;
    }
    for hsp in hits {
        shift_hit(&mut hsp.hit, queries, subjects);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1275-1283
/// The block header of one query in outfmt 7: the program, the query, the subjects, and,
/// for a searched query (`num_hits` is `Some`), the fields when there are hits and the
/// count of HSPs.
pub(crate) fn write_outfmt7_query_header<W: Write>(
    writer: &mut W,
    query_title: &[u8],
    database: &[u8],
    num_hits: Option<usize>,
) -> io::Result<()> {
    writeln!(writer, "# TBLASTX {}", crate::wording::REPORT_VERSION)?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1305-1308
    // The title's bytes (`write_outfmt7_query_line`).
    crate::report::outfmt6::write_outfmt7_query_line(writer, query_title)?;
    writer.write_all(b"# Database: ")?;
    writer.write_all(database)?;
    writer.write_all(b"\n")?;
    let Some(num_hits) = num_hits else {
        return Ok(());
    };
    if num_hits > 0 {
        writeln!(
            writer,
            "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
        )?;
    }
    writeln!(writer, "# {num_hits} hits found")
}

/// Writes the outfmt 6 or 7 rows of the final hit list (in its order); outfmt 7 adds the
/// block header of every query and, with the epilog (`Some` of the queries that earlier calls
/// wrote, `report_stream`), the count of queries.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tabular<W: crate::report::ostream::ReportWrite>(
    hits: &[TblastxHsp],
    writer: &mut W,
    comments: bool,
    query_titles: &[&[u8]],
    database: &[u8],
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    loads: &crate::report::db_subject::HitLoads,
    unsearched: &[bool],
    epilog: Option<usize>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = OutputConfig::ncbi_compat();
    let mut by_query: Vec<Vec<(HspIndex, &TblastxHsp)>> = vec![Vec::new(); query_titles.len()];
    for (hsp_index, hsp) in hits.iter().enumerate() {
        by_query[hsp.hit.q_idx as usize].push((hsp_index, hsp));
    }
    for (q_idx, query_hits) in by_query.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        if comments {
            write_outfmt7_query_header(
                writer,
                query_titles[q_idx],
                database,
                (!unsearched[q_idx]).then_some(query_hits.len()),
            )?;
        }
        let mut previous: Option<u32> = None;
        // The query's formatter is new: its rows keep no subject before the first one that
        // loads it (`TabularRows`).
        let mut rows = crate::report::tabular::TabularRows::default();
        for &(hsp_index, hsp) in query_hits {
            // The subject's load throws where the query's batch has its prefetch off
            // (`report::db_subject::HitLoads`): each row is written as the catch leaves it
            // (`TabularRows::write_thrown_default`).
            let thrown = loads.thrown(hsp.hit.s_idx as usize).is_some();
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
            // A subject's first row loads it where the query's batch did not prefetch it
            // (`QueryWarnings::subject_read`; a subject's HSPs follow each other).
            if previous != Some(hsp.hit.s_idx) {
                previous = Some(hsp.hit.s_idx);
                if let (Some(warnings), false) = (warnings.as_deref_mut(), thrown) {
                    warnings.subject_read(q_idx, hsp.hit.s_idx, &mut *writer)?;
                }
            }
            // One printed row is one HSP; the probe marks it without changing it.
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.begin(hsp_index);
            }
            if thrown {
                let query_id = query_ids
                    .get(hsp.hit.q_idx as usize)
                    .map_or(&b""[..], |id| id.as_ref());
                rows.write_thrown_default(
                    writer,
                    query_id,
                    |kept| crate::algorithm::common::tabular::printed_id(subject_ids, kept),
                    hsp.hit.e_value,
                    hsp.hit.bit_score,
                    &config,
                )?;
            } else {
                rows.found(hsp.hit.s_idx as usize);
                write_tabular_row(writer, hsp, query_ids, subject_ids, &config)?;
            }
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.end(hsp_index);
            }
        }
        // `~CBlastTabularInfo` of the query's report flushes the stream.
        writer.ncbi_flush()?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2238
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1322-1325
    if let (true, Some(queries_before)) = (comments, epilog) {
        writeln!(
            writer,
            "# BLAST processed {} queries",
            queries_before + query_titles.len()
        )?;
        // The epilog's `CBlastTabularInfo` flushes the stream in its destructor.
        writer.ncbi_flush()?;
    }
    Ok(())
}

/// The outfmt 6 / 7 row of one HSP (`write_tabular`).
pub(crate) fn write_tabular_row<W: Write>(
    writer: &mut W,
    hsp: &TblastxHsp,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    config: &OutputConfig,
) -> io::Result<()> {
    let hit = &hsp.hit;
    // The query ID of the record (`FastaRecord::shown_id`) and the subject's
    // `GetSeqIdList` ID (`tabular_subject_id`), as bytes.
    let query_id = query_ids
        .get(hit.q_idx as usize)
        .map_or(&b""[..], |id| id.as_ref());
    // A failed lookup's rows are refused before (`check_failed_lookups`).
    let subject_id = crate::algorithm::common::tabular::printed_id(subject_ids, hit.s_idx as usize);
    write_hit_fields_bytes(
        writer,
        query_id,
        subject_id,
        hit.identity,
        hit.num_ident,
        hit.length,
        hit.mismatch,
        hit.gapopen,
        hit.q_start,
        hit.q_end,
        hit.s_start,
        hit.s_end,
        hit.e_value,
        hit.bit_score,
        config,
    )
}

/// The TBLASTX output formats (`-outfmt` 0 to 16 and 20).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TblastxOutputFormat {
    /// The pairwise report (outfmt 0) and the query-anchored formats (1 to 4), which are outfmt
    /// 0's report with another alignment section (`PairwiseView`).
    Pairwise(crate::report::pairwise::PairwiseView),
    Tabular,
    TabularWithComments,
    /// The tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
    /// default `std` fields with a tab) and outfmt 10 and 20 (`report::tabular`, through
    /// `algorithm::tblastx::tabular`). The default outfmt 6 and 7 keep `Tabular` and
    /// `TabularWithComments`, whose writers (`write_tabular`, the report pool) this report
    /// never enters.
    CustomTabular(&'static crate::report::tabular::TabularSpec),
    /// A Seq-annot of Seq-aligns (outfmt 8, 9, 12).
    Seqalign(crate::report::serial::SerialFormat),
    /// BLAST XML (outfmt 5, `tblastx/blast_xml.rs`).
    Xml,
    /// The XML2 and JSON documents (outfmt 16, 15, 14 and 13), written from the same records as
    /// outfmt 5 (`tblastx/blast_xml.rs`); 14 and 13 as side files with this output as the master.
    Xml2(crate::report::blast_xml2::Xml2Kind),
    /// The BLAST archive (outfmt 11, `tblastx/archive.rs`): one document per query batch.
    Archive,
    /// SAM (outfmt 17, `tblastx/sam.rs`), which only BLASTN's application writes: a TBLASTX search
    /// strategy that it runs (a TBLASTX run itself fails at its formatting handler, `app::
    /// formatting_handler_check`). The document is written at the end of the run.
    Sam(crate::report::sam::SamSpec),
    /// A value that TBLASTX does not write (a value that NCBI cannot read), which the
    /// command line, the library and LOSAT Web refuse before the search
    /// (`ReportFormat::writes`, `app::report_format`): its writer refuses it, no report is
    /// written for it.
    Unsupported,
}

impl TblastxOutputFormat {
    /// Whether the format is the pairwise report or a query-anchored format (outfmt 0 to 4).
    pub(crate) fn is_pairwise(self) -> bool {
        matches!(self, Self::Pairwise(_))
    }

    /// Whether the format is a table (outfmt 6, 7, 10 and 20), which reads the tabular subject
    /// ids (`report::tabular::writes_rows`: not the organism report).
    pub(crate) fn is_tabular(self) -> bool {
        match self {
            Self::Tabular | Self::TabularWithComments => true,
            Self::CustomTabular(spec) => spec.writes_rows(),
            _ => false,
        }
    }

    /// Whether the format is a table (`is_tabular`) whose fields read the sequence that NCBI's
    /// formatter finds by a subject's Seq-align id (`report::tabular::reads_subject`: the default
    /// fields do).
    pub(crate) fn tabular_reads_subject(self) -> bool {
        self.is_tabular()
            && crate::report::tabular::reads_subject(
                self.tabular_spec(),
                crate::blastinput::app::Program::Tblastx,
            )
    }

    /// The field list of the custom tabular report (`CustomTabular`).
    pub(crate) fn tabular_spec(self) -> Option<&'static crate::report::tabular::TabularSpec> {
        match self {
            Self::CustomTabular(spec) => Some(spec),
            _ => None,
        }
    }

    /// Whether the format writes the displayed identity counts of the rows
    /// (`set_displayed_identities`): not the structured formats, which read the engine's counts
    /// (outfmt 8, 9, 12) or count the rows again themselves (outfmt 5, 15, 16), and not the
    /// custom tabular reports, which make the rows they show themselves (`tabular`).
    pub(crate) fn reads_displayed_counts(self) -> bool {
        !matches!(
            self,
            Self::Seqalign(_)
                | Self::Xml
                | Self::Xml2(_)
                | Self::Archive
                | Self::Sam(_)
                | Self::CustomTabular(_)
                | Self::Unsupported
        )
    }

    /// Whether the format is BLAST XML (5) or one of its XML2 / JSON siblings (15, 16), which
    /// are made from the same records and refuse the same names.
    pub(crate) fn is_structured_xml(self) -> bool {
        matches!(self, Self::Xml | Self::Xml2(_))
    }

    /// Whether the format writes the subjects' taxids or their names: the XML2 and JSON
    /// documents and a table of taxonomy fields (`TabularNeeds::taxonomy`).
    pub(crate) fn writes_taxonomy(self) -> bool {
        match self {
            Self::Xml2(_) => true,
            Self::CustomTabular(spec) => spec.needs().taxonomy,
            _ => false,
        }
    }

    /// The subject ids that the format prints (`PrintedIds`): none but in a table.
    pub(crate) fn printed_ids(self) -> crate::report::tabular::PrintedIds {
        match self {
            Self::Tabular | Self::TabularWithComments => {
                crate::report::tabular::PrintedIds::of(None)
            }
            Self::CustomTabular(spec) => crate::report::tabular::PrintedIds::of(Some(spec)),
            _ => crate::report::tabular::PrintedIds::default(),
        }
    }

    /// Whether the format is a table that counts following subjects whose Seq-align ids match as
    /// one (`TabularNeeds::joins_subjects`: `qcovs`, `qcovus`).
    pub(crate) fn tabular_joins_subjects(self) -> bool {
        matches!(self, Self::CustomTabular(spec) if spec.needs().joins_subjects())
    }

    /// Whether NCBI prefetches the subjects that each query batch reports before the batch's
    /// reports (outfmt 0 to 7 and 13 to 18, 10 and 20 as 6 and 7; with
    /// `PRE_FETCH_SEQS_LIMIT=2147483647` also 8, 9 and 12, `database::prefetches_subjects`): as the
    /// call's first report, the posts of that prefetch go with the query warnings.
    pub(crate) fn prefetches_subjects(self) -> bool {
        // The archive fetches nothing; the XML2 and JSON `Hit`s write only the texts of their own
        // lookups (`database::RowLookups::Taxids`).
        crate::algorithm::common::database::prefetches_subjects(Some(self.prefetch_format()))
    }

    /// The format as the first format of a TBLASTX run that NCBI's prefetch of each query batch
    /// asks (`common::database::PrefetchFormat`, the one owner of whether a batch prefetches).
    pub(crate) fn prefetch_format(self) -> crate::algorithm::common::database::PrefetchFormat {
        use crate::blastinput::app::output_format::{ARCHIVE_FORMAT, END_VALUE, JSON, SAM, XML};
        crate::algorithm::common::database::PrefetchFormat::of(
            match self {
                // `s_IsPrefetchFormat` takes 1 to 4 as 0 (DECISIONS T3-4).
                Self::Pairwise(view) => view.outfmt(),
                Self::Xml => XML,
                Self::Tabular => 6,
                Self::TabularWithComments => 7,
                // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (6, 7, 10 and
                // 20): they read the subjects as 6 and 7 do, and show no taxonomy name that the
                // prefetch would look up (`report::tabular` writes `N/A`; a name of the taxonomy
                // database is LOSAT's refusal, `tabular`).
                Self::CustomTabular(spec) => match spec.kind {
                    crate::report::tabular::TabularKind::Commented => 7,
                    _ => 6,
                },
                Self::Seqalign(serial) => serial.outfmt(),
                Self::Sam(_) => SAM,
                // The XML2 and JSON documents prefetch as each of 13 to 16 does.
                Self::Xml2(_) => JSON,
                Self::Archive => ARCHIVE_FORMAT,
                // Refused before the search: no report, no prefetch.
                Self::Unsupported => END_VALUE,
            },
            crate::report::archive::options::RequestKind::Tblastx,
        )
    }
}

/// The format of a validated `-outfmt` value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
///
/// A tabular format with a field list or a delimiter, and every outfmt 10 and 20, is the custom
/// tabular report (`report::tabular::TabularSpec::of_app`); the default outfmt 6 and 7 (the `std`
/// fields with a tab) are `Tabular` and `TabularWithComments`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
pub(crate) fn output_format(outfmt: impl AsRef<[u8]>) -> TblastxOutputFormat {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt);
    // SAM reads its custom specification (`SQ`, `SR`) and ignores a delimiter.
    if let Ok(choice) = &choice {
        if choice.number == crate::blastinput::app::output_format::SAM {
            return TblastxOutputFormat::Sam(crate::report::sam::SamSpec::parse(&choice.spec));
        }
    }
    // Only a value with a custom specification, or outfmt 10 or 20, reads the field list;
    // outfmt 18 (the organism report) is written from the same hit list.
    if let Ok(choice) = &choice {
        if choice.custom() || matches!(choice.number, 10 | 18 | 20) {
            if let Some(spec) = crate::blastinput::app::ReportFormat::from_number(choice.number)
                .and_then(|format| {
                    crate::report::tabular::TabularSpec::of_app(
                        format,
                        choice,
                        crate::report::archive::options::RequestKind::Tblastx,
                    )
                })
            {
                return TblastxOutputFormat::CustomTabular(spec.intern());
            }
        }
    }
    let number = choice.map(|choice| choice.number);
    // Outfmt 0 to 4: the pairwise report and its query-anchored views (line RG T3c).
    if let Some(view) = number
        .as_ref()
        .ok()
        .and_then(|&number| crate::blastinput::app::ReportFormat::from_number(number))
        .and_then(crate::report::pairwise::PairwiseView::of)
    {
        return TblastxOutputFormat::Pairwise(view);
    }
    match number {
        Ok(6) => TblastxOutputFormat::Tabular,
        Ok(7) => TblastxOutputFormat::TabularWithComments,
        Ok(8) => TblastxOutputFormat::Seqalign(crate::report::serial::SerialFormat::AsnText),
        Ok(9) => TblastxOutputFormat::Seqalign(crate::report::serial::SerialFormat::AsnBinary),
        Ok(12) => TblastxOutputFormat::Seqalign(crate::report::serial::SerialFormat::Json),
        Ok(5) => TblastxOutputFormat::Xml,
        Ok(11) => TblastxOutputFormat::Archive,
        Ok(16) => TblastxOutputFormat::Xml2(crate::report::blast_xml2::Xml2Kind::Single16),
        Ok(15) => TblastxOutputFormat::Xml2(crate::report::blast_xml2::Xml2Kind::Single15),
        Ok(14) => TblastxOutputFormat::Xml2(crate::report::blast_xml2::Xml2Kind::Files14),
        Ok(13) => TblastxOutputFormat::Xml2(crate::report::blast_xml2::Xml2Kind::Files13),
        // Never a table by default: a value the callers refuse first (`Unsupported`).
        _ => TblastxOutputFormat::Unsupported,
    }
}

/// The filters of NCBI's `s_FilterBlastResults` that TBLASTX has: `-max_hsps` (`max_hsps`,
/// 0 when omitted) and `-qcov_hsp_perc` (`query_cov_hsp_perc`, 0 when omitted), with the
/// nucleotide lengths of the queries (`query_lengths`, read for the coverage only).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1521-1538
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1766
/// TBLASTX is ungapped, so the preliminary stage removes the HSPs below the query coverage
/// after the e-value reap, before the hit list (or the culling writer) receives a subject's
/// HSPs, and a subject left without HSPs never reaches it. After the traceback (and the
/// culling pipe) `s_FilterBlastResults` cuts each list to `-max_hsps` in its order (score
/// order, or e-value order once the hit list has filled), without setting its best e-value
/// again; its coverage reap then finds nothing to remove, because NCBI never splits the query
/// of an ungapped search (split_query_cxx.cpp:60-61), so both reaps read the same context
/// lengths.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct HitFilters<'a> {
    pub max_hsps: usize,
    pub query_cov_hsp_perc: f64,
    pub query_lengths: &'a [usize],
    pub subject_besthit: bool,
}

impl HitFilters<'_> {
    /// The lengths of the contexts (frames) of query `q_idx`, read only by the coverage.
    pub(super) fn context_lengths(&self, q_idx: u32) -> [i32; 6] {
        if self.query_cov_hsp_perc > 0.0 {
            super::hsp_culling::tblastx_context_lengths(self.query_lengths[q_idx as usize])
        } else {
            [0; 6]
        }
    }

    /// `Blast_HSPListReapByQueryCoverage` on one subject's HSPs; returns whether none is
    /// left (NCBI's collector and culling writer skip an empty list).
    pub(super) fn reap_by_query_coverage(
        &self,
        hsps: &mut Vec<TblastxHsp>,
        context_lengths: &[i32; 6],
    ) -> bool {
        if self.query_cov_hsp_perc <= 0.0 {
            return false;
        }
        crate::core::blast_hits::reap_by_query_coverage(hsps, self.query_cov_hsp_perc, |hsp| {
            let context = super::hsp_culling::tblastx_context(hsp.hit.query_frame);
            (
                hsp.hit.sort_query_offset as i64,
                hsp.hit.sort_query_end as i64,
                i64::from(context_lengths[context]),
            )
        });
        hsps.is_empty()
    }

    /// `-subject_besthit` on one final list of a query (`s_FilterBlastResults`, after `-max_hsps`
    /// and the coverage, which the preliminary reap already applied to these HSPs): the contexts
    /// are the query's frames. The best e-value is not set again, as NCBI's is not. Every HSP here
    /// has a positive length (the one-hit word finder's HSPs of length 0 or less are dropped
    /// before the re-evaluation, after its fetch of the subject, as NCBI's re-evaluation deletes
    /// them: `subject_link::drop_unread_int4_inverted`), so the `usize` offsets are NCBI's `Int4`
    /// ones.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:858-861
    pub(super) fn subject_besthit(&self, hsps: &mut Vec<TblastxHsp>) {
        if !self.subject_besthit {
            return;
        }
        crate::core::blast_hits::subject_besthit(hsps, false, |hsp| {
            crate::core::blast_hits::BestHitExtent {
                context: super::hsp_culling::tblastx_context(hsp.hit.query_frame) as i64,
                query_frame: hsp.hit.query_frame,
                offset: hsp.hit.sort_query_offset as i64,
                end: hsp.hit.sort_query_end as i64,
                query_length: 0,
            }
        });
    }
}

/// The final HSP list of a run in NCBI's order, which every format prints: per query (in
/// input order) the subjects that the hit list keeps, best first, and each subject's
/// HSPs in the order of its Seq-align.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1803
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:143-150
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1772
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
/// The subjects reach the hit list in oid order, each with its HSPs in score order (the
/// state after linking). A query with more subjects than the hit list size sorts every list
/// by e-value when the list fills (`Blast_HitListUpdate`), so the comparison of its subjects
/// reads the first HSP in e-value order; the others read the highest-scoring HSP. The
/// preliminary hit list of tblastx has the final size (`GetPrelimHitlistSize` without
/// gapping or composition-based statistics, blast_hits.c:44-71), and the traceback of an
/// ungapped search keeps it (`s_BlastPruneExtraHits` removes nothing).
///
/// `filters` (`-qcov_hsp_perc`, `-max_hsps`) apply where NCBI applies them
/// (`HitFilters`).
///
/// The search of a query batch runs the same steps on its queries' hit lists as the subjects'
/// lists come (`hit_collector::QueryLists`, in OID order); this is the reference of the tests,
/// on the HSPs of a run at once.
#[cfg(test)]
pub(crate) fn final_hit_order(
    hits: Vec<TblastxHsp>,
    hitlist_size: usize,
    filters: &HitFilters<'_>,
) -> Vec<TblastxHsp> {
    let mut pairs: std::collections::BTreeMap<(u32, u32), Vec<TblastxHsp>> =
        std::collections::BTreeMap::new();
    for hsp in hits {
        pairs
            .entry((hsp.hit.q_idx, hsp.hit.s_idx))
            .or_default()
            .push(hsp);
    }
    let mut lists = QueryLists::new(hitlist_size);
    for ((q_idx, oid), hsps) in pairs {
        if let Some(hsps) = super::hit_collector::subject_list(q_idx, hsps, filters) {
            lists.add(q_idx, oid, hsps);
        }
    }
    lists.finish(filters, true)
}

/// The final HSP list of a run with `-culling_limit` (`culling_limit` > 0) in NCBI's order:
/// `final_hit_order` with NCBI's culling writer in the preliminary stage and its culling
/// pipe after the traceback (`hsp_culling.rs`). `query_lengths` are the queries'
/// nucleotide lengths (the lengths of their frames are the culling trees' ranges).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1341
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:330-341
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:386-394
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:699-730
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1720-1721
/// The culling writer replaces the preliminary hit list (the collector): each subject's
/// HSPs (after the `-evalue` reap, in score order) go into the trees, subject by subject
/// in OID order, with the merit `culling_max` (+3 above 1; in `i64`, where NCBI's `Int4`
/// overflows near `INT4_MAX`, NCBI parity DECISIONS 17). The trees' HSPs then reach the hit list of the traceback in OID order, where
/// the `-max_target_seqs` cut is made (`Blast_HitListUpdate`), and the pipe culls the kept
/// lists again, sorted by e-value, with the merit `culling_limit`. The pipe's lists come out in
/// the rip order of its trees (hspfilter_culling.c:501-580); they are sorted by e-value only where
/// NCBI sorts the results (`sorts_results`: a positive `BlastSeqSrcGetTotLen`,
/// blast_traceback.c:1768-1772), then cut to `-max_target_seqs` (blast_traceback.c:1777), and each
/// list's HSPs sorted by e-value. The trees belong to the query's contexts, so the queries are
/// culled one after another.
pub(crate) fn culled_hit_order(
    hits: Vec<TblastxHsp>,
    hitlist_size: usize,
    culling_limit: i32,
    query_lengths: &[usize],
    filters: &HitFilters<'_>,
    sorts_results: bool,
) -> Vec<TblastxHsp> {
    use super::hsp_culling::{tblastx_context_lengths, CullingWriter};
    use crate::algorithm::blastn::hsp::{HitList, HitListEntry};
    use crate::core::hspfilter_culling::CullingForest;
    type Subjects = std::collections::BTreeMap<u32, Vec<TblastxHsp>>;
    let mut lists: std::collections::BTreeMap<u32, Subjects> = std::collections::BTreeMap::new();
    for hsp in hits {
        lists
            .entry(hsp.hit.q_idx)
            .or_default()
            .entry(hsp.hit.s_idx)
            .or_default()
            .push(hsp);
    }
    let prelim_max = CullingForest::<TblastxHsp>::prelim_max(culling_limit);
    let mut ordered = Vec::new();
    for (q_idx, subjects) in lists {
        let context_lengths = tblastx_context_lengths(query_lengths[q_idx as usize]);
        let mut writer = CullingWriter::new(prelim_max, context_lengths);
        for (oid, mut hsps) in subjects {
            hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
            if filters.reap_by_query_coverage(&mut hsps, &context_lengths) {
                continue;
            }
            writer.run(hsps, oid);
        }
        let mut culled = writer.finalize();
        culled.sort_by_key(|(oid, _)| *oid);
        let mut hit_list: HitList<TblastxHspList> = HitList::new(hitlist_size);
        for (oid, hsps) in culled {
            hit_list.update(TblastxHspList {
                oid,
                hsps,
                best_evalue: 0.0,
            });
        }
        let mut kept = hit_list.hsplist_array;
        for list in &mut kept {
            HitListEntry::sort_by_evalue(list);
            list.best_evalue = list.hsps[0].hit.e_value;
        }
        let mut pipe_lists: HitList<TblastxHspList> = HitList::new(hitlist_size);
        pipe_lists.hsplist_count = kept.len();
        pipe_lists.hsplist_array = kept;
        pipe_lists.sort_by_evalue();
        let mut pipe = CullingWriter::new(i64::from(culling_limit), context_lengths);
        for list in pipe_lists.hsplist_array {
            pipe.run(list.hsps, list.oid);
        }
        let mut final_lists: HitList<TblastxHspList> = HitList::new(hitlist_size);
        for (oid, hsps) in pipe.finalize() {
            let mut list = TblastxHspList {
                oid,
                hsps,
                best_evalue: 0.0,
            };
            list.update_best_evalue();
            trim_by_max_hsps(&mut list.hsps, filters.max_hsps);
            filters.subject_besthit(&mut list.hsps);
            final_lists.hsplist_array.push(list);
        }
        final_lists.hsplist_count = final_lists.hsplist_array.len();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1772 (cited at
        // `DbStatistics::sorts_results`)
        if sorts_results {
            final_lists.sort_by_evalue();
        }
        final_lists.prune_by_size(hitlist_size);
        for mut list in final_lists.hsplist_array {
            HitListEntry::sort_by_evalue(&mut list);
            ordered.extend(list.hsps);
        }
    }
    ordered
}

#[cfg(test)]
#[path = "report_tests.rs"]
mod tests;

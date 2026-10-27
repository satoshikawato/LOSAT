//! The TBLASTX reports from the final HSP list of a run: the displayed rows of every HSP
//! (outfmt 0, the hit records, and the identity counts of outfmt 6 and 7), the per-query
//! blocks of outfmt 7, and the data of the pairwise report.

use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

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
use super::scoring::{Blosum62Score, TbxScore};

/// The version that NCBI's TBLASTX reports print.
///
/// NCBI reference: c++/src/algo/blast/format/blastfmtutil.cpp:60-64
/// ```c++
/// string CBlastFormatUtil::BlastGetVersion(const string program)
/// {
///     string program_uc = program;
///     return NStr::ToUpper(program_uc) + " " + blast::CBlastVersion().Print();
/// }
/// ```
pub(crate) const NCBI_TBLASTX_VERSION: &str = "2.17.0+";

/// The 4-bit IUPAC code of a nucleotide as the display's translation table reads it
/// (`A` 1, `C` 2, `G` 4, `T` 8, ambiguity letters their union, `U` as `T`, `X` as `N`);
/// another byte is a gap (0), which translates to `X`.
///
/// NCBI reference: c++/src/objects/seqfeat/Genetic_code_table.cpp:94,97-112
/// ```c++
///     static char  charToBase [17] = "-ACMGRSVTWYHKDBN";
///     ...
///     // illegal characters map to 0
///     for (i = 0; i < 256; i++) {
///         sm_BaseToIdx [i] = 0;
///     }
///
///     // map iupacna alphabet to EBaseCode
///     for (i = eBase_gap; i <= eBase_N; i++) {
///         ch = charToBase [i];
///         sm_BaseToIdx [(int) ch] = i;
///         ch = (unsigned char)tolower (ch);
///         sm_BaseToIdx [(int) ch] = i;
///     }
///     sm_BaseToIdx [(int) 'U'] = eBase_T;
///     sm_BaseToIdx [(int) 'u'] = eBase_T;
///     sm_BaseToIdx [(int) 'X'] = eBase_N;
///     sm_BaseToIdx [(int) 'x'] = eBase_N;
/// ```
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
// display table (`LOSAT_X_DISPTABLE`). NCBI's `sm_BaseToIdx` (above) is the plus-strand table.
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

// No NCBI counterpart: reads the LOSAT_X_DISPTABLE variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_DISPTABLE): a switch of the strict set (`utils/xswitch.rs`): on unless it
/// is `0`.
fn x_display_table() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DISPTABLE"))
}

// NCBI reference: c++/src/objects/seqfeat/Genetic_code_table.cpp:155-160,228-231 (quoted at
// `blastx::report::display_codon`): NCBI's display reads the residue of a codon's state from a
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
            crate::algorithm::blastx::report::display_codon(masks, code)
        }))
    })
}

/// A genetic code as the display translates with it (`display_residue`): with
/// `LOSAT_X_DISPTABLE` the code's table (`display_table`), looked up once for many residues.
#[derive(Clone, Copy)]
pub(crate) struct DisplayCode<'a> {
    code: &'a GeneticCode,
    table: Option<&'static [u8; 4096]>,
}

impl<'a> DisplayCode<'a> {
    pub(crate) fn new(code: &'a GeneticCode) -> Self {
        Self {
            code,
            table: x_display_table().then(|| display_table(code)),
        }
    }

    /// `display_residue` with this code.
    pub(crate) fn residue(&self, sequence: &[u8], frame: i8, offset: usize) -> u8 {
        let Some(table) = self.table else {
            return display_residue_scan(sequence, frame, offset, self.code);
        };
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
}

// NCBI reference: c++/src/objtools/alnmgr/alnvec.hpp:292-302
// ```c++
//     if (GetWidth(row) == 3) {
//         string buff;
//         buffer.erase();
//         if (IsPositiveStrand(row)) {
//             x_GetSeqVector(row).GetSeqData(seq_from, seq_to + 1, buff);
//         } else {
//             CSeqVector& seq_vec = x_GetSeqVector(row);
//             TSeqPos size = seq_vec.size();
//             seq_vec.GetSeqData(size - seq_to - 1, size - seq_from, buff);
//         }
//         TranslateNAToAA(buff, buffer, GetGenCode(row));
// ```
// NCBI reference: c++/src/objtools/alnmgr/alnvec.cpp:116-121
// ```c++
//         CBioseq_Handle h = GetBioseqHandle(row);
//         CSeqVector vec = h.GetSeqVector
//             (CBioseq_Handle::eCoding_Iupac,
//              IsPositiveStrand(row) ?
//              CBioseq_Handle::eStrand_Plus :
//              CBioseq_Handle::eStrand_Minus);
// ```
// NCBI reference: c++/src/objtools/alnmgr/alnvec.cpp:911-918
// ```c++
//     int state = 0;
//     size_t aa_i = 0;
//     for (size_t na_i = 0; na_i < na_size; ) {
//         for (size_t i = 0; i < 3; i++) {
//             state = tbl.NextCodonState(state, na[na_i++]);
//         }
//         aa[aa_i++] = tbl.GetCodonResidue(state);
//     }
// ```
/// The displayed residue of codon `offset` of `frame` (+1..+3, -1..-3) of `sequence`: the
/// codon's bases read on the frame's strand (the IUPAC complement on the minus strand) and
/// translated with the display's table (`display_codon`, which merges D/N, E/Q and I/L).
pub(crate) fn display_residue(sequence: &[u8], frame: i8, offset: usize, code: &GeneticCode) -> u8 {
    DisplayCode::new(code).residue(sequence, frame, offset)
}

/// `display_residue` without the display table (`LOSAT_X_DISPTABLE=0`): each base's code found
/// in the IUPAC letters and the codon translated by `display_codon`.
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
    crate::algorithm::blastx::report::display_codon(masks, code)
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

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1948-1953
// ```c++
//             CRef<CAlnVec> avRef = x_GetAlnVecForSeqalign(**iter);
//     ...
//                 avRef->SetGenCode(m_SlaveGeneticCode);
//                 avRef->SetGenCode(m_MasterGeneticCode, 0);
// ```
// NCBI reference: c++/src/objtools/align_format/tabular.cpp:982-986
// ```c++
//         alnVec->SetGapChar('-');
//         alnVec->SetGenCode(m_QueryGeneticCode, 0);
//         alnVec->SetGenCode(m_DbGeneticCode, 1);
//         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
//         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
// ```
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
    /// The reported sequences of a `-db` search by rank (`ordinals`, the ordinal id of each
    /// rank, into `stores`), or (LOSAT_X_TBXSTORE) the packed records of `-subject` (ordinals
    /// `0..n`): a row reads the letters of its HSP's codons only, made from the stored codes
    /// (`codon_window`).
    Stored {
        stores: &'a [PackedNucleotide],
        ordinals: &'a [u32],
    },
}

impl<'a> SubjectLetters<'a> {
    /// The number of bases of subject `s_idx`.
    pub(crate) fn length(&self, s_idx: usize) -> usize {
        match self {
            Self::Records(records) => records[s_idx].seq().len(),
            Self::Stored { stores, ordinals } => stores[ordinals[s_idx] as usize].len(),
        }
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
                let store = &stores[ordinals[s_idx] as usize];
                let (from, to) = (hsp.hit.sort_subject_offset, hsp.hit.sort_subject_end);
                let (start, length) = codon_window(store.len(), hsp.subject_frame, from, to);
                store.letters_into(start, start + length, window);
                (window, from)
            }
        }
    }

    /// The engine counts of `hsp` (`hsp_engine_counts`) against its subject; `window` is the
    /// buffer of a stored subject's letters.
    fn engine_counts<S: TbxScore>(
        &self,
        hsp: &TblastxHsp,
        query: &[u8],
        query_code: &GeneticCode,
        db_code: &GeneticCode,
        window: &mut Vec<u8>,
        sc: S,
    ) -> (usize, usize) {
        let (subject, subject_first) = self.letters(hsp, window);
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

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2131-2150
// ```c++
//         if(sequence_standard[i]==sequence[i]){
//     ...
//             match ++;
//         } else {
//             if ((m_AlignType&eProt)
//                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//                 positive ++;
// ```
// NCBI reference: c++/src/objtools/align_format/tabular.cpp:1003-1021
// ```c++
//             for (unsigned int i = 0;
//                  i < min(m_QuerySeq.size(), m_SubjectSeq.size());
//                  ++i) {
//                 if (m_QuerySeq[i] == m_SubjectSeq[i]) {
//                     ++num_ident;
//                     ++num_positives;
//                     ++num_matches;
//                 } else {
//     ...
//                     if (matrix && !matrix->GetData().empty() &&
//                            (*matrix)(m_QuerySeq[i], m_SubjectSeq[i]) > 0) {
//                         ++num_positives;
//                     }
// ```
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
/// NCBI reference: c++/src/objtools/align_format/tabular.cpp:1460-1465
/// ```c++
///     case ePercentIdentical:
///         x_PrintPercentIdentical(); break;
///     case eNumIdentical:
///         x_PrintNumIdentical(); break;
///     case eMismatches:
///         x_PrintMismatches(); break;
/// ```
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
        hit.mismatch = hit.length - identities;
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
/// ```c
///       if (!delete_hsp) {
///           const Uint1* query_nomask = query_blk->sequence_nomask +
///               query_info->contexts[context].query_offset;
///           Int4 align_length = 0;
///           Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
///                    							    subject_start,
///                                					hsp,
///                                					score_params->options,
///                                					&align_length,
///                                					sbp);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:770-792
/// ```c
///    if(NULL != sbp)
///    {
/// 	   if(sbp->protein_alphabet)
/// 		   matrix = sbp->matrix->data;
///    }
///
///    if (!hsp->gap_info) {
///       /* Ungapped case. Check that lengths are the same in query and subject,
///          then count number of matches. */
///       if (q_length != s_length)
///          return -1;
///       align_length = q_length;
///       for (i=0; i<align_length; i++) {
///          if (*q == *s)
///             num_ident++;
///          else if (NULL != matrix) {
///         	 if (matrix[*q][*s] > 0)
///         		 num_pos ++;
///             }
///          q++;
///          s++;
///       }
///    	}
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:831-832
/// ```c
///    if(NULL != matrix)
/// 	   *num_pos_ptr = num_pos + num_ident;
/// ```
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
/// (`hsp_engine_counts`: the same letters give the same identities, checked in debug builds).
/// Made only when a Seq-align format is requested.
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
/// NCBI reference: c++/src/algo/blast/core/blast_filter.c:1256
/// ```c
///     BlastSeqLocCombine(filter_out, 0);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_filter.c:908-931
/// ```c
///    for (index=0; index < query_info->num_queries; ++index)
///    {
///        Int4 frame_start = index*NUM_FRAMES;
///        Int4 frame_index;
///        Int4 dna_length = BlastQueryInfoGetQueryLength(query_info,
///                                                       eBlastTypeBlastx,
///                                                       index);
///    ...
///                if (frame < 0) {
///                    to = dna_length - CODON_LENGTH*seq_range->left + frame;
///                    from = dna_length - CODON_LENGTH*seq_range->right + frame + 1;
///                } else {
///                    from = CODON_LENGTH*seq_range->left + frame - 1;
///                    to = CODON_LENGTH*seq_range->right + frame - 1;
///                }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_query_info.c:143-160
/// ```c
///     Int4 start_context = NUM_FRAMES*query_index;
///     Int4 dna_length = 2;
///     ...
///     if (query_info->contexts[start_context].query_length == 0)
///         start_context += 3;
///
///     for (index = start_context; index < start_context + 3; ++index)
///         dna_length += query_info->contexts[index].query_length;
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_util.c:923-929
/// ```c
/// BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
/// {
///     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
///         return 0;
///     }
///     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
/// }
/// ```
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
            let combined = crate::algorithm::blastx::query_setup::combine_masks(ranges);
            (
                *frame,
                crate::algorithm::blastx::query_setup::protein_to_dna_masks(
                    &combined, dna_length, *frame,
                ),
            )
        })
        .collect()
}

// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:936-951
// ```c++
//                 TSeqRange masked_range(loc->ssr->left, loc->ssr->right);
//                 TSeqRange range(Map(kTarget, masked_range));
//                 if (range.NotEmpty() && range != kTarget) {
//                     int frame = BLAST_ContextToFrame(program, index);
//     ...
//                     CRef<CSeqLocInfo> seqloc_info
//                         (new CSeqLocInfo(seqint, frame));
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:4305-4326
// ```c++
//             if(interval.GetId().Match(m_AV->GetSeqId(i)) &&
//                m_AV->GetSeqRange(i).IntersectingWith(loc_range)){
//                 int actualAlnStart = 0, actualAlnStop = 0;
//                 if(m_AV->IsPositiveStrand(i)){
//                     actualAlnStart =
//                         m_AV->GetAlnPosFromSeqPos(i,
//                                                   interval.GetFrom(),
//                                                           CAlnMap::eBackwards, true);
//                     actualAlnStop =
//                         m_AV->GetAlnPosFromSeqPos(i,
//                                                   interval.GetTo(),
//                                                   CAlnMap::eBackwards, true);
//                 } else {
//                     actualAlnStart =
//                         m_AV->GetAlnPosFromSeqPos(i,
//                                                   interval.GetTo(),
//                                                   CAlnMap::eBackwards, true);
//                     actualAlnStop =
//                         m_AV->GetAlnPosFromSeqPos(i,
//                                                   interval.GetFrom(),
//                                                   CAlnMap::eBackwards, true);
//                 }
// ```
// NCBI reference: c++/src/objtools/alnmgr/alnmap.cpp:551-558,593-595
// ```c++
//         if ((plus ? seq_pos < start : seq_pos > stop)) {
//             return GetAlnStart(seg.GetAlnSeg());
//         }
//         if ((plus ? seq_pos > stop : seq_pos < start)) {
//             return GetAlnStop(seg.GetAlnSeg());
//         }
//     ...
//     TSeqPos delta = (seq_pos - start) / GetWidth(row);
//     return m_AlnStarts[seg.GetAlnSeg()]
//         + (plus ? delta : m_Lens[raw_seg] - 1 - delta);
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2500-2521
// ```c++
//             int locFrame = (*iter)->seqloc->GetFrame();
//             if(id.Match((*iter)->seqloc->GetInterval().GetId())
//                && locFrame == frame){
//     ...
//                     } else if (m_SeqLocChar==eLowerCase){
//                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
// ```
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
    let target = (offset as i32, (offset + query_length) as i32 - 1);
    masks
        .into_iter()
        .map(|(frame, ranges)| {
            let shown = ranges
                .into_iter()
                .filter_map(|(from, to)| {
                    if from > to || from > target.1 || from + target.0 > target.1 {
                        return None;
                    }
                    let mapped = (target.0 + from, (target.0 + to).min(target.1));
                    (mapped != target).then_some(mapped)
                })
                .collect();
            (frame, shown)
        })
        .collect()
}

/// The masked regions of a query as NCBI hands them to the BLAST XML writer
/// (`Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959): the SEG intervals of each frame in
/// nucleotides, joined as `BlastSeqLocCombine(filter_out, 0)` joins them (only overlapping ones,
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

/// Lowercases the residues of a displayed query row that a shown SEG mask (record
/// coordinates, `shown_query_masks`) labelled with the row's frame covers. `[low, high]` is
/// the 0-based nucleotide range of the row on the plus strand of the record, and `frame`
/// the row's frame in the record (`s_GetStdsegMasterFrame`). With `-query_loc` a mask's
/// label is the frame of the searched letters, which can differ from the frame the record
/// gives the same letters: NCBI then shows a mask of another frame on the row (`locFrame ==
/// frame`), and so does LOSAT.
///
/// With `seek` (`LOSAT_X_TBXMASKSEEK`) the scan of a frame's ranges starts at the first range
/// whose end, or the end of a range before it, reaches `low`, and stops at the first range
/// that starts after `high`; the ranges it leaves out are the ones the full scan skips.
fn lowercase_query_row(
    row: &mut [u8],
    frame: i8,
    low: i32,
    high: i32,
    masks: &[FrameMasks],
    seek: bool,
) {
    let columns = row.len() as i32;
    let column = |nt: i32| -> i32 {
        let column = if frame > 0 {
            if nt < low {
                0
            } else if nt > high {
                columns - 1
            } else {
                (nt - low) / 3
            }
        } else if nt > high {
            0
        } else if nt < low {
            columns - 1
        } else {
            columns - 1 - (nt - low) / 3
        };
        column.clamp(0, columns - 1)
    };
    for frame_masks in masks.iter().filter(|masks| masks.frame == frame) {
        let first_range = if seek {
            frame_masks.ends_max.partition_point(|&end| end < low)
        } else {
            0
        };
        for &(from, to) in &frame_masks.ranges[first_range..] {
            if seek && from > high {
                break;
            }
            if to < low || from > high {
                continue;
            }
            let (first, last) = if frame > 0 {
                (column(from), column(to))
            } else {
                (column(to), column(from))
            };
            for residue in &mut row[first as usize..=last as usize] {
                *residue = residue.to_ascii_lowercase();
            }
        }
    }
}

/// Every HSP of the final hit list (in the order of all formats) with its displayed rows,
/// frames, counts, linked-set size and subject.
// Kept out of line: this runs only for outfmt 0 or hit records.
#[inline(never)]
///
/// With `-query_loc` and `-subject_loc` the rows are read from the searched letters, and
/// the hits move into the records (`shift_to_records`) with the records' frames and
/// lengths. The subjects' titles come from `subjects`, their letters and lengths from
/// `subject_letters`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn pairwise_hits(
    hits: &[TblastxHsp],
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    query_stats: &[TblastxQueryStats],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    query_placements: &Placements,
    subject_placements: &Placements,
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
    let mut window = Vec::new();
    hits.iter().map(|hsp| rows.hit(hsp, &mut window)).collect()
}

/// The positives of the outfmt 0 rows `hits` by `matrix`, a matrix other than BLOSUM62, counted
/// again from the rows (`row_counts`) after `PairwiseRows::hit` counted BLOSUM62's: the query
/// row's masked letters are lower case there, and the letters of the rows are upper case. One
/// `hit` for every matrix keeps the helpers that it calls once (`SubjectLetters::rows`,
/// `lowercase_query_row`) in its body.
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

// No NCBI counterpart: the masks of `shown_query_masks` arranged for `lowercase_query_row`. The
// lowercase letters of a row do not depend on the order of the ranges, so a frame's ranges are
// sorted by start, and the running maximum of their ends lets the row's scan start at the first
// range that can reach it (`LOSAT_X_TBXMASKSEEK`). Before, every row scanned all ranges of its
// frame: in TBLASTX of two ~650 kb genomes that scan was most of the row's time.
/// The shown SEG masks of one frame of a query, sorted by start.
struct FrameMasks {
    frame: i8,
    ranges: Vec<(i32, i32)>,
    /// `ends_max[i]`: the largest end of `ranges[..=i]` (non-decreasing).
    ends_max: Vec<i32>,
}

impl FrameMasks {
    fn new(frame: i8, mut ranges: Vec<(i32, i32)>) -> Self {
        ranges.sort_unstable();
        let ends_max = ranges
            .iter()
            .scan(i32::MIN, |max, &(_, to)| {
                *max = (*max).max(to);
                Some(*max)
            })
            .collect();
        Self {
            frame,
            ranges,
            ends_max,
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
    mask_seek: bool,
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
            mask_seek: crate::utils::xswitch::on("LOSAT_X_TBXMASKSEEK"),
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
    /// subject); `window` is the buffer of a stored subject's letters.
    #[cfg_attr(target_arch = "wasm32", inline(never))]
    pub(crate) fn hit(&self, hsp: &TblastxHsp, window: &mut Vec<u8>) -> PairwiseHit {
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
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:430-440
        // ```c++
        // static int s_GetStdsegMasterFrame(const CStd_seg& ss, CScope& scope)
        // {
        //     const CRef<CSeq_loc> slc = ss.GetLoc().front();
        //     ENa_strand strand = GetStrand(*slc);
        //     int frame = s_GetFrame(strand ==  eNa_strand_plus ?
        //                            GetStart(*slc, &scope) : GetStop(*slc, &scope),
        //                            strand ==  eNa_strand_plus ?
        //                            eNa_strand_plus : eNa_strand_minus,
        //                            *(ss.GetIds().front()), scope);
        //     return frame;
        // }
        // ```
        // The frames that the report prints come from the records' coordinates and
        // lengths (`seq_range::record_frame`); without ranges they are the search's.
        let query_frame =
            record_frame(hit.query_frame > 0, hit.q_start as usize - 1, query_length) as i8;
        let subject_frame = record_frame(
            hsp.subject_frame > 0,
            hit.s_start as usize - 1,
            subject_length,
        ) as i8;
        lowercase_query_row(
            &mut query_row,
            query_frame,
            low,
            high,
            &self.query_masks[q_idx],
            self.mask_seek,
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
            // The title after its first word (the `bio` description of the inputs
            // that `bio` read alike), as BLASTN's; the reports take the subject's
            // title bytes by subject index (`write_tblastx_pairwise_report`).
            subject_title: subject
                .title
                .iter()
                .position(|&byte| byte == b' ')
                .map(|space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned()),
            comp_adjust_method: None,
            sum_n: Some(hsp.num),
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
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
/// ```c++
/// 	if (seqinfo_src->CanReturnPartialSequence() == true)
/// 	{
///         	CConstRef<CSeq_loc> subj_loc = seqinfo_src->GetSeqLoc(kOid);
///         	NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
///              	   RemapToQueryLoc(*iter, query_loc);
///              	   if ( !is_ooframe )
///                    	s_RemapToSubjectLoc(*iter, *subj_loc);
/// ```
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

// NCBI reference: c++/src/objtools/align_format/tabular.cpp:1275-1283
// ```c++
//     x_PrintQueryAndDbNames(program_version, bioseq, dbname, rid, iteration, subj_bioseq);
//     // Print number of alignments found, but only if it has been set.
//     if (align_set) {
//        int num_hits = align_set->Get().size();
//        if (num_hits != 0) {
//            PrintFieldNames(is_csv);
//        }
//        m_Ostream << "# " << num_hits << " hits found" << "\n";
//     }
// ```
/// The block header of one query in outfmt 7: the program, the query, the subjects, and,
/// for a searched query (`num_hits` is `Some`), the fields when there are hits and the
/// count of HSPs.
pub(crate) fn write_outfmt7_query_header<W: Write>(
    writer: &mut W,
    query_title: &[u8],
    database: &str,
    num_hits: Option<usize>,
) -> io::Result<()> {
    writeln!(writer, "# TBLASTX {NCBI_TBLASTX_VERSION}")?;
    // NCBI reference: c++/src/objtools/align_format/tabular.cpp:1305-1308
    // ```c
    //     CAlignFormatUtil::AcknowledgeBlastQuery(bioseq, kLineLength, m_Ostream,
    //                                             m_ParseLocalIds, kHtmlFormat,
    //                                             kTabularFormat, rid);
    // ```
    // The title's bytes (`write_outfmt7_query_line`).
    crate::report::outfmt6::write_outfmt7_query_line(writer, query_title)?;
    writeln!(writer, "# Database: {database}")?;
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
/// block header of every query and, after the last batch, the count of queries.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
/// ```c
/// void CBlastTabularInfo::Print()
/// {
///     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
///         if (iter != m_FieldsToShow.begin())
///             m_Ostream << m_FieldDelimiter;
///         x_PrintField(*iter);
///     }
///     m_Ostream << "\n";
/// }
/// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tabular<W: Write>(
    hits: &[TblastxHsp],
    writer: &mut W,
    comments: bool,
    query_titles: &[&[u8]],
    database: &str,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    unsearched: &[bool],
    epilog: bool,
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
        for &(hsp_index, hsp) in query_hits {
            // One printed row is one HSP; the probe marks it without changing it.
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.begin(hsp_index);
            }
            write_tabular_row(writer, hsp, query_ids, subject_ids, &config)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.end(hsp_index);
            }
        }
    }
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2233-2238
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
    //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
    //         tabinfo.PrintNumProcessed(m_QueriesFormatted);
    // ```
    // NCBI reference: c++/src/objtools/align_format/tabular.cpp:1322-1325
    // ```c++
    // void CBlastTabularInfo::PrintNumProcessed(int num_queries)
    // {
    //     m_Ostream << "# BLAST processed " << num_queries << " queries\n";
    // }
    // ```
    if comments && epilog {
        writeln!(writer, "# BLAST processed {} queries", query_titles.len())?;
    }
    Ok(())
}

/// The outfmt 6 / 7 row of one HSP (`write_tabular`).
pub(crate) fn write_tabular_row<W: Write>(
    writer: &mut W,
    hsp: &TblastxHsp,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    config: &OutputConfig,
) -> io::Result<()> {
    let hit = &hsp.hit;
    // The query ID of the record (`FastaRecord::shown_id`) and the subject's
    // `GetSeqIdList` ID (`tabular_subject_id`), as bytes.
    let query_id = query_ids
        .get(hit.q_idx as usize)
        .map_or(&b""[..], |id| id.as_ref());
    let subject_id = subject_ids
        .get(hit.s_idx as usize)
        .map_or(&b""[..], |id| id.as_ref());
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

/// The TBLASTX output formats (`-outfmt` 0, 5, 6, 7, 8, 9, 11, 12, 15 and 16, `tblastx_outfmt`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TblastxOutputFormat {
    Pairwise,
    Tabular,
    TabularWithComments,
    /// A Seq-annot of Seq-aligns (outfmt 8, 9, 12).
    Seqalign(crate::report::serial::SerialFormat),
    /// BLAST XML (outfmt 5, `tblastx/blast_xml.rs`).
    Xml,
    /// The XML2 and JSON documents (outfmt 16, 15, 14 and 13), written from the same records as
    /// outfmt 5 (`tblastx/blast_xml.rs`); 14 and 13 as side files with this output as the master.
    Xml2(crate::report::blast_xml2::Xml2Kind),
    /// The BLAST archive (outfmt 11, `tblastx/archive.rs`): one document per query batch.
    Archive,
}

impl TblastxOutputFormat {
    /// Whether the format is a table (outfmt 6, 7).
    pub(crate) fn is_tabular(self) -> bool {
        matches!(self, Self::Tabular | Self::TabularWithComments)
    }

    /// Whether the format writes the displayed identity counts of the rows
    /// (`set_displayed_identities`): not the structured formats, which read the engine's counts
    /// (outfmt 8, 9, 12) or count the rows again themselves (outfmt 5, 15, 16).
    pub(crate) fn reads_displayed_counts(self) -> bool {
        !matches!(
            self,
            Self::Seqalign(_) | Self::Xml | Self::Xml2(_) | Self::Archive
        )
    }

    /// Whether the format is BLAST XML (5) or one of its XML2 / JSON siblings (15, 16), which
    /// are made from the same records and refuse the same names.
    pub(crate) fn is_structured_xml(self) -> bool {
        matches!(self, Self::Xml | Self::Xml2(_))
    }

    /// Whether the format writes a query's own warnings to standard error
    /// (`ReportFormat::writes_query_warnings`: not the structured formats).
    pub(crate) fn writes_query_warnings(self) -> bool {
        self.reads_displayed_counts()
    }
}

/// The format of a validated `-outfmt` value (`tblastx_outfmt`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
/// ```c
///     if (args[kArgOutputFormat]) {
///         string fmt_choice =
///             NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
/// ```
pub(crate) fn output_format(outfmt: &str) -> TblastxOutputFormat {
    match crate::blastinput::app::parse_formatting_string(outfmt).map(|choice| choice.number) {
        Ok(0) => TblastxOutputFormat::Pairwise,
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
        _ => TblastxOutputFormat::Tabular,
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
//                       Set to 0 if not applicable */
//    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
//    Int4 hspcnt; /**< Number of HSPs saved */
//    ...
//    double best_evalue; /**< Smallest e-value for HSPs in this list. Filled after
//                           e-values are calculated. Necessary because HSPs are
//                           sorted by score, but highest scoring HSP may not have
//                           the lowest e-value if sum statistics is used. */
// } BlastHSPList;
// ```
/// The final HSPs of one query and one subject, as the hit list holds them.
struct TblastxHspList {
    oid: u32,
    hsps: Vec<TblastxHsp>,
    best_evalue: f64,
}

impl crate::algorithm::blastn::hsp::HitListEntry for TblastxHspList {
    fn oid(&self) -> u32 {
        self.oid
    }

    fn hsp_count(&self) -> usize {
        self.hsps.len()
    }

    fn best_evalue(&self) -> f64 {
        self.best_evalue
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1749
    // ```c
    // static double
    // s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
    // {
    //     int index = 0;
    //     double best_evalue = (double) INT4_MAX;
    //
    //     for (index=0; index<hsp_list->hspcnt; index++)
    //        best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
    //
    //     return best_evalue;
    // }
    // ```
    fn update_best_evalue(&mut self) {
        let mut best = i32::MAX as f64;
        for hsp in &self.hsps {
            if hsp.hit.e_value < best {
                best = hsp.hit.e_value;
            }
        }
        self.best_evalue = best;
    }

    fn first_score(&self) -> Option<i32> {
        self.hsps.first().map(|hsp| hsp.hit.raw_score)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1452-1454
    // ```c
    //       qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //             s_EvalueCompareHSPs);
    // ```
    fn sort_by_evalue(&mut self) {
        self.hsps
            .sort_by(|a, b| crate::common::evalue_compare_hsps(&a.hit, &b.hit));
    }
}

/// The filters of NCBI's `s_FilterBlastResults` that TBLASTX has: `-max_hsps` (`max_hsps`,
/// 0 when omitted) and `-qcov_hsp_perc` (`query_cov_hsp_perc`, 0 when omitted), with the
/// nucleotide lengths of the queries (`query_lengths`, read for the coverage only).
///
/// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1533-1538
/// ```c
///                } else {
///         	   status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
///                }
///
///                Blast_HSPListReapByQueryCoverage(hsp_list,hit_params->options, query_info, program_number);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_traceback.c:1763-1766
/// ```c
///     if(hit_params->options->query_cov_hsp_perc > 0 || hit_params->options->max_hsps_per_subject > 0 ||
///        (hit_params->options->hsp_filt_opt != NULL && hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
///     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
///     }
/// ```
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
}

impl HitFilters<'_> {
    /// The lengths of the contexts (frames) of query `q_idx`, read only by the coverage.
    fn context_lengths(&self, q_idx: u32) -> [i32; 6] {
        if self.query_cov_hsp_perc > 0.0 {
            super::hsp_culling::tblastx_context_lengths(self.query_lengths[q_idx as usize])
        } else {
            [0; 6]
        }
    }

    /// `Blast_HSPListReapByQueryCoverage` on one subject's HSPs; returns whether none is
    /// left (NCBI's collector and culling writer skip an empty list).
    fn reap_by_query_coverage(
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
}

/// The final HSP list of a run in NCBI's order, which every format prints: per query (in
/// input order) the subjects that the hit list keeps, best first, and each subject's
/// HSPs in the order of its Seq-align.
///
/// NCBI reference: c++/src/algo/blast/core/link_hsps.c:1802-1803
/// ```c
///     /* Sort the HSP array by score */
///     Blast_HSPListSortByScore(hsp_list);
/// ```
/// NCBI reference: c++/src/algo/blast/core/hspfilter_collector.c:143-150
/// ```c
///       for (index = 0; index < results->num_queries; index++) {
///          if (hsp_list_array[index]) {
///             if (!results->hitlist_array[index]) {
///                results->hitlist_array[index] =
///                   Blast_HitListNew(params->prelim_hitlist_size);
///             }
///             Blast_HitListUpdate(results->hitlist_array[index],
///                                 hsp_list_array[index]);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_traceback.c:1768-1772
/// ```c
///     /* Re-sort the hit lists according to their best e-values, because they
///        could have changed. Only do this for a database search. */
///     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
///         Blast_HSPResultsSortByEvalue(results);
///     }
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
/// ```c
///         // Sort HSPs with e-values as first priority and scores as
///         // tie-breakers, since that is the order we want to see them in
///         // in Seq-aligns.
///         Blast_HSPListSortByEvalue(hsp_list);
/// ```
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
pub(crate) fn final_hit_order(
    hits: Vec<TblastxHsp>,
    hitlist_size: usize,
    filters: &HitFilters<'_>,
    // EXPERIMENT (LOSAT_X_TBXRPTPAR, x_rpt_par.rs item 4): the pool the lists' sorts run on, `None`
    // for this thread.
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Vec<TblastxHsp> {
    use crate::algorithm::blastn::hsp::HitList;
    let mut lists: std::collections::BTreeMap<(u32, u32), Vec<TblastxHsp>> =
        std::collections::BTreeMap::new();
    for hsp in hits {
        lists
            .entry((hsp.hit.q_idx, hsp.hit.s_idx))
            .or_default()
            .push(hsp);
    }
    let mut lists: Vec<((u32, u32), Vec<TblastxHsp>)> = lists.into_iter().collect();
    // Each list's sort by score depends on that list alone: on the pool or here, the same order.
    super::x_rpt_par::for_each_list(x_pool, &mut lists, |(_, hsps)| {
        hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
    });
    let mut ordered = Vec::new();
    let mut lists = lists.into_iter().peekable();
    while let Some(&((q_idx, _), _)) = lists.peek() {
        let context_lengths = filters.context_lengths(q_idx);
        let mut hit_list: HitList<TblastxHspList> = HitList::new(hitlist_size);
        while let Some(((_, oid), mut hsps)) = lists.next_if(|((q, _), _)| *q == q_idx) {
            if filters.reap_by_query_coverage(&mut hsps, &context_lengths) {
                continue;
            }
            hit_list.update(TblastxHspList {
                oid,
                hsps,
                best_evalue: 0.0,
            });
        }
        if filters.max_hsps > 0 {
            for list in &mut hit_list.hsplist_array {
                trim_by_max_hsps(&mut list.hsps, filters.max_hsps);
            }
        }
        hit_list.sort_by_evalue();
        // Each list's sort by e-value depends on that list alone (as the sort by score above).
        super::x_rpt_par::for_each_list(x_pool, &mut hit_list.hsplist_array, |list| {
            crate::algorithm::blastn::hsp::HitListEntry::sort_by_evalue(list);
        });
        for list in hit_list.hsplist_array {
            ordered.extend(list.hsps);
        }
    }
    ordered
}

/// The final HSP list of a run with `-culling_limit` (`culling_limit` > 0) in NCBI's order:
/// `final_hit_order` with NCBI's culling writer in the preliminary stage and its culling
/// pipe after the traceback (`hsp_culling.rs`). `query_lengths` are the queries'
/// nucleotide lengths (the lengths of their frames are the culling trees' ranges).
///
/// NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1341
/// ```c
///     if (s <= 0) {
///         return;
///     }
/// ...
///     if (m_HitSaveOpts->hsp_filt_opt->culling_opts == NULL) {
///         BlastHSPCullingOptions* culling = BlastHSPCullingOptionsNew(s);
///         BlastHSPFilteringOptions_AddCulling(m_HitSaveOpts->hsp_filt_opt,
///                                             &culling,
///                                             eBoth);
/// ```
/// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:330-341
/// ```c
///         else if (filt_opts->culling_opts &&
///                  (filt_opts->culling_stage & ePrelimSearch))
///         {
///             BlastHSPCullingParams* params =
///                 BlastHSPCullingParamsNew(opts_memento->m_HitSaveOpts,
///                      filt_opts->culling_opts,
///                      opts_memento->m_ExtnOpts->compositionBasedStats,
///                      opts_memento->m_ScoringOpts->gapped_calculation);
///             if(params->culling_max > 1){
///             	params->culling_max += 3;
///             }
///             writer_info = BlastHSPCullingInfoNew(params);
///         }
/// ```
/// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:386-394
/// ```c
///         } else if (filt_opts->culling_opts &&
///                    (filt_opts->culling_stage & eTracebackSearch)) {
///             BlastHSPCullingParams* params =
///                 BlastHSPCullingParamsNew(opts_memento->m_HitSaveOpts,
///                      filt_opts->culling_opts,
///                      opts_memento->m_ExtnOpts->compositionBasedStats,
///                      opts_memento->m_ScoringOpts->gapped_calculation);
///             BlastHSPPipeInfo_Add(&pipe_info,
///                                  BlastHSPCullingPipeInfoNew(params));
/// ```
/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:699-730
/// ```c
///    s_BlastHSPCullingInit(data, results);
///    for (qid = 0; qid < results->num_queries; ++qid) {
///          if (!(results->hitlist_array[qid])) continue;
///          num_list = results->hitlist_array[qid]->hsplist_count;
///          for (sid = 0; sid < num_list; ++sid) {
///         	 hsp_list = results->hitlist_array[qid]->hsplist_array[sid];
///         	 Blast_HSPListSortByEvalue(hsp_list);
///         	 hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
///          }
///          Blast_HitListSortByEvalue(results->hitlist_array[qid]);
///    }
///
///    for (qid = 0; qid < results->num_queries; ++qid) {
/// ...
///       for (sid = 0; sid < num_list; ++sid) {
///          s_BlastHSPCullingRun(data,
///                    results->hitlist_array[qid]->hsplist_array[sid]);
/// ...
///    s_BlastHSPCullingFinal(data, results);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_traceback.c:1720-1721
/// ```c
///         /* post-traceback pipes */
///         BlastHSPStreamTBackClose(hsp_stream, results);
/// ```
/// The culling writer replaces the preliminary hit list (the collector): each subject's
/// HSPs (after the `-evalue` reap, in score order) go into the trees, subject by subject
/// in OID order, with the merit `culling_max` (+3 above 1; in `i64`, where NCBI's `Int4`
/// overflows near `INT4_MAX`, NCBI parity DECISIONS 17). The trees' HSPs then reach the hit list of the traceback in OID order, where
/// the `-max_target_seqs` cut is made (`Blast_HitListUpdate`), and the pipe culls the kept
/// lists again, sorted by e-value, with the merit `culling_limit`. The rest is
/// `final_hit_order`'s: the lists sorted by e-value (their HSPs in score order) and each
/// list's HSPs sorted by e-value. The trees belong to the query's contexts, so the
/// queries are culled one after another.
pub(crate) fn culled_hit_order(
    hits: Vec<TblastxHsp>,
    hitlist_size: usize,
    culling_limit: i32,
    query_lengths: &[usize],
    filters: &HitFilters<'_>,
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
            final_lists.hsplist_array.push(list);
        }
        final_lists.hsplist_count = final_lists.hsplist_array.len();
        final_lists.sort_by_evalue();
        final_lists.prune_by_size(hitlist_size);
        for mut list in final_lists.hsplist_array {
            HitListEntry::sort_by_evalue(&mut list);
            ordered.extend(list.hsps);
        }
    }
    ordered
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Hit;

    fn hsp(q_idx: u32, s_idx: u32, raw_score: i32, e_value: f64, offset: usize) -> TblastxHsp {
        TblastxHsp {
            hit: Hit {
                identity: 0.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 30,
                s_start: 1,
                s_end: 30,
                e_value,
                bit_score: raw_score as f64 / 2.0,
                num_ident: 0,
                query_frame: 1,
                query_length: 0,
                q_idx,
                s_idx,
                raw_score,
                sort_query_offset: offset,
                sort_query_end: offset + 10,
                sort_subject_offset: offset,
                sort_subject_end: offset + 10,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 0,
            },
            subject_frame: 1,
            num: 1,
        }
    }

    // NCBI c++/src/objects/seqfeat/Genetic_code_table.cpp:94-112: the display reads the
    // IUPAC letters of both cases, U as T and X as N; another byte is a gap.
    #[test]
    fn display_bases_follow_the_translation_table() {
        assert_eq!(display_base(b'a'), 1);
        assert_eq!(display_base(b'T'), 8);
        assert_eq!(display_base(b'u'), 8);
        assert_eq!(display_base(b'X'), 15);
        assert_eq!(display_base(b'n'), 15);
        assert_eq!(display_base(b'R'), 5);
        assert_eq!(display_base(b'*'), 0);
    }

    // The display table (LOSAT_X_DISPTABLE) holds the scan's answers: each byte's code on both
    // strands, each codon of 4-bit codes of each genetic code, and so each residue of a row.
    #[test]
    fn display_table_gives_the_residues_of_the_scan() {
        let reversed = |mask: u8| {
            ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
        };
        for byte in 0..=255u8 {
            assert_eq!(DISPLAY_BASE_PLUS[byte as usize], display_base(byte));
            assert_eq!(
                DISPLAY_BASE_MINUS[byte as usize],
                reversed(display_base(byte))
            );
        }
        let letters = b"ACGTUNXRYMKSWBDHVacgtunxry-*";
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let sequence: Vec<u8> = (0..61)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                letters[(state % letters.len() as u64) as usize]
            })
            .collect();
        let mut codes = 0;
        for id in 0..=255u8 {
            let Ok(code) = GeneticCode::try_from_id(id) else {
                continue;
            };
            codes += 1;
            let table = display_table(&code);
            for (codon, &residue) in table.iter().enumerate() {
                let masks = [
                    (codon >> 8) as u8,
                    ((codon >> 4) & 15) as u8,
                    (codon & 15) as u8,
                ];
                assert_eq!(
                    residue,
                    crate::algorithm::blastx::report::display_codon(masks, &code)
                );
            }
            let display = DisplayCode {
                code: &code,
                table: Some(table),
            };
            for frame in [1i8, 2, 3, -1, -2, -3] {
                for offset in 0..(sequence.len() + 1 - frame.unsigned_abs() as usize) / 3 {
                    assert_eq!(
                        display.residue(&sequence, frame, offset),
                        display_residue_scan(&sequence, frame, offset, &code),
                        "code {id} frame {frame} offset {offset}"
                    );
                }
            }
        }
        assert!(codes > 20);
    }

    // NCBI c++/src/objects/seqfeat/Genetic_code_table.cpp:191-207 merges D/N to B, E/Q
    // to Z and I/L to J; the minus strand is the IUPAC complement in reverse.
    #[test]
    fn displayed_residues_merge_ambiguous_codons_on_both_strands() {
        let code = GeneticCode::from_id(1);
        assert_eq!(display_residue(b"GCN", 1, 0, &code), b'A');
        assert_eq!(display_residue(b"RAY", 1, 0, &code), b'B');
        assert_eq!(display_residue(b"SAR", 1, 0, &code), b'Z');
        assert_eq!(display_residue(b"MTA", 1, 0, &code), b'J');
        assert_eq!(display_residue(b"NNN", 1, 0, &code), b'X');
        assert_eq!(display_residue(b"TAR", 1, 0, &code), b'*');
        // Frame -1 of ATG reads CAT; frame +2 of AATGC reads ATG.
        assert_eq!(display_residue(b"ATG", -1, 0, &code), b'H');
        assert_eq!(display_residue(b"AATGC", 2, 0, &code), b'M');
        // The complement of RTY is RAY.
        assert_eq!(display_residue(b"RTY", -1, 0, &code), b'B');
    }

    // NCBI c++/src/objtools/align_format/tabular.cpp:1003-1021 and showalign.cpp:2131-2150:
    // equal letters are identities (X/X too); B/D is positive, J/I is not (J is not in the
    // display matrix).
    #[test]
    fn row_counts_use_the_displayed_letters() {
        assert_eq!(row_counts(Blosum62Score, b"XAB*J", b"XADKI"), (2, 3));
    }

    // The identities pass counts without the rows: the counts of a residue-by-residue loop
    // over `display_residue` (the query with code 1, the subject with code 4), on every frame
    // pair, for HSPs at the start, inside and at the end of the subject's frame, for a record and
    // a stored subject with ambiguous bases; and the same as `row_counts` of the rows.
    #[test]
    fn displayed_counts_equal_the_counts_of_the_rows() {
        let (query_code, db_code) = (GeneticCode::from_id(1), GeneticCode::from_id(4));
        let subject: Vec<u8> = (0..120u32)
            .map(|i| b"ACGTACGTRYNACGT"[(i * 7 + i / 5) as usize % 15])
            .collect();
        let query: Vec<u8> = (0..120u32)
            .map(|i| b"ACGTACGTN"[(i * 11 + i / 4) as usize % 9])
            .collect();
        let records = [FastaRecord::new("s", b"s", &subject)];
        let stores = [PackedNucleotide::from_subject_letters(&subject).expect("letters")];
        let ordinals = [0u32];
        let forms = [
            SubjectLetters::Records(&records),
            SubjectLetters::Stored {
                stores: &stores,
                ordinals: &ordinals,
            },
        ];
        let frame_length = |frame: i8| (120 - usize::from(frame.unsigned_abs()) + 1) / 3;
        let mut window = Vec::new();
        let (mut mismatched, mut positive) = (0, 0);
        for subjects in forms {
            for query_frame in [1i8, 2, 3, -1, -2, -3] {
                for subject_frame in [1i8, 2, 3, -1, -2, -3] {
                    let last = frame_length(subject_frame) - 10;
                    for (query_offset, subject_offset) in [(0, 0), (3, 7), (25, last)] {
                        let mut hsp = hsp(0, 0, 10, 1.0, 0);
                        hsp.hit.query_frame = query_frame as i32;
                        hsp.hit.sort_query_offset = query_offset;
                        hsp.hit.sort_query_end = query_offset + 10;
                        hsp.hit.sort_subject_offset = subject_offset;
                        hsp.hit.sort_subject_end = subject_offset + 10;
                        hsp.subject_frame = subject_frame;
                        let mut expected = (0, 0);
                        for i in 0..10 {
                            let q =
                                display_residue(&query, query_frame, query_offset + i, &query_code);
                            let s = display_residue(
                                &subject,
                                subject_frame,
                                subject_offset + i,
                                &db_code,
                            );
                            if q == s {
                                expected.0 += 1;
                                expected.1 += 1;
                            } else if protein_display_score(ScoringMatrix::Blosum62, q, s) > 0 {
                                expected.1 += 1;
                            }
                        }
                        let counts = subjects.counts(
                            &hsp,
                            &query,
                            &query_code,
                            &db_code,
                            &mut window,
                            Blosum62Score,
                        );
                        assert_eq!(
                            counts, expected,
                            "frames {query_frame} {subject_frame} at {subject_offset}"
                        );
                        let (query_row, subject_row) =
                            subjects.rows(&hsp, &query, &query_code, &db_code, &mut window);
                        assert_eq!(
                            row_counts(Blosum62Score, &query_row, &subject_row),
                            expected
                        );
                        mismatched += usize::from(expected.0 < 10);
                        positive += usize::from(expected.1 > expected.0);
                    }
                }
            }
        }
        // Not all identities, and some positives that are not identities.
        assert!(mismatched > 0 && positive > 0, "{mismatched} {positive}");
    }

    // The rows are ASCII letters (uppercase and the lowercase of the masks): the String keeps
    // the row; any other byte gives the text of `from_utf8_lossy`.
    #[test]
    fn row_strings_are_the_lossy_text_of_the_rows() {
        for row in [&b"ARNDxyz*-"[..], b"", &[b'A', 0xff, b'C', 0xc3]] {
            assert_eq!(
                row_string(row.to_vec()),
                String::from_utf8_lossy(row).into_owned()
            );
        }
        let row = b"MKV".to_vec();
        let pointer = row.as_ptr();
        assert_eq!(row_string(row).as_ptr(), pointer);
    }

    // NCBI c++/src/algo/blast/core/blast_filter.c:925-931: a minus-frame run starts one
    // codon too high, so its last residue is not masked.
    #[test]
    fn query_masks_convert_with_ncbis_frame_formulas() {
        let stats = TblastxQueryStats {
            seg_masks: vec![(1, vec![(2, 5)]), (-1, vec![(2, 5)])],
            ..TblastxQueryStats::default()
        };
        // A 30-nt query: frame lengths 10, 9, 9, so the DNA length is 30.
        let masks = query_dna_masks(&stats, 30);
        assert_eq!(masks[0], (1, vec![(6, 12)]));
        assert_eq!(masks[1], (-1, vec![(18, 23)]));
        let masks: Vec<FrameMasks> = shown_query_masks(masks, 30, 0)
            .into_iter()
            .map(|(frame, ranges)| FrameMasks::new(frame, ranges))
            .collect();
        for seek in [false, true] {
            let mut plus = b"AAAAAAAAAA".to_vec();
            lowercase_query_row(&mut plus, 1, 0, 29, &masks, seek);
            assert_eq!(&plus, b"AAaaaAAAAA");
            // Frame -1 residue i covers 0-based nucleotides 27-3i..29-3i: the mask 18-23 covers
            // residues 2 and 3 only.
            let mut minus = b"AAAAAAAAAA".to_vec();
            lowercase_query_row(&mut minus, -1, 0, 29, &masks, seek);
            assert_eq!(&minus, b"AAaaAAAAAA");
        }
    }

    #[test]
    fn mask_seek_lowercases_the_letters_of_the_full_scan() {
        // Ranges in any order, overlapping, empty (start after end) or outside the row, on rows
        // of both frame signs.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = |bound: i32| -> i32 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as i32
        };
        for _ in 0..2000 {
            let ranges: Vec<(i32, i32)> = (0..next(9))
                .map(|_| {
                    let from = next(70) - 5;
                    (from, from + next(16) - 3)
                })
                .collect();
            let masks = [
                FrameMasks::new(1, ranges.clone()),
                FrameMasks::new(-1, ranges),
            ];
            let low = next(55);
            let high = low + 2 + next(58 - low);
            for frame in [1, -1] {
                let mut rows = [false, true].map(|seek| {
                    let mut row = vec![b'A'; ((high - low + 1) / 3) as usize];
                    lowercase_query_row(&mut row, frame, low, high, &masks, seek);
                    row
                });
                let seek = rows[1].split_off(0);
                assert_eq!(seek, rows[0], "frame {frame} low {low} high {high}");
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_filter.c:993 (BlastSeqLocCombine(filter_out, 0)) joins
    // only intervals that overlap (adjacent ones stay apart), per frame, before the conversion
    // to nucleotides; the XML writer reads that list (blast_aux.cpp:930-959), in context order,
    // each mask with its context's frame, moved into the record by the interval's start.
    #[test]
    fn xml_masks_join_overlapping_intervals_only_and_keep_the_frames() {
        let stats = TblastxQueryStats {
            // Frame 1: residues 2-4 and 5-7 (adjacent, apart), frame 2: 2-6 and 4-7 (joined),
            // frame -1: one run.
            seg_masks: vec![
                (1, vec![(5, 8), (2, 5)]),
                (2, vec![(2, 7), (4, 8)]),
                (-1, vec![(2, 5)]),
            ],
            ..TblastxQueryStats::default()
        };
        let mut out = Vec::new();
        xml_query_masks(&stats, 30, 0, &mut out);
        let got: Vec<(i32, i32, i32)> = out.iter().map(|m| (m.from, m.to, m.frame)).collect();
        assert_eq!(got, vec![(6, 12, 1), (15, 21, 1), (7, 22, 2), (18, 23, -1)]);
        // `-query_loc`: moved by the interval's start, buffer cleared first.
        xml_query_masks(&stats, 30, 100, &mut out);
        assert_eq!(out.len(), 4);
        assert_eq!((out[0].from, out[0].to), (106, 112));
        // No SEG masks, no masks.
        xml_query_masks(&TblastxQueryStats::default(), 30, 0, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn the_structured_formats_read_no_displayed_counts_and_write_no_query_warnings() {
        assert_eq!(output_format("5"), TblastxOutputFormat::Xml);
        // The XML2 and JSON documents of 16 and 15 take the path of outfmt 5 (`is_structured_xml`).
        use crate::report::blast_xml2::Xml2Kind;
        assert_eq!(
            output_format(" 16 std"),
            TblastxOutputFormat::Xml2(Xml2Kind::Single16)
        );
        assert_eq!(
            output_format("15"),
            TblastxOutputFormat::Xml2(Xml2Kind::Single15)
        );
        for (outfmt, structured) in [
            ("5", true),
            ("13", true),
            ("14", true),
            ("15", true),
            ("16", true),
            ("0", false),
            ("6", false),
            ("7", false),
            ("8", false),
            ("9", false),
            ("12", false),
        ] {
            assert_eq!(
                output_format(outfmt).is_structured_xml(),
                structured,
                "{outfmt}"
            );
        }
        for (outfmt, counts) in [
            ("0", true),
            ("5", false),
            ("6", true),
            ("7", true),
            ("8", false),
            ("9", false),
            ("12", false),
            ("15", false),
            ("16", false),
        ] {
            let format = output_format(outfmt);
            assert_eq!(format.reads_displayed_counts(), counts, "{outfmt}");
            assert_eq!(format.writes_query_warnings(), counts, "{outfmt}");
        }
    }

    // NCBI c++/src/algo/blast/api/blast_aux.cpp:825-842 (`Map`): with `-query_loc` a mask
    // moves by the interval's start and is cut at its end; one that covers the whole
    // interval, or starts past it, is not shown.
    #[test]
    fn shown_query_masks_follow_ncbis_map() {
        let masks = vec![(2, vec![(0, 9), (3, 20), (0, 29), (31, 40), (5, 4)])];
        assert_eq!(
            shown_query_masks(masks, 30, 100),
            vec![(2, vec![(100, 109), (103, 120)])]
        );
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:3243-3299 (Blast_HitListUpdate) and
    // 3078-3107 (s_EvalueCompareHSPLists): equal best e-values are ordered by the score of
    // the first HSP, which is the highest score until the hit list overflows and its lists
    // are sorted by e-value.
    #[test]
    fn hit_list_keeps_ncbis_subjects_on_ties() {
        let hits = vec![
            hsp(0, 0, 120, 1e-3, 0),
            hsp(0, 0, 100, 1e-5, 20),
            hsp(0, 1, 110, 1e-5, 0),
        ];
        let order = |size| -> Vec<(u32, i32)> {
            final_hit_order(hits.clone(), size, &HitFilters::default(), None)
                .iter()
                .map(|hsp| (hsp.hit.s_idx, hsp.hit.raw_score))
                .collect()
        };
        // Not full: subject 0 compares with its highest score (120 > 110).
        assert_eq!(order(2), vec![(0, 100), (0, 120), (1, 110)]);
        // Full at 1: the lists are sorted by e-value (100 < 110), so subject 1 is kept.
        assert_eq!(order(1), vec![(1, 110)]);
    }

    #[test]
    fn outfmt7_headers_follow_print_header() {
        let header = |hits| {
            let mut out = Vec::new();
            write_outfmt7_query_header(&mut out, b"q1 title", "db", hits).unwrap();
            String::from_utf8(out).unwrap()
        };
        assert_eq!(
            header(None),
            "# TBLASTX 2.17.0+\n# Query: q1 title\n# Database: db\n"
        );
        assert!(header(Some(0)).ends_with("# Database: db\n# 0 hits found\n"));
        assert!(header(Some(2)).contains("# Fields: query acc.ver, subject acc.ver"));
        assert!(header(Some(2)).ends_with("# 2 hits found\n"));
    }

    #[test]
    fn output_formats_are_read_after_trimming() {
        assert_eq!(output_format(" 0"), TblastxOutputFormat::Pairwise);
        assert_eq!(output_format("6"), TblastxOutputFormat::Tabular);
        assert_eq!(
            output_format("7 "),
            TblastxOutputFormat::TabularWithComments
        );
    }
}

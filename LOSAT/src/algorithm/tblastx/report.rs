//! The TBLASTX reports from the final HSP list of a run: the displayed rows of every HSP
//! (outfmt 0, the hit records, and the identity counts of outfmt 6 and 7), the per-query
//! blocks of outfmt 7, and the data of the pairwise report.

use std::io::{self, Write};
use std::sync::Arc;

use bio::io::fasta;

use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::seq_range::{record_frame, Placements};
use crate::config::ScoringMatrix;
use crate::report::outfmt6::{write_hit_fields, OutputConfig};
use crate::report::pairwise::PairwiseHit;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::protein_display_score;

use super::blast_engine::{TblastxHsp, TblastxQueryStats};

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
/// The displayed query and subject rows of an HSP, in uppercase: the HSP's codons of each
/// sequence translated again from the nucleotides, the query with `-query_gencode` and the
/// subject with `-db_gencode` (`query_code`, `db_code`). The HSP's residue offsets within
/// its frames are those of the search (`Hit::sort_*`).
pub(crate) fn displayed_rows(
    hsp: &TblastxHsp,
    query: &[u8],
    subject: &[u8],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) -> (Vec<u8>, Vec<u8>) {
    let hit = &hsp.hit;
    let query_frame = hit.query_frame as i8;
    let query_row = (hit.sort_query_offset..hit.sort_query_end)
        .map(|offset| display_residue(query, query_frame, offset, query_code))
        .collect();
    let subject_row = (hit.sort_subject_offset..hit.sort_subject_end)
        .map(|offset| display_residue(subject, hsp.subject_frame, offset, db_code))
        .collect();
    (query_row, subject_row)
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
/// The identities and the positives (identities included) of two displayed rows.
fn row_counts(query_row: &[u8], subject_row: &[u8]) -> (usize, usize) {
    let mut identities = 0;
    let mut positives = 0;
    for (&q, &s) in query_row.iter().zip(subject_row) {
        if q == s {
            identities += 1;
            positives += 1;
        } else if protein_display_score(ScoringMatrix::Blosum62, q, s) > 0 {
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
pub(crate) fn set_displayed_identities(
    hits: &mut [TblastxHsp],
    queries: &[fasta::Record],
    subjects: &[fasta::Record],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) {
    for hsp in hits {
        let (query_row, subject_row) = displayed_rows(
            hsp,
            queries[hsp.hit.q_idx as usize].seq(),
            subjects[hsp.hit.s_idx as usize].seq(),
            query_code,
            db_code,
        );
        let (identities, positives) = row_counts(&query_row, &subject_row);
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
fn query_dna_masks(stats: &TblastxQueryStats, query_length: usize) -> Vec<(i8, Vec<(i32, i32)>)> {
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

/// Lowercases the residues of a displayed query row that a shown SEG mask (record
/// coordinates, `shown_query_masks`) labelled with the row's frame covers. `[low, high]` is
/// the 0-based nucleotide range of the row on the plus strand of the record, and `frame`
/// the row's frame in the record (`s_GetStdsegMasterFrame`). With `-query_loc` a mask's
/// label is the frame of the searched letters, which can differ from the frame the record
/// gives the same letters: NCBI then shows a mask of another frame on the row (`locFrame ==
/// frame`), and so does LOSAT.
fn lowercase_query_row(
    row: &mut [u8],
    frame: i8,
    low: i32,
    high: i32,
    masks: &[(i8, Vec<(i32, i32)>)],
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
    for (_, frame_masks) in masks.iter().filter(|(mask_frame, _)| *mask_frame == frame) {
        for &(from, to) in frame_masks {
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
/// lengths.
#[allow(clippy::too_many_arguments)]
pub(crate) fn pairwise_hits(
    hits: &[TblastxHsp],
    queries: &[fasta::Record],
    subjects: &[fasta::Record],
    query_stats: &[TblastxQueryStats],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    query_placements: &Placements,
    subject_placements: &Placements,
) -> Vec<PairwiseHit> {
    let query_masks: Vec<Vec<(i8, Vec<(i32, i32)>)>> = query_stats
        .iter()
        .zip(queries)
        .enumerate()
        .map(|(q_idx, (stats, query))| {
            shown_query_masks(
                query_dna_masks(stats, query.seq().len()),
                query.seq().len(),
                query_placements.offset(q_idx),
            )
        })
        .collect();
    hits.iter()
        .map(|hsp| {
            let query = &queries[hsp.hit.q_idx as usize];
            let subject = &subjects[hsp.hit.s_idx as usize];
            let (mut query_row, subject_row) =
                displayed_rows(hsp, query.seq(), subject.seq(), query_code, db_code);
            let (identities, positives) = row_counts(&query_row, &subject_row);
            let (q_idx, s_idx) = (hsp.hit.q_idx as usize, hsp.hit.s_idx as usize);
            let mut hit = hsp.hit.clone();
            shift_hit(&mut hit, query_placements, subject_placements);
            let query_length = query_placements.length(q_idx, query.seq().len());
            let subject_length = subject_placements.length(s_idx, subject.seq().len());
            let (low, high) = (
                hit.q_start.min(hit.q_end) as i32 - 1,
                hit.q_start.max(hit.q_end) as i32 - 1,
            );
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:427-438
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
            lowercase_query_row(&mut query_row, query_frame, low, high, &query_masks[q_idx]);
            hit.num_ident = identities;
            hit.num_positives = positives;
            PairwiseHit {
                hit,
                query_seq: Some(String::from_utf8_lossy(&query_row).into_owned()),
                subject_seq: Some(String::from_utf8_lossy(&subject_row).into_owned()),
                query_frame: Some(query_frame),
                subject_frame: Some(subject_frame),
                positives: Some(positives),
                gaps: Some(0),
                subject_length: Some(subject_length),
                subject_title: subject.desc().map(str::to_string),
                comp_adjust_method: None,
                sum_n: Some(hsp.num),
            }
        })
        .collect()
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

/// The `# Query:` and `Query=` text of a query: its FASTA defline.
///
/// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:729-734
/// ```c++
///         out << label << "= ";
///     ...
///         string all_id_str = GetSeqIdString(cbs, believe_query);
///         all_id_str += " ";
///         all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
/// ```
pub(crate) fn fasta_defline(record: &fasta::Record) -> String {
    match record.desc() {
        Some(desc) => format!("{} {desc}", record.id()),
        None => record.id().to_string(),
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
fn write_outfmt7_query_header<W: Write>(
    writer: &mut W,
    query_title: &str,
    database: &str,
    num_hits: Option<usize>,
) -> io::Result<()> {
    writeln!(writer, "# TBLASTX {NCBI_TBLASTX_VERSION}")?;
    writeln!(writer, "# Query: {query_title}")?;
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
    query_titles: &[String],
    database: &str,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
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
                &query_titles[q_idx],
                database,
                (!unsearched[q_idx]).then_some(query_hits.len()),
            )?;
        }
        for &(hsp_index, hsp) in query_hits {
            let hit = &hsp.hit;
            let (query_id, subject_id) = hit.resolve_ids(query_ids, subject_ids);
            // One printed row is one HSP; the probe marks it without changing it.
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.begin(hsp_index);
            }
            write_hit_fields(
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
                &config,
            )?;
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

/// The TBLASTX output formats (`-outfmt` 0, 6 and 7, `tblastx_outfmt`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TblastxOutputFormat {
    Pairwise,
    Tabular,
    TabularWithComments,
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
pub(crate) fn final_hit_order(hits: Vec<TblastxHsp>, hitlist_size: usize) -> Vec<TblastxHsp> {
    use crate::algorithm::blastn::hsp::HitList;
    let mut lists: std::collections::BTreeMap<(u32, u32), Vec<TblastxHsp>> =
        std::collections::BTreeMap::new();
    for hsp in hits {
        lists
            .entry((hsp.hit.q_idx, hsp.hit.s_idx))
            .or_default()
            .push(hsp);
    }
    let mut ordered = Vec::new();
    let mut lists = lists.into_iter().peekable();
    while let Some(&((q_idx, _), _)) = lists.peek() {
        let mut hit_list: HitList<TblastxHspList> = HitList::new(hitlist_size);
        while let Some(((_, oid), mut hsps)) = lists.next_if(|((q, _), _)| *q == q_idx) {
            hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
            hit_list.update(TblastxHspList {
                oid,
                hsps,
                best_evalue: 0.0,
            });
        }
        hit_list.sort_by_evalue();
        for mut list in hit_list.hsplist_array {
            crate::algorithm::blastn::hsp::HitListEntry::sort_by_evalue(&mut list);
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
/// in OID order, with the merit `culling_max` (+3 above 1, an `Int4` that NCBI's binary
/// wraps). The trees' HSPs then reach the hit list of the traceback in OID order, where
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
) -> Vec<TblastxHsp> {
    use super::hsp_culling::{tblastx_context_lengths, CullingWriter};
    use crate::algorithm::blastn::hsp::{HitList, HitListEntry};
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
    let prelim_max = if culling_limit > 1 {
        culling_limit.wrapping_add(3)
    } else {
        culling_limit
    };
    let mut ordered = Vec::new();
    for (q_idx, subjects) in lists {
        let context_lengths = tblastx_context_lengths(query_lengths[q_idx as usize]);
        let mut writer = CullingWriter::new(prelim_max, context_lengths);
        for (oid, mut hsps) in subjects {
            hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
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
        let mut pipe = CullingWriter::new(culling_limit, context_lengths);
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
        assert_eq!(row_counts(b"XAB*J", b"XADKI"), (2, 3));
    }

    // NCBI c++/src/algo/blast/core/blast_filter.c:925-931: a minus-frame run starts one
    // codon too high, so its last residue is not masked.
    #[test]
    fn query_masks_convert_with_ncbis_frame_formulas() {
        let stats = TblastxQueryStats {
            karlin: None,
            eff_searchsp: 0,
            seg_masks: vec![(1, vec![(2, 5)]), (-1, vec![(2, 5)])],
        };
        // A 30-nt query: frame lengths 10, 9, 9, so the DNA length is 30.
        let masks = query_dna_masks(&stats, 30);
        assert_eq!(masks[0], (1, vec![(6, 12)]));
        assert_eq!(masks[1], (-1, vec![(18, 23)]));
        let masks = shown_query_masks(masks, 30, 0);
        let mut plus = b"AAAAAAAAAA".to_vec();
        lowercase_query_row(&mut plus, 1, 0, 29, &masks);
        assert_eq!(&plus, b"AAaaaAAAAA");
        // Frame -1 residue i covers 0-based nucleotides 27-3i..29-3i: the mask 18-23 covers
        // residues 2 and 3 only.
        let mut minus = b"AAAAAAAAAA".to_vec();
        lowercase_query_row(&mut minus, -1, 0, 29, &masks);
        assert_eq!(&minus, b"AAaaAAAAAA");
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
            final_hit_order(hits.clone(), size)
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
            write_outfmt7_query_header(&mut out, "q1 title", "db", hits).unwrap();
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

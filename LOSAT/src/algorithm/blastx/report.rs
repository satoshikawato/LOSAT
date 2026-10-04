//! BLASTX formatting from real final native runtime results.
use super::{args::ResolvedOptions, input::FastaRecord, results::Hsp, runtime::BatchResults};
use crate::{
    common::GapEditOp,
    report::outfmt6::{
        format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
    },
};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
// ```c++
//     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
//
//     SNCBIFullScoreMatrix mtx;
//     NCBISM_Unpack(packed_mtx, &mtx);
//
//     for(int i = 0; i < ePMatrixSize; ++i){
//         for(int j = 0; j < ePMatrixSize; ++j){
//             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
//                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
//         }
//     }
//     for(int i = 0; i < ePMatrixSize; ++i) {
//         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
//     }
//     retval('*', '*') = 1;
//     // this is to count Selenocysteine to Cysteine matches as positive
//     retval('U', 'U') = retval('C', 'C');
//     retval('U', 'C') = retval('C', 'C');
//     retval('C', 'U') = retval('C', 'C');
// ```
use crate::utils::matrix::protein_display_score;
use anyhow::{ensure, Result};
use std::{io::Write, path::Path};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-987
// ```c++
//     if (x_IsFieldRequested(eQuerySeq) ||
//         x_IsFieldRequested(eSubjectSeq) ||
//         x_IsFieldRequested(ePositives) ||
//         x_IsFieldRequested(ePercentPositives) ||
//         x_IsFieldRequested(eBTOP) ||
//         (x_IsFieldRequested(eNumIdentical) && !kNoFetchSequence) ||
//         (x_IsFieldRequested(eMismatches) && !kNoFetchSequence) ||
//         (x_IsFieldRequested(ePercentIdentical) && !kNoFetchSequence)) {
//
//         alnVec->SetGapChar('-');
//         alnVec->SetGenCode(m_QueryGeneticCode, 0);
//         alnVec->SetGenCode(m_DbGeneticCode, 1);
//         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
//         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
//
// ```
#[derive(Debug)]
pub(crate) struct DisplayAlignment {
    pub qseq: String,
    pub sseq: String,
    pub masked_qseq: String,
    pub btop: String,
    pub identity: usize,
    pub positive: usize,
    pub gaps: usize,
    pub gap_opens: usize,
    pub qstart: i32,
    pub qend: i32,
    pub sstart: i32,
    pub send: i32,
    pub frame: i8,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1348-1361
// ```c++
//     if (hsp->query.frame == 0) {
//        query_loc->SetInt().SetFrom(hsp->query.offset);
//        query_loc->SetInt().SetTo(hsp->query.end - 1);
//     } else if (hsp->query.frame > 0) {
//        query_loc->SetInt().SetFrom(CODON_LENGTH*(hsp->query.offset) +
//                                    hsp->query.frame - 1);
//        query_loc->SetInt().SetTo(CODON_LENGTH*(hsp->query.end) +
//                                  hsp->query.frame - 2);
//     } else {
//        query_loc->SetInt().SetFrom(query_length -
//            CODON_LENGTH*(hsp->query.end) + hsp->query.frame + 1);
//        query_loc->SetInt().SetTo(query_length - CODON_LENGTH*hsp->query.offset
//                                  + hsp->query.frame);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1034-1043
// ```c++
//         // For translated search, for a negative query frame, reverse its start
//         // and end offsets.
//         if (kTranslated && ds.GetSeqStrand(kQueryRow) == eNa_strand_minus) {
//             q_start = alnVec->GetSeqStop(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStart(kQueryRow) + 1;
//         } else {
//             q_start = alnVec->GetSeqStart(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStop(kQueryRow) + 1;
//         }
//     }
// ```
pub(crate) fn query_endpoints(length: usize, frame: i8, start: i32, end: i32) -> (i32, i32) {
    let f = i32::from(frame);
    if f > 0 {
        (3 * start + f, 3 * end + f - 1)
    } else {
        (
            length as i32 - 3 * start + f + 1,
            length as i32 - 3 * end + f + 2,
        )
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
// ```c++
//     if (hsp->score == 0) {
//         return CRef<CSeq_align>();
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1484-1490
// ```c++
//         } else {
//             seqalign =
//                 s_BlastHSP2SeqAlign(program, hsp, query_id, subject_id,
//                                     query_length, subject_length);
//         }
//
//         if (seqalign.Empty()) continue;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1606-1626
// ```c++
//         GetFilteredRedundantSeqids(*seqinfo_src, hsp_list->oid, seqid_list, subject_id->IsGi());
//         // stores a CSeq_align for each matching sequence
//         vector<CRef<CSeq_align > > hit_align;
//         if (is_gapped) {
//                 BLASTHspListToSeqAlign(prog,
//                                        hsp_list,
//                                        query_id,
//                                        subject_id,
//                                        query_length,
//                                        subj_length,
//                                        is_ooframe,
//                                        seqid_list,
//                                        hit_align);
//         } else {
//                 BLASTUngappedHspListToSeqAlign(prog,
//                                                hsp_list,
//                                                query_id,
//                                                subject_id,
//                                                query_length,
//                                                subj_length,
//                                                seqid_list,
// ```
// This is Seq-align construction, after runtime search/pruning and e-value
fn has_seqalign(hsp: &Hsp, gapped: bool) -> bool {
    !gapped || hsp.hsp.score != 0
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:94-112
// ```c++
//     static char  charToBase [17] = "-ACMGRSVTWYHKDBN";
//     static char  baseToComp [17] = "-TGKCYSBAWRDMHVN";
//
//     // illegal characters map to 0
//     for (i = 0; i < 256; i++) {
//         sm_BaseToIdx [i] = 0;
//     }
//
//     // map iupacna alphabet to EBaseCode
//     for (i = eBase_gap; i <= eBase_N; i++) {
//         ch = charToBase [i];
//         sm_BaseToIdx [(int) ch] = i;
//         ch = (unsigned char)tolower (ch);
//         sm_BaseToIdx [(int) ch] = i;
//     }
//     sm_BaseToIdx [(int) 'U'] = eBase_T;
//     sm_BaseToIdx [(int) 'u'] = eBase_T;
//     sm_BaseToIdx [(int) 'X'] = eBase_N;
//     sm_BaseToIdx [(int) 'x'] = eBase_N;
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:145-148
// ```c++
//     static int  expansions [4] = {eBase_A, eBase_C, eBase_G, eBase_T};
//                                 // T = 0, C = 1, A = 2, G = 3
//     static int  codonIdx [9] = {0, 2, 1, 0, 3, 0, 0, 0, 0};
//
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:170-209
// ```c++
//                 // expand ambiguous IJK nucleotide symbols into component bases XYZ
//                 for (p = 0; p < 4 && go_on; p++) {
//                     x = expansions [p];
//                     if ((x & i) != 0) {
//                         for (q = 0; q < 4 && go_on; q++) {
//                             y = expansions [q];
//                             if ((y & j) != 0) {
//                                 for (r = 0; r < 4 && go_on; r++) {
//                                     z = expansions [r];
//                                     if ((z & k) != 0) {
//
//                                         // calculate offset in genetic code string
//
//                                         // the T = 0, C = 1, A = 2, G = 3 order is
//                                         // necessary because the genetic code strings
//                                         // are presented in TCAG order in printed tables
//                                         // and in the genetic code strings
//                                         cd = 16 * codonIdx [x] + 4 * codonIdx [y] + codonIdx [z];
//
//                                         // lookup amino acid for codon XYZ
//                                         ch = (*ncbieaa) [cd];
//                                         if (aa == '\0') {
//                                             aa = ch;
//                                         } else if (aa != ch) {
//                                             // allow Asx (Asp or Asn) and Glx (Glu or Gln)
//                                             if ((aa == 'B' || aa == 'D' || aa == 'N') &&
//                                                 (ch == 'D' || ch == 'N')) {
//                                                 aa = 'B';
//                                             } else if ((aa == 'Z' || aa == 'E' || aa == 'Q') &&
//                                                        (ch == 'E' || ch == 'Q')) {
//                                                 aa = 'Z';
//                                             } else if ((aa == 'J' || aa == 'I' || aa == 'L') &&
//                                                 (ch == 'I' || ch == 'L')) {
//                                                 aa = 'J';
//                                             } else {
//                                                 aa = 'X';
//                                             }
//                                         }
//
//                                         // lookup translation start flag
// ```

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:155-160
// ```c++
//     // ambiguous codons map to unknown amino acid or not start
//     for (i = 0; i <= 4096; i++) {
//         m_AminoAcid [i] = 'X';
//         m_OrfStart [i] = '-';
//         m_OrfStop [i] = '-';
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:228-231
// ```c++
//                 // assign amino acid
//                 if (aa != '\0') {
//                     m_AminoAcid [st] = aa;
//                 }
// ```
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
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:162-176
// ```c++
//
//         if (chunk->GetType() & fSeq) {
//             // add the sequence string
//             if (IsPositiveStrand(row)) {
//                 seq_vec.GetSeqData(chunk->GetRange().GetFrom(),
//                                    chunk->GetRange().GetTo() + 1,
//                                    buff);
//             } else {
//                 seq_vec.GetSeqData(seq_vec_size - chunk->GetRange().GetTo() - 1,
//                                    seq_vec_size - chunk->GetRange().GetFrom(),
//                                    buff);
//             }
//             if (GetWidth(row) == 3) {
//                 TranslateNAToAA(buff, buff, GetGenCode(row));
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:903-918
// ```c++
//     const CTrans_table& tbl = CGen_code_table::GetTransTable(gencode);
//
//     size_t na_size = na.size();
//
//     if (&aa != &na) {
//         aa.resize(na_size / 3);
//     }
//
//     int state = 0;
//     size_t aa_i = 0;
//     for (size_t na_i = 0; na_i < na_size; ) {
//         for (size_t i = 0; i < 3; i++) {
//             state = tbl.NextCodonState(state, na[na_i++]);
//         }
//         aa[aa_i++] = tbl.GetCodonResidue(state);
//     }
// ```

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-125
// ```c++
//         CBioseq_Handle h = GetBioseqHandle(row);
//         CSeqVector vec = h.GetSeqVector
//             (CBioseq_Handle::eCoding_Iupac,
//              IsPositiveStrand(row) ?
//              CBioseq_Handle::eStrand_Plus :
//              CBioseq_Handle::eStrand_Minus);
//         seq_vec.Reset(new CSeqVector(vec));
//         m_SeqVectorCache[row] = seq_vec;
//     }
//     if ( seq_vec->IsNucleotide() ) {
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:125-133
// ```c++
//     for (i = eBase_gap, st = 1; i <= eBase_N; i++) {
//         for (j = eBase_gap, nx = 1; j <= eBase_N; j++) {
//             for (k = eBase_gap; k <= eBase_N; k++, st++, nx += 16) {
//                 sm_NextState [st] = nx;
//                 p = sm_BaseToIdx [(int) (Uint1) baseToComp [k]];
//                 q = sm_BaseToIdx [(int) (Uint1) baseToComp [j]];
//                 r = sm_BaseToIdx [(int) (Uint1) baseToComp [i]];
//                 sm_RvCmpState [st] = 256 * p + 16 * q + r + 1;
//             }
// ```
fn display_query_residue(
    sequence: &[u8],
    frame: i8,
    offset: usize,
    code: &crate::utils::genetic_code::GeneticCode,
) -> u8 {
    let first = 3 * offset + frame.unsigned_abs() as usize - 1;
    let masks = std::array::from_fn(|i| {
        let position = if frame > 0 {
            first + i
        } else {
            sequence.len() - first - i - 1
        };
        let mask = b"-ACMGRSVTWYHKDBN"
            .iter()
            .position(|&c| c == sequence[position])
            .expect("validated IUPAC query residue") as u8;
        if frame > 0 {
            mask
        } else {
            ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
        }
    });
    display_codon(masks, code)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:200-244
// ```c++
//         case eGapAlignDecline:
//         case eGapAlignSub:
//             m_start =
//                 s_GetAlignmentStart(start1, esp->num[esp_index], m_strand, translate1, length1,
//                                     query_length, hsp->query.frame);
//
//             s_start =
//                 s_GetAlignmentStart(start2, esp->num[esp_index], s_strand, translate2, length2,
//                                     subject_length, hsp->subject.frame);
//
//             strands.push_back(m_strand);
//             strands.push_back(s_strand);
//             starts.push_back(m_start);
//             starts.push_back(s_start);
//             break;
//
//         // Insertion on the master sequence (gap on slave)
//         case eGapAlignIns:
//             m_start =
//                 s_GetAlignmentStart(start1, esp->num[esp_index], m_strand, translate1, length1,
//                                     query_length, hsp->query.frame);
//
//             s_start = GAP_VALUE;
//
//             strands.push_back(m_strand);
//             strands.push_back(esp_index == 0 ? eNa_strand_unknown : s_strand);
//             starts.push_back(m_start);
//             starts.push_back(s_start);
//             break;
//
//         // Deletion on master sequence (gap; insertion on slave)
//         case eGapAlignDel:
//             m_start = GAP_VALUE;
//
//             s_start =
//                 s_GetAlignmentStart(start2, esp->num[esp_index], s_strand, translate2, length2,
//                                     subject_length, hsp->subject.frame);
//
//             strands.push_back(esp_index == 0 ? eNa_strand_unknown : m_strand);
//             strands.push_back(s_strand);
//             starts.push_back(m_start);
//             starts.push_back(s_start);
//             break;
//
//         default:
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:996-1028
// ```c++
//             int num_matches = 0;
//             num_ident = 0;
//             // The query and subject sequence strings must be the same size in a correct
//             // alignment, but if alignment extends beyond the end of sequence because of
//             // a bug, one of the sequence strings may be truncated, hence it is
//             // necessary to take a minimum here.
//             /// @todo FIXME: Should an exception be thrown instead?
//             for (unsigned int i = 0;
//                  i < min(m_QuerySeq.size(), m_SubjectSeq.size());
//                  ++i) {
//                 if (m_QuerySeq[i] == m_SubjectSeq[i]) {
//                     ++num_ident;
//                     ++num_positives;
//                     ++num_matches;
//                 } else {
//                     if(num_matches > 0) {
//                         btop_string +=  NStr::Int8ToString(num_matches);
//                         num_matches=0;
//                     }
//                     btop_string += m_QuerySeq[i];
//                     btop_string += m_SubjectSeq[i];
//                     if (matrix && !matrix->GetData().empty() &&
//                            (*matrix)(m_QuerySeq[i], m_SubjectSeq[i]) > 0) {
//                         ++num_positives;
//                     }
//                 }
//             }
//
//             if (num_matches > 0) {
//                 btop_string +=  NStr::Int8ToString(num_matches);
//             }
//             SetBTOP(btop_string);
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2494-2529
// ```c++
//
//     if(id.Which() != CSeq_id::e_not_set){
//         /*only do this for sequence but not for others like middle line,
//           features*/
//         ITERATE(TSAlnSeqlocInfoList, iter, loc_list) {
//             int from=(*iter)->aln_range.GetFrom();
//             int to=(*iter)->aln_range.GetTo();
//             int locFrame = (*iter)->seqloc->GetFrame();
//             if(id.Match((*iter)->seqloc->GetInterval().GetId())
//                && locFrame == frame){
//                 bool isFirstChar = true;
//                 CRange<int> eachSeqloc(0, 0);
//                 //go through each residule and mask it
//                 for (int i=max<int>(from, start);
//                      i<=min<int>(to, start+len -1); i++){
//                     //store seqloc start for font tag below
//                     if ((m_AlignOption & eHtml) && isFirstChar){
//                         isFirstChar = false;
//                         eachSeqloc.Set(i, eachSeqloc.GetTo());
//                     }
//                     if (m_SeqLocChar==eX){
//                         if(isalpha((unsigned char) actualSeq[i-start])){
//                             actualSeq[i-start]='X';
//                         }
//                     } else if (m_SeqLocChar==eN){
//                         actualSeq[i-start]='n';
//                     } else if (m_SeqLocChar==eLowerCase){
//                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
//                     }
//                     //store seqloc start for font tag below
//                     if ((m_AlignOption & eHtml)
//                         && i == min<int>(to, start+len)){
//                         eachSeqloc.Set(eachSeqloc.GetFrom(), i);
//                     }
//                 }
//                 if(!(eachSeqloc.GetFrom()==0&&eachSeqloc.GetTo()==0)){
// ```
pub(crate) fn alignment(
    batch: &BatchResults,
    hsp: &Hsp,
    subject: &FastaRecord,
    query_record: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
) -> Result<DisplayAlignment> {
    let raw = &hsp.hsp;
    let context = &batch.prepared.contexts[raw.context];
    let query_length = query_record.sequence.len();
    let mut qi = usize::try_from(raw.q_start)?;
    let mut si = usize::try_from(raw.s_start)?;
    let mut qseq = Vec::new();
    let mut sseq = Vec::new();
    let mut masked = Vec::new();
    let mut query_columns = Vec::new();
    let mut gap_opens = 0;
    let script = if hsp.edit_script.is_empty() {
        vec![GapEditOp::Sub(u32::try_from(raw.q_end - raw.q_start)?)]
    } else {
        hsp.edit_script.clone()
    };
    for op in script {
        let (count, consume_q, consume_s) = match op {
            GapEditOp::Sub(n) => (n, true, true),
            GapEditOp::Ins(n) => (n, true, false),
            GapEditOp::Del(n) => (n, false, true),
        };
        if !consume_q || !consume_s {
            gap_opens += 1;
        }
        for _ in 0..count {
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:162-176
            // ```c++
            //
            //         if (chunk->GetType() & fSeq) {
            //             // add the sequence string
            //             if (IsPositiveStrand(row)) {
            //                 seq_vec.GetSeqData(chunk->GetRange().GetFrom(),
            //                                    chunk->GetRange().GetTo() + 1,
            //                                    buff);
            //             } else {
            //                 seq_vec.GetSeqData(seq_vec_size - chunk->GetRange().GetTo() - 1,
            //                                    seq_vec_size - chunk->GetRange().GetFrom(),
            //                                    buff);
            //             }
            //             if (GetWidth(row) == 3) {
            //                 TranslateNAToAA(buff, buff, GetGenCode(row));
            //             }
            // ```
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
            // ```c++
            //         alnVec->SetGapChar('-');
            //         alnVec->SetGenCode(m_QueryGeneticCode, 0);
            //         alnVec->SetGenCode(m_DbGeneticCode, 1);
            //         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
            //         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
            // ```
            let qc = if consume_q {
                display_query_residue(&query_record.sequence, raw.frame, qi, code)
            } else {
                b'-'
            };
            let sc = if consume_s {
                subject.sequence[si]
            } else {
                b'-'
            };
            if consume_q {
                query_columns.push(qseq.len());
            }
            qseq.push(qc);
            sseq.push(sc);
            masked.push(qc);
            qi += usize::from(consume_q);
            si += usize::from(consume_s);
        }
    }
    ensure!(
        qi == raw.q_end as usize && si == raw.s_end as usize,
        "BLASTX edit script does not reach final endpoints"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:930-959
    // ```c++
    //                                 (*query_interval)->GetTo());
    //         TMaskedQueryRegions query_masks;
    //         for (unsigned int index = 0; index < kNumContexts; index++) {
    //
    //             BlastSeqLoc* loc = mask->seqloc_array[qindex*kNumContexts+index];
    //             for ( ; loc; loc = loc->next) {
    //                 TSeqRange masked_range(loc->ssr->left, loc->ssr->right);
    //                 TSeqRange range(Map(kTarget, masked_range));
    //                 if (range.NotEmpty() && range != kTarget) {
    //                     int frame = BLAST_ContextToFrame(program, index);
    //                     if (frame == INT1_MAX) {
    //                         string msg("Conversion from context to frame failed ");
    //                         msg += "for '" + Blast_ProgramNameFromType(program)
    //                             + "'";
    //                         NCBI_THROW(CBlastException, eCoreBlastError, msg);
    //                     }
    //                     CRef<CSeq_interval> seqint(new CSeq_interval);
    //                     seqint->SetId().Assign((*query_interval)->GetId());
    //                     seqint->SetFrom(range.GetFrom());
    //                     seqint->SetTo(range.GetTo());
    //                     CRef<CSeqLocInfo> seqloc_info
    //                         (new CSeqLocInfo(seqint, frame));
    //                     query_masks.push_back(seqloc_info);
    //                 }
    //             }
    //         }
    //         mask_v.push_back(query_masks);
    //         qindex++;
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4305-4326
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
    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:551-558
    // ```c++
    //         // then return the edge alnpos
    //         if ((plus ? seq_pos < start : seq_pos > stop)) {
    //             return GetAlnStart(seg.GetAlnSeg());
    //         }
    //         if ((plus ? seq_pos > stop : seq_pos < start)) {
    //             return GetAlnStop(seg.GetAlnSeg());
    //         }
    //
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:593-595
    // ```c++
    //     TSeqPos delta = (seq_pos - start) / GetWidth(row);
    //     return m_AlnStarts[seg.GetAlnSeg()]
    //         + (plus ? delta : m_Lens[raw_seg] - 1 - delta);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4303-4306
    // ```c++
    //         	const CSeq_interval& interval = (*iter)->GetInterval();
    //             TSeqRange loc_range(interval.GetFrom(), interval.GetTo());
    //             if(interval.GetId().Match(m_AV->GetSeqId(i)) &&
    //                m_AV->GetSeqRange(i).IntersectingWith(loc_range)){
    // ```
    // Every consumed query residue has width three; query gaps retain columns
    // between the mapped mask endpoints. The masks here are the returned DNA
    // intervals, after the core protein-to-DNA conversion, not lookup AA masks.
    let (qstart, qend) = query_endpoints(query_length, raw.frame, raw.q_start, raw.q_end);
    for &(from, to) in &context.dna_masks {
        if from > to
            || (from == 0 && to == query_length as i32 - 1)
            || to < qstart.min(qend) - 1
            || from > qstart.max(qend) - 1
        {
            continue;
        }
        let map_nt = |nt: i32| {
            let aa = if raw.frame > 0 {
                (nt - i32::from(raw.frame) + 1).div_euclid(3)
            } else {
                (query_length as i32 + i32::from(raw.frame) - nt).div_euclid(3)
            };
            query_columns[(aa.clamp(raw.q_start, raw.q_end - 1) - raw.q_start) as usize]
        };
        let (first, last) = if raw.frame > 0 {
            (map_nt(from), map_nt(to))
        } else {
            (map_nt(to), map_nt(from))
        };
        for c in &mut masked[first..=last] {
            *c = c.to_ascii_lowercase();
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
    // ```c++
    //     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
    //
    //     SNCBIFullScoreMatrix mtx;
    //     NCBISM_Unpack(packed_mtx, &mtx);
    //
    //     for(int i = 0; i < ePMatrixSize; ++i){
    //         for(int j = 0; j < ePMatrixSize; ++j){
    //             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
    //                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
    //         }
    //     }
    //     for(int i = 0; i < ePMatrixSize; ++i) {
    //         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
    //     }
    //     retval('*', '*') = 1;
    //     // this is to count Selenocysteine to Cysteine matches as positive
    //     retval('U', 'U') = retval('C', 'C');
    //     retval('U', 'C') = retval('C', 'C');
    //     retval('C', 'U') = retval('C', 'C');
    // ```
    let mut btop = String::new();
    let mut run = 0;
    let mut identity = 0;
    let mut positive = 0;
    let mut gaps = 0;
    for (&qc, &sc) in qseq.iter().zip(&sseq) {
        if qc == sc {
            identity += 1;
            positive += 1;
            run += 1;
        } else {
            if run > 0 {
                btop.push_str(&run.to_string());
                run = 0;
            }
            btop.push(qc as char);
            btop.push(sc as char);
            if qc != b'-'
                && sc != b'-'
                && protein_display_score(crate::config::ScoringMatrix::Blosum62, qc, sc) > 0
            {
                positive += 1;
            }
        }
        gaps += usize::from(qc == b'-' || sc == b'-');
    }
    if run > 0 {
        btop.push_str(&run.to_string());
    }
    Ok(DisplayAlignment {
        qseq: String::from_utf8(qseq)?,
        sseq: String::from_utf8(sseq)?,
        masked_qseq: String::from_utf8(masked)?,
        btop,
        identity,
        positive,
        gaps,
        gap_opens,
        qstart,
        qend,
        sstart: raw.s_start + 1,
        send: raw.s_end,
        frame: raw.frame,
    })
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
// ```c++
// CRef<CSeq_id> s_ReplaceLocalId(const CBioseq_Handle& bh, CConstRef<CSeq_id> sid_in, bool parse_local)
// {
//
//         CRef<CSeq_id> retval(new CSeq_id());
//
//         // Local ids are usually fake. If a title exists, use the first token
//         // of the title instead of the local id. If no title or if the local id
//         // should be parsed, use the local id, but without the "lcl|" prefix.
//         if (sid_in->IsLocal()) {
//             string id_token;
//             vector<string> title_tokens;
//             title_tokens =
//                 NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
//             if(title_tokens.empty()){
//                 id_token = NcbiEmptyString;
//             } else {
//                 id_token = title_tokens[0];
//             }
//
//             if (id_token == NcbiEmptyString || parse_local) {
//                 const CObject_id& obj_id = sid_in->GetLocal();
//                 if (obj_id.IsStr())
//                     id_token = obj_id.GetStr();
//                 else
//                     id_token = NStr::IntToString(obj_id.GetId());
//             }
//             CObject_id* obj_id = new CObject_id();
//             obj_id->SetStr(id_token);
//             retval->SetLocal(*obj_id);
//         } else {
//             retval->Assign(*sid_in);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:175-187
// ```c++
//          if (found_gi)
//              id_string += "|";
//
//          if (best_id->IsLocal()) {
//              string id_token;
//              best_id->GetLabel(&id_token, CSeq_id::eContent, 0);
//              id_string += id_token;
//          }
//          else
//              id_string += best_id->AsFastaString();
//     }
//
//     return id_string;
// ```
fn local_id(record: &FastaRecord) -> &str {
    record
        .title
        .split(' ')
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or(&record.internal_id)
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:169-191
// ```c++
//     string id_str = NcbiEmptyString;
//
//     switch (id_type) {
//     case CBlastTabularInfo::eFullId:
//         id_str = CShowBlastDefline::GetSeqIdListString(id, true);
//         break;
//     case CBlastTabularInfo::eAccession:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, 0);
//         break;
//     }
//     case CBlastTabularInfo::eAccVersion:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, CSeq_id::fLabel_Version);
//         break;
//     }
//     case CBlastTabularInfo::eGi:
//         id_str = NStr::NumericToString(FindGi(id));
//         break;
//     default: break;
//     }
// ```
fn id_field(record: &FastaRecord, field: &str) -> String {
    // Local FASTA uses generated local Seq-ids, replaced by title tokens.
    // Keep all three source field paths distinct; local content has no accession version object.
    match field {
        "qseqid" | "sseqid" => local_id(record).to_owned(),
        "qacc" | "sacc" => local_id(record).to_owned(),
        "qaccver" | "saccver" => local_id(record).to_owned(),
        _ => unreachable!("validated identifier field"),
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1128-1196
// ```c++
//         case eQuerySeqId:
//             m_Ostream << "query id"; break;
//         case eQueryGi:
//             m_Ostream << "query gi"; break;
//         case eQueryAccession:
//             m_Ostream << "query acc."; break;
//         case eQueryAccessionVersion:
//             m_Ostream << "query acc.ver"; break;
//         case eQueryLength:
//             m_Ostream << "query length"; break;
//         case eSubjectSeqId:
//             m_Ostream << "subject id"; break;
//         case eSubjectAllSeqIds:
//             m_Ostream << "subject ids"; break;
//         case eSubjectGi:
//             m_Ostream << "subject gi"; break;
//         case eSubjectAllGis:
//             m_Ostream << "subject gis"; break;
//         case eSubjectAccession:
//             m_Ostream << "subject acc."; break;
//         case eSubjAccessionVersion:
//             m_Ostream << "subject acc.ver"; break;
//         case eSubjectAllAccessions:
//             m_Ostream << "subject accs."; break;
//         case eSubjectLength:
//             m_Ostream << "subject length"; break;
//         case eQueryStart:
//             m_Ostream << "q. start"; break;
//         case eQueryEnd:
//             m_Ostream << "q. end"; break;
//         case eSubjectStart:
//             m_Ostream << "s. start"; break;
//         case eSubjectEnd:
//             m_Ostream << "s. end"; break;
//         case eQuerySeq:
//             m_Ostream << "query seq"; break;
//         case eSubjectSeq:
//             m_Ostream << "subject seq"; break;
//         case eEvalue:
//             m_Ostream << "evalue"; break;
//         case eBitScore:
//             m_Ostream << "bit score"; break;
//         case eScore:
//             m_Ostream << "score"; break;
//         case eAlignmentLength:
//             m_Ostream << "alignment length"; break;
//         case ePercentIdentical:
//             m_Ostream << "% identity"; break;
//         case eNumIdentical:
//             m_Ostream << "identical"; break;
//         case eMismatches:
//             m_Ostream << "mismatches"; break;
//         case ePositives:
//             m_Ostream << "positives"; break;
//         case eGapOpenings:
//             m_Ostream << "gap opens"; break;
//         case eGaps:
//             m_Ostream << "gaps"; break;
//         case ePercentPositives:
//             m_Ostream << "% positives"; break;
//         case eFrames:
//             m_Ostream << "query/sbjct frames"; break;
//         case eQueryFrame:
//             m_Ostream << "query frame"; break;
//         case eSubjFrame:
//             m_Ostream << "sbjct frame"; break;
//         case eBTOP:
//             m_Ostream << "BTOP"; break;
//         case eSubjectTaxIds:
// ```
fn field_name(field: &str) -> &str {
    match field {
        "qseqid" => "query id",
        "qacc" => "query acc.",
        "qaccver" => "query acc.ver",
        "qlen" => "query length",
        "sseqid" => "subject id",
        "sacc" => "subject acc.",
        "saccver" => "subject acc.ver",
        "slen" => "subject length",
        "qstart" => "q. start",
        "qend" => "q. end",
        "sstart" => "s. start",
        "send" => "s. end",
        "qseq" => "query seq",
        "sseq" => "subject seq",
        "evalue" => "evalue",
        "bitscore" => "bit score",
        "score" => "score",
        "length" => "alignment length",
        "pident" => "% identity",
        "nident" => "identical",
        "mismatch" => "mismatches",
        "positive" => "positives",
        "gapopen" => "gap opens",
        "gaps" => "gaps",
        "ppos" => "% positives",
        "qframe" => "query frame",
        "sframe" => "sbjct frame",
        "frames" => "query/sbjct frames",
        "btop" => "BTOP",
        "stitle" => "subject title",
        _ => unreachable!("validated BLASTX field"),
    }
}
// NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:437-479
// ```c++
// inline void CBlastTabularInfo::x_PrintPercentIdentical(void)
// {
//     double perc_ident =
//         (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
//     m_Ostream << NStr::DoubleToString(perc_ident, 3);
// }
//
// inline void CBlastTabularInfo::x_PrintPercentPositives(void)
// {
//     double perc_positives =
//         (m_AlignLength > 0 ? ((double)m_NumPositives)/m_AlignLength * 100 : 0);
//     m_Ostream << NStr::DoubleToString(perc_positives, 2);
// }
//
// inline void CBlastTabularInfo::x_PrintFrames(void)
// {
//     m_Ostream << m_QueryFrame << "/" << m_SubjectFrame;
// }
//
// inline void CBlastTabularInfo::x_PrintQueryFrame(void)
// {
//     m_Ostream << m_QueryFrame;
// }
//
// inline void CBlastTabularInfo::x_PrintSubjectFrame(void)
// {
//     m_Ostream << m_SubjectFrame;
// }
//
// inline void CBlastTabularInfo::x_PrintBTOP(void)
// {
//     m_Ostream << m_BTOP;
// }
//
// inline void CBlastTabularInfo::x_PrintNumIdentical(void)
// {
//     m_Ostream << m_NumIdent;
// }
//
// inline void CBlastTabularInfo::x_PrintMismatches(void)
// {
//     int num_mismatches = m_AlignLength - m_NumIdent - m_NumGaps;
//     m_Ostream << num_mismatches;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:415-436
// ```c++
// void CBlastTabularInfo::x_PrintSubjectTitle()
// {
// 	if(m_SubjectDefline.NotEmpty() && m_SubjectDefline->CanGet() &&
// 	   m_SubjectDefline->IsSet() && !m_SubjectDefline->Get().empty())
// 	{
// 		const list<CRef<CBlast_def_line> > & defline = m_SubjectDefline->Get();
//
// 		if(defline.empty())
// 			m_Ostream << NA;
// 		else
// 		{
// 			if(defline.front()->IsSetTitle())
// 			{
// 				if(defline.front()->GetTitle().empty())
// 					m_Ostream << NA;
// 				else
// 					m_Ostream << defline.front()->GetTitle();
// 			}
// 			else
// 				m_Ostream << NA;
// 		}
// 	}
// ```
fn field_value(
    field: &str,
    query: &FastaRecord,
    subject: &FastaRecord,
    hsp: &Hsp,
    a: &DisplayAlignment,
) -> String {
    match field {
        "qseqid" | "qacc" | "qaccver" => id_field(query, field),
        "sseqid" | "sacc" | "saccver" => id_field(subject, field),
        "qlen" => query.sequence.len().to_string(),
        "slen" => subject.sequence.len().to_string(),
        "qstart" => a.qstart.to_string(),
        "qend" => a.qend.to_string(),
        "sstart" => a.sstart.to_string(),
        "send" => a.send.to_string(),
        "qseq" => a.qseq.clone(),
        "sseq" => a.sseq.clone(),
        "evalue" => format_evalue_ncbi_tabular(hsp.evalue),
        "bitscore" => format_bitscore_ncbi(hsp.bit_score),
        "score" => hsp.hsp.score.to_string(),
        "length" => a.qseq.len().to_string(),
        "pident" => format_percent_identity_ncbi(a.identity, a.qseq.len(), 3),
        "nident" => a.identity.to_string(),
        "mismatch" => (a.qseq.len() - a.identity - a.gaps).to_string(),
        "positive" => a.positive.to_string(),
        "gapopen" => a.gap_opens.to_string(),
        "gaps" => a.gaps.to_string(),
        "ppos" => format_percent_identity_ncbi(a.positive, a.qseq.len(), 2),
        "qframe" => a.frame.to_string(),
        "sframe" => "0".into(),
        "frames" => format!("{}/0", a.frame),
        "btop" => a.btop.clone(),
        // Local-subject FASTA has no Blast-def-line-set; the source title field is NA.
        "stitle" => "N/A".into(),
        _ => unreachable!("validated BLASTX field"),
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1578
// ```c++
//     for (int index = 0; index < hit_list->hsplist_count; index++) {
//         BlastHSPList* hsp_list = hit_list->hsplist_array[index];
//         if (!hsp_list)
//             continue;
//
//         // Sort HSPs with e-values as first priority and scores as
//         // tie-breakers, since that is the order we want to see them in
//         // in Seq-aligns.
//         Blast_HSPListSortByEvalue(hsp_list);
//
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1265-1283
// ```c++
// void
// CBlastTabularInfo::PrintHeader(const string& program_version,
//        const CBioseq& bioseq,
//        const string& dbname,
//        const string& rid /* = kEmptyStr */,
//        unsigned int iteration /* = numeric_limits<unsigned int>::max() */,
//        const CSeq_align_set* align_set /* = 0 */,
//        CConstRef<CBioseq> subj_bioseq /* = CConstRef<CBioseq>() */,
//        bool is_csv /* = false */)
// {
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
pub fn render(
    writer: &mut impl Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    subject_path: &Path,
    batches: &mut [BatchResults],
) -> Result<()> {
    write_prolog(writer, options, subjects, subject_path)?;
    render_queries(
        writer,
        &mut std::io::stderr().lock(),
        options,
        queries,
        subjects,
        subject_path,
        batches,
    )?;
    write_epilog(writer, options, subjects, subject_path, queries.len())
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
// ```c++
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// ```c++
//         formatter.PrintEpilog(opt);
// ```
pub fn render_queries(
    writer: &mut impl Write,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    subject_path: &Path,
    batches: &mut [BatchResults],
) -> Result<()> {
    for batch in batches.iter_mut() {
        for lists in &mut batch.queries {
            for list in lists {
                list.sort_evalue();
            }
        }
    }
    if options.outfmt == 0 {
        return write_pairwise(
            writer,
            stderr,
            options,
            queries,
            subjects,
            subject_path,
            batches,
        );
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
    // ```c++
    //         alnVec->SetGapChar('-');
    //         alnVec->SetGenCode(m_QueryGeneticCode, 0);
    //         alnVec->SetGenCode(m_DbGeneticCode, 1);
    //         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
    //         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
    // ```
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    let comments = options.outfmt == 7;
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            // ```c++
            //     if (results.HasWarnings()) {
            //         ERR_POST(Warning << results.GetWarningStrings());
            //     }
            // ```
            stderr.write_all(&query_warning(query, batch.search_skipped))?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:779-783
            // ```c++
            //               m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader)
            //              ? CBlastTabularInfo::eComma : CBlastTabularInfo::eTab);
            //
            //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
            //         if(!m_CustomDelim.empty()) {
            // ```
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
            // ```c++
            // CBlastTabularInfo::~CBlastTabularInfo()
            // {
            //     m_Ostream.flush();
            // }
            // ```
            // Tabular stream failures unwind through this local formatter.
            // Its destructor flushes the already-failed stream and terminates,
            // including failures before the final explicit flush is reached.
            (|| -> Result<()> {
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:762-765
                // ```c++
                //     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
                //     if (m_IsUngappedSearch && results.HasAlignments()) {
                //         aln_set.Reset(CDisplaySeqalign::PrepareBlastUngappedSeqalign(*aln_set));
                //     }
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3162-3185
                // ```c++
                // CRef<CSeq_align_set>
                // CDisplaySeqalign::PrepareBlastUngappedSeqalign(const CSeq_align_set& alnset)
                // {
                //     CRef<CSeq_align_set> alnSetRef(new CSeq_align_set);
                //
                //     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()){
                //         const CSeq_align::TSegs& seg = (*iter)->GetSegs();
                //         if(seg.Which() == CSeq_align::C_Segs::e_Std){
                //             if(seg.GetStd().size() > 1){
                //                 //has more than one stdseg. Need to seperate as each
                //                 //is a distinct HSP
                //                 ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
                //                     CRef<CSeq_align> aln(new CSeq_align);
                //                     if((*iterStdseg)->IsSetScores()){
                //                         aln->SetScore() = (*iterStdseg)->GetScores();
                //                     }
                //                     aln->SetSegs().SetStd().push_back(*iterStdseg);
                //                     alnSetRef->Set().push_back(aln);
                //                 }
                //
                //             } else {
                //                 alnSetRef->Set().push_back(*iter);
                //             }
                //         } else if(seg.Which() == CSeq_align::C_Segs::e_Dendiag){
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
                // ```c++
                //     if (align_set) {
                //        int num_hits = align_set->Get().size();
                // ```
                let count: usize = lists
                    .iter()
                    .map(|l| {
                        l.hsps
                            .iter()
                            .filter(|h| has_seqalign(h, options.gapped))
                            .count()
                    })
                    .sum();
                if comments {
                    writeln!(writer, "# BLASTX 2.17.0+")?;
                    writeln!(
                        writer,
                        "# Query: {}",
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:630-637
                        // ```c++
                        // CAlignFormatUtil::GetSeqIdString(const list<CRef<CSeq_id> > & ids, bool believe_local_id)
                        // {
                        //     string all_id_str = NcbiEmptyString;
                        //     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
                        //
                        //     if (wid && (wid->Which()!= CSeq_id::e_Local || believe_local_id)){
                        //         TGi gi = FindGi(ids);
                        //
                        // ```
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-739
                        // ```c++
                        //     string all_id_str = GetSeqIdString(cbs, believe_query);
                        //     all_id_str += " ";
                        //     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
                        //
                        //     // For tabular output, there is no limit on the line length.
                        //     // There is also no extra line with the sequence length.
                        //     if (tabular) {
                        //         out << all_id_str;
                        // ```
                        &query.title
                    )?;
                    writeln!(
                        writer,
                        "# Database: User specified sequence set (Input: {})",
                        subject_path.display()
                    )?;
                    if count > 0 {
                        writeln!(
                            writer,
                            "# Fields: {}",
                            options
                                .fields
                                .iter()
                                .map(|f| field_name(f))
                                .collect::<Vec<_>>()
                                .join(", ")
                        )?;
                    }
                    if !batch.search_skipped {
                        writeln!(writer, "# {count} hits found")?;
                    }
                }
                for list in lists {
                    let subject = &subjects[list.oid];
                    for hsp in &list.hsps {
                        if !has_seqalign(hsp, options.gapped) {
                            continue;
                        }
                        let a = alignment(batch, hsp, subject, query, &code)?;
                        let values: Vec<_> = options
                            .fields
                            .iter()
                            .map(|f| field_value(f, query, subject, hsp, &a))
                            .collect();
                        writeln!(writer, "{}", values.join("	"))?;
                    }
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
                // ```c++
                // CBlastTabularInfo::~CBlastTabularInfo()
                // {
                //     m_Ostream.flush();
                // }
                // ```
                writer.flush().map_err(TabularFlushError)?;
                Ok(())
            })()
            .map_err(tabular_io_error)?;
        }
    }
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1420
// ```c++
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
//                         CConstRef<blast::CBlastQueryVector> queries,
//                         unsigned int itr_num
//                         /* = numeric_limits<unsigned int>::max() */,
//                         blast::CPsiBlastIterationState::TSeqIds prev_seqids
//                         /* = CPsiBlastIterationState::TSeqIds() */,
//                         bool is_deltablast_domain_result /* = false */)
// {
//     // For remote searches, we don't retrieve the sequence data for the query
//     // sequence when initially sending the request to the BLAST server (if it's
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1578
// ```c++
//     for (int index = 0; index < hit_list->hsplist_count; index++) {
//         BlastHSPList* hsp_list = hit_list->hsplist_array[index];
//         if (!hsp_list)
//             continue;
//
//         // Sort HSPs with e-values as first priority and scores as
//         // tie-breakers, since that is the order we want to see them in
//         // in Seq-aligns.
//         Blast_HSPListSortByEvalue(hsp_list);
//
// ```
fn write_pairwise(
    writer: &mut impl Write,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    subject_path: &Path,
    batches: &[BatchResults],
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
    // ```c++
    //         alnVec->SetGapChar('-');
    //         alnVec->SetGenCode(m_QueryGeneticCode, 0);
    //         alnVec->SetGenCode(m_DbGeneticCode, 1);
    //         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
    //         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
    // ```
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    use crate::report::pairwise::{
        write_blastx_pairwise_report, BlastpPairwiseQuery, BlastpPairwiseReport,
        BlastxPairwiseOptions, PairwiseConfig, PairwiseHit,
    };
    use crate::{common::Hit, config::ScoringMatrix};
    let mut hits = Vec::new();
    let mut qinfo = Vec::new();
    let mut valid = Vec::new();
    let mut skipped = Vec::new();
    let mut warnings = Vec::new();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let absolute = batch.query_ordinal + qidx;
            let query = &queries[absolute];
            let first = (qidx * 6..qidx * 6 + 6).find(|&c| batch.parameters[c].valid);
            let stat = first.map(|c| &batch.parameters[c]);
            qinfo.push(BlastpPairwiseQuery {
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:630-637
                // ```c++
                // CAlignFormatUtil::GetSeqIdString(const list<CRef<CSeq_id> > & ids, bool believe_local_id)
                // {
                //     string all_id_str = NcbiEmptyString;
                //     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
                //
                //     if (wid && (wid->Which()!= CSeq_id::e_Local || believe_local_id)){
                //         TGi gi = FindGi(ids);
                //
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-743
                // ```c++
                //     string all_id_str = GetSeqIdString(cbs, believe_query);
                //     all_id_str += " ";
                //     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
                //
                //     // For tabular output, there is no limit on the line length.
                //     // There is also no extra line with the sequence length.
                //     if (tabular) {
                //         out << all_id_str;
                //     } else {
                //         x_WrapOutputLine(all_id_str, line_len, out, html);
                //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
                //             out << "\nLength=";
                // ```
                query_name: query.title.clone(),
                query_length: query.sequence.len(),
                ungapped_karlin: stat.map_or(batch.parameters[qidx * 6].ungapped, |p| p.ungapped),
                effective_search_space: stat.map_or(0, |p| p.search_space),
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            // ```c++
            //     if (results.HasWarnings()) {
            //         ERR_POST(Warning << results.GetWarningStrings());
            //     }
            // ```
            warnings.push(query_warning(query, batch.search_skipped));
            valid.push(first.is_some());
            skipped.push(batch.search_skipped);
            for list in lists {
                let subject = &subjects[list.oid];
                for hsp in &list.hsps {
                    if !has_seqalign(hsp, options.gapped) {
                        continue;
                    }
                    let a = alignment(batch, hsp, subject, query, &code)?;
                    let raw = &hsp.hsp;
                    let len = a.qseq.len();
                    let hit = Hit {
                        identity: 100.0 * a.identity as f64 / len as f64,
                        length: len,
                        mismatch: len - a.identity - a.gaps,
                        gapopen: a.gap_opens,
                        q_start: usize::try_from(a.qstart)?,
                        q_end: usize::try_from(a.qend)?,
                        s_start: usize::try_from(a.sstart)?,
                        s_end: usize::try_from(a.send)?,
                        e_value: hsp.evalue,
                        bit_score: hsp.bit_score,
                        num_ident: a.identity,
                        query_frame: i32::from(raw.frame),
                        query_length: query.sequence.len(),
                        q_idx: u32::try_from(absolute)?,
                        s_idx: u32::try_from(list.oid)?,
                        raw_score: raw.score,
                        sort_query_offset: usize::try_from(raw.q_start)?,
                        sort_query_end: usize::try_from(raw.q_end)?,
                        sort_subject_offset: usize::try_from(raw.s_start)?,
                        sort_subject_end: usize::try_from(raw.s_end)?,
                        has_sort_offsets: true,
                        gap_info: Some(hsp.edit_script.clone()),
                        num_positives: a.positive,
                    };
                    hits.push(PairwiseHit {
                        hit,
                        query_seq: Some(a.masked_qseq),
                        subject_seq: Some(a.sseq),
                        query_frame: Some(raw.frame),
                        subject_frame: None,
                        positives: Some(a.positive),
                        gaps: Some(a.gaps),
                        subject_length: Some(subject.sequence.len()),
                        subject_title: subject.title.split_once(' ').map(|(_, t)| t.to_owned()),
                        comp_adjust_method: Some(u8::try_from(hsp.composition_method)?),
                        sum_n: Some(hsp.num),
                    });
                }
            }
        }
    }
    let report = pairwise_metadata(options, subjects, subject_path);
    let ids: Vec<std::sync::Arc<str>> = subjects
        .iter()
        .map(|s| std::sync::Arc::from(local_id(s)))
        .collect();
    let config = PairwiseConfig {
        program: "blastx".into(),
        protein_matrix: ScoringMatrix::Blosum62,
        ..PairwiseConfig::default()
    };
    write_blastx_pairwise_report(
        &hits,
        writer,
        stderr,
        &warnings,
        &config,
        &qinfo,
        &valid,
        &skipped,
        &ids,
        &report,
        &BlastxPairwiseOptions {
            word_threshold: options.threshold,
            gapped: options.gapped,
            sum_stats: options.sum_stats,
            num_descriptions: options.num_descriptions as usize,
            num_alignments: options.num_alignments as usize,
        },
    )?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:226-244
// ```c++
//         CBlastFormat formatter(opt, *db_adapter,
//                                fmt_args->GetFormattedOutputChoice(),
//                                query_opts->GetParseDeflines(),
//                                m_CmdLineArgs->GetOutputStream(),
//                                fmt_args->GetNumDescriptions(),
//                                fmt_args->GetNumAlignments(),
//                                *scope,
//                                opt.GetMatrixName(),
//                                fmt_args->ShowGis(),
//                                fmt_args->DisplayHtmlOutput(),
//                                opt.GetQueryGeneticCode(),
//                                opt.GetDbGeneticCode(),
//                                opt.GetSumStatisticsMode(),
//                                m_CmdLineArgs->ExecuteRemotely(),
//                                db_adapter->GetFilteringAlgorithm(),
//                                fmt_args->GetCustomOutputFormatSpec(),
//                                false, false, NULL, NULL,
//                                GetCmdlineArgs(GetArguments()),
// 				GetSubjectFile(args));
// ```
fn pairwise_metadata(
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    subject_path: &Path,
) -> crate::report::pairwise::BlastpPairwiseReport {
    use crate::{
        config::{ProteinScoringSpec, ScoringMatrix},
        report::pairwise::BlastpPairwiseReport,
        stats::{spouge::lookup_protein_gumbel_params, tables::lookup_protein_params},
    };
    let scoring = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    };
    let total: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    BlastpPairwiseReport {
        version: "2.17.0+".into(),
        database_name: format!(
            "User specified sequence set (Input: {})",
            subject_path.display()
        ),
        database_num_sequences: subjects.len(),
        database_total_letters: total,
        matrix_name: options.matrix.clone(),
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
        word_threshold: 0.0,
        window_size: options.window_size,
        gapped_karlin: lookup_protein_params(&scoring),
        gumbel: lookup_protein_gumbel_params(&scoring, total as i64)
            .expect("pinned BLOSUM62 gap11/1 Gumbel block"),
        // BLASTX writes its own description table and alignments (`write_blastx_report`).
        num_descriptions: usize::MAX,
        num_alignments: usize::MAX,
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
pub fn write_prolog(
    writer: &mut impl Write,
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    subject_path: &Path,
) -> Result<()> {
    if options.outfmt == 0 {
        crate::report::pairwise::write_blastx_prolog(
            writer,
            &pairwise_metadata(options, subjects, subject_path),
        )?;
    }
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2223-2245
// ```c++
//     	if (m_FormatType == CFormattingArgs::eXml2
//     		|| m_FormatType == CFormattingArgs::eXml2_S) {
//     		x_GenerateXML2MasterFile();
//     	}
//     	else {
//     		x_GenerateJSONMasterFile();
//     	}
//     	return;
//     }
//
//     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
//         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
//         tabinfo.PrintNumProcessed(m_QueriesFormatted);
//         return;
//     } else if (m_FormatType >= CFormattingArgs::eTabular)
//         return;  // No footer for these.
//
//     // Most of XML is printed as it's finished.
//     // the epilog closes the report.
//     if (m_FormatType == CFormattingArgs::eXml) {
//         m_Outfile << m_BlastXMLIncremental->m_SerialXmlEnd << endl;
//         m_AccumulatedResults.clear();
//         m_AccumulatedQueries->clear();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
// ```c++
//                         options.GetMismatchPenalty() << "\n";
//     }
//     else {
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
//
//     if (options.GetGappedMode() == true) {
//         double gap_extension = (double) options.GetGapExtensionCost();
//         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
//         { // Formula from PMID 10890397 applies if both gap values are zero.
//                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
//                gap_extension /= 2.0;
//         }
//         m_Outfile << "Gap Penalties: Existence: "
//                 << options.GetGapOpeningCost() << ", Extension: "
//                 << gap_extension << "\n";
//     }
//     if (options.GetWordThreshold()) {
//         m_Outfile << "Neighboring words threshold: " <<
//                         options.GetWordThreshold() << "\n";
//     }
//     if (options.GetWindowSize()) {
//         m_Outfile << "Window for multiple hits: " <<
//                         options.GetWindowSize() << "\n";
//     }
// ```
pub fn write_epilog(
    writer: &mut impl Write,
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    subject_path: &Path,
    queries: usize,
) -> Result<()> {
    if options.outfmt == 0 {
        crate::report::pairwise::write_blastx_epilog(
            writer,
            &pairwise_metadata(options, subjects, subject_path),
            &crate::report::pairwise::BlastxPairwiseOptions {
                word_threshold: options.threshold,
                gapped: options.gapped,
                sum_stats: options.sum_stats,
                num_descriptions: options.num_descriptions as usize,
                num_alignments: options.num_alignments as usize,
            },
        )?;
    } else if options.outfmt == 7 {
        writeln!(writer, "# BLAST processed {queries} queries")?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
        // ```c++
        //     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
        //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
        //         tabinfo.PrintNumProcessed(m_QueriesFormatted);
        //         return;
        // ```
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
        // ```c++
        // CBlastTabularInfo::~CBlastTabularInfo()
        // {
        //     m_Ostream.flush();
        // }
        // ```
        writer.flush().map_err(TabularFlushError)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:534-543
// ```c++
//                 const string kTitle = queries.GetTitle(index);
//                 string query_id = id->GetSeqIdString();
//                 if (kTitle != kEmptyStr) {
//                     query_id += " " + kTitle;
//                 }
//                  if(query_id.size() > 35) {
//                 	 query_id = query_id.substr(0, 25) + ".. ";
//                  }
//
//                 messages[index].SetQueryId(query_id);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:282-292
// ```c++
//
//     string retval(m_Errors.GetQueryId());
//     if ( !retval.empty() ) {    // in case the query id is not known
//         retval += ": ";
//     }
//     ITERATE(TQueryMessages, iter, m_Errors) {
//         if ((**iter).GetSeverity() == eBlastSevWarning) {
//             retval += (*iter)->GetMessage(false) + " ";
//         }
//     }
//     return retval;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:37-40
// ```c++
// const char* kBlastErrMsg_CantCalculateUngappedKAParams
//     = "Could not calculate ungapped Karlin-Altschul parameters due "
//       "to an invalid query sequence or its translation. Please verify the "
//       "query sequence(s) and/or filtering options";
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2815-2822
// ```c++
//    if (valid_context == FALSE)
//    {   /* No valid contexts were found. */
//        /* Message for non-translated search issued above. */
//        if (Blast_QueryIsTranslated(program) ) {
//             Blast_MessageWrite(blast_message, eBlastSevWarning, kBlastMessageNoContext,
//             kBlastErrMsg_CantCalculateUngappedKAParams);
//        }
//        status = 1;  /* Not a single context was valid. */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:94-101
// ```c++
//         } else {
//             // applies to all queries
//             CRef<CSearchMessage> sm(new CSearchMessage(blmsg->severity,
//                                                        kBlastMessageNoContext,
//                                                        msg));
//             NON_CONST_ITERATE(TSearchMessages, query_messages, messages) {
//                 query_messages->push_back(sm);
//             }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:634-639
// ```c++
//             // to determine whether the message should contain a warning or
//             // error?
//             CRef<CSearchMessage> m
//                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
//             messages[index].push_back(m);
//             s_InvalidateQueryContexts(qinfo, index);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1042-1052
// ```c++
// void
// TSearchMessages::RemoveDuplicates()
// {
//     NON_CONST_ITERATE(TSearchMessages, sm, (*this)) {
//         if (sm->empty()) {
//             continue;
//         }
//         sort(sm->begin(), sm->end(), TQueryMessagesLessComparator());
//         TQueryMessages::iterator new_end =
//             unique(sm->begin(), sm->end(), TQueryMessagesEqualComparator());
//         sm->erase(new_end, sm->end());
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:294-303
// ```c++
// inline bool
// CSearchMessage::operator<(const CSearchMessage& rhs) const
// {
//     if (m_ErrorId < rhs.m_ErrorId ||
//         m_Severity < rhs.m_Severity ||
//         m_Message < rhs.m_Message) {
//         return true;
//     } else {
//         return false;
//     }
// ```
// Global errorId0 "Could..." precedes per-query "Sequence..." after Combine.
fn query_warning(query: &FastaRecord, search_skipped: bool) -> Vec<u8> {
    if !query.sequence.is_empty() && !search_skipped {
        return Vec::new();
    }
    let mut id = query.internal_id.as_bytes().to_vec();
    if !query.title.is_empty() {
        id.push(b' ');
        id.extend_from_slice(query.title.as_bytes());
    }
    if id.len() > 35 {
        id.truncate(25);
        id.extend_from_slice(b".. ");
    }
    let mut bytes = b"Warning: [blastx] ".to_vec();
    bytes.extend_from_slice(&id);
    bytes.extend_from_slice(b": ");
    if search_skipped {
        bytes.extend_from_slice(b"Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options ");
    }
    if query.sequence.is_empty() {
        bytes.extend_from_slice(b"Sequence contains no data ");
    }
    bytes.push(b'\n');
    bytes
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c++
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
fn tabular_io_error(error: anyhow::Error) -> anyhow::Error {
    match error.downcast::<std::io::Error>() {
        Ok(error) => TabularFlushError(error).into(),
        Err(error) => error,
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c++
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
#[derive(Debug)]
pub(crate) struct TabularFlushError(pub std::io::Error);
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
impl std::fmt::Display for TabularFlushError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
impl std::error::Error for TabularFlushError {}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:123-135
// ```c++
//         if (translate)
//             retval = original_length -
//                 CODON_LENGTH*(s_GetCurrPos(curr_pos, num) + num)
//                 + frame + 1;
//         else
//             retval = length - s_GetCurrPos(curr_pos, num) - num;
//
//     } else {
//
//         if (translate)
//             retval = frame - 1 + CODON_LENGTH*s_GetCurrPos(curr_pos, num);
//         else
//             retval = s_GetCurrPos(curr_pos, num);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1036-1042
// ```c++
//         if (kTranslated && ds.GetSeqStrand(kQueryRow) == eNa_strand_minus) {
//             q_start = alnVec->GetSeqStop(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStart(kQueryRow) + 1;
//         } else {
//             q_start = alnVec->GetSeqStart(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStop(kQueryRow) + 1;
//         }
// ```
#[cfg(test)]
mod tests {
    use super::query_endpoints;
    // NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:195-209
    // ```c++
    //                                             if ((aa == 'B' || aa == 'D' || aa == 'N') &&
    //                                                 (ch == 'D' || ch == 'N')) {
    //                                                 aa = 'B';
    //                                             } else if ((aa == 'Z' || aa == 'E' || aa == 'Q') &&
    //                                                        (ch == 'E' || ch == 'Q')) {
    //                                                 aa = 'Z';
    //                                             } else if ((aa == 'J' || aa == 'I' || aa == 'L') &&
    //                                                 (ch == 'I' || ch == 'L')) {
    //                                                 aa = 'J';
    //                                             } else {
    //                                                 aa = 'X';
    //                                             }
    //                                         }
    //
    //                                         // lookup translation start flag
    // ```
    // Expected values are the verbatim C++ function over all 26 BLASTX
    // genetic-code tables and 4,096 ncbi4na triplets, not Rust output.
    #[test]
    fn display_codons_match_pinned_cpp_all_26_codes() {
        let expected =
            include_str!("../../../tests/unit/blastx_stage_e_display_codon_expected.tsv");
        let mut last_code = 0;
        let mut code = crate::utils::genetic_code::GeneticCode::from_id(1);
        let mut count = 0;
        for row in expected.lines() {
            let values: Vec<_> = row.split('\t').collect();
            let id: u8 = values[0].parse().unwrap();
            if id != last_code {
                code = crate::utils::genetic_code::GeneticCode::from_id(id);
                last_code = id;
            }
            let masks = [
                values[1].parse().unwrap(),
                values[2].parse().unwrap(),
                values[3].parse().unwrap(),
            ];
            assert_eq!(
                super::display_codon(masks, &code),
                values[4].as_bytes()[0],
                "{row}"
            );
            count += 1;
        }
        assert_eq!(count, 106496);
    }
    #[test]
    fn translated_endpoints_match_independent_ncbi_cpp_boundaries() {
        for row in
            include_str!("../../../tests/unit/blastx_stage_e_coordinates_expected.tsv").lines()
        {
            let f: Vec<i32> = row.split('\t').map(|x| x.parse().unwrap()).collect();
            assert_eq!(
                query_endpoints(f[0] as usize, f[1] as i8, f[2], f[3]),
                (f[4], f[5]),
                "{row}"
            );
        }
    }
}

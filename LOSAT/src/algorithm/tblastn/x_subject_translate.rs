//! EXPERIMENT (LOSAT_X_TBNXLATE / LOSAT_X_TBNXLATESHADOW): the six translated frames of a
//! TBLASTN subject read from tables, as NCBI makes them from the compressed subject.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1071-1074,1080-1090
//! ```c
//!    } else {
//!       translation_table = s_BlastGetTranslationTable(genetic_code, FALSE);
//!       translation_table_rc = s_BlastGetTranslationTable(genetic_code, TRUE);
//!    }
//! ...
//!    for (context = 0; context < NUM_FRAMES; ++context) {
//!       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//!       if (encoding == eBlastEncodingNcbi2na) {
//!          if (frame > 0) {
//!             length =
//!                BLAST_TranslateCompressedSequence(translation_table,
//!                   nucl_length, nucl_seq, frame, translation_buffer+offset);
//!          } else {
//!             length =
//!                BLAST_TranslateCompressedSequence(translation_table_rc,
//!                   nucl_length, nucl_seq, frame, translation_buffer+offset);
//! ```
//! The preliminary search reads the subject in ncbi2na, so `BLAST_GetAllTranslations` takes this
//! branch: one table of 64 residues per strand and one table read per codon. The reference path
//! of LOSAT (`resolve_local_subject_ncbi2na`, then `tblastx::translation::generate_frames`)
//! reaches the same residues through letters: the subject as ncbi4na, as ncbi2na, as "ACGT", its
//! reverse complement, and `GeneticCode::get` per codon. With `LOSAT_X_TBNXLATE` the subject is
//! turned into ncbi2na in one pass and the frames are read from the two tables. The tables are
//! filled from the reference path's own functions, so every residue is the reference's by
//! construction; the ambiguous letters draw from the same generator with the same seed in the
//! same order.
//!
//! `LOSAT_X_TBNXLATE` is a switch of the strict set (`utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. `LOSAT_X_TBNXLATESHADOW` is not: it counts when set to
//! any value and wins over the switch. They act where the subject plans of the unit path are made
//! (`search_seed::x_scan_plan`, `LOSAT_X_TBNPAR`, on by default); the reference loop that
//! searches a subject without that switch translates it itself and reads neither variable.
//! `LOSAT_X_TBNXLATESHADOW` makes both translations for every such plan and stops the search
//! if a frame, a unit or an error differs.

use super::search_seed::{translated_chunk_scan_ranges, XSubjectPlan};
use crate::algorithm::tblastx::constants::SENTINEL_BYTE;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::core::blast_encoding::{iupacna_to_ncbi4na, NcbiRandom};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbistdaa;
use anyhow::{bail, Result};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

/// EXPERIMENT (LOSAT_X_TBNXLATE / LOSAT_X_TBNXLATESHADOW): which translation makes a subject
/// plan.
// No NCBI counterpart: names the LOSAT code path; it does not change any value NCBI computes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XTranslate {
    /// The reference path.
    Reference,
    /// The tables.
    Tables,
    /// Both; the search stops if they differ.
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_TBNXLATE / LOSAT_X_TBNXLATESHADOW variables once; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNXLATE / LOSAT_X_TBNXLATESHADOW): `LOSAT_X_TBNXLATE` is a strict-set
/// switch (on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set); the shadow counts when
/// set to any value and wins.
pub(super) fn x_translate() -> XTranslate {
    static MODE: OnceLock<XTranslate> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBNXLATESHADOW").is_some() {
            XTranslate::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_TBNXLATE") {
            XTranslate::Tables
        } else {
            XTranslate::Reference
        }
    })
}

// No NCBI counterpart: counters for the shadow summary and LOSAT_X_STATS that tell which
// translation ran; they do not change any value NCBI computes.
/// Subject plans made from the tables (in shadow mode the table half of each comparison).
static TABLE_PLANS: AtomicU64 = AtomicU64::new(0);
/// Shadow mode: the subjects whose two plans were compared, the translated residues in them
/// (sentinels included), and the subjects that both translations refused with the same error.
static SHADOW_PLANS: AtomicU64 = AtomicU64::new(0);
static SHADOW_RESIDUES: AtomicU64 = AtomicU64::new(0);
static SHADOW_ERRORS: AtomicU64 = AtomicU64::new(0);

/// The summary lines (printed by `main` at exit): the shadow summary, and with LOSAT_X_STATS
/// the number of subject plans made from the tables.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if x_translate() == XTranslate::Shadow {
        let plans = SHADOW_PLANS.load(Ordering::Relaxed);
        let errors = SHADOW_ERRORS.load(Ordering::Relaxed);
        eprintln!(
            "[X_SHADOW] TBNXLATE plans_compared={plans} residues_compared={} errors_compared={errors} ({})",
            SHADOW_RESIDUES.load(Ordering::Relaxed),
            if plans + errors == 0 {
                "nothing compared: no subject plan of the unit path was made"
            } else {
                "frames, units and errors identical"
            }
        );
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBNXLATE table_plans={}",
            TABLE_PLANS.load(Ordering::Relaxed)
        );
    }
}

// No NCBI counterpart: the two marks of `x_letter_table` below, which joins LOSAT's read of a
// letter to NCBI's `ctable` read; they are no value NCBI computes.
/// Entry of `x_letter_table` for a byte that is no subject letter.
const X_NOT_A_LETTER: u8 = 0xFF;
/// Flag of an entry of `x_letter_table` that holds the ncbi4na value of a gap or an ambiguous
/// letter.
const X_AMBIGUOUS: u8 = 0x10;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:424-425,441-445
// ```c++
// static unsigned char ctable[16] = {0xFF, 0x00, 0x01, 0xFF, 0x02, 0xFF, 0xFF, 0xFF,
// 		                           0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
// ...
//         unsigned char c = ctable[b];
//
//         if (c  != 0xFF) {
//             // No ambiguities, so we can do this the easy way.
//         	ncbi2na[i] = c;
// ```
// NCBI reads the ncbi2na value of an unambiguous ncbi4na value from `ctable`. This table joins
// that read to the one that precedes it in LOSAT, letter to ncbi4na: per byte, the ncbi2na
// value of `A`, `C`, `G`, `T` (either case), the flagged ncbi4na value of a gap or an ambiguous
// letter, or `X_NOT_A_LETTER`. It is filled from the reference path's own
// `iupacna_to_ncbi4na`, with the reference path's rule for `U`.
/// EXPERIMENT (LOSAT_X_TBNXLATE): what a subject byte is, in one read.
fn x_letter_table() -> &'static [u8; 256] {
    static TABLE: OnceLock<[u8; 256]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [X_NOT_A_LETTER; 256];
        for byte in 0..=u8::MAX {
            // TBLASTN reads the IUPAC letters without `U`.
            let ncbi4na = iupacna_to_ncbi4na(byte).filter(|_| !byte.eq_ignore_ascii_case(&b'U'));
            if let Some(mask) = ncbi4na {
                table[byte as usize] = match mask {
                    1 => 0,
                    2 => 1,
                    4 => 2,
                    8 => 3,
                    _ => X_AMBIGUOUS | mask,
                };
            }
        }
        table
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:448-455,461-475
// ```c++
//             if (b == 0 || b == 0x0F) {
//             	//gap or N
//                 ncbi2na[i] = random.GetRand() & 0x3;
//             }
//             else {
//
//             	int bitcount = ((b & 1) + ((b >> 1) & 1) +
//             	               ((b >> 2) & 1) + ((b >> 3) & 1));
// ...
//             	int pick = random.GetRand() % bitcount;
//
//             	for(int j = 0; j < 4; j++) {
//             		// skip 0 bits in input.
//             	    if ((b & (1 << j)) == 0)
//             	    	continue;
//
//             	    // If the bitcount is zero, this is the bit we want.
//             	    if (! pick) {
//             	    	ncbi2na[i] = j;
//             	    	break;
//             	    }
//             	    // Else, decrement.
//             	    pick--;
//             	}
// ```
// One draw per gap or ambiguous letter: two bits of it for a gap or `N`, otherwise the set bit
// of the ncbi4na value that the draw picks.
/// EXPERIMENT (LOSAT_X_TBNXLATE): the ncbi2na value drawn for a gap or an ambiguous letter
/// with the ncbi4na value `mask`.
fn x_draw_base(mask: u8, random: &mut NcbiRandom) -> u8 {
    if mask == 0 || mask == 0x0F {
        return (random.get_rand() & 0x3) as u8;
    }
    let mut pick = random.get_rand() % mask.count_ones();
    for bit in 0..4u8 {
        if mask & (1 << bit) == 0 {
            continue;
        }
        if pick == 0 {
            return bit;
        }
        pick -= 1;
    }
    0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-445
// ```c++
// void s_Ncbi4naToNcbi2na(const string & ncbi4na, int base_length,
//                         unsigned char * ncbi2na)
// {
//     int inp_bytes   = base_length;
//     CRandom random(base_length);
//
//     for(int i = 0; i < inp_bytes; i++) {
//         // one input byte
//         unsigned char inp = ncbi4na[i];
//
//         // represents 2 bases
//         unsigned char b = inp & 0xF;
//
//         // compress each to 2 bits
//         unsigned char c = ctable[b];
//
//         if (c  != 0xFF) {
//             // No ambiguities, so we can do this the easy way.
//         	ncbi2na[i] = c;
// ```
// The subject in ncbi2na, one value per byte, in one pass over its letters. The generator is
// seeded with the length of the subject, as `CRandom random(base_length)`, and is made when the
// first gap or ambiguous letter asks for a draw: its values depend on the seed and on the
// number of draws before, not on when it is made. The length check and the error of a byte
// that is no subject letter are those of `resolve_local_subject_ncbi2na`. That function checks
// every letter before its first draw; this loop may have drawn for the letters before the
// byte it stops at. The error is the same one, for the first such byte of the subject, and
// the draws made before it are dropped with the generator, which belongs to this call.
/// EXPERIMENT (LOSAT_X_TBNXLATE): `resolve_local_subject_ncbi2na` as ncbi2na values (0 A, 1 C,
/// 2 G, 3 T) instead of letters.
fn x_subject_ncbi2na(subject: &[u8]) -> Result<Vec<u8>> {
    let length = u32::try_from(subject.len())?;
    let table = x_letter_table();
    let mut random = None;
    let mut ncbi2na = vec![0u8; subject.len()];
    for (code, &letter) in ncbi2na.iter_mut().zip(subject) {
        let entry = table[letter as usize];
        if entry < 4 {
            *code = entry;
        } else if entry == X_NOT_A_LETTER {
            bail!(
                "TBLASTN subject base '{}' is not implemented",
                letter as char
            );
        } else {
            let random = random.get_or_insert_with(|| NcbiRandom::new(length));
            *code = x_draw_base(entry & 0x0F, random);
        }
    }
    Ok(ncbi2na)
}

/// EXPERIMENT (LOSAT_X_TBNXLATE): the residue (ncbistdaa) of every codon of one genetic code,
/// by the three ncbi2na values of the plus strand in sequence order.
// No NCBI counterpart for the type: the two arrays are NCBI's `translation_table` and
// `translation_table_rc` (see `x_codon_tables`).
struct XCodonTables {
    /// The codon read on the plus strand.
    forward: [u8; 64],
    /// The codon read on the minus strand at the same three bases.
    reverse: [u8; 64],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:999-1002,1022-1037
// ```c
//   	static Uint1 mapping[4] = {2, /* A in ncbi2na */
//                               1, /* C in ncbi2na. */
//                               3, /* G in ncbi2na. */
//                               0 /* T in ncbi2na. */ };
// ...
// 			   if (reverse_complement)
//             {
//                bp1 = 3 - index1;
//                bp2 = 3 - index2;
//                bp3 = 3 - index3;
// 			   	codon = (mapping[bp1]<<4) + (mapping[bp2]<<2) + (mapping[bp3]);
// 			   	translation[(index3<<4) + (index2<<2) + index1] =
//                   genetic_code[codon];
// 			   }
// 			   else
// 			   {
// 			   	codon = (mapping[index1]<<4) + (mapping[index2]<<2) +
//                   (mapping[index3]);
// 			   	translation[(index1<<4) + (index2<<2) + index3] =
//                   genetic_code[codon];
// 			   }
// ```
// `s_BlastGetTranslationTable`: the forward table is indexed by the three ncbi2na values of a
// codon, and the reverse-complement table by the three values as they stand on the plus strand
// (the codon of the minus strand is their complements, 3 - value, in the opposite order). NCBI
// reads the residue from the genetic-code string; here each entry is what the reference path
// gives for that codon as letters (`GeneticCode::get`, then `aa_char_to_ncbistdaa`, the two
// calls of `tblastx::translation::codon_to_aa_ncbistdaa`), so a table read returns the
// reference path's residue. The tables of a genetic code are made once per process.
/// EXPERIMENT (LOSAT_X_TBNXLATE): the two codon tables of the genetic code `db_gencode`.
fn x_codon_tables(db_gencode: u8) -> Result<&'static XCodonTables> {
    static TABLES: [OnceLock<XCodonTables>; 256] = [const { OnceLock::new() }; 256];
    let cell = &TABLES[db_gencode as usize];
    if let Some(tables) = cell.get() {
        return Ok(tables);
    }
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    Ok(cell.get_or_init(|| {
        const LETTERS: [u8; 4] = *b"ACGT";
        let residue = |codon: [u8; 3]| {
            aa_char_to_ncbistdaa(code.get(&[
                LETTERS[codon[0] as usize],
                LETTERS[codon[1] as usize],
                LETTERS[codon[2] as usize],
            ]))
        };
        let mut tables = XCodonTables {
            forward: [0; 64],
            reverse: [0; 64],
        };
        for index1 in 0..4u8 {
            for index2 in 0..4u8 {
                for index3 in 0..4u8 {
                    let index = x_codon_index(&[index1, index2, index3]);
                    tables.forward[index] = residue([index1, index2, index3]);
                    tables.reverse[index] = residue([3 - index3, 3 - index2, 3 - index1]);
                }
            }
        }
        tables
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1035
// ```c
// 			   	translation[(index1<<4) + (index2<<2) + index3] =
// ```
// The index of three ncbi2na values in a codon table. The values are below 4, so the mask
// changes nothing; it tells the compiler that the index is below 64.
/// EXPERIMENT (LOSAT_X_TBNXLATE): where a codon is in a table of `XCodonTables`.
#[inline(always)]
fn x_codon_index(codon: &[u8]) -> usize {
    (((codon[0] as usize) << 4) | ((codon[1] as usize) << 2) | (codon[2] as usize)) & 63
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:521-526,800-802
// ```c
//    if (nt_seq == NULL || prot_seq == NULL ||
//        (length-ABS(frame)+1) < CODON_LENGTH)
//       return prot_length;
//
//    *prot_seq = NULLB;
//    prot_seq++;
// ...
//    *prot_seq = NULLB;
//
//    return (Int4)(prot_seq - prot_seq_start);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:544-547,704-707
// ```c
//          case 0:
//             codon = (byte_value >> 2);
//             *prot_seq = translation[codon];
//             prot_seq++;
// ...
//          case 3:
//             codon = (byte_value & 63);
//             *prot_seq = translation[codon];
//             prot_seq++;
// ```
// `BLAST_TranslateCompressedSequence` once per frame: a sentinel, one table read per codon, a
// sentinel. A plus frame reads the codons from its first base on; a minus frame reads the
// reverse-complement table from the last base of the frame back. NCBI takes the six bits of a
// codon out of the packed bytes; here the three values are separate bytes. The frames, their
// order and every field are those of `tblastx::translation::generate_frames` on the resolved
// subject and its reverse complement: frame `i + 1` holds the codons of `seq[i..]`, frame
// `-(i + 1)` those of `revcomp(seq)[i..]`, which are the plus-strand triplets of
// `seq[..len - i]` counted from its end.
/// EXPERIMENT (LOSAT_X_TBNXLATE): `generate_frames` from ncbi2na values and the codon tables.
fn x_generate_frames(ncbi2na: &[u8], tables: &XCodonTables) -> Vec<QueryFrame> {
    let seq_len = ncbi2na.len();
    let mut frames = Vec::with_capacity(6);
    // Forward (Frames 1,2,3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            let codons = ncbi2na[i..].chunks_exact(3);
            frames.push(x_frame(
                (i as i8) + 1,
                seq_len,
                codons.len(),
                codons.map(|codon| tables.forward[x_codon_index(codon)]),
            ));
        }
    }
    // Reverse (Frames -1,-2,-3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            let codons = ncbi2na[..seq_len - i].rchunks_exact(3);
            frames.push(x_frame(
                -((i as i8) + 1),
                seq_len,
                codons.len(),
                codons.map(|codon| tables.reverse[x_codon_index(codon)]),
            ));
        }
    }
    frames
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:525-526,800-802
// ```c
//    *prot_seq = NULLB;
//    prot_seq++;
// ...
//    *prot_seq = NULLB;
//
//    return (Int4)(prot_seq - prot_seq_start);
// ```
// One translated frame: the sentinel, the residues, the sentinel, and the fields that
// `generate_frames` gives a frame (`translate_sequence` for the bytes and the length).
/// EXPERIMENT (LOSAT_X_TBNXLATE): a frame of `aa_len` residues.
#[inline(always)]
fn x_frame(
    frame: i8,
    seq_len: usize,
    aa_len: usize,
    residues: impl Iterator<Item = u8>,
) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(aa_len + 2);
    aa_seq.push(SENTINEL_BYTE);
    aa_seq.extend(residues);
    aa_seq.push(SENTINEL_BYTE);
    QueryFrame {
        frame,
        aa_seq,
        aa_seq_nomask: None,
        aa_len,
        orig_len: seq_len,
        seg_masks: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
// The subject plan as the reference closure of `x_scan_plan` makes it, with the frames from the
// tables. The chunks and scan ranges of the frames come from the same function as there
// (`translated_chunk_scan_ranges`), in the same order.
/// EXPERIMENT (LOSAT_X_TBNXLATE): the plan of one subject from the tables.
fn x_table_plan(subject: &[u8], db_gencode: u8, mask_lowercase: bool) -> Result<XSubjectPlan> {
    let tables = x_codon_tables(db_gencode)?;
    let ncbi2na = x_subject_ncbi2na(subject)?;
    let frames = x_generate_frames(&ncbi2na, tables);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    let mut units = Vec::new();
    for (index, frame) in frames.iter().enumerate() {
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-500:
        // SUBJECT_SPLIT_NO_RANGE skips WordFinder for this chunk.
        for (chunk, ranges) in translated_chunk_scan_ranges(
            subject,
            frame.frame,
            frame.aa_len,
            negative_first_length,
            mask_lowercase,
        )? {
            units.push((index, chunk, ranges));
        }
    }
    TABLE_PLANS.fetch_add(1, Ordering::Relaxed);
    Ok(XSubjectPlan { frames, units })
}

// No NCBI counterpart: the comparison of the shadow mode; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNXLATESHADOW): what differs between the plan from the tables and the
/// reference plan of one subject, `None` when nothing does. Two errors are the same when their
/// texts are.
fn x_plan_difference(
    tables: &Result<XSubjectPlan>,
    reference: &Result<Arc<XSubjectPlan>>,
) -> Option<String> {
    let (tables, reference) = match (tables, reference) {
        (Ok(tables), Ok(reference)) => (tables, reference),
        (Err(tables), Err(reference)) => {
            let (tables, reference) = (format!("{tables:#}"), format!("{reference:#}"));
            return (tables != reference)
                .then(|| format!("error \"{tables}\" for the reference's \"{reference}\""));
        }
        (Ok(_), Err(reference)) => {
            return Some(format!(
                "a plan for the reference's error \"{reference:#}\""
            ))
        }
        (Err(tables), Ok(_)) => {
            return Some(format!("error \"{tables:#}\" for a plan of the reference"))
        }
    };
    // Every field by name, so that a field added to a plan or to a frame has to be compared here.
    let XSubjectPlan { frames, units } = tables;
    let XSubjectPlan {
        frames: reference_frames,
        units: reference_units,
    } = &**reference;
    if frames.len() != reference_frames.len() {
        return Some(format!(
            "{} frames for the reference's {}",
            frames.len(),
            reference_frames.len()
        ));
    }
    for (ours, theirs) in frames.iter().zip(reference_frames) {
        let QueryFrame {
            frame,
            aa_seq,
            aa_seq_nomask,
            aa_len,
            orig_len,
            seg_masks,
        } = ours;
        let QueryFrame {
            frame: reference_frame,
            aa_seq: reference_aa_seq,
            aa_seq_nomask: reference_aa_seq_nomask,
            aa_len: reference_aa_len,
            orig_len: reference_orig_len,
            seg_masks: reference_seg_masks,
        } = theirs;
        if frame != reference_frame
            || aa_len != reference_aa_len
            || orig_len != reference_orig_len
            || aa_seq_nomask != reference_aa_seq_nomask
            || seg_masks != reference_seg_masks
        {
            return Some(format!(
                "frame {frame} (length {aa_len}, subject length {orig_len}) for the reference's frame {reference_frame} (length {reference_aa_len}, subject length {reference_orig_len})"
            ));
        }
        if aa_seq != reference_aa_seq {
            let at = aa_seq
                .iter()
                .zip(reference_aa_seq)
                .position(|(ours, theirs)| ours != theirs)
                .unwrap_or_else(|| aa_seq.len().min(reference_aa_seq.len()));
            return Some(format!(
                "frame {frame}: the translation differs at byte {at} (lengths {} and {})",
                aa_seq.len(),
                reference_aa_seq.len()
            ));
        }
    }
    (units != reference_units).then(|| "the units of the frames differ".to_string())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
// Dispatch point: how the translation of one subject is made. `Reference`: the closure of
// `x_scan_plan`, as before. `Tables`: `x_table_plan`. `Shadow`: both, and the search stops if a
// frame, a unit or an error differs; the plan from the tables is the one searched.
/// EXPERIMENT (LOSAT_X_TBNXLATE / LOSAT_X_TBNXLATESHADOW): the plan of one subject.
pub(super) fn x_subject_plan(
    mode: XTranslate,
    subject: &[u8],
    db_gencode: u8,
    mask_lowercase: bool,
    reference: impl Fn() -> Result<Arc<XSubjectPlan>>,
) -> Result<Arc<XSubjectPlan>> {
    match mode {
        XTranslate::Reference => reference(),
        XTranslate::Tables => x_table_plan(subject, db_gencode, mask_lowercase).map(Arc::new),
        XTranslate::Shadow => {
            let tables = x_table_plan(subject, db_gencode, mask_lowercase);
            let reference = reference();
            if let Some(difference) = x_plan_difference(&tables, &reference) {
                panic!(
                    "LOSAT_X_TBNXLATESHADOW: the subject translation from the tables differs from the reference one (subject of {} bases, genetic code {db_gencode}): {difference}",
                    subject.len()
                );
            }
            match &tables {
                Ok(plan) => {
                    SHADOW_PLANS.fetch_add(1, Ordering::Relaxed);
                    let residues: usize = plan.frames.iter().map(|frame| frame.aa_seq.len()).sum();
                    SHADOW_RESIDUES.fetch_add(residues as u64, Ordering::Relaxed);
                }
                Err(_) => {
                    SHADOW_ERRORS.fetch_add(1, Ordering::Relaxed);
                }
            }
            tables.map(Arc::new)
        }
    }
}

/// Bases whose codes `x_packed_frames` unpacks at a time: a multiple of three, so that every
/// block starts a codon of frame 1.
const X_PACKED_BLOCK: usize = 3 << 14;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1058,1063-1064
// ```c
//    Uint4 buffer_length =2*(nucl_length+1)+2;
// ...
//    if ((translation_buffer =
//         (Uint1*) malloc(buffer_length)) == NULL)
// ```
// No NCBI counterpart for the size: NCBI holds the six frames of a subject in one buffer of
// about two bytes per base, which glibc maps for a chromosome and unmaps when NCBI frees it.
// Here each frame is an allocation of a third of a byte per base. glibc maps an allocation of
// at least its mmap threshold, but it raises the threshold to the size of each mapped block
// that is freed (up to 32 MiB), so after the first chromosome the frames of the next ones
// (13 to 26 MB each on the zebrafish genome) come from a thread's arena, and the freed frames
// stay resident there (zebrafish genome, 61 proteins, 8 threads: 1,612 MiB, 1,366 MiB with a
// fixed threshold). A frame of at least `X_MAPPED_FRAME` bytes is given a capacity above that
// ceiling, so that it is always mapped and returned when dropped; the pages beyond its length
// are never touched and are not resident. Other allocators return large blocks on their own.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the frames from which a frame is mapped (Linux with glibc).
const X_MAPPED_FRAME: usize = 1 << 20;

/// EXPERIMENT (LOSAT_X_TBNSTORE): the capacity of a frame of `length` bytes (see
/// `X_MAPPED_FRAME`).
fn x_frame_capacity(length: usize) -> usize {
    if cfg!(all(target_os = "linux", target_env = "gnu")) && length >= X_MAPPED_FRAME {
        // glibc's largest dynamic mmap threshold (DEFAULT_MMAP_THRESHOLD_MAX on 64-bit), and
        // one page more
        length.max((32 << 20) + 4096)
    } else {
        length
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1090
// ```c
//    for (context = 0; context < NUM_FRAMES; ++context) {
//       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//       if (encoding == eBlastEncodingNcbi2na) {
//          if (frame > 0) {
//             length =
//                BLAST_TranslateCompressedSequence(translation_table,
//                   nucl_length, nucl_seq, frame, translation_buffer+offset);
//          } else {
//             length =
//                BLAST_TranslateCompressedSequence(translation_table_rc,
//                   nucl_length, nucl_seq, frame, translation_buffer+offset);
// ```
// NCBI translates the six frames from the packed subject, holding nothing but the frames.
// `x_generate_frames` on the ncbi2na values of the store, with the values unpacked one block
// at a time instead of for the whole subject: the codon that starts at base `p` is residue
// `p / 3` of frame `p % 3 + 1`, and, read on the minus strand, residue `q / 3` of frame
// `-(q % 3 + 1)` where `q = length - 3 - p` (the codons of `seq[..length - i]` counted from its
// end). Each frame is filled at those places, so every residue is the one `x_generate_frames`
// puts there.
/// EXPERIMENT (LOSAT_X_TBNSTORE): `x_generate_frames` of a packed subject's codes.
fn x_packed_frames(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    tables: &XCodonTables,
) -> Vec<QueryFrame> {
    let seq_len = store.len();
    let frame_count = (0..3).filter(|&i| i + 3 <= seq_len).count();
    let empty_frame = |i: usize| {
        let aa_len = (seq_len - i) / 3;
        let mut aa_seq = Vec::with_capacity(x_frame_capacity(aa_len + 2));
        aa_seq.resize(aa_len + 2, 0u8);
        aa_seq[0] = SENTINEL_BYTE;
        aa_seq[aa_len + 1] = SENTINEL_BYTE;
        aa_seq
    };
    let mut forward: Vec<Vec<u8>> = (0..frame_count).map(empty_frame).collect();
    let mut reverse: Vec<Vec<u8>> = (0..frame_count).map(empty_frame).collect();
    let mut codes = vec![0u8; X_PACKED_BLOCK + 2];
    let mut start = 0;
    // the codons that start in [start, end), with the two bases after the last one
    while start + 3 <= seq_len {
        let end = (start + X_PACKED_BLOCK).min(seq_len - 2);
        let codes = &mut codes[..end - start + 2];
        store.codes_into(start, codes);
        for phase in 0..3 {
            let first = start + phase;
            if first >= end {
                break;
            }
            let minus = seq_len - 3 - first;
            let (plus_frame, minus_frame) = (&mut forward[phase], &mut reverse[minus % 3]);
            let mut minus_residue = 1 + minus / 3;
            for (residue, codon) in plus_frame[1 + first / 3..]
                .iter_mut()
                .zip(codes[phase..].chunks_exact(3))
            {
                let index = x_codon_index(codon);
                *residue = tables.forward[index];
                minus_frame[minus_residue] = tables.reverse[index];
                minus_residue -= 1;
            }
        }
        start = end;
    }
    let frame = |frame: i8, aa_seq: Vec<u8>| QueryFrame {
        frame,
        aa_len: aa_seq.len() - 2,
        aa_seq,
        aa_seq_nomask: None,
        orig_len: seq_len,
        seg_masks: Vec::new(),
    };
    let mut frames = Vec::with_capacity(2 * frame_count);
    frames.extend((1..).zip(forward).map(|(i, aa_seq)| frame(i, aa_seq)));
    frames.extend((1..).zip(reverse).map(|(i, aa_seq)| frame(-i, aa_seq)));
    frames
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
// The plan of a packed subject as `x_table_plan` makes it from the subject's letters. Only a
// stage that does not read the case of the letters takes the store (`-lcase_masking` gives the
// letters), so every frame has one unit per chunk with one range, as without lowercase masking.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the plan of one packed subject, from its codes.
fn x_packed_plan(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
) -> Result<XSubjectPlan> {
    let tables = x_codon_tables(db_gencode)?;
    let frames = x_packed_frames(store, tables);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    let mut units = Vec::new();
    for (index, frame) in frames.iter().enumerate() {
        // Without lowercase masking the ranges read no letter of the subject.
        for (chunk, ranges) in translated_chunk_scan_ranges(
            &[],
            frame.frame,
            frame.aa_len,
            negative_first_length,
            false,
        )? {
            units.push((index, chunk, ranges));
        }
    }
    super::x_subject_store::x_count_packed_plan();
    Ok(XSubjectPlan { frames, units })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
// Dispatch point: a packed subject is translated from its codes (`x_packed_plan`) whichever
// translation the variables select, since the reference path reads letters. With
// `LOSAT_X_TBNXLATESHADOW` the letters the reference path reads (the resolved bases) are made
// too and the search stops if the two plans differ.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the plan of one packed subject; `reference` makes the plan of
/// the letters it is given.
pub(super) fn x_packed_subject_plan(
    mode: XTranslate,
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
    reference: impl Fn(&[u8]) -> Result<Arc<XSubjectPlan>>,
) -> Result<Arc<XSubjectPlan>> {
    let packed = x_packed_plan(store, db_gencode);
    if mode == XTranslate::Shadow {
        let mut letters = Vec::new();
        store.resolved_letters_into(&mut letters);
        if let Some(difference) = x_plan_difference(&packed, &reference(&letters)) {
            panic!(
                "LOSAT_X_TBNXLATESHADOW: the subject translation from the packed subject differs from the reference one (subject of {} bases, genetic code {db_gencode}): {difference}",
                store.len()
            );
        }
        match &packed {
            Ok(plan) => {
                SHADOW_PLANS.fetch_add(1, Ordering::Relaxed);
                let residues: usize = plan.frames.iter().map(|frame| frame.aa_seq.len()).sum();
                SHADOW_RESIDUES.fetch_add(residues as u64, Ordering::Relaxed);
            }
            Err(_) => {
                SHADOW_ERRORS.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
    packed.map(Arc::new)
}

#[cfg(test)]
mod tests {
    use super::super::search_seed::{resolve_local_subject_ncbi2na, x_scan_plan};
    use super::*;
    use crate::algorithm::tblastx::translation::generate_frames;
    use crate::config::ScoringMatrix;

    // No NCBI counterpart: the three helpers below prepare the inputs of the tests; they compute
    // no value NCBI computes.
    /// The genetic codes of NCBI's `gc.prt` that LOSAT knows.
    fn genetic_codes() -> Vec<u8> {
        let codes: Vec<u8> = (0..=u8::MAX)
            .filter(|&id| GeneticCode::try_from_id(id).is_ok())
            .collect();
        assert_eq!(codes.len(), 27);
        codes
    }

    /// The reference plan of one subject: a copy of the reference closure of `x_scan_plan`.
    /// `x_table_plans_are_the_reference_plans` compares it with what `x_scan_plan` itself
    /// gives for every subject it takes, so the copy cannot drift from that closure.
    fn reference_plan(
        subject: &[u8],
        db_gencode: u8,
        mask_lowercase: bool,
    ) -> Result<Arc<XSubjectPlan>> {
        let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
        let resolved_subject = resolve_local_subject_ncbi2na(subject)?;
        let frames = generate_frames(&resolved_subject, &code);
        let negative_first_length = frames
            .iter()
            .find(|frame| frame.frame == -1)
            .map(|frame| frame.aa_len)
            .unwrap_or(0);
        let mut units = Vec::new();
        for (index, frame) in frames.iter().enumerate() {
            for (chunk, ranges) in translated_chunk_scan_ranges(
                subject,
                frame.frame,
                frame.aa_len,
                negative_first_length,
                mask_lowercase,
            )? {
                units.push((index, chunk, ranges));
            }
        }
        Ok(Arc::new(XSubjectPlan { frames, units }))
    }

    /// `length` letters of `alphabet`, drawn by a fixed linear congruential sequence.
    fn letters(alphabet: &[u8], length: usize, seed: u64) -> Vec<u8> {
        let mut state = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                alphabet[(state >> 33) as usize % alphabet.len()]
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:441-450
    // ```c++
    //         unsigned char c = ctable[b];
    //
    //         if (c  != 0xFF) {
    //             // No ambiguities, so we can do this the easy way.
    //         	ncbi2na[i] = c;
    //
    //         } else {
    //             if (b == 0 || b == 0x0F) {
    //             	//gap or N
    //                 ncbi2na[i] = random.GetRand() & 0x3;
    // ```
    // Every byte as a subject of one letter: the ncbi2na value is the letter the reference path
    // resolves it to, and a byte that is no subject letter (`U` among them) gives the reference
    // path's error text.
    #[test]
    fn x_every_byte_resolves_as_in_the_reference_path() {
        let mut letters = 0;
        for byte in 0..=u8::MAX {
            let subject = [byte];
            match (
                x_subject_ncbi2na(&subject),
                resolve_local_subject_ncbi2na(&subject),
            ) {
                (Ok(ncbi2na), Ok(resolved)) => {
                    assert_eq!(ncbi2na.len(), 1);
                    assert_eq!(b"ACGT"[ncbi2na[0] as usize], resolved[0], "byte {byte}");
                    letters += 1;
                }
                (Err(ours), Err(theirs)) => {
                    assert_eq!(ours.to_string(), theirs.to_string(), "byte {byte}");
                }
                (ours, theirs) => panic!("byte {byte}: {ours:?} for the reference's {theirs:?}"),
            }
        }
        // Fifteen IUPAC letters in either case and the gap.
        assert_eq!(letters, 31);
        assert!(x_subject_ncbi2na(b"ACGU").is_err());
        assert!(x_subject_ncbi2na(b"acgu").is_err());
        // The error names the first byte of the subject that is no letter, also after letters
        // that draw and before other such bytes.
        for (subject, named) in [
            (&b"A!A?"[..], '!'),
            (b"?A!", '?'),
            (b"NNRY!ACGT?", '!'),
            (b"ACGTNu*", 'u'),
            (b"-N-\xE9!", '\u{e9}'),
        ] {
            let expected = format!("TBLASTN subject base '{named}' is not implemented");
            let ours = x_subject_ncbi2na(subject).err().expect("no letter");
            let theirs = resolve_local_subject_ncbi2na(subject)
                .err()
                .expect("no letter");
            assert_eq!(ours.to_string(), expected);
            assert_eq!(theirs.to_string(), expected);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:999-1002,1022-1037
    // ```c
    //   	static Uint1 mapping[4] = {2, /* A in ncbi2na */
    //                               1, /* C in ncbi2na. */
    //                               3, /* G in ncbi2na. */
    //                               0 /* T in ncbi2na. */ };
    // ...
    // 			   if (reverse_complement)
    //             {
    //                bp1 = 3 - index1;
    //                bp2 = 3 - index2;
    //                bp3 = 3 - index3;
    // 			   	codon = (mapping[bp1]<<4) + (mapping[bp2]<<2) + (mapping[bp3]);
    // 			   	translation[(index3<<4) + (index2<<2) + index1] =
    //                   genetic_code[codon];
    // 			   }
    // 			   else
    // 			   {
    // 			   	codon = (mapping[index1]<<4) + (mapping[index2]<<2) +
    //                   (mapping[index3]);
    // 			   	translation[(index1<<4) + (index2<<2) + index3] =
    //                   genetic_code[codon];
    // 			   }
    // ```
    // The two tables of every genetic code against this loop over the genetic-code string
    // (`GeneticCode::table`, in ncbieaa), and against frames 1 and -1 of the reference path's
    // `generate_frames` on the codon as letters.
    #[test]
    fn x_codon_tables_are_ncbis_translation_tables_for_every_genetic_code() {
        let mapping = [2usize, 1, 3, 0];
        for db_gencode in genetic_codes() {
            let code = GeneticCode::try_from_id(db_gencode).unwrap();
            let tables = x_codon_tables(db_gencode).unwrap();
            let mut forward = [0u8; 64];
            let mut reverse = [0u8; 64];
            for index1 in 0..4usize {
                for index2 in 0..4usize {
                    for index3 in 0..4usize {
                        let (bp1, bp2, bp3) = (3 - index1, 3 - index2, 3 - index3);
                        let codon = (mapping[bp1] << 4) + (mapping[bp2] << 2) + mapping[bp3];
                        reverse[(index3 << 4) + (index2 << 2) + index1] =
                            aa_char_to_ncbistdaa(code.table[codon]);
                        let codon =
                            (mapping[index1] << 4) + (mapping[index2] << 2) + mapping[index3];
                        forward[(index1 << 4) + (index2 << 2) + index3] =
                            aa_char_to_ncbistdaa(code.table[codon]);
                    }
                }
            }
            assert_eq!(tables.forward, forward, "genetic code {db_gencode}");
            assert_eq!(tables.reverse, reverse, "genetic code {db_gencode}");
            for index in 0..64usize {
                let codon = [
                    b"ACGT"[index >> 4],
                    b"ACGT"[(index >> 2) & 3],
                    b"ACGT"[index & 3],
                ];
                let frames = generate_frames(&codon, &code);
                assert_eq!(frames.len(), 2);
                assert_eq!(
                    (frames[0].frame, frames[0].aa_seq[1]),
                    (1, tables.forward[index]),
                    "genetic code {db_gencode}, codon {index}"
                );
                assert_eq!(
                    (frames[1].frame, frames[1].aa_seq[1]),
                    (-1, tables.reverse[index]),
                    "genetic code {db_gencode}, codon {index}"
                );
            }
        }
        assert_eq!(
            x_codon_tables(0).err().map(|error| error.to_string()),
            GeneticCode::try_from_id(0).err()
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1091,1100-1102
    // ```c
    //    for (context = 0; context < NUM_FRAMES; ++context) {
    //       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //       if (encoding == eBlastEncodingNcbi2na) {
    //          if (frame > 0) {
    //             length =
    //                BLAST_TranslateCompressedSequence(translation_table,
    //                   nucl_length, nucl_seq, frame, translation_buffer+offset);
    //          } else {
    //             length =
    //                BLAST_TranslateCompressedSequence(translation_table_rc,
    //                   nucl_length, nucl_seq, frame, translation_buffer+offset);
    //          }
    // ...
    //       offset += length + 1;
    //       frame_offsets[context+1] = offset;
    //    }
    // ```
    // The plan from the tables against the reference plan: every length from 0 to 40 bases (the
    // subjects too short for a frame, and every remainder of the three frames on both strands)
    // and longer ones, of the four bases, of the bases with lowercase runs, of every IUPAC
    // letter with gaps (one draw each, in subject order), and of `N` only; with and without
    // lowercase masking, for every genetic code on the long subjects.
    #[test]
    fn x_table_plans_are_the_reference_plans() {
        let alphabets: [&[u8]; 4] = [
            b"ACGT",
            b"ACGTACGTACGTacgtn",
            b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn",
            b"N",
        ];
        let mut compared = 0usize;
        let mut residues = 0usize;
        let mut scan_plans = 0usize;
        // One lookup table for every call of `x_scan_plan` (as `LOSAT_X_TBNQSIDE` shares it).
        let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
        let lookup = OnceLock::new();
        for (which, alphabet) in alphabets.iter().enumerate() {
            let lengths = (0..=40).chain([99, 100, 101, 1_000, 4_999, 20_003]);
            for length in lengths {
                let subject = letters(alphabet, length, (which * 100_000 + length) as u64);
                let codes = if length >= 1_000 {
                    genetic_codes()
                } else {
                    vec![1, 4, 11]
                };
                for db_gencode in codes {
                    for mask_lowercase in [false, true] {
                        let tables = x_table_plan(&subject, db_gencode, mask_lowercase);
                        let reference = reference_plan(&subject, db_gencode, mask_lowercase);
                        assert_eq!(
                            x_plan_difference(&tables, &reference),
                            None,
                            "alphabet {which}, {length} bases, genetic code {db_gencode}, lowercase masking {mask_lowercase}"
                        );
                        let plan = tables.unwrap();
                        assert_eq!(plan.frames.len(), 2 * length.saturating_sub(2).min(3));
                        residues += plan.frames.iter().map(|f| f.aa_len).sum::<usize>();
                        compared += 1;
                        // The plan `x_scan_plan` makes (by its reference closure unless the
                        // test run has one of the variables set) is the copy's too.
                        let scan = x_scan_plan(
                            &[query],
                            crate::sequence::packed_subject::SubjectNt::Letters(&subject),
                            db_gencode,
                            None,
                            11,
                            mask_lowercase,
                            ScoringMatrix::Blosum62,
                            3,
                            Some(&lookup),
                            None,
                        )
                        .unwrap();
                        assert_eq!(scan.is_some(), length >= 11);
                        if let Some(scan) = scan {
                            let copy = reference
                                .map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                            assert_eq!(
                                x_plan_difference(&copy, &Ok(scan.subject)),
                                None,
                                "x_scan_plan: alphabet {which}, {length} bases, genetic code {db_gencode}, lowercase masking {mask_lowercase}"
                            );
                            scan_plans += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(compared, 4 * (44 * 3 + 3 * 27) * 2);
        assert_eq!(scan_plans, 4 * (33 * 3 + 3 * 27) * 2);
        assert!(residues > 1_000_000);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:431,461
    // ```c++
    //     CRandom random(base_length);
    // ...
    //             	int pick = random.GetRand() % bitcount;
    // ```
    // The generator is seeded with the length of the subject and drawn from once per gap or
    // ambiguous letter, in subject order. A subject whose only such letter is its last one, and
    // subjects that differ in length only, draw as the reference path does.
    #[test]
    fn x_draws_follow_the_seed_and_the_order_of_the_reference_path() {
        for length in [1usize, 2, 3, 33, 34, 257, 1_000] {
            for last in *b"NRYKMSWBDHV-" {
                let mut subject = vec![b'A'; length];
                subject[length - 1] = last;
                let ncbi2na = x_subject_ncbi2na(&subject).unwrap();
                let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
                let ours: Vec<u8> = ncbi2na.iter().map(|&code| b"ACGT"[code as usize]).collect();
                assert_eq!(ours, resolved, "{length} bases ending in {}", last as char);
            }
        }
        let subject = letters(b"NRYKMSWBDHV-", 50_000, 7);
        let ncbi2na = x_subject_ncbi2na(&subject).unwrap();
        let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
        let ours: Vec<u8> = ncbi2na.iter().map(|&code| b"ACGT"[code as usize]).collect();
        assert_eq!(ours, resolved);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
    // ```c
    //             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
    //                                      backup.full_range.right,
    //                                      subject->gen_code_string, &translation_buffer,
    //                                      &frame_offsets, NULL);
    // ```
    // The three ways of `x_subject_plan` give one plan, and `x_scan_plan` (whichever way the
    // environment of the test run selects) gives it too. A subject with a byte that is no
    // letter gives the reference path's error in all three.
    #[test]
    fn x_subject_plans_agree_in_every_mode_and_with_the_scan_plan() {
        let subject = letters(b"ACGTACGTACGTacgtnRY", 3_001, 11);
        for (db_gencode, mask_lowercase) in [(1, false), (1, true), (4, true), (11, false)] {
            let reference = || reference_plan(&subject, db_gencode, mask_lowercase);
            let expected = reference();
            for mode in [
                XTranslate::Reference,
                XTranslate::Tables,
                XTranslate::Shadow,
            ] {
                let plan = x_subject_plan(mode, &subject, db_gencode, mask_lowercase, reference)
                    .map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                assert_eq!(x_plan_difference(&plan, &expected), None, "{mode:?}");
            }
            let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
            let scan = x_scan_plan(
                &[query],
                crate::sequence::packed_subject::SubjectNt::Letters(&subject),
                db_gencode,
                None,
                11,
                mask_lowercase,
                ScoringMatrix::Blosum62,
                3,
                None,
                None,
            )
            .unwrap()
            .expect("a subject long enough to search");
            let tables = x_table_plan(&subject, db_gencode, mask_lowercase);
            assert_eq!(x_plan_difference(&tables, &Ok(scan.subject)), None);
        }
        let mut bad = subject.clone();
        bad[1_500] = b'u';
        let reference = || reference_plan(&bad, 1, false);
        let expected = reference().err().expect("no letter").to_string();
        assert_eq!(expected, "TBLASTN subject base 'u' is not implemented");
        for mode in [
            XTranslate::Reference,
            XTranslate::Tables,
            XTranslate::Shadow,
        ] {
            let error = x_subject_plan(mode, &bad, 1, false, reference)
                .err()
                .expect("no letter");
            assert_eq!(error.to_string(), expected, "{mode:?}");
        }
    }

    // No NCBI counterpart: the comparison of the shadow mode. It must see a changed residue, a
    // changed field of a frame, a missing frame, changed units and a different error.
    #[test]
    fn x_plan_difference_sees_every_part_of_a_plan() {
        let subject = letters(b"ACGTacgt", 300, 3);
        let reference = reference_plan(&subject, 1, true);
        let table_plan = || x_table_plan(&subject, 1, true).unwrap();
        assert_eq!(x_plan_difference(&Ok(table_plan()), &reference), None);

        let mut plan = table_plan();
        plan.frames[4].aa_seq[17] ^= 1;
        let difference = x_plan_difference(&Ok(plan), &reference).expect("a residue");
        assert!(difference.contains("frame -2") && difference.contains("byte 17"));

        let mut plan = table_plan();
        plan.frames[2].aa_len += 1;
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.frames[0].orig_len += 1;
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.frames[5].frame = 3;
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.frames[1].aa_seq_nomask = Some(Vec::new());
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.frames[1].seg_masks.push((0, 1));
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.frames.pop();
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        assert!(
            plan.units[3].2.len() > 1,
            "lowercase runs split the scan ranges"
        );
        plan.units[3].2[0].1 += 1;
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());
        let mut plan = table_plan();
        plan.units.pop();
        assert!(x_plan_difference(&Ok(plan), &reference).is_some());

        let failed = || -> Result<XSubjectPlan> { bail!("one") };
        assert!(x_plan_difference(&failed(), &reference).is_some());
        assert!(x_plan_difference(&Ok(table_plan()), &Err(anyhow::anyhow!("one"))).is_some());
        assert!(x_plan_difference(&failed(), &Err(anyhow::anyhow!("two"))).is_some());
        assert_eq!(
            x_plan_difference(&failed(), &Err(anyhow::anyhow!("one"))),
            None
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1090
    // ```c
    //    for (context = 0; context < NUM_FRAMES; ++context) {
    //       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //       if (encoding == eBlastEncodingNcbi2na) {
    //          if (frame > 0) {
    //             length =
    //                BLAST_TranslateCompressedSequence(translation_table,
    //                   nucl_length, nucl_seq, frame, translation_buffer+offset);
    // ```
    // The plan of a packed subject against the reference plan of its letters without lowercase
    // masking: every length from 0 to 40 bases, lengths around the blocks of `x_packed_frames`
    // (every remainder of a block and of a frame), of the four bases, with lowercase runs, of
    // every IUPAC letter with gaps, and of `N` only; for every genetic code on the long subjects.
    // In every mode of the variables, and as `x_scan_plan` makes it from the store.
    #[test]
    fn x_packed_plans_are_the_reference_plans() {
        use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};
        let alphabets: [&[u8]; 4] = [
            b"ACGT",
            b"ACGTACGTACGTacgtn",
            b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn",
            b"N",
        ];
        let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
        let lookup = OnceLock::new();
        let mut compared = 0usize;
        for (which, alphabet) in alphabets.iter().enumerate() {
            let around_blocks = [1, 2, 3]
                .into_iter()
                .flat_map(|blocks| (0..7).map(move |extra| blocks * X_PACKED_BLOCK + extra - 3));
            let lengths = (0..=40)
                .chain([99, 100, 101, 1_000, 4_999])
                .chain(around_blocks);
            for length in lengths {
                let subject = letters(alphabet, length, (which * 100_000 + length) as u64);
                let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
                let codes = if (1_000..10_000).contains(&length) {
                    genetic_codes()
                } else {
                    vec![1, 4, 11]
                };
                for db_gencode in codes {
                    let reference = reference_plan(&subject, db_gencode, false);
                    let packed = x_packed_plan(&store, db_gencode);
                    assert_eq!(
                        x_plan_difference(&packed, &reference),
                        None,
                        "alphabet {which}, {length} bases, genetic code {db_gencode}"
                    );
                    compared += 1;
                    if db_gencode != 11 {
                        continue;
                    }
                    for mode in [
                        XTranslate::Reference,
                        XTranslate::Tables,
                        XTranslate::Shadow,
                    ] {
                        let plan = x_packed_subject_plan(mode, &store, db_gencode, |letters| {
                            reference_plan(letters, db_gencode, false)
                        })
                        .map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                        assert_eq!(x_plan_difference(&plan, &reference), None, "{mode:?}");
                    }
                    let scan = x_scan_plan(
                        &[query],
                        SubjectNt::Packed(&store),
                        db_gencode,
                        None,
                        11,
                        false,
                        ScoringMatrix::Blosum62,
                        3,
                        Some(&lookup),
                        None,
                    )
                    .unwrap();
                    assert_eq!(scan.is_some(), length >= 11);
                    if let Some(scan) = scan {
                        let reference =
                            reference.map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                        assert_eq!(x_plan_difference(&reference, &Ok(scan.subject)), None);
                    }
                }
            }
        }
        assert_eq!(compared, 4 * (44 * 3 + 2 * 27 + 21 * 3));
        // a packed subject is never searched with lowercase masking (the letters are given)
        let store = PackedNucleotide::from_subject_letters(&letters(b"ACGT", 300, 1)).unwrap();
        let masked = x_scan_plan(
            &[query],
            SubjectNt::Packed(&store),
            1,
            None,
            11,
            true,
            ScoringMatrix::Blosum62,
            3,
            None,
            None,
        );
        assert!(masked.is_err());
    }

    // No NCBI counterpart: the shadow mode stops the search when the plan of a packed subject
    // differs from the reference plan of its letters.
    #[test]
    #[should_panic(
        expected = "LOSAT_X_TBNXLATESHADOW: the subject translation from the packed subject differs"
    )]
    fn x_shadow_mode_stops_on_a_difference_of_a_packed_plan() {
        let subject = letters(b"ACGT", 300, 5);
        let other = letters(b"ACGT", 300, 6);
        let store =
            crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&subject)
                .unwrap();
        let _ = x_packed_subject_plan(XTranslate::Shadow, &store, 1, |_| {
            reference_plan(&other, 1, false)
        });
    }

    // No NCBI counterpart: the shadow mode stops the search when the two plans differ.
    #[test]
    #[should_panic(
        expected = "LOSAT_X_TBNXLATESHADOW: the subject translation from the tables differs"
    )]
    fn x_shadow_mode_stops_on_a_difference() {
        let subject = letters(b"ACGT", 300, 5);
        let other = letters(b"ACGT", 300, 6);
        let _ = x_subject_plan(XTranslate::Shadow, &subject, 1, false, || {
            reference_plan(&other, 1, false)
        });
    }
}

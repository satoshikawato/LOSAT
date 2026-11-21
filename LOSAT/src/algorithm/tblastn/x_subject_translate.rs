//! The six translated frames of a TBLASTN subject read from tables, as NCBI makes them from the
//! compressed subject (graduated 2026-10-22: `LOSAT_X_TBNXLATE` and `LOSAT_X_TBNXLATESHADOW` are
//! gone, this is the only translation).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1071-1074,1080-1090
//! The preliminary search reads the subject in ncbi2na, so `BLAST_GetAllTranslations` takes this
//! branch: one table of 64 residues per strand and one table read per codon. The subject is turned
//! into ncbi2na in one pass and the frames are read from the two tables. The tables are filled
//! from `GeneticCode::get` and `aa_char_to_ncbistdaa`, so every residue is the one the letter-wise
//! translation of `tblastx::translation::generate_frames` gives; the ambiguous letters draw from
//! the same generator with the same seed in the same order (`resolve_local_subject_ncbi2na`, kept in
//! this module's tests as the reference).
//!
//! A packed subject (`LOSAT_X_TBNSTORE`, every `-db` sequence) is translated from its packed
//! bytes (`x_packed_frames`), 32 codons of a phase at a time with AVX2 (`LOSAT_X_TBNXLATESIMD`,
//! in the strict set: on unless it is `0`), 16 at a time with 128 bits where AVX2 does not run
//! (`LOSAT_X_TBNXLATE128`: wasm32 simd128, and x86-64 with SSSE3 and SSE4.1; written in the form
//! V8 compiles to the fewest instructions, `XW128`), elsewhere eight at a time from one word.
//! With `LOSAT_X_STATS` set, `main` prints `[X_STATS] TBNXLATESIMD codons_vector= codons_scalar=`
//! at exit, when a packed subject was translated.
//!
//! A long packed subject whose plan is made for its use alone (the search keeps no plan of it:
//! fewer than four query batches, or a plan over the limit of kept plans, `XSubjectSide::keeps`)
//! gets a plan without residues (`x_packed_lean_plan`;
//! `LOSAT_X_TBNUNITXLATE`, in the strict set: on unless it is `0`): its frames have their lengths
//! only, and each (frame, chunk) unit translates the residues it reads into a buffer the threads
//! of a query batch share (`x_unit_residues`; `search_gapped::XUnitResiduesPool`, freed at the
//! end of the batch's preliminary stage), one strand of one phase from the packed bytes, by the
//! same tables and loops (32 codons with AVX2, else eight from one word; counted in `[X_STATS]
//! TBNXLATESIMD`). With `LOSAT_X_STATS` set, `main` prints `[X_STATS] TBNUNITXLATE plans= units=
//! residues= buffers= freed=` when a plan without residues was made.

use super::search_seed::XSubjectPlan;
use crate::algorithm::common::translated_chunks::{
    translated_chunk_scan_ranges, translated_db_mask_chunks, TranslatedChunk,
};
use crate::algorithm::common::x_db_stream::x_frames_pooled;
use crate::algorithm::tblastx::constants::SENTINEL_BYTE;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::core::blast_encoding::{iupacna_to_ncbi4na, NcbiRandom};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbistdaa;
use anyhow::{bail, Result};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

// No NCBI counterpart: counter for LOSAT_X_STATS that tells how the packed subjects were
// translated; it does not change any value NCBI computes.
/// The codons of packed subjects read by the vector loop (`[0]`) and by the word loop or one at
/// a time (`[1]`).
static PACKED_CODONS: [AtomicU64; 2] = [AtomicU64::new(0), AtomicU64::new(0)];

// No NCBI counterpart: counters for LOSAT_X_STATS; they do not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the plans without residues, the units that translated their
/// residues, those residues, the units' buffers made and those freed (with LOSAT_X_STATS).
static UNIT_RESIDUES: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];

// No NCBI counterpart: reads LOSAT_X_STATS once; it does not change any value NCBI computes.
fn x_stats() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_STATS").is_some())
}

/// The summary lines (printed by `main` at exit) with LOSAT_X_STATS: the codons of packed subjects
/// read by the vector loop and the units that translated their residues.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let [vector, scalar] = PACKED_CODONS
            .each_ref()
            .map(|codons| codons.load(Ordering::Relaxed));
        if vector + scalar > 0 {
            eprintln!("[X_STATS] TBNXLATESIMD codons_vector={vector} codons_scalar={scalar}");
        }
        let [plans, units, residues, buffers, freed] = UNIT_RESIDUES
            .each_ref()
            .map(|count| count.load(Ordering::Relaxed));
        if plans > 0 {
            eprintln!(
                "[X_STATS] TBNUNITXLATE plans={plans} units={units} residues={residues} buffers={buffers} freed={freed}"
            );
        }
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
// NCBI reads the ncbi2na value of an unambiguous ncbi4na value from `ctable`. This table joins
// that read to the one that precedes it in LOSAT, letter to ncbi4na: per byte, the ncbi2na
// value of `A`, `C`, `G`, `T` (either case), the flagged ncbi4na value of a gap or an ambiguous
// letter, or `X_NOT_A_LETTER`. It is filled from the reference path's own
// `iupacna_to_ncbi4na`, with the reference path's rule for `U`.
/// What a subject byte is, in one read.
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
// One draw per gap or ambiguous letter: two bits of it for a gap or `N`, otherwise the set bit
// of the ncbi4na value that the draw picks.
/// The ncbi2na value drawn for a gap or an ambiguous letter
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
// The subject in ncbi2na, one value per byte, in one pass over its letters. The generator is
// seeded with the length of the subject, as `CRandom random(base_length)`, and is made when the
// first gap or ambiguous letter asks for a draw: its values depend on the seed and on the
// number of draws before, not on when it is made. The length check and the error of a byte
// that is no subject letter are those of `resolve_local_subject_ncbi2na`. That function checks
// every letter before its first draw; this loop may have drawn for the letters before the
// byte it stops at. The error is the same one, for the first such byte of the subject, and
// the draws made before it are dropped with the generator, which belongs to this call.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:428-473 resolves the ambiguities during
// the ncbi2na encoding; c++/src/algo/blast/api/blast_setup_cxx.cpp:742-745,830-853 keeps the
// ncbi4na sequence for the later reevaluation.
/// `resolve_local_subject_ncbi2na` as ncbi2na values (0 A, 1 C,
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

/// The residue (ncbistdaa) of every codon of one genetic code,
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
// `s_BlastGetTranslationTable`: the forward table is indexed by the three ncbi2na values of a
// codon, and the reverse-complement table by the three values as they stand on the plus strand
// (the codon of the minus strand is their complements, 3 - value, in the opposite order). NCBI
// reads the residue from the genetic-code string; here each entry is what the reference path
// gives for that codon as letters (`GeneticCode::get`, then `aa_char_to_ncbistdaa`, the two
// calls of `tblastx::translation::codon_to_aa_ncbistdaa`), so a table read returns the
// reference path's residue. The tables of a genetic code are made once per process.
/// The two codon tables of the genetic code `db_gencode`.
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
// The index of three ncbi2na values in a codon table. The values are below 4, so the mask
// changes nothing; it tells the compiler that the index is below 64.
/// Where a codon is in a table of `XCodonTables`.
#[inline(always)]
fn x_codon_index(codon: &[u8]) -> usize {
    (((codon[0] as usize) << 4) | ((codon[1] as usize) << 2) | (codon[2] as usize)) & 63
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:521-526,800-802
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:544-547,704-707
// `BLAST_TranslateCompressedSequence` once per frame: a sentinel, one table read per codon, a
// sentinel. A plus frame reads the codons from its first base on; a minus frame reads the
// reverse-complement table from the last base of the frame back. NCBI takes the six bits of a
// codon out of the packed bytes; here the three values are separate bytes. The frames, their
// order and every field are those of `tblastx::translation::generate_frames` on the resolved
// subject and its reverse complement: frame `i + 1` holds the codons of `seq[i..]`, frame
// `-(i + 1)` those of `revcomp(seq)[i..]`, which are the plus-strand triplets of
// `seq[..len - i]` counted from its end.
/// The frames `generate_frames` gives, made from ncbi2na values and the codon tables.
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
// One translated frame: the sentinel, the residues, the sentinel, and the fields that
// `generate_frames` gives a frame (`translate_sequence` for the bytes and the length).
/// A frame of `aa_len` residues.
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
// The subject plan as `x_scan_plan` makes it, with the frames from the tables. The chunks and
// scan ranges of the frames come from the same function as there
// (`translated_chunk_scan_ranges`), in the same order.
/// The plan of one subject from the tables.
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-500:
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
    Ok(XSubjectPlan { frames, units })
}

// No NCBI counterpart: the comparison of two plans of one subject for the tests; it does not
// change any value NCBI computes.
/// What differs between two plans of one subject, `None` when nothing does.
#[cfg(test)]
pub(super) fn x_plan_fields_difference(
    tables: &XSubjectPlan,
    reference: &XSubjectPlan,
) -> Option<String> {
    // Every field by name, so that a field added to a plan or to a frame has to be compared here.
    let XSubjectPlan { frames, units } = tables;
    let XSubjectPlan {
        frames: reference_frames,
        units: reference_units,
    } = reference;
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
// The translation of one subject is made from the tables (`x_table_plan`).
/// The plan of one subject.
pub(super) fn x_subject_plan(
    subject: &[u8],
    db_gencode: u8,
    mask_lowercase: bool,
) -> Result<Arc<XSubjectPlan>> {
    x_table_plan(subject, db_gencode, mask_lowercase).map(Arc::new)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1058,1063-1064
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

// No NCBI counterpart for the reuse (graduated 2026-10-27; the switch was
// `LOSAT_X_TBNFRAMEREUSE`): a thread keeps the frames of a dropped subject plan (`x_reuse_frame`)
// for the next plan it translates from a packed subject (`x_frame_buffer`).
// No NCBI counterpart for the reuse: NCBI allocates a subject's translation buffer
// (`BLAST_GetAllTranslations`) and frees it after the subject (blast_engine.c:627-631); here a
// thread gives the frames of a dropped plan to its next plan. The residues written are the same.
/// The largest frame capacity a thread keeps for its next
/// plan, and how many frames it keeps (those of one plan). A thread holds at most 1.5 MiB so.
const X_REUSED_FRAME: usize = 256 << 10;
const X_REUSED_FRAMES: usize = 6;

thread_local! {
    /// The frames this thread keeps between subject plans.
    static X_FRAMES: std::cell::RefCell<Vec<Vec<u8>>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Frees the frames this thread keeps for its next plan, at the end of its loop over the chunks
/// of the subject flow (`x_oid_flow::x_loop_free`, LOSAT_X_TBNLOOPFREE).
#[cfg(feature = "parallel")]
pub(super) fn x_free_kept_frames() {
    let _ = X_FRAMES.try_with(|frames| {
        if let Ok(mut frames) = frames.try_borrow_mut() {
            *frames = Vec::new();
        }
    });
}

/// The number of frames this thread keeps for its next plan.
#[cfg(all(test, feature = "parallel"))]
pub(super) fn x_kept_frames() -> usize {
    X_FRAMES.with(|frames| frames.borrow().len())
}

/// A frame of `length` bytes with a sentinel at both ends,
/// whose other bytes the caller writes, every one of them. A frame of at most
/// `X_REUSED_FRAME` bytes is one this thread kept (`x_reuse_frame`), its bytes not cleared, which
/// saves the allocation and the zero fill of a new one; otherwise a new zeroed frame of
/// `x_frame_capacity`. While a `-db` stream is open, a frame of
/// `X_MAPPED_FRAME` bytes or more is an idle buffer of `utils::x_mapped::X_BYTES` that touched at
/// most twice its bytes, if it holds one (its pages resident), zeroed to `length`.
fn x_frame_buffer(length: usize) -> Vec<u8> {
    let kept = if length <= X_REUSED_FRAME {
        X_FRAMES.with(|frames| frames.borrow_mut().pop())
    } else if length >= X_MAPPED_FRAME && x_frames_pooled() {
        crate::utils::x_mapped::X_BYTES.x_take_near(length)
    } else {
        None
    };
    let mut frame = match kept {
        Some(mut frame) => {
            if frame.len() >= length {
                frame.truncate(length);
            } else {
                // Grown to the length only, so that it stays a frame the thread keeps.
                frame.reserve_exact(length - frame.len());
                frame.resize(length, 0);
            }
            frame
        }
        None => {
            let mut frame = Vec::with_capacity(x_frame_capacity(length));
            frame.resize(length, 0);
            frame
        }
    };
    frame[0] = SENTINEL_BYTE;
    frame[length - 1] = SENTINEL_BYTE;
    frame
}

/// Keeps `frame`, a frame of a dropped subject plan, for
/// this thread's next `x_frame_buffer`, if its capacity is at most
/// `X_REUSED_FRAME` and the thread keeps fewer than `X_REUSED_FRAMES`; drops it otherwise.
/// While a `-db` stream is open, a frame of a capacity of
/// `X_MAPPED_FRAME` bytes or more goes to `utils::x_mapped::X_BYTES` instead (which frees the idle
/// buffers beyond its bound).
pub(super) fn x_reuse_frame(mut frame: Vec<u8>) {
    if frame.capacity() == 0 || frame.capacity() > X_REUSED_FRAME {
        if frame.capacity() >= X_MAPPED_FRAME && x_frames_pooled() {
            crate::utils::x_mapped::X_BYTES.give(&mut frame);
        }
        return;
    }
    // A thread that is ending keeps nothing.
    let _ = X_FRAMES.try_with(|frames| {
        let mut frames = frames.borrow_mut();
        if frames.len() < X_REUSED_FRAMES {
            frames.push(frame);
        }
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1090
// NCBI translates the six frames from the packed subject, holding nothing but the frames.
// `x_generate_frames` on the ncbi2na values of the store, read from the packed bytes: the codon
// that starts at base `p` is residue `p / 3` of frame `p % 3 + 1`, and, read on the minus strand,
// residue `q / 3` of frame `-(q % 3 + 1)` where `q = length - 3 - p` (the codons of
// `seq[..length - i]` counted from its end). So the codons of one phase (`p % 3`) fill one plus
// frame forward and one minus frame backward, and every residue is the one `x_generate_frames`
// puts there.
/// EXPERIMENT (LOSAT_X_TBNSTORE): `x_generate_frames` of a packed subject's codes, from its packed
/// bytes (`PackedNucleotide::packed`) and its length `seq_len`.
fn x_packed_frames(packed: &[u8], seq_len: usize, tables: &XCodonTables) -> Vec<QueryFrame> {
    x_packed_frames_with(packed, seq_len, tables, x_xlate_simd())
}

/// EXPERIMENT (LOSAT_X_TBNSTORE / LOSAT_X_TBNXLATESIMD): [`x_packed_frames`], with the vector
/// loop when `simd` (and the CPU has it).
fn x_packed_frames_with(
    packed: &[u8],
    seq_len: usize,
    tables: &XCodonTables,
    simd: bool,
) -> Vec<QueryFrame> {
    let frame_count = (0..3).filter(|&i| i + 3 <= seq_len).count();
    // x_phase_residues writes every residue (positions 1 to aa_len) of the frames.
    let empty_frame = |i: usize| x_frame_buffer((seq_len - i) / 3 + 2);
    let mut forward: Vec<Vec<u8>> = (0..frame_count).map(empty_frame).collect();
    let mut reverse: Vec<Vec<u8>> = (0..frame_count).map(empty_frame).collect();
    // codons read by the vector loop, and in all (for `[X_STATS] TBNXLATESIMD`)
    let (mut vector, mut all) = (0, 0);
    for (phase, plus) in forward.iter_mut().enumerate() {
        // the codons at `phase`, `phase + 3`, ...: residues 1.. of plus frame `phase + 1` and,
        // from the last one down, of minus frame `-(r + 1)`, which has as many
        let codons = (seq_len - phase) / 3;
        let minus = &mut reverse[(seq_len - phase) % 3];
        vector += x_phase_residues(packed, phase, codons, tables, plus, minus, simd);
        all += codons;
    }
    if all > 0 && x_stats() {
        PACKED_CODONS[0].fetch_add(vector as u64, Ordering::Relaxed);
        PACKED_CODONS[1].fetch_add((all - vector) as u64, Ordering::Relaxed);
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:544-553
// The six bits of the codon at base `p`, from the packed bytes (4 bases per byte from the high
// bits): bits `2 * (p % 4)` on of byte `p / 4`, running into the next byte when `p % 4 >= 2`
// (which then holds base `p + 2`).
/// EXPERIMENT (LOSAT_X_TBNSTORE): the codon-table index of the codon at base `p`
/// (`p + 3 <= length`).
#[inline(always)]
fn x_packed_codon(packed: &[u8], p: usize) -> usize {
    let at = p >> 2;
    let pair = (u16::from(packed[at]) << 8) | u16::from(packed.get(at + 1).copied().unwrap_or(0));
    usize::from((pair >> (10 - 2 * (p & 3))) & 63)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596
// NCBI's fast loop reads three bytes for four codons. Here one phase's codons `j..j + 8` are
// the 48 bits after bit `2 * phase` of the eight bytes from byte `6 * j / 8` (24 bases per six
// bytes), read as one big-endian word; the vector loop (`x_phase_residues_avx2`) does 32 at a
// time first, the codons after the last whole group come one at a time.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the residues of the `codons` codons at bases `phase`,
/// `phase + 3`, ...: codon `j` at `plus[1 + j]` (forward table) and `minus[codons - j]` (reverse
/// table). Returns the codons the vector loop read.
fn x_phase_residues(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    plus: &mut [u8],
    minus: &mut [u8],
    simd: bool,
) -> usize {
    assert!(phase < 3 && plus.len() == codons + 2 && minus.len() == codons + 2);
    let mut done = 0;
    #[cfg(target_arch = "x86_64")]
    match x_vector_loop(simd) {
        XVectorLoop::Avx2 => {
            // SAFETY: the CPU has AVX2; the bounds are checked inside.
            done = unsafe { x_phase_residues_avx2(packed, phase, codons, tables, plus, minus) };
        }
        XVectorLoop::W128 => {
            // SAFETY: the CPU has SSSE3 and SSE4.1; the bounds are checked inside.
            done = unsafe { x_phase_residues_w128(packed, phase, codons, tables, plus, minus) };
        }
        XVectorLoop::Words => {}
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    if x_vector_loop(simd) == XVectorLoop::W128 {
        // SAFETY: simd128 is part of the build; the bounds are checked inside.
        done = unsafe { x_phase_residues_w128(packed, phase, codons, tables, plus, minus) };
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = simd;
    let vector = done;
    done = match phase {
        0 => x_phase_words::<0>(packed, codons, tables, plus, minus, done),
        1 => x_phase_words::<1>(packed, codons, tables, plus, minus, done),
        _ => x_phase_words::<2>(packed, codons, tables, plus, minus, done),
    };
    for j in done..codons {
        let codon = x_packed_codon(packed, phase + 3 * j);
        plus[1 + j] = tables.forward[codon];
        minus[codons - j] = tables.reverse[codon];
    }
    vector
}

/// EXPERIMENT (LOSAT_X_TBNSTORE): [`x_phase_residues`] by groups of eight codons from codon
/// `from` (a multiple of 8) while the group's eight bytes are in `packed`; returns the codons done.
#[inline(always)]
fn x_phase_words<const PHASE: usize>(
    packed: &[u8],
    codons: usize,
    tables: &XCodonTables,
    plus: &mut [u8],
    minus: &mut [u8],
    from: usize,
) -> usize {
    let mut j = from;
    while j + 8 <= codons && j / 8 * 6 + 8 <= packed.len() {
        let at = j / 8 * 6;
        let word = u64::from_be_bytes(packed[at..at + 8].try_into().unwrap());
        let plus8: &mut [u8; 8] = (&mut plus[1 + j..9 + j]).try_into().unwrap();
        let minus8: &mut [u8; 8] = (&mut minus[codons - j - 7..=codons - j])
            .try_into()
            .unwrap();
        for i in 0..8 {
            let codon = ((word >> (58 - 2 * PHASE - 6 * i)) & 63) as usize;
            plus8[i] = tables.forward[codon];
            minus8[7 - i] = tables.reverse[codon];
        }
        j += 8;
    }
    j
}

// No NCBI counterpart: reads LOSAT_X_TBNXLATESIMD once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNXLATESIMD): whether a packed subject's frames are read with the
/// vector loop (AVX2 or the 128-bit loop, `x_vector_loop`; where neither runs the variable
/// changes nothing); a switch of the strict set, on unless it is `0`.
fn x_xlate_simd() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNXLATESIMD"))
}

// No NCBI counterpart: whether the CPU has AVX2; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn x_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

// No NCBI counterpart: reads LOSAT_X_TBNXLATE128 once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNXLATE128): 0 = off, 1 = the 128-bit loop where the AVX2 loop does not
/// run (wasm32 simd128; x86-64 with SSSE3 and SSE4.1 but without AVX2), 2 = the 128-bit loop on
/// any x86-64 CPU with SSSE3 and SSE4.1, AVX2 or not (diagnostic: the tests run it end to end). A
/// switch of the strict set, on unless it is `0`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_xlate128() -> usize {
    static ON: OnceLock<usize> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_TBNXLATE128", 1..=2, 1).unwrap_or(0))
}

/// EXPERIMENT (LOSAT_X_TBNXLATESIMD, LOSAT_X_TBNXLATE128): the loop that reads a phase's codons
/// before the word loop (`x_phase_words` and the tail).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum XVectorLoop {
    /// 32 codons a step (`x_phase_residues_avx2`, `x_phase_one_avx2`).
    #[cfg(target_arch = "x86_64")]
    Avx2,
    /// 16 codons a step (`x_phase_residues_w128`, `x_phase_one_w128`).
    W128,
    /// None: the word loop from the first codon.
    Words,
}

// No NCBI counterpart: chooses LOSAT's loop; every loop gives the same residues.
/// EXPERIMENT (LOSAT_X_TBNXLATE128): the vector loop for `simd` (LOSAT_X_TBNXLATESIMD) on this
/// CPU: AVX2 where the CPU has it (unless LOSAT_X_TBNXLATE128 is `2`), else the 128-bit loop where
/// LOSAT_X_TBNXLATE128 is on and the CPU has SSSE3 and SSE4.1 (always on wasm32 simd128).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_vector_loop(simd: bool) -> XVectorLoop {
    if !simd {
        return XVectorLoop::Words;
    }
    let w128 = x_xlate128();
    #[cfg(target_arch = "x86_64")]
    {
        if w128 != 2 && x_avx2_available() {
            return XVectorLoop::Avx2;
        }
        if w128 != 0
            && std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("sse4.1")
        {
            return XVectorLoop::W128;
        }
        XVectorLoop::Words
    }
    #[cfg(not(target_arch = "x86_64"))]
    if w128 != 0 {
        XVectorLoop::W128
    } else {
        XVectorLoop::Words
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
// `x_phase_words`)
// The same residues, 32 codons at a time: the 24 bytes after bit `2 * phase` (made from two
// loads one byte apart), each three bytes spread to four bytes and split into four six-bit
// codons (bits 18-23, 12-17, 6-11, 0-5 of the three bytes: the shuffle puts the bytes in two
// 16-bit halves and two multiplications move each field to its byte), and each table of 64
// residues read as four tables of 16 (`vpshufb`) chosen by bits 4 and 5 of the index. The minus
// residues are stored in reverse order. Same reads, same tables as the scalar loop.
/// EXPERIMENT (LOSAT_X_TBNXLATESIMD): [`x_phase_residues`] for the groups of 32 codons whose
/// bytes can be read with 29 bytes from the group's first; returns the codons done (a multiple
/// of 32).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_phase_residues_avx2(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    plus: &mut [u8],
    minus: &mut [u8],
) -> usize {
    use std::arch::x86_64::*;
    let table = |values: &[u8; 64], k: usize| {
        // SAFETY: 16 bytes of the 64.
        _mm256_broadcastsi128_si256(unsafe { _mm_loadu_si128(values.as_ptr().add(16 * k).cast()) })
    };
    let forward = [0, 1, 2, 3].map(|k| table(&tables.forward, k));
    let reverse = [0, 1, 2, 3].map(|k| table(&tables.reverse, k));
    let spread = _mm256_setr_epi8(
        1, 0, 2, 1, 4, 3, 5, 4, 7, 6, 8, 7, 10, 9, 11, 10, 1, 0, 2, 1, 4, 3, 5, 4, 7, 6, 8, 7, 10,
        9, 11, 10,
    );
    let backward = _mm256_setr_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6,
        5, 4, 3, 2, 1, 0,
    );
    let left = _mm_cvtsi32_si128(2 * phase as i32);
    let right = _mm_cvtsi32_si128(8 - 2 * phase as i32);
    let keep_left = _mm256_set1_epi8((0xFFu8 << (2 * phase)) as i8);
    let keep_right = _mm256_set1_epi8(((1u16 << (2 * phase)) - 1) as u8 as i8);
    let look = |values: &[__m256i; 4], index: __m256i, bit4: __m256i, bit5: __m256i| {
        let low = _mm256_blendv_epi8(
            _mm256_shuffle_epi8(values[0], index),
            _mm256_shuffle_epi8(values[1], index),
            bit4,
        );
        let high = _mm256_blendv_epi8(
            _mm256_shuffle_epi8(values[2], index),
            _mm256_shuffle_epi8(values[3], index),
            bit4,
        );
        _mm256_blendv_epi8(low, high, bit5)
    };
    let mut j = 0;
    while j + 32 <= codons && 3 * j / 4 + 29 <= packed.len() {
        // SAFETY: bytes `3 * j / 4 ..+ 29` are in `packed` (checked above).
        let at = unsafe { packed.as_ptr().add(3 * j / 4) };
        let (bytes, next) = unsafe {
            (
                _mm256_loadu2_m128i(at.add(12).cast(), at.cast()),
                _mm256_loadu2_m128i(at.add(13).cast(), at.add(1).cast()),
            )
        };
        // byte k of each half: bits 2 * phase .. 2 * phase + 8 of bytes k, k + 1
        let shifted = _mm256_or_si256(
            _mm256_and_si256(_mm256_sll_epi16(bytes, left), keep_left),
            _mm256_and_si256(_mm256_srl_epi16(next, right), keep_right),
        );
        let spread = _mm256_shuffle_epi8(shifted, spread);
        let index = _mm256_or_si256(
            _mm256_mulhi_epu16(
                _mm256_and_si256(spread, _mm256_set1_epi32(0x0FC0_FC00)),
                _mm256_set1_epi32(0x0400_0040),
            ),
            _mm256_mullo_epi16(
                _mm256_and_si256(spread, _mm256_set1_epi32(0x003F_03F0)),
                _mm256_set1_epi32(0x0100_0010),
            ),
        );
        // bits 4 and 5 of each index in the byte's top bit (a 16-bit shift: the bit that
        // crosses into the next byte lands below its top bit)
        let bit4 = _mm256_slli_epi16(index, 3);
        let bit5 = _mm256_slli_epi16(index, 2);
        let plus_residues = look(&forward, index, bit4, bit5);
        let minus_residues = _mm256_permute4x64_epi64(
            _mm256_shuffle_epi8(look(&reverse, index, bit4, bit5), backward),
            0x4E,
        );
        // SAFETY: `plus[1 + j..1 + j + 32]` and `minus[codons - j - 31..=codons - j]` lie in the
        // frames (`j + 32 <= codons`, length `codons + 2`).
        unsafe {
            _mm256_storeu_si256(plus.as_mut_ptr().add(1 + j).cast(), plus_residues);
            _mm256_storeu_si256(
                minus.as_mut_ptr().add(codons - j - 31).cast(),
                minus_residues,
            );
        }
        j += 32;
    }
    j
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// The plan of a packed subject as `x_table_plan` makes it from the subject's letters. Only a
// stage that does not read the case of the letters, or a subject without lower-case letters,
// takes the store (`-lcase_masking` gives the letters of the others), so every frame has one unit
// per chunk with one range, as without lowercase masking.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the plan of one packed subject, from its codes.
fn x_packed_plan(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
) -> Result<XSubjectPlan> {
    x_bytes_plan(store.packed(), store.len(), db_gencode, None)
}

/// The database masking of a subject's plan (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b
/// PR 2): the sequence's ranges in bases (`XDbStream::mask_ranges`, the complement of its stored
/// masks) and whether the masking is hard.
#[derive(Clone, Copy)]
pub(super) struct XDbMaskRanges<'r> {
    pub(super) ranges: &'r [(i32, i32)],
    pub(super) hard: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-830
// In a masked search the frames' ranges are mapped from the sequence's ranges in bases, once per
// strand, and each frame's chunks and scan ranges follow from them as for a protein subject
// (`translated_db_mask_chunks`).
/// EXPERIMENT (LOSAT_X_TBNSTORE): [`x_packed_plan`] of a subject's packed bytes and length, with
/// the units of its database masking `db_mask` (`None`: a search without the options).
fn x_bytes_plan(
    packed: &[u8],
    length: usize,
    db_gencode: u8,
    db_mask: Option<XDbMaskRanges<'_>>,
) -> Result<XSubjectPlan> {
    let tables = x_codon_tables(db_gencode)?;
    let frames = x_packed_frames(packed, length, tables);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    let mut units = Vec::new();
    let mut split = crate::core::subject_split::HardSplitState::default();
    for (index, frame) in frames.iter().enumerate() {
        // Without lowercase masking the ranges read no letter of the subject.
        let calls = match db_mask {
            None => translated_chunk_scan_ranges(
                &[],
                frame.frame,
                frame.aa_len,
                negative_first_length,
                false,
            )?,
            Some(mask) => translated_db_mask_chunks(
                mask.ranges,
                mask.hard,
                frame.frame,
                frame.aa_len,
                negative_first_length,
                &mut split,
            )?,
        };
        for (chunk, ranges) in calls {
            units.push((index, chunk, ranges));
        }
    }
    super::x_subject_store::x_count_packed_plan();
    Ok(XSubjectPlan { frames, units })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// A packed subject is translated from its codes (`x_packed_plan`).
/// EXPERIMENT (LOSAT_X_TBNSTORE): the plan of one packed subject.
pub(super) fn x_packed_subject_plan(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
) -> Result<Arc<XSubjectPlan>> {
    x_packed_plan(store, db_gencode).map(Arc::new)
}

/// [`x_packed_subject_plan`] of a database sequence's packed
/// bytes and length, read in a block of the sequence file (`XDbStream::with_block`), without its
/// store: the plan reads nothing else of a packed subject.
pub(super) fn x_packed_bytes_plan(
    packed: &[u8],
    length: usize,
    db_gencode: u8,
    db_mask: Option<XDbMaskRanges<'_>>,
) -> Result<Arc<XSubjectPlan>> {
    x_bytes_plan(packed, length, db_gencode, db_mask).map(Arc::new)
}

// No NCBI counterpart: reads LOSAT_X_TBNUNITXLATE once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): whether the units of a long packed subject translate their
/// own residues (`x_packed_lean_plan`, `x_unit_residues`) instead of reading the frames of the
/// plan. A switch of the strict set (`utils/xswitch.rs`), on unless it is `0`.
fn x_units_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNUNITXLATE"))
}

// No NCBI counterpart for the rule: NCBI translates the six frames of every subject into one
// buffer (blast_engine.c:772-775). A frame of `X_MAPPED_FRAME` bytes or more is mapped when it is
// made and unmapped when its plan is dropped (`x_frame_capacity`); on a genome database that is
// six fresh frames of up to 83 MB per chromosome, made by one thread before the units of the
// chromosome start, and their pages are faulted in and returned again (human genome, one protein,
// 16 threads: 5.2 GiB at most, most of it the frames of the chromosomes searched at the same
// time). A unit reads its chunk and nothing else.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): whether the plan of a packed subject of `length` bases
/// made for its use alone (a search that keeps no plan of it: no `XSubjectSide`, or one whose
/// limit of kept plans refuses it, `XSubjectSide::keeps`) is made
/// without residues (`x_packed_lean_plan`): with the switch, when its first frame has
/// `X_MAPPED_FRAME` bytes or more (the frames that are mapped on Linux with glibc; the rule is the
/// same everywhere).
pub(super) fn x_units_translate(length: usize) -> bool {
    x_units_on() && length / 3 + 2 >= X_MAPPED_FRAME
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// The frames and units of `x_packed_plan`, without residues: a frame's length is that of the
// frame `x_packed_frames` makes, `(length - i) / 3` for frames `i + 1` and `-(i + 1)`, and the
// units follow from the lengths alone (no lowercase masking).
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the plan of one packed subject whose frames have an empty
/// `aa_seq`; each unit makes its residues with `x_unit_residues`.
pub(super) fn x_packed_lean_plan(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
) -> Result<Arc<XSubjectPlan>> {
    // An unknown genetic code is refused here, as by x_packed_plan.
    x_codon_tables(db_gencode)?;
    let seq_len = store.len();
    let mut frames = Vec::with_capacity(6);
    for sign in [1i8, -1] {
        for i in 0..3usize {
            if i + 3 <= seq_len {
                frames.push(QueryFrame {
                    frame: sign * (i as i8 + 1),
                    aa_seq: Vec::new(),
                    aa_seq_nomask: None,
                    aa_len: (seq_len - i) / 3,
                    orig_len: seq_len,
                    seg_masks: Vec::new(),
                });
            }
        }
    }
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    let mut units = Vec::new();
    for (index, frame) in frames.iter().enumerate() {
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
    if x_stats() {
        UNIT_RESIDUES[0].fetch_add(1, Ordering::Relaxed);
    }
    Ok(Arc::new(XSubjectPlan { frames, units }))
}

// No NCBI counterpart for the margin: NCBI's scan reads the letters of its chunk
// (aa_ungapped.c:496-501), as the scan loops here do (their blocks follow the scan range); one of
// them (x_scan_byte) also checks that its last word fits in the slice it is given. With a margin
// past the chunk, or the frame's end sentinel when it is nearer, a unit's slice is longer than
// every read, as the frame's slice is; 4,096 residues are 0.08 % of a chunk of 5,000,000.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the residues a unit's translation runs past its chunk.
const X_UNIT_MARGIN: usize = 4096;

// No NCBI counterpart: the range of a unit's translation (`X_UNIT_MARGIN`); it computes no value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): positions `from..to` of the frame's `aa_seq` that the unit
/// of `chunk` reads: from the chunk's first residue to `X_UNIT_MARGIN` past its end, or to the
/// end sentinel.
fn x_unit_range(aa_len: usize, chunk: TranslatedChunk) -> (usize, usize) {
    let from = 1 + chunk.offset;
    (from, (aa_len + 2).min(from + chunk.length + X_UNIT_MARGIN))
}

// No NCBI counterpart: counters for LOSAT_X_STATS; they do not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): counts `made` buffers for the units' residues and `freed`
/// ones (`search_gapped::XUnitResiduesPool`), with LOSAT_X_STATS.
pub(super) fn x_count_unit_buffers(made: u64, freed: u64) {
    if x_stats() {
        UNIT_RESIDUES[3].fetch_add(made, Ordering::Relaxed);
        UNIT_RESIDUES[4].fetch_add(freed, Ordering::Relaxed);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812
// The part of `subject->sequence` a unit reads, made when the unit runs.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the residues of `frame` (of a plan without residues) that
/// the unit of `chunk` reads, made into `buffer`; returns the range of `buffer` that holds them,
/// from the chunk's first residue on (see `x_unit_range`).
pub(super) fn x_unit_residues(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    db_gencode: u8,
    frame: &QueryFrame,
    chunk: TranslatedChunk,
    buffer: &mut Vec<u8>,
) -> Result<std::ops::Range<usize>> {
    let tables = x_codon_tables(db_gencode)?;
    let (from, to) = x_unit_range(frame.aa_len, chunk);
    let (range, codons, vector) = x_frame_residues(
        store,
        tables,
        frame.frame,
        frame.aa_len,
        from,
        to,
        buffer,
        x_xlate_simd(),
    );
    if x_stats() {
        UNIT_RESIDUES[1].fetch_add(1, Ordering::Relaxed);
        UNIT_RESIDUES[2].fetch_add(range.len() as u64, Ordering::Relaxed);
        PACKED_CODONS[0].fetch_add(vector as u64, Ordering::Relaxed);
        PACKED_CODONS[1].fetch_add((codons - vector) as u64, Ordering::Relaxed);
    }
    Ok(range)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1090 (see
// `x_packed_frames`)
// The residues at positions `from..to` of the `aa_seq` that `x_packed_frames` makes for frame
// `frame`: residue `i` of plus frame `r + 1` is the codon at base `r + 3 * i`; residue `i` of
// minus frame `-(r + 1)` is codon `aa_len - 1 - i` of the phase `p` with `(length - p) % 3 == r`,
// read from the reverse table. The codons are read from the packed bytes that start with the one
// holding the first codon's first base; when that base is the byte's last (a phase of 3 within
// the byte, which the loops do not take) the codon before it is read too and left out.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): positions `from..to` (`1 <= from <= to <= aa_len + 2`) of
/// frame `frame` of `store` (`aa_len` residues), made into `buffer` by one strand of one phase;
/// returns the range of `buffer` that holds them, the codons read and those the vector loop read.
/// `simd`: the vector loop where the CPU has it. The buffer only grows (its other bytes are left
/// as they are), so a buffer is filled with zeros once; from 1 MiB on its capacity is one glibc
/// maps (`utils::x_mapped::x_reserve_mapped`), so that its pages are returned when it is freed.
#[allow(clippy::too_many_arguments)]
fn x_frame_residues(
    store: &crate::sequence::packed_subject::PackedNucleotide,
    tables: &XCodonTables,
    frame: i8,
    aa_len: usize,
    from: usize,
    to: usize,
    buffer: &mut Vec<u8>,
    simd: bool,
) -> (std::ops::Range<usize>, usize, usize) {
    let seq_len = store.len();
    assert!(1 <= from && from <= to && to <= aa_len + 2);
    // residues `i0..i1`; position `aa_len + 1` is the end sentinel
    let (i0, i1) = (from - 1, (to - 1).min(aa_len));
    let sentinel = from < to && to == aa_len + 2;
    if i0 >= i1 {
        if !sentinel {
            return (0..0, 0, 0);
        }
        if buffer.is_empty() {
            buffer.push(0);
        }
        buffer[0] = SENTINEL_BYTE;
        return (0..1, 0, 0);
    }
    let r = usize::from(frame.unsigned_abs()) - 1;
    let phase = if frame > 0 { r } else { (seq_len - r) % 3 };
    let codons = (seq_len - phase) / 3;
    assert_eq!(codons, aa_len, "frame {frame} of {seq_len} bases");
    // the first codon read, in the phase's order (a minus frame's residues run backwards)
    let mut first_codon = if frame > 0 { i0 } else { codons - i1 };
    let mut back = 0;
    if (phase + 3 * first_codon) % 4 == 3 {
        first_codon -= 1;
        back = 1;
    }
    let count = i1 - i0 + back;
    let first_base = phase + 3 * first_codon;
    let byte = first_base / 4;
    let packed = store.packed();
    // the bytes of the codons, and 32 more for the groups of the loops
    let end = packed
        .len()
        .min(byte + (first_base % 4 + 3 * count).div_ceil(4) + 32);
    if buffer.len() < count + 2 {
        crate::utils::x_mapped::x_reserve_mapped(buffer, count + 2 - buffer.len());
        buffer.resize(count + 2, 0);
    }
    let bytes = &packed[byte..end];
    let out = &mut buffer[..count + 2];
    let vector = if frame > 0 {
        x_phase_one::<false>(bytes, first_base % 4, count, tables, out, simd)
    } else {
        x_phase_one::<true>(bytes, first_base % 4, count, tables, out, simd)
    };
    // codon `j` read is at `buffer[1 + j]` (plus) or `buffer[count - j]` (minus)
    let start = if frame > 0 { 1 + back } else { 1 };
    let mut length = i1 - i0;
    if sentinel {
        buffer[start + length] = SENTINEL_BYTE;
        length += 1;
    }
    (start..start + length, count, vector)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
// `x_phase_words`)
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): [`x_phase_residues`] for one strand: codon `j` at
/// `out[1 + j]` (forward table) or, with `MINUS`, at `out[codons - j]` (reverse table). Returns
/// the codons the vector loop read.
fn x_phase_one<const MINUS: bool>(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    out: &mut [u8],
    simd: bool,
) -> usize {
    assert!(phase < 3 && out.len() == codons + 2);
    let mut j = 0;
    #[cfg(target_arch = "x86_64")]
    match x_vector_loop(simd) {
        XVectorLoop::Avx2 => {
            // SAFETY: the CPU has AVX2; the bounds are checked inside.
            j = unsafe { x_phase_one_avx2::<MINUS>(packed, phase, codons, tables, out) };
        }
        XVectorLoop::W128 => {
            // SAFETY: the CPU has SSSE3 and SSE4.1; the bounds are checked inside.
            j = unsafe { x_phase_one_w128::<MINUS>(packed, phase, codons, tables, out) };
        }
        XVectorLoop::Words => {}
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    if x_vector_loop(simd) == XVectorLoop::W128 {
        // SAFETY: simd128 is part of the build; the bounds are checked inside.
        j = unsafe { x_phase_one_w128::<MINUS>(packed, phase, codons, tables, out) };
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = simd;
    let vector = j;
    let table = if MINUS {
        &tables.reverse
    } else {
        &tables.forward
    };
    let shift = 58 - 2 * phase;
    while j + 8 <= codons && j / 8 * 6 + 8 <= packed.len() {
        let at = j / 8 * 6;
        let word = u64::from_be_bytes(packed[at..at + 8].try_into().unwrap());
        if MINUS {
            let out8: &mut [u8; 8] = (&mut out[codons - j - 7..=codons - j]).try_into().unwrap();
            for i in 0..8 {
                out8[7 - i] = table[((word >> (shift - 6 * i)) & 63) as usize];
            }
        } else {
            let out8: &mut [u8; 8] = (&mut out[1 + j..9 + j]).try_into().unwrap();
            for i in 0..8 {
                out8[i] = table[((word >> (shift - 6 * i)) & 63) as usize];
            }
        }
        j += 8;
    }
    for j in j..codons {
        let codon = x_packed_codon(packed, phase + 3 * j);
        if MINUS {
            out[codons - j] = table[codon];
        } else {
            out[1 + j] = table[codon];
        }
    }
    vector
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
// `x_phase_residues_avx2`)
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): [`x_phase_residues_avx2`] for one strand; returns the
/// codons done (a multiple of 32).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_phase_one_avx2<const MINUS: bool>(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    out: &mut [u8],
) -> usize {
    use std::arch::x86_64::*;
    let values = if MINUS {
        &tables.reverse
    } else {
        &tables.forward
    };
    // SAFETY: 16 bytes of the 64.
    let table = |k: usize| {
        _mm256_broadcastsi128_si256(unsafe { _mm_loadu_si128(values.as_ptr().add(16 * k).cast()) })
    };
    let t = [table(0), table(1), table(2), table(3)];
    let spread = _mm256_setr_epi8(
        1, 0, 2, 1, 4, 3, 5, 4, 7, 6, 8, 7, 10, 9, 11, 10, 1, 0, 2, 1, 4, 3, 5, 4, 7, 6, 8, 7, 10,
        9, 11, 10,
    );
    let backward = _mm256_setr_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6,
        5, 4, 3, 2, 1, 0,
    );
    let left = _mm_cvtsi32_si128(2 * phase as i32);
    let right = _mm_cvtsi32_si128(8 - 2 * phase as i32);
    let keep_left = _mm256_set1_epi8((0xFFu8 << (2 * phase)) as i8);
    let keep_right = _mm256_set1_epi8(((1u16 << (2 * phase)) - 1) as u8 as i8);
    let mut j = 0;
    while j + 32 <= codons && 3 * j / 4 + 29 <= packed.len() {
        // SAFETY: bytes `3 * j / 4 ..+ 29` are in `packed` (checked above).
        let at = unsafe { packed.as_ptr().add(3 * j / 4) };
        let (bytes, next) = unsafe {
            (
                _mm256_loadu2_m128i(at.add(12).cast(), at.cast()),
                _mm256_loadu2_m128i(at.add(13).cast(), at.add(1).cast()),
            )
        };
        // the codons of the group as in `x_phase_residues_avx2`
        let shifted = _mm256_or_si256(
            _mm256_and_si256(_mm256_sll_epi16(bytes, left), keep_left),
            _mm256_and_si256(_mm256_srl_epi16(next, right), keep_right),
        );
        let spread = _mm256_shuffle_epi8(shifted, spread);
        let index = _mm256_or_si256(
            _mm256_mulhi_epu16(
                _mm256_and_si256(spread, _mm256_set1_epi32(0x0FC0_FC00)),
                _mm256_set1_epi32(0x0400_0040),
            ),
            _mm256_mullo_epi16(
                _mm256_and_si256(spread, _mm256_set1_epi32(0x003F_03F0)),
                _mm256_set1_epi32(0x0100_0010),
            ),
        );
        let bit4 = _mm256_slli_epi16(index, 3);
        let bit5 = _mm256_slli_epi16(index, 2);
        let low = _mm256_blendv_epi8(
            _mm256_shuffle_epi8(t[0], index),
            _mm256_shuffle_epi8(t[1], index),
            bit4,
        );
        let high = _mm256_blendv_epi8(
            _mm256_shuffle_epi8(t[2], index),
            _mm256_shuffle_epi8(t[3], index),
            bit4,
        );
        let residues = _mm256_blendv_epi8(low, high, bit5);
        // SAFETY: `out[1 + j..1 + j + 32]` and `out[codons - j - 31..=codons - j]` lie in `out`
        // (`j + 32 <= codons`, length `codons + 2`).
        unsafe {
            if MINUS {
                _mm256_storeu_si256(
                    out.as_mut_ptr().add(codons - j - 31).cast(),
                    _mm256_permute4x64_epi64(_mm256_shuffle_epi8(residues, backward), 0x4E),
                );
            } else {
                _mm256_storeu_si256(out.as_mut_ptr().add(1 + j).cast(), residues);
            }
        }
        j += 32;
    }
    j
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
// `x_phase_words`)
// The loop of `x_phase_residues_avx2` with one 128-bit register: 16 codons a step, written once
// over the 128-bit layer of `utils/xdrop_simd.rs` for x86-64 (SSSE3 and SSE4.1) and wasm32
// (simd128), in the form V8 compiles to the fewest instructions. The loop starts at codon `phase`
// of the phase, the first whose first base is the first of a byte (base `4 * phase`, the first of
// byte `phase`), so that each step reads the 12 bytes from a byte with one load; codons
// `0..phase` are read one at a time. Each three bytes `b0 b1 b2` (four codons) are spread to the 16-bit halves `b0 b1` and
// `b1 b2`, the second moved up four bits by a multiplication, so that every half holds its two
// codons at bits 10-15 and 4-9 and two shifts give the residue index of every byte. Each table of
// 64 residues is read as four tables of 16 at the index minus 0, 16, 32 and 48: an index byte
// outside a table's 16 entries gives 0 (`pickz`, the rule of wasm32's `i8x16.swizzle`), so the
// four reads are or-ed. The constants are made once before the loop (`hoisted`): in a loop LLVM
// makes a wasm32 vector constant again at every use, which V8 builds with three or four
// instructions, and turns the multiplication by powers of two into shifts by lane, which wasm32
// lacks (eight lane moves). Same reads, same tables as the scalar loop.

/// EXPERIMENT (LOSAT_X_TBNXLATE128): the constants of the 128-bit loops, made once.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[derive(Clone, Copy)]
struct XW128 {
    /// The bytes `b1 b0 b2 b1` of every three bytes `b0 b1 b2`: the 16-bit halves `b0 b1`, `b1 b2`.
    spread: crate::utils::xdrop_simd::arch::V,
    /// 1 in the low 16-bit half of every 32-bit lane, 16 in the high one.
    up: crate::utils::xdrop_simd::arch::V,
    /// Bits 8-13 of every 16-bit half.
    high: crate::utils::xdrop_simd::arch::V,
    /// -16 in every byte: an index minus 16, 32, 48 is that of the second, third, fourth 16
    /// entries.
    minus16: crate::utils::xdrop_simd::arch::V,
    /// The byte order that reverses a register's 16 bytes.
    backward: crate::utils::xdrop_simd::arch::V,
}

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
impl XW128 {
    // No NCBI counterpart: the constants of LOSAT's vector loop; they are no value NCBI computes.
    /// The constants (each `hoisted`).
    #[inline(always)]
    unsafe fn new() -> Self {
        use crate::utils::xdrop_simd::arch::*;
        const SPREAD: [u8; 16] = [1, 0, 2, 1, 4, 3, 5, 4, 7, 6, 8, 7, 10, 9, 11, 10];
        const BACKWARD: [u8; 16] = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0];
        // SAFETY (the block): loads of the two 16-byte constants; register operations.
        unsafe {
            Self {
                spread: hoisted(load(SPREAD.as_ptr())),
                up: hoisted(splat_i32(0x0010_0001)),
                high: hoisted(splat16(0x3F00)),
                minus16: hoisted(splat8(-16)),
                backward: hoisted(load(BACKWARD.as_ptr())),
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
    // `x_phase_words`)
    /// One step's residue indexes: the 16 codons of the 12 bytes from `at`, the first codon's
    /// first base the first of `at`'s byte (`at + 16` readable). Codon 0 of four (bytes `b0 b1
    /// b2`) is bits 10-15 of the low 16-bit half (`b0 b1`), codon 1 bits 4-9, codon 2 bits 10-15
    /// of the high half (`b1 b2`, moved up four bits), codon 3 bits 4-9.
    #[inline(always)]
    unsafe fn index(&self, at: *const u8) -> crate::utils::xdrop_simd::arch::V {
        use crate::utils::xdrop_simd::arch::*;
        // SAFETY (the block): the caller makes bytes `at ..+ 16` readable.
        unsafe {
            let halves = mullo16(pick(load(at), self.spread), self.up);
            or(shr16n(halves, 10), and(shl16n(halves, 4), self.high))
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:595 (see `x_phase_words`)
    /// Entry `index` (0..64) of the 64-byte table `table`, held as four registers of 16 entries,
    /// for every byte of `index`.
    #[inline(always)]
    unsafe fn look(
        &self,
        table: &[crate::utils::xdrop_simd::arch::V; 4],
        index: crate::utils::xdrop_simd::arch::V,
    ) -> crate::utils::xdrop_simd::arch::V {
        use crate::utils::xdrop_simd::arch::*;
        // SAFETY (the block): register operations of the build's 128-bit instructions.
        unsafe {
            let second = add8(index, self.minus16);
            let third = add8(second, self.minus16);
            let fourth = add8(third, self.minus16);
            or(
                or(pickz(table[0], index), pickz(table[1], second)),
                or(pickz(table[2], third), pickz(table[3], fourth)),
            )
        }
    }
}

/// EXPERIMENT (LOSAT_X_TBNXLATE128): the 16-entry quarters of a 64-residue codon table.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn x_table_w128(values: &[u8; 64]) -> [crate::utils::xdrop_simd::arch::V; 4] {
    // SAFETY: bytes `16 * k ..+ 16` of the 64.
    [0, 1, 2, 3]
        .map(|k| unsafe { crate::utils::xdrop_simd::arch::load(values.as_ptr().add(16 * k)) })
}

// No NCBI counterpart: the bounds of LOSAT's vector loop; they are no value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNXLATE128): the steps of the 128-bit loops of a phase of `codons`
/// codons whose first base is base `phase` of `bytes` bytes: step `s` reads codons `phase + 16 *
/// s ..+ 16` from bytes `phase + 12 * s ..+ 16`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
fn x_steps_w128(bytes: usize, phase: usize, codons: usize) -> usize {
    let by_codons = codons.saturating_sub(phase) / 16;
    let by_bytes = bytes
        .checked_sub(phase + 16)
        .map_or(0, |spare| spare / 12 + 1);
    by_codons.min(by_bytes)
}

/// EXPERIMENT (LOSAT_X_TBNXLATE128): [`x_phase_residues`] for the codons of the groups of 16
/// codons from codon `phase` that can be read with 16 bytes from the group's first byte, and the
/// `phase` codons before them; returns `16 *` the groups (every codon below it is written; the
/// loop may also write the `phase` codons after it). The 128-bit form of
/// [`x_phase_residues_avx2`].
///
/// # Safety
/// On x86-64 the CPU has SSSE3 and SSE4.1.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1"))]
unsafe fn x_phase_residues_w128(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    plus: &mut [u8],
    minus: &mut [u8],
) -> usize {
    use crate::utils::xdrop_simd::arch::*;
    assert!(phase < 3 && plus.len() == codons + 2 && minus.len() == codons + 2);
    let steps = x_steps_w128(packed.len(), phase, codons);
    // SAFETY (the block): step `s` reads `packed[phase + 12 * s ..+ 16]` and the 64-byte tables
    // and writes `plus[1 + j ..+ 16]` and `minus[codons - j - 15 ..= codons - j]` for `j = phase +
    // 16 * s` (`x_steps_w128`: `j + 16 <= codons`, length `codons + 2`).
    unsafe {
        let k = XW128::new();
        let forward = x_table_w128(&tables.forward);
        let reverse = x_table_w128(&tables.reverse);
        for step in 0..steps {
            let j = phase + 16 * step;
            let index = k.index(packed.as_ptr().add(phase + 12 * step));
            store(plus.as_mut_ptr().add(1 + j), k.look(&forward, index));
            store(
                minus.as_mut_ptr().add(codons - j - 15),
                pick(k.look(&reverse, index), k.backward),
            );
        }
    }
    if steps > 0 {
        for j in 0..phase {
            let codon = x_packed_codon(packed, phase + 3 * j);
            plus[1 + j] = tables.forward[codon];
            minus[codons - j] = tables.reverse[codon];
        }
    }
    16 * steps
}

/// EXPERIMENT (LOSAT_X_TBNXLATE128): [`x_phase_one`]'s vector loop with 16 codons a step, as
/// [`x_phase_residues_w128`] for one strand; returns `16 *` the groups. The 128-bit form of
/// [`x_phase_one_avx2`].
///
/// # Safety
/// On x86-64 the CPU has SSSE3 and SSE4.1.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1"))]
unsafe fn x_phase_one_w128<const MINUS: bool>(
    packed: &[u8],
    phase: usize,
    codons: usize,
    tables: &XCodonTables,
    out: &mut [u8],
) -> usize {
    use crate::utils::xdrop_simd::arch::*;
    assert!(phase < 3 && out.len() == codons + 2);
    let values = if MINUS {
        &tables.reverse
    } else {
        &tables.forward
    };
    let steps = x_steps_w128(packed.len(), phase, codons);
    // SAFETY (the block): step `s` reads `packed[phase + 12 * s ..+ 16]` and the 64-byte table
    // and writes `out[1 + j ..+ 16]` or `out[codons - j - 15 ..= codons - j]` for `j = phase + 16
    // * s` (`x_steps_w128`: `j + 16 <= codons`, length `codons + 2`).
    unsafe {
        let k = XW128::new();
        let table = x_table_w128(values);
        for step in 0..steps {
            let j = phase + 16 * step;
            let residues = k.look(&table, k.index(packed.as_ptr().add(phase + 12 * step)));
            if MINUS {
                store(
                    out.as_mut_ptr().add(codons - j - 15),
                    pick(residues, k.backward),
                );
            } else {
                store(out.as_mut_ptr().add(1 + j), residues);
            }
        }
    }
    if steps > 0 {
        for j in 0..phase {
            let residue = values[x_packed_codon(packed, phase + 3 * j)];
            if MINUS {
                out[codons - j] = residue;
            } else {
                out[1 + j] = residue;
            }
        }
    }
    16 * steps
}

#[cfg(test)]
#[path = "x_subject_translate_tests.rs"]
mod tests;

use super::*;
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject_one_range;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;

struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// The frame tblastn's `encode_tblastn_lookup_query(query, None, false)` gives a query of the
// 20 common letters (blastp's `encode_protein_query_frame_with_seg` without SEG), built here
// so that the tests stay below tblastn (the layering of `tests/check_program_edges.py`).
fn protein_frame(query: &[u8]) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(query.len() + 2);
    aa_seq.push(0);
    aa_seq.extend(
        query
            .iter()
            .map(|&a| crate::utils::matrix::aa_char_to_ncbistdaa(a)),
    );
    aa_seq.push(0);
    QueryFrame {
        frame: 0,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: query.len(),
        orig_len: query.len(),
        seg_masks: Vec::new(),
    }
}

fn lookup_for(queries: &[Vec<u8>], threshold: i32) -> BlastAaLookupTable {
    let frames: Vec<_> = queries
        .iter()
        .map(|query| vec![protein_frame(query)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    build_ncbi_lookup_for_profile(
        &frames,
        threshold,
        &karlin,
        false,
        ScoringMatrix::Blosum62,
        3,
    )
    .0
}

/// The scan calls of one range as `x_scan_unit` makes them: (offset pairs, s_range[1]) after
/// each call; `avx2` = `None` for NCBI's loop. `RESUMED` counts the calls that resume the
/// block held by the call before: a block is held and was kept for this call's subject and
/// `s_range` (`resume`: one `XScanHeld` for the range, as `LOSAT_X_TBNSCANRESUME` keeps it;
/// else a new one per call).
fn calls(
    lookup: &BlastAaLookupTable,
    presence: &XPresence,
    subject: &[u8],
    range: (i32, i32),
    array_size: usize,
    form: Option<(bool, bool)>,
) -> Vec<(Vec<(u32, u32)>, i32)> {
    calls_held(lookup, presence, subject, range, array_size, form, false)
}

thread_local! {
    static RESUMED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// `calls`, with one `XScanHeld` for the range when `resume`. `form` = `scan`'s `(avx2, run)`;
/// `None` = NCBI's scan.
fn calls_held(
    lookup: &BlastAaLookupTable,
    presence: &XPresence,
    subject: &[u8],
    range: (i32, i32),
    array_size: usize,
    form: Option<(bool, bool)>,
    resume: bool,
) -> Vec<(Vec<(u32, u32)>, i32)> {
    let mut held = XScanHeld::new();
    // Guard entries after the slack: no form may write past `array_size + X_SCAN_BYTE_SLACK`.
    let guard = BlastOffsetPair {
        q_off: 0xdead_beef,
        s_off: 0xfeed_f00d,
    };
    let end = array_size + X_SCAN_BYTE_SLACK;
    let mut pairs = vec![BlastOffsetPair::default(); end];
    pairs.resize(end + 2 * X_SCAN_BYTE_SLACK, guard);
    let mut scan_range = [0, range.0, range.1];
    let mut out = Vec::new();
    while scan_range[1] <= scan_range[2] {
        let hits = match form {
            Some((avx2, run)) => {
                let at = (subject.as_ptr() as usize, subject.len());
                if !resume {
                    held = XScanHeld::new();
                } else if held.next < held.n
                    && held.at == (at.0, at.1, scan_range[1], scan_range[2])
                {
                    RESUMED.with(|r| r.set(r.get() + 1));
                }
                scan(
                    avx2,
                    run,
                    lookup,
                    presence,
                    subject,
                    &mut pairs,
                    array_size as i32,
                    &mut scan_range,
                    &mut held,
                )
            }
            None => s_blast_aa_scan_subject_one_range(
                lookup,
                subject,
                &mut pairs,
                array_size as i32,
                &mut scan_range,
            ),
        };
        let got = pairs[..hits as usize]
            .iter()
            .map(|pair| (pair.q_off, pair.s_off))
            .collect();
        out.push((got, scan_range[1]));
        assert!(
            pairs[end..]
                .iter()
                .all(|pair| (pair.q_off, pair.s_off) == (guard.q_off, guard.s_off)),
            "write past the slack (array_size {array_size}, form {form:?})"
        );
        assert!(out.len() < 1_000_000, "scan does not advance");
    }
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127 (the function
// compared against, as `s_blast_aa_scan_subject_one_range`). Not a port of an NCBI unit case:
// the pairs and the `s_range[1]` of every call must equal the reference's, with offset arrays
// that stop the scan after a few words (as aascan_unit_test.cpp's ScanOffsetTest scans), cells
// past `AA_HITS_PER_CELL` hits, every letter value of ncbistdaa, and block ends at every
// alignment.
#[test]
fn x_scan_byte_equals_ncbi_scan_on_random_subjects() {
    RESUMED.with(|r| r.set(0));
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    // (avx2, run): the plain loops, the AVX2 forms, and those with the chain runs
    // (LOSAT_X_CHAINRUN).
    let forms: Vec<(bool, bool)> = if avx2_detected() {
        vec![(false, false), (true, false), (true, true)]
    } else {
        vec![(false, false)]
    };
    let mut compared = 0usize;
    let mut stops = 0usize;
    let mut big_cells = 0usize;
    // Chains longer than a run, and chains too near the overflow array's end for one.
    let mut past_run = 0usize;
    let mut near_end = 0usize;
    for case in 0..24 {
        let num_queries = [1, 1, 3, 12][case % 4];
        let queries: Vec<Vec<u8>> = (0..num_queries)
            .map(|_| {
                let len = 20 + rng.below(380);
                (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
            })
            .collect();
        let threshold = [11, 13, 15][case % 3];
        let lookup = lookup_for(&queries, threshold);
        assert!(fits(&lookup), "tblastn's table must fit");
        let presence = build_presence(&lookup);
        big_cells += lookup
            .backbone
            .iter()
            .filter(|cell| cell.num_used as usize > AA_HITS_PER_CELL)
            .count();
        for cell in &lookup.backbone {
            if cell.num_used as usize > AA_HITS_PER_CELL {
                past_run += usize::from(cell.num_used as usize > X_CHAIN_RUN);
                near_end += usize::from(cell.cursor() + X_CHAIN_RUN > lookup.overflow.len());
            }
        }
        let longest = lookup.longest_chain.max(1) as usize;
        // A frame in ncbistdaa (0..=27; a translated subject has 1..=25), sentinel 0 after it.
        let len = [3, 4, 255, 256, 257, 271, 1_000, 20_000][case % 8];
        let mut subject: Vec<u8> = (0..len).map(|_| rng.below(28) as u8).collect();
        subject.push(0);
        let last = len - 3;
        let ranges = [
            (0, last.max(0)),
            (1.min(last.max(0)), last.max(0)),
            (len / 3, last.max(len / 3)),
            (5.min(last.max(0)), (last - 7).max(5.min(last.max(0)))),
        ];
        for &range in &ranges {
            for array_size in [
                longest,
                longest + 1,
                longest + 2,
                3 * longest,
                3 * longest + 1,
                4096 + longest,
                32768 + longest,
                longest * 1024,
            ] {
                let want = calls(&lookup, &presence, &subject, range, array_size, None);
                for &form in &forms {
                    for resume in [false, true] {
                        let got = calls_held(
                            &lookup,
                            &presence,
                            &subject,
                            range,
                            array_size,
                            Some(form),
                            resume,
                        );
                        assert_eq!(
                            got, want,
                            "case {case} range {range:?} array_size {array_size} form {form:?} \
                                 resume {resume}"
                        );
                    }
                }
                compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                stops += want.len().saturating_sub(1);
            }
        }
    }
    // The cases above must reach the capacity stop (each one resumed in its block with
    // `resume`), many hits and cells of the overflow.
    assert!(compared > 100_000, "{compared} pairs compared");
    assert!(stops > 100, "{stops} capacity stops");
    assert_eq!(
        RESUMED.with(|r| r.get()),
        stops * forms.len(),
        "every stop resumed in its block"
    );
    assert!(big_cells > 0, "no cell past {AA_HITS_PER_CELL} hits");
    assert!(past_run > 0, "no chain past a run of {X_CHAIN_RUN} hits");
    assert!(
        near_end > 0,
        "no chain too near the overflow array's end for a run"
    );
}

// No NCBI counterpart: LOSAT_X_TBXSCANBYTE decides from `pv` (`x_table_takes`) what the kept
// bytes decide (`x_kept_presence` of `x_table_presence`), on sparse and dense tables, in every
// mode, and the ones counted from `pv` are the ones of the bytes.
#[test]
fn table_takes_equals_kept_presence() {
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    let modes = [
        XScanByte::Off,
        XScanByte::On { avx2: false },
        XScanByte::On { avx2: true },
        XScanByte::Sparse,
    ];
    let mut sides = [0usize; 2];
    for case in 0..24 {
        let queries: Vec<Vec<u8>> = (0..1 + case % 4)
            .map(|_| {
                let len = 8 + rng.below([40, 400, 4_000][case % 3]);
                (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
            })
            .collect();
        let lookup = lookup_for(&queries, [11, 13, 15][case % 3]);
        let kept = x_table_presence(&lookup);
        let (bytes, _) = kept.as_ref().expect("a tblastn-shaped table fits");
        let ones = bytes.iter().filter(|&&b| b == 1).count();
        assert_eq!(presence_ones(&lookup), ones, "case {case}");
        for mode in modes {
            let want = x_kept_presence(mode, &kept).is_some();
            assert_eq!(x_table_takes(mode, &lookup), want, "case {case} {mode:?}");
            if mode == XScanByte::Sparse {
                sides[usize::from(want)] += 1;
            }
        }
        assert_eq!(*x_build_presence(&lookup), **bytes, "case {case}");
    }
    assert!(
        sides[0] > 0 && sides[1] > 0,
        "both sides of the sparse rule: {sides:?}"
    );
}

// No NCBI counterpart: the index made from the letters is the rolled one
// (`ComputeTableIndexIncremental` from `ComputeTableIndex`) for every byte value, and the
// presence bytes are the bits of `pv`.
#[test]
fn indexes_and_presence_equal_the_rolled_index_and_pv() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let letters: Vec<u8> = (0..X_BLOCK + 2).map(|_| rng.next() as u8).collect();
    // Test protein: the first residues of E. coli AsnA (RefSeq NP_418200.1).
    let lookup = lookup_for(&[b"MKTAYIAKQRQISFVKSHFSRQ".to_vec()], 11);
    let presence = build_presence(&lookup);
    assert_eq!(presence.len(), X_MASK + 1);
    for ix in 0..=X_MASK {
        let bit = lookup
            .pv
            .get(ix >> PV_ARRAY_BTS)
            .is_some_and(|w| (w >> (ix & PV_ARRAY_MASK)) & 1 == 1);
        assert_eq!(presence[ix] == 1, bit, "index {ix}");
        assert!(
            presence[ix] == 0 || ix < lookup.backbone.len(),
            "index {ix}"
        );
    }
    for words in [1, 15, 16, 17, 31, 32, 33, 100, X_BLOCK] {
        let block = &letters[..words + 2];
        let mut forms = vec![false];
        if avx2_detected() {
            forms.push(true);
        }
        for avx2 in forms {
            let mut index = [0u16; X_BLOCK];
            let mut cand = [MaybeUninit::<u16>::uninit(); X_BLOCK];
            let n = run_candidates(avx2, block, &presence, &mut index, &mut cand);
            let mut rolled = ((block[0] as usize) << 5 | block[1] as usize) & X_MASK;
            let mut want = Vec::new();
            for k in 0..words {
                rolled = ((rolled << 5) | block[k + 2] as usize) & X_MASK;
                assert_eq!(
                    index[k] as usize, rolled,
                    "word {k} of {words}, avx2 {avx2}"
                );
                if presence[rolled] == 1 {
                    want.push(k as u16);
                }
            }
            // SAFETY: the first `n` entries were written.
            let got: Vec<u16> = cand[..n]
                .iter()
                .map(|c| unsafe { c.assume_init() })
                .collect();
            assert_eq!(got, want, "{words} words, avx2 {avx2}");
        }
    }
}

/// `candidates_avx2` when `avx2` (detected by the caller), else `candidates`.
fn run_candidates(
    avx2: bool,
    letters: &[u8],
    presence: &[u8; X_MASK + 1],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: AVX2 was detected.
        return unsafe { candidates_avx2(letters, presence, index, cand) };
    }
    assert!(!avx2, "AVX2 off x86-64");
    candidates(letters, presence, index, cand)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115 (the copy of an
// overflow chain).
// Not a port of an NCBI unit case: the chain that ends the overflow array, of a length that is
// not a multiple of four, is copied exactly (the AVX2 form then writes its last pairs one by
// one instead of reading past the array), with offset arrays from the longest chain up (as
// `x_scan_unit`'s: a smaller one would never hold the longest cell, and NCBI's loop would stop
// at that word on every call).
#[test]
fn the_last_overflow_chain_is_copied_exactly() {
    let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
    let mut tails = 0usize;
    for case in 0..64 {
        let queries: Vec<Vec<u8>> = (0..1 + case % 5)
            .map(|_| {
                let len = 30 + rng.below(600);
                (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
            })
            .collect();
        let lookup = lookup_for(&queries, [11, 12, 13][case % 3]);
        let Some((ix, cell)) = lookup
            .backbone
            .iter()
            .enumerate()
            .filter(|(_, cell)| cell.num_used as usize > AA_HITS_PER_CELL)
            .max_by_key(|(_, cell)| cell.cursor())
        else {
            continue;
        };
        let hits = cell.num_used as usize;
        let cursor = cell.cursor();
        if cursor + hits != lookup.overflow.len() || hits.is_multiple_of(4) {
            continue;
        }
        tails += 1;
        let presence = build_presence(&lookup);
        // The chain's word three times, then the sentinel.
        let word = [(ix >> 10) as u8, ((ix >> 5) & 31) as u8, (ix & 31) as u8];
        let mut subject: Vec<u8> = word.repeat(3);
        subject.push(0);
        let longest = lookup.longest_chain.max(1) as usize;
        for array_size in [
            longest,
            longest + 1,
            longest + hits,
            3 * longest,
            4096 + longest,
            32768 + longest,
            longest * 1024,
        ] {
            let range = (0, subject.len() as i32 - 4);
            let want = calls(&lookup, &presence, &subject, range, array_size, None);
            assert!(
                want[0].0.len() >= hits,
                "case {case}: the chain was not copied"
            );
            let mut forms = vec![(false, false)];
            if avx2_detected() {
                forms.extend([(true, false), (true, true)]);
            }
            for form in forms {
                for resume in [false, true] {
                    let got = calls_held(
                        &lookup,
                        &presence,
                        &subject,
                        range,
                        array_size,
                        Some(form),
                        resume,
                    );
                    assert_eq!(
                        got, want,
                        "case {case} array_size {array_size} form {form:?} resume {resume}"
                    );
                }
            }
        }
    }
    assert!(tails >= 4, "{tails} tables whose last chain ends the array");
}

/// Where `copy_block` leaves the cells of up to `AA_HITS_PER_CELL` hits for a candidate of
/// `cand[first..n]` (the room of `room_for_rest`, computed at the start and again after every
/// larger cell or tested candidate): `[ran out, larger cell fit and room held, larger cell fit
/// and room lost, larger cell did not fit]`. A statement of the rule, made from the
/// candidates' hit counts, to show which exits a case reaches.
fn copy_exits(hits: &[i32], first: usize, totalhits: i32, array_size: i32, exits: &mut [usize; 4]) {
    let left = |t: i32, i: usize| t as usize + AA_HITS_PER_CELL * (hits.len() - i);
    let (mut i, mut t) = (first, totalhits);
    while i < hits.len() {
        if left(t, i) <= array_size as usize {
            while i < hits.len() && hits[i] as usize <= AA_HITS_PER_CELL {
                t += hits[i];
                i += 1;
            }
            if i == hits.len() {
                exits[0] += 1;
                return;
            }
            if hits[i] > array_size - t {
                exits[3] += 1;
                return;
            }
            t += hits[i];
            i += 1;
            if i < hits.len() {
                exits[if left(t, i) <= array_size as usize {
                    1
                } else {
                    2
                }] += 1;
            }
        } else {
            if hits[i] > array_size - t {
                return;
            }
            t += hits[i];
            i += 1;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:91-123 (the statements
// compared against, cited at `copy_block`). Not a port of an NCBI unit case: the copy stage of
// every form, called on hand-made tables and blocks (cells of 1-3 hits, chains of 4-40 hits,
// the last chain ending the overflow array, blocks of 1-256 words, offset arrays that stop the
// copy, resumed blocks), returns NCBI's total and stopping word and writes NCBI's pairs, and
// nothing before `totalhits` or past the slack; the cases reach every exit of the loop of the
// small cells (`copy_exits`), chains past a run and chains too near the overflow array's end.
#[test]
fn the_copy_stage_equals_ncbis_statements_at_every_exit() {
    let mut rng = XorShift(0x4f1b_bcdc_bfa5_3e0b);
    let forms: Vec<(bool, bool)> = if avx2_detected() {
        vec![(false, false), (true, false), (true, true)]
    } else {
        vec![(false, false)]
    };
    let mut exits = [0usize; 4];
    let (mut stops, mut chains, mut past_run, mut near_end, mut full_blocks) = (0, 0, 0, 0, 0);
    for case in 0..3000 {
        // A table of 48 cells: 1-3 hits, or a chain of 4-40 hits in the overflow array (the
        // last one ends it in every other case).
        let small_tenths = [9, 7, 5, 2][case % 4];
        let mut backbone = vec![BackboneCell::default(); 48];
        let mut overflow: Vec<i32> = Vec::new();
        for cell in &mut backbone {
            if rng.below(10) < small_tenths {
                cell.num_used = 1 + rng.below(AA_HITS_PER_CELL) as i32;
                for e in &mut cell.entries {
                    *e = rng.below(100_000) as i32;
                }
            } else {
                cell.num_used = 4 + rng.below(37) as i32;
                cell.set_cursor(overflow.len());
                overflow.extend((0..cell.num_used).map(|_| rng.below(100_000) as i32));
            }
        }
        if case % 2 == 0 {
            overflow.extend(std::iter::repeat_n(-1, rng.below(8)));
        }
        let table = Table {
            backbone: backbone.as_ptr(),
            overflow: overflow.as_ptr(),
            overflow_len: overflow.len(),
        };
        // A block: the cell of each word, and the words that are candidates, ascending.
        let words = if case % 5 == 0 {
            X_BLOCK
        } else {
            1 + rng.below(X_BLOCK)
        };
        full_blocks += usize::from(words == X_BLOCK);
        let mut held = XScanHeld::new();
        for k in 0..words {
            held.index[k] = rng.below(backbone.len()) as u16;
        }
        let density = 1 + rng.below(100);
        let mut n = 0;
        for k in 0..words {
            if rng.below(100) < density || k + 1 == words {
                held.cand[n] = MaybeUninit::new(k as u16);
                n += 1;
            }
        }
        held.n = n;
        let cand: Vec<usize> = held.cand[..n]
            .iter()
            .map(|c| unsafe { c.assume_init() } as usize)
            .collect();
        let cells: Vec<&BackboneCell> = cand
            .iter()
            .map(|&k| &backbone[held.index[k] as usize])
            .collect();
        let hits: Vec<i32> = cells.iter().map(|cell| cell.num_used).collect();
        let longest = *hits.iter().max().unwrap() as usize;
        // Offset arrays from none to more than the block needs; in every fifth case room for
        // three hits of each of the last few candidates and a few more (a larger cell may then
        // not fit).
        let mut first = if case % 3 == 0 { rng.below(n) } else { 0 };
        let (array_size, totalhits) = if case % 5 == 1 {
            first = n - 1 - rng.below(n.min(8));
            let totalhits = rng.below(64);
            let room = AA_HITS_PER_CELL * (n - first) + rng.below(24);
            ((totalhits + room) as i32, totalhits as i32)
        } else {
            let array_size = match rng.below(4) {
                0 => rng.below(4),
                1 => rng.below(64),
                2 => 64 + rng.below(1024),
                _ => 3 * X_BLOCK + longest + rng.below(4096),
            };
            (array_size as i32, rng.below(array_size + 1) as i32)
        };
        let block_first = rng.below(1 << 20);
        // NCBI's statements (c++/src/algo/blast/core/blast_aascan.c:91-123) from candidate `first`.
        let (mut want, mut total, mut stop) = (Vec::new(), totalhits, None);
        for (i, cell) in cells.iter().enumerate().skip(first) {
            let s_off = (block_first + cand[i]) as u32;
            let numhits = cell.num_used;
            if numhits > array_size - total {
                stop = Some(s_off);
                break;
            }
            for j in 0..numhits as usize {
                let q_off = if numhits as usize <= AA_HITS_PER_CELL {
                    cell.entries[j]
                } else {
                    overflow[cell.cursor() + j]
                };
                want.push((q_off as u32, s_off));
            }
            if numhits as usize > AA_HITS_PER_CELL {
                chains += 1;
                past_run += usize::from(numhits as usize > X_CHAIN_RUN);
                near_end += usize::from(cell.cursor() + X_CHAIN_RUN > overflow.len());
            }
            total += numhits;
        }
        stops += usize::from(stop.is_some());
        copy_exits(&hits, first, totalhits, array_size, &mut exits);
        let before = BlastOffsetPair {
            q_off: 0x0bad_cafe,
            s_off: 0x0dd_ba11,
        };
        let guard = BlastOffsetPair {
            q_off: 0xdead_beef,
            s_off: 0xfeed_f00d,
        };
        let end = array_size as usize + X_SCAN_BYTE_SLACK;
        for &(avx2, run) in &forms {
            let mut pairs = vec![before; end];
            pairs.resize(end + 2 * X_SCAN_BYTE_SLACK, guard);
            // SAFETY: every candidate's cell is in `backbone`; `pairs` has `array_size +
            // X_SCAN_BYTE_SLACK` entries; AVX2 only where it was detected.
            let got = unsafe {
                copy_any(
                    avx2,
                    run,
                    &table,
                    pairs.as_mut_ptr(),
                    &held,
                    first,
                    block_first,
                    totalhits,
                    array_size,
                )
            };
            let what = format!("case {case} form {:?}", (avx2, run));
            assert_eq!(got, (total, stop), "{what}");
            let written: Vec<(u32, u32)> = pairs[totalhits as usize..total as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off))
                .collect();
            assert_eq!(written, want, "{what}");
            let is =
                |p: &BlastOffsetPair, q: BlastOffsetPair| (p.q_off, p.s_off) == (q.q_off, q.s_off);
            assert!(
                pairs[..totalhits as usize].iter().all(|p| is(p, before)),
                "{what}"
            );
            assert!(pairs[end..].iter().all(|p| is(p, guard)), "{what}");
        }
    }
    // Every exit of the loop of the small cells, stops, and chains of every kind.
    assert!(exits.iter().all(|&e| e >= 100), "exits {exits:?}");
    assert!(stops >= 500, "{stops} stops");
    assert!(chains >= 10_000, "{chains} chains");
    assert!(past_run >= 2_000, "{past_run} chains past a run");
    assert!(near_end >= 300, "{near_end} chains near the overflow end");
    assert!(full_blocks >= 300, "{full_blocks} blocks of X_BLOCK words");
}

// No NCBI counterpart: the kept presence bytes follow the table (a new table, also one at a
// freed table's address, gets its own).
#[test]
fn presence_bytes_are_rebuilt_for_another_table() {
    let one = Arc::new(lookup_for(&[b"MKTAYIAKQRQISFVKSHFSRQ".to_vec()], 11));
    let two = Arc::new(lookup_for(&[b"WWWWCCCCHHHHMMMM".to_vec()], 11));
    assert_eq!(*presence_bytes(&one).0, *build_presence(&one));
    assert_eq!(*presence_bytes(&two).0, *build_presence(&two));
    assert_eq!(*presence_bytes(&one).0, *build_presence(&one));
    drop(one);
    let three = Arc::new(lookup_for(&[b"PPPPGGGGAAAA".to_vec()], 11));
    assert_eq!(*presence_bytes(&three).0, *build_presence(&three));
}

// No NCBI counterpart: a table is sparse when at most three tenths of its backbone cells have
// hits; the table of one random protein of 400 letters is (0.21 at threshold 13), that of 64 is
// not (0.54), and `XScanByte::Sparse` scans only the sparse one this way.
#[test]
fn one_protein_makes_a_sparse_table_and_many_do_not() {
    let mut rng = XorShift(0x0bad_5eed_1234_5678);
    let mut protein = || -> Vec<u8> { (0..400).map(|_| AMINO[rng.below(AMINO.len())]).collect() };
    let one = Arc::new(lookup_for(&[protein()], 13));
    let many = Arc::new(lookup_for(
        &(0..64).map(|_| protein()).collect::<Vec<_>>(),
        11,
    ));
    for (table, want) in [(&one, true), (&many, false)] {
        let bytes = build_presence(table);
        let ones = bytes.iter().filter(|&&b| b == 1).count();
        let cells = table.backbone.len();
        assert_eq!(sparse(&bytes, table), ones * 10 <= cells * X_SPARSE_TENTHS);
        assert_eq!(presence_bytes(table).1, want, "{ones} of {cells} cells");
        assert_eq!(x_presence(XScanByte::Sparse, table).is_some(), want);
        let mode = XScanByte::On { avx2: false };
        assert!(x_presence(mode, table).is_some(), "{mode:?}");
        assert!(x_presence(XScanByte::Off, table).is_none());
    }
}

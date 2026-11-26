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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203 (cited at
// `seqalign_identities`): each row keeps its HSP's Seq-align `num_ident`, the engine's
// identities that the HSP holds before `set_displayed_identities` replaces them, while
// `hit.num_ident` is the count of the displayed rows (here 10 of 10 equal residues).
#[test]
fn rows_keep_the_seqalign_identities_of_their_hsps() {
    let code = GeneticCode::from_id(1);
    let letters = b"ATGGCTAAACTGGAAGGTCCATTTAGCCGA";
    let queries = [FastaRecord::new("q", b"q", letters)];
    let subjects = [FastaRecord::new("s", b"s", letters)];
    let mut hits = [hsp(0, 0, 50, 1e-5, 0), hsp(0, 0, 40, 1e-4, 0)];
    (hits[0].hit.num_ident, hits[1].hit.num_ident) = (7, 9);
    let identities = seqalign_identities(&hits);
    assert_eq!(identities, [7, 9]);
    let placements = Placements::default();
    let rows = pairwise_hits(
        &hits,
        Some(&identities),
        &queries,
        &subjects,
        SubjectLetters::Records(&subjects),
        &[TblastxQueryStats::default()],
        &code,
        &code,
        &placements,
        &placements,
        None,
    );
    let counts: Vec<(usize, usize)> = rows
        .iter()
        .map(|row| (row.hit.num_ident, row.seqalign_num_ident))
        .collect();
    assert_eq!(counts, [(10, 7), (10, 9)]);
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:94-112: the display reads the
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

// The display table (graduated `LOSAT_X_DISPTABLE`) holds the scan's answers: each byte's code
// on both strands, each codon of 4-bit codes of each genetic code, and so each residue of a row.
#[test]
fn display_table_gives_the_residues_of_the_scan() {
    let reversed =
        |mask: u8| ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3);
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
            assert_eq!(residue, display_codon(masks, &code));
        }
        let display = DisplayCode { table };
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

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:191-207 merges D/N to B, E/Q
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1003-1021 and c++/src/objtools/align_format/showalign.cpp:2131-2150:
// equal letters are identities (X/X too); B/D is positive, J/I is not (J is not in the
// display matrix).
#[test]
fn row_counts_use_the_displayed_letters() {
    assert_eq!(row_counts(Blosum62Score, b"XAB*J", b"XADKI"), (2, 3));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2717 and 746-792
// Review F1 of M4b PR 2: with `-db_hard_mask` the engine's counts (the Seq-aligns'
// `num_ident`, `num_positives`) read a database sequence as the re-evaluation read it, `B` over
// the masked bases, while the rows read its stored letters. Frame +2's first codon (bases 1-3,
// `GCT`, A) has its first base masked: `BCT` is X, which scores -1 against every letter (the
// X row of NCBI's matrices), so NCBI's re-evaluation trims it from an HSP's end (rows
// db.tmask.tx.hard.te.11.*: no HSP keeps such a residue). The HSP here is placed over it to
// show which letters each count reads: X for the engine's counts (neither an identity nor a
// positive), the stored A for the rows' (an identity).
#[test]
fn engine_counts_read_the_hard_masked_letters() {
    use crate::algorithm::common::db_mask_search::mask_nucleotide_store;
    use crate::algorithm::tblastx::database::ReportedStore;
    let code = GeneticCode::from_id(1);
    let mut subject = b"T".to_vec();
    subject.extend(b"GCT".repeat(10));
    subject.extend(b"AA");
    let query = b"GCT".repeat(10);
    let shown = Arc::new(PackedNucleotide::from_subject_letters(&subject).expect("letters"));
    let engine = Arc::new(mask_nucleotide_store(&shown, &[(0, 2)]));
    let stores = [ReportedStore { shown, engine }];
    let lengths = [subject.len() as u32];
    let subjects = SubjectLetters::Database {
        lengths: &lengths,
        stores: &stores,
    };
    let records = [FastaRecord::new("q", b"q", &query)];
    let mut hsp = hsp(0, 0, 10, 1.0, 0);
    hsp.subject_frame = 2;
    hsp.hit.num_ident = 9;
    let counts = engine_counts(
        std::slice::from_ref(&hsp),
        &records,
        subjects,
        &code,
        &code,
        Blosum62Score,
    );
    assert_eq!(counts, [(9, 9)]);
    let mut window = Vec::new();
    assert_eq!(
        subjects.engine_identities(&hsp, &query, &code, &code, &mut window),
        9
    );
    assert_eq!(
        subjects.counts(&hsp, &query, &code, &code, &mut window, Blosum62Score),
        (10, 10)
    );
    let (_, subject_row) = subjects.rows(&hsp, &query, &code, &code, &mut window);
    assert_eq!(subject_row, b"AAAAAAAAAA");
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
                        let q = display_residue(&query, query_frame, query_offset + i, &query_code);
                        let s =
                            display_residue(&subject, subject_frame, subject_offset + i, &db_code);
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
    let text = row_string(row);
    assert_eq!(text.as_ptr(), pointer);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:925-931: a minus-frame run starts one
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
    // One `HspMasks` for the three rows, as the rows' builders keep one.
    let mut row_masks = crate::report::query_anchored::HspMasks::default();
    let mut plus = b"AAAAAAAAAA".to_vec();
    assert!(query_row_masks(&mut plus, 1, 0, 29, &masks, &mut row_masks));
    assert_eq!(&plus, b"AAaaaAAAAA");
    // Frame -1 residue i covers 0-based nucleotides 27-3i..29-3i: the mask 18-23 covers
    // residues 2 and 3 only.
    let mut minus = b"AAAAAAAAAA".to_vec();
    assert!(query_row_masks(
        &mut minus,
        -1,
        0,
        29,
        &masks,
        &mut row_masks
    ));
    assert_eq!(&minus, b"AAaaAAAAAA");
    // A row of another frame shows none.
    let mut other = b"AAAAAAAAA".to_vec();
    assert!(!query_row_masks(
        &mut other,
        2,
        1,
        27,
        &masks,
        &mut row_masks
    ));
    assert_eq!(row_masks.row(0), None);
    assert_eq!(&other, b"AAAAAAAAA");
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4305-4326 (cited at
// `shown_query_masks`): `x_FillLocList` hands the masks to `x_OutputSeq` in the order of the
// shown list, which a minus frame's conversion makes descending (`protein_to_dna_masks`); the
// row's scan reads them sorted by start (`FrameMasks`) and gives them back in the list's order.
#[test]
fn row_masks_keep_the_order_of_the_shown_masks() {
    let masks = [FrameMasks::new(-1, vec![(18, 23), (6, 11)])];
    let mut row = b"AAAAAAAAAA".to_vec();
    let mut row_masks = crate::report::query_anchored::HspMasks::default();
    assert!(
        query_row_masks(&mut row, -1, 0, 29, &masks, &mut row_masks),
        "masks meet the row"
    );
    assert_eq!(row_masks.row(0), Some((&[(2, 3, -1), (6, 7, -1)][..], -1)));
    assert_eq!(row_masks.row(1), None);
    assert_eq!(&row, b"AAaaAAaaAA");
}

// The scan of all the ranges of a frame, which every row made before `LOSAT_X_TBXMASKSEEK`
// graduated: the reference of `mask_seek_lowercases_the_letters_of_the_full_scan`.
fn lowercase_query_row_full_scan(
    row: &mut [u8],
    frame: i8,
    low: i32,
    high: i32,
    masks: &[FrameMasks],
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
        for &(from, to, _) in frame_masks.seek.ranges() {
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

#[test]
fn mask_seek_lowercases_the_letters_of_the_full_scan() {
    // Ranges in any order, overlapping, empty (start after end) or outside the row, on rows
    // of both frame signs, of a residue for each codon of their range (TBLASTX's rows).
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = |bound: i32| -> i32 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % bound as u64) as i32
    };
    // One `HspMasks` for every row, as the rows' builders keep one.
    let mut row_masks = crate::report::query_anchored::HspMasks::default();
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
        let high = low + 3 * (1 + next((60 - low) / 3)) - 1;
        for frame in [1, -1] {
            let mut full = vec![b'A'; ((high - low + 1) / 3) as usize];
            let mut seek = full.clone();
            lowercase_query_row_full_scan(&mut full, frame, low, high, &masks);
            let met = query_row_masks(&mut seek, frame, low, high, &masks, &mut row_masks);
            assert_eq!(seek, full, "frame {frame} low {low} high {high}");
            // A lower-case letter has its mask (an empty range meets a row and shows none).
            assert!(met || !full.iter().any(u8::is_ascii_lowercase));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:993
// `BlastSeqLocCombine` (with `link_value` 0) joins
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
    }
    // The custom tabular reports make their own rows.
    for outfmt in ["6 qseqid", "10", "20 std"] {
        let format = output_format(outfmt);
        assert!(!format.reads_displayed_counts(), "{outfmt}");
        assert!(format.is_tabular(), "{outfmt}");
    }
    // DBID item T: the formats whose batches' prefetch writes the taxonomy database's
    // messages (not the archive; the Seq-annots only with `PRE_FETCH_SEQS_LIMIT=2147483647`, DBID
    // P3, which the test process may have set; the XML2 and JSON documents too, line DBID
    // taxpf-1).
    let seq_annots = |number| {
        crate::algorithm::common::database::seq_annots_prefetch(Some(
            crate::algorithm::common::database::PrefetchFormat::of(
                number,
                crate::report::archive::options::RequestKind::Tblastx,
            ),
        ))
    };
    for (outfmt, prefetches) in [
        ("0", true),
        ("5", true),
        ("6", true),
        ("7 qseqid", true),
        ("10", true),
        ("20 qseqid", true),
        ("8", seq_annots(8)),
        ("9", seq_annots(9)),
        ("11", false),
        ("12", seq_annots(12)),
        ("13", true),
        ("15", true),
        ("16", true),
    ] {
        assert_eq!(
            output_format(outfmt).prefetches_subjects(),
            prefetches,
            "{outfmt}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:825-842 (`Map`): with `-query_loc` a mask
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3299 (Blast_HitListUpdate) and
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
        final_hit_order(hits.clone(), size, &HitFilters::default())
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
        write_outfmt7_query_header(&mut out, b"q1 title", b"db", hits).unwrap();
        String::from_utf8(out).unwrap()
    };
    assert_eq!(
        header(None),
        concat!(
            "# TBLASTX (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: q1 title\n# Database: db\n"
        )
    );
    assert!(header(Some(0)).ends_with("# Database: db\n# 0 hits found\n"));
    assert!(header(Some(2)).contains("# Fields: query acc.ver, subject acc.ver"));
    assert!(header(Some(2)).ends_with("# 2 hits found\n"));
}

#[test]
fn output_formats_are_read_after_trimming() {
    assert_eq!(
        output_format(" 0"),
        TblastxOutputFormat::Pairwise(crate::report::pairwise::PairwiseView::Pairwise)
    );
    assert_eq!(output_format("6"), TblastxOutputFormat::Tabular);
    assert_eq!(
        output_format("7 "),
        TblastxOutputFormat::TabularWithComments
    );
}

// Outfmt 1 to 4 are the pairwise report's query-anchored views (line RG T3c); a value that
// TBLASTX does not write is never read as a table.
#[test]
fn the_query_anchored_formats_are_pairwise_views_and_others_unsupported() {
    use crate::report::pairwise::PairwiseView;
    for (outfmt, identity, anchored) in [
        ("1", true, true),
        ("2", false, true),
        ("3", true, false),
        (" 4 qseqid", false, false),
    ] {
        assert_eq!(
            output_format(outfmt),
            TblastxOutputFormat::Pairwise(PairwiseView::QueryAnchored { identity, anchored }),
            "{outfmt}"
        );
    }
    assert_eq!(output_format("17"), TblastxOutputFormat::Unsupported);
    assert_eq!(output_format("99"), TblastxOutputFormat::Unsupported);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786: every tabular kind goes to `CBlastTabularInfo`; LOSAT's
// default 6 and 7 (the `std` fields with a tab) keep `write_tabular` and the report pool.
#[test]
fn only_a_field_list_a_delimiter_and_10_and_20_take_the_shared_writer() {
    use crate::report::tabular::TabularKind;
    for (outfmt, format) in [
        ("6", TblastxOutputFormat::Tabular),
        ("6 std", TblastxOutputFormat::Tabular),
        (" 7 std ", TblastxOutputFormat::TabularWithComments),
        ("6 -qseqid", TblastxOutputFormat::Tabular),
        ("6 delim=", TblastxOutputFormat::Tabular),
        (
            "0 qseqid",
            TblastxOutputFormat::Pairwise(crate::report::pairwise::PairwiseView::Pairwise),
        ),
        ("11 qseqid", TblastxOutputFormat::Archive),
    ] {
        assert_eq!(output_format(outfmt), format, "{outfmt}");
    }
    for (outfmt, kind, delimiter) in [
        ("6 qseqid", TabularKind::Tabular, "\t"),
        ("6 delim=, std", TabularKind::Tabular, ","),
        ("7 qaccver foo", TabularKind::Commented, "\t"),
        ("10", TabularKind::Csv, ","),
        ("10 delim=; std", TabularKind::Csv, ";"),
        ("20 qseqid", TabularKind::CsvWithHeader, ","),
    ] {
        match output_format(outfmt) {
            TblastxOutputFormat::CustomTabular(spec) => {
                assert_eq!(
                    (spec.kind, spec.delimiter.as_slice()),
                    (kind, delimiter.as_bytes())
                )
            }
            other => panic!("{outfmt}: {other:?}"),
        }
    }
}

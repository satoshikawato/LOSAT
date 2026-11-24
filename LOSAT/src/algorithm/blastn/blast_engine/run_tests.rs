use super::*;

// No NCBI counterpart (a check of the scans' reads): every megablast scan
// (`scan_subject_kmers_range_mb_*`, which read their three or four bytes with one test of the
// length, `x_be_bytes`, or none in the rounds `x_whole_rounds` counts) gives, for random packed subjects (with and without a byte past the
// last base) and ranges from every offset modulo 4, the words at `start`, `start + scan_step`,
// ... up to `end`, each the word `packed_kmer_at` reads there; the scan of any stride also for
// a masked subject (no aligned read), and its aligned read from aligned starts. The stride-2
// scans of 9 and 10 bases take even starts only (NCBI's s_MBScanSubject_9_2 and _10_2 enter at
// an offset of 0 or 2 modulo 4, blast_nascan.c:1733, 1863).
#[test]
fn the_megablast_scans_give_the_word_at_each_step_of_their_range() {
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    macro_rules! got {
        ($scan:expr) => {{
            let mut got = Vec::new();
            $scan(&mut |pos: usize, kmer: u64| {
                got.push((pos, kmer));
                true
            });
            got
        }};
    }
    let mut checked = 0usize;
    for round in 0..4000 {
        let len = 12 + next() % 150;
        let packed: Vec<u8> = (0..len.div_ceil(COMPRESSION_RATIO) + round % 2)
            .map(|_| next() as u8)
            .collect();
        let p = packed.as_slice();
        // `m`: the starts a scan takes are multiples of `m` (NCBI's stride-2 scans of 9 and 10
        // bases start on an even offset, its aligned read on a byte).
        let range = |lut: usize, m: usize, next: &mut dyn FnMut() -> usize| {
            let last = len - lut;
            let mut start = next() % (last + 1);
            start -= start % m;
            (start, start + next() % (last - start + 1))
        };
        let want = |lut: usize, step: usize, (start, end): (usize, usize)| -> Vec<(usize, u64)> {
            (start..=end)
                .step_by(step)
                .map(|pos| (pos, packed_kmer_at(p, pos, lut)))
                .collect()
        };
        type Scan = fn(&[u8], usize, usize, usize, usize, &mut dyn FnMut(usize, u64) -> bool);
        let fixed: [(usize, usize, Scan); 9] = [
            (9, 1, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_9_1(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (9, 2, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_9_2(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (10, 1, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_10_1(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (10, 2, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_10_2(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (10, 3, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_10_3(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (11, 1, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_11_1mod4(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (11, 2, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_11_2mod4(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (11, 3, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_11_3mod4(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
            (11, 1, |p, l, s, a, b, f| {
                scan_subject_kmers_range_mb_11_1mod4(p, l, s, a, b, &mut |x, y| f(x, y))
            }),
        ];
        for (i, &(lut, base_step, scan)) in fixed.iter().enumerate() {
            // The 11-base scans take any stride of their class modulo 4 (the last row: 17).
            let step = if lut == 11 {
                base_step + 4 * (next() % 5)
            } else {
                base_step
            };
            let step = if i == 8 { 17 } else { step };
            let m = if base_step == 2 && lut < 11 { 2 } else { 1 };
            // A random range, and the whole subject from each offset modulo 4 that the scan
            // takes (the rounds of `x_whole_rounds` up to the last byte, then the tail).
            let r = range(lut, m, &mut next);
            let whole = (0..COMPRESSION_RATIO).step_by(m).map(|s| (s, len - lut));
            for r in std::iter::once(r).chain(whole.filter(|&(s, last)| s <= last)) {
                let got =
                    got!(|f: &mut dyn FnMut(usize, u64) -> bool| scan(p, len, step, r.0, r.1, f));
                assert_eq!(
                    got,
                    want(lut, step, r),
                    "lut {lut} step {step} range {r:?} len {len}"
                );
                checked += got.len();
            }
        }
        for lut in 9..=12 {
            for step in [1 + next() % 24, 4 * (1 + next() % 6)] {
                for masked in [false, true] {
                    let aligned = step % COMPRESSION_RATIO == 0 && !masked;
                    let r = range(lut, if aligned { COMPRESSION_RATIO } else { 1 }, &mut next);
                    let got = got!(|f: &mut dyn FnMut(usize, u64) -> bool| {
                        scan_subject_kmers_range_mb_any(
                            p,
                            len,
                            lut,
                            step,
                            masked,
                            r.0,
                            r.1,
                            &mut |x, y| f(x, y),
                        )
                    });
                    assert_eq!(
                        got,
                        want(lut, step, r),
                        "any: lut {lut} step {step} masked {masked} range {r:?} len {len}"
                    );
                    checked += got.len();
                }
            }
        }
    }
    assert!(checked > 400_000, "{checked} words");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1668 (cited at
// `scan_subject_kmers_with_ranges`); SURVEY §10 of M4a-2, review 1 T3 of M4a-2a. With
// `-db_soft_mask` every database sequence is a masked subject, those without masks too
// (decision D3: `subject_masked` is the option, seqsrc_seqdb.cpp:399-405): its one range
// `[0, length)` is scanned from `word_length - lut_word_length` with the scanner of any
// stride, where an unmasked subject is scanned from 0 with the scanner of the lookup width
// and stride (the aligned one of a stride of four bases). Both visit every `scan_step`-th
// offset to `length - lut_word_length`.
#[test]
fn a_masked_subject_without_masks_starts_its_scan_after_the_word_s_head() {
    let letters: Vec<u8> = (0..403u32)
        .map(|i| b"ACGT"[((i * 7 + i / 3 + i / 11) % 4) as usize])
        .collect();
    let packed = crate::core::blast_encoding::encode_subject_ncbi2na_packed(&letters);
    let length = letters.len();
    let whole = [(0i32, length as i32)];
    // (word_length, lut_word_length, scan_step): megablast with lookup widths 11 and 12,
    // blastn (8-base lookup, aligned scan), blastn-short (7 and 7).
    for (word_length, lut_word_length, scan_step) in
        [(28, 11, 18), (28, 12, 17), (11, 8, 4), (7, 7, 1)]
    {
        let starts = |masked: bool| {
            let mut starts = Vec::new();
            scan_subject_kmers_with_ranges(
                &packed,
                length,
                word_length,
                lut_word_length,
                scan_step,
                &whole,
                masked,
                select_mb_scan_kind(lut_word_length, scan_step, masked),
                |start, _, _| starts.push(start),
            );
            starts
        };
        let last = length - lut_word_length;
        let case = (word_length, lut_word_length, scan_step);
        assert_eq!(
            starts(false),
            (0..=last).step_by(scan_step).collect::<Vec<_>>(),
            "unmasked {case:?}"
        );
        assert_eq!(
            starts(true),
            (word_length - lut_word_length..=last)
                .step_by(scan_step)
                .collect::<Vec<_>>(),
            "masked {case:?}"
        );
    }
}

// `-soft_masking false`: each context's masks become Ns in the buffer (the minus context's
// are its own, already reversed), the returned copy keeps the letters, and the contexts no
// longer hold their masks; without any mask nothing is copied and nothing changes.
#[test]
fn the_masks_are_written_into_the_searched_query() {
    let context = |offset: i32, len: usize, masks: Vec<MaskedInterval>| QueryContext {
        query_idx: 0,
        frame: if offset == 0 { 1 } else { -1 },
        query_offset: offset,
        len,
        masks,
    };
    let plus = vec![MaskedInterval::new(1, 3)];
    let minus = reverse_mask_intervals(&plus, 5);
    let buffer = vec![15, 0, 1, 2, 3, 0, 15, 3, 2, 1, 0, 3, 15];
    let mut searched = buffer.clone();
    let mut contexts = vec![context(0, 5, plus.clone()), context(6, 5, minus.clone())];
    let nomask = mask_query_letters(&mut searched, &mut contexts, &[plus, minus]);
    assert_eq!(nomask.as_deref(), Some(buffer.as_slice()));
    assert_eq!(searched, [15, 0, 14, 14, 3, 0, 15, 3, 2, 14, 14, 3, 15]);
    assert!(contexts.iter().all(|context| context.masks.is_empty()));

    let mut searched = buffer.clone();
    let mut contexts = vec![context(0, 5, Vec::new()), context(6, 5, Vec::new())];
    let nomask = mask_query_letters(&mut searched, &mut contexts, &[Vec::new(), Vec::new()]);
    assert!(nomask.is_none());
    assert_eq!(searched, buffer);
}

// QueryContextIndex answers the value of `direct_map` at every offset, through whichever form
// the layout allows (LOSAT_X_CTXO1: division, or the block table, here with the
// runs for its crowded blocks in the second case: a query of 3 letters among 500 of 300).
#[test]
fn the_context_index_answers_the_table_at_every_offset() {
    let mut crowded = vec![300usize; 1000];
    crowded[200] = 3;
    crowded[201] = 3;
    crowded.extend([0, 0]);
    for (case, lengths) in [
        vec![24usize; 40],
        crowded,
        vec![0, 0, 18, 18, 30, 30, 25, 25, 19, 19],
        vec![5000, 5000],
    ]
    .into_iter()
    .enumerate()
    {
        let mut contexts = Vec::new();
        let mut offset = 0usize;
        for (i, &len) in lengths.iter().enumerate() {
            contexts.push(QueryContext {
                query_idx: (i / 2) as u32,
                frame: if i % 2 == 0 { 1 } else { -1 },
                query_offset: offset as i32,
                len,
                masks: Vec::new(),
            });
            offset += if len == 0 { 0 } else { len + 1 };
        }
        let total = contexts
            .iter()
            .map(|c| c.query_offset as usize + c.len)
            .max()
            .unwrap();
        let mut table = vec![0usize; total];
        for (idx, c) in contexts.iter().enumerate() {
            let start = c.query_offset as usize;
            table[start..start + c.len].fill(idx);
        }
        let index = QueryContextIndex::new(&contexts, false);
        if case == 1 {
            if let Some(o1) = &index.x_o1 {
                assert!(o1.needs_runs() && !index.x_ranges.is_empty());
            }
        }
        for (n, &want) in table.iter().enumerate() {
            assert_eq!(
                index.context_for_offset(n),
                want,
                "offset {n} of {lengths:?}"
            );
        }
        // past the table: NCBI's search, as before the forms
        for n in total..total + 40 {
            assert_eq!(index.context_for_offset(n), bsearch_context_info(n, &index));
        }
    }
}

/// The records of a FASTA text as NCBI's reader reads them (`config`).
fn reader_records(fasta: &str, config: ReaderConfig) -> Vec<FastaRecord> {
    crate::blastinput::fasta_reader::read_all(
        &mut FastaInputSource::from_bytes(fasta.as_bytes(), config),
        &mut |_| Ok(()),
    )
    .unwrap()
}

/// The query buffer as it was built from per-context IUPAC copies: every context encoded
/// on its own, then copied into the block and into the block with sentinels.
fn query_buffers_from_context_copies(
    queries: &[&[u8]],
    offsets: &[i32],
    concat_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut logical = vec![NUCL_SENTINEL; concat_len];
    let mut with_sentinels = vec![NUCL_SENTINEL; concat_len + 2];
    for (q, seq) in queries.iter().enumerate() {
        let contexts = [
            encode_iupac_to_blastna(seq),
            encode_iupac_to_blastna(&super::super::super::lookup::reverse_complement(seq)),
        ];
        for (encoded, &offset) in contexts.iter().zip(&offsets[2 * q..2 * q + 2]) {
            let start = offset.max(0) as usize;
            let end = start.saturating_add(encoded.len());
            if start >= end || end > logical.len() {
                continue;
            }
            logical[start..end].copy_from_slice(encoded);
            with_sentinels[start + 1..end + 1].copy_from_slice(encoded);
        }
    }
    (logical, with_sentinels)
}

// One buffer written context by context equals the two buffers of the IUPAC copies,
// with an empty query (no room) between queries and letters of both cases.
#[test]
fn kept_query_buffer_equals_a_new_one() {
    // (letters, offsets, concatenated length) of a batch of `queries`, laid out as above.
    let layout = |queries: &[&[u8]]| {
        let mut offsets = Vec::new();
        let mut at = 0usize;
        let room = |length: usize| if length == 0 { 0 } else { length + 1 };
        for seq in queries {
            offsets.push(at as i32);
            at += room(seq.len());
            offsets.push(at as i32);
            at += room(seq.len());
        }
        (offsets, at - 1)
    };
    let (long, mid) = (vec![b'G'; 3 << 20], vec![b'c'; 1 << 20]);
    let huge = vec![b'T'; 17 << 20];
    let batches: [Vec<&[u8]>; 8] = [
        vec![&b"ACGTTGCAAC"[..], &b"acgtnRYKMSWBDHVNu-acgt"[..]],
        vec![&b"A"[..]],
        vec![&long[..], &b"TTTT"[..]],
        vec![&mid[..], &b"GA"[..]],
        vec![&b"GGGGCCCCTTTTAAAAccgg*0"[..], &b""[..], &b"CA"[..]],
        vec![&long[..]],
        vec![&huge[..], &b"G"[..]],
        vec![&b"ACGT"[..]],
    ];
    let mut keep = Vec::new();
    let mut blocks = Vec::new();
    for queries in &batches {
        let (offsets, concat_len) = layout(&queries[..]);
        fill_query_blastna_buffer(&queries[..], &offsets, concat_len, &mut keep);
        assert_eq!(
            keep,
            build_query_blastna_buffer(&queries[..], &offsets, concat_len)
        );
        blocks.push((keep.as_ptr(), keep.capacity()));
        for strand in [QueryStrand::Plus, QueryStrand::Minus] {
            let mut kept = keep.clone();
            kept.extend_from_slice(b"stale bytes past the length");
            fill_strand_query_blastna_buffer(&queries[..], &offsets, concat_len, strand, &mut kept);
            assert_eq!(
                kept,
                build_strand_query_blastna_buffer(&queries[..], &offsets, concat_len, strand)
            );
        }
    }
    if crate::utils::x_mapped::X_MAPPED_ALLOCATOR {
        // the 3 MiB batch made a mapped buffer; the shorter one of the same class, the batch
        // below 1 MiB and the next 3 MiB batch wrote the same block again (it is not mapped
        // again); the 34 MiB batch, past the ceiling, took a block of its own length, and the
        // short batch after it did not keep that block
        assert!(blocks[2].1 > 32 << 20);
        assert_eq!(blocks[3], blocks[2]);
        assert_eq!(blocks[4], blocks[2]);
        assert_eq!(blocks[5], blocks[2]);
        assert!(blocks[6].1 > blocks[2].1);
        assert!(blocks[7].1 < 1 << 20);
    }
}

#[test]
fn query_buffer_equals_the_buffers_of_the_context_copies() {
    let queries: [&[u8]; 5] = [
        b"ACGTTGCAAC",
        b"",
        b"acgtnRYKMSWBDHVNu-acgt",
        b"A",
        b"GGGGCCCCTTTTAAAAccgg*0",
    ];
    let mut offsets = Vec::new();
    let mut at = 0usize;
    let room = |length: usize| if length == 0 { 0 } else { length + 1 };
    for seq in &queries {
        offsets.push(at as i32);
        at += room(seq.len());
        offsets.push(at as i32);
        at += room(seq.len());
    }
    let concat_len = at - 1;
    let (logical, with_sentinels) =
        query_buffers_from_context_copies(&queries, &offsets, concat_len);
    let buffer = build_query_blastna_buffer(&queries[..], &offsets, concat_len);
    assert_eq!(buffer, with_sentinels);
    assert_eq!(&buffer[1..concat_len + 1], logical.as_slice());
    assert_eq!(
        (buffer[0], buffer[concat_len + 1]),
        (NUCL_SENTINEL, NUCL_SENTINEL)
    );
}

/// The contexts of `lengths` (query letters) searched on `strand`: offset and length of
/// each of the two contexts of a query, laid out as NCBI's `s_QueryInfo_SetContext` does.
fn strand_contexts(lengths: &[usize], strand: QueryStrand) -> Vec<QueryContext> {
    let (plus_on, minus_on) = strand.searched();
    let room = |length: usize| if length == 0 { 0 } else { length + 1 };
    let mut contexts = Vec::new();
    let mut at = 0usize;
    for (q, &length) in lengths.iter().enumerate() {
        for (frame, on) in [(1, plus_on), (-1, minus_on)] {
            let len = if on { length } else { 0 };
            contexts.push(QueryContext {
                query_idx: q as u32,
                frame,
                query_offset: at as i32,
                len,
                masks: Vec::new(),
            });
            at += room(len);
        }
    }
    contexts
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81,229-243
// The buffer of a search of one strand holds only the contexts of that strand, each
// followed by a sentinel, and equals the one laid out naively; queries of 10, 3, 10
// letters and an empty query.
#[test]
fn a_single_strand_query_buffer_holds_only_the_strand_it_searches() {
    let queries: [&[u8]; 5] = [b"ACGTTGCAAC", b"ACN", b"", b"GGRYCCTTAA", b"T"];
    for strand in [QueryStrand::Plus, QueryStrand::Minus] {
        let lengths: Vec<usize> = queries.iter().map(|q| q.len()).collect();
        let contexts = strand_contexts(&lengths, strand);
        let offsets: Vec<i32> = contexts.iter().map(|c| c.query_offset).collect();
        let last = contexts.last().unwrap();
        let concat_len =
            (last.query_offset as usize + last.len + if last.len == 0 { 1 } else { 2 })
                .saturating_sub(2);
        // the naive buffer: a sentinel, then each searched context and a sentinel
        let mut naive = vec![NUCL_SENTINEL];
        for letters in &queries {
            if letters.is_empty() {
                continue;
            }
            let encoded = if strand == QueryStrand::Plus {
                encode_iupac_to_blastna(letters)
            } else {
                encode_iupac_to_blastna(&super::super::super::lookup::reverse_complement(letters))
            };
            naive.extend(encoded);
            naive.push(NUCL_SENTINEL);
        }
        assert_eq!(naive.len(), concat_len + 2, "{strand:?}");
        let buffer = build_strand_query_blastna_buffer(&queries[..], &offsets, concat_len, strand);
        assert_eq!(buffer, naive, "{strand:?}");
    }
    // an empty query alone, and a single query, build a block of sentinels / its letters
    let buffer =
        build_strand_query_blastna_buffer(&[&b"AC"[..]][..], &[0, 3], 2, QueryStrand::Plus);
    assert_eq!(buffer, [NUCL_SENTINEL, 0, 1, NUCL_SENTINEL]);
    let buffer =
        build_strand_query_blastna_buffer(&[&b"AC"[..]][..], &[0, 0], 2, QueryStrand::Minus);
    assert_eq!(buffer, [NUCL_SENTINEL, 2, 3, NUCL_SENTINEL]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
// With one strand searched the context of every offset inside a valid context is that
// context, through the index (whichever form the environment chooses), the table, the
// runs, the constant-time form and `bsearch_context_info` (with `min_length` 0).
#[test]
fn a_single_strand_layout_maps_every_offset_to_its_context() {
    for lengths in [
        vec![10usize, 3, 10],
        vec![10, 3, 0, 10, 1],
        vec![0, 7, 7, 7],
        vec![5000, 17, 5000, 2],
        vec![1],
    ] {
        for strand in [QueryStrand::Plus, QueryStrand::Minus, QueryStrand::Both] {
            let contexts = strand_contexts(&lengths, strand);
            let index = QueryContextIndex::new(&contexts, strand != QueryStrand::Both);
            let total = contexts
                .iter()
                .map(|c| c.query_offset as usize + c.len)
                .max()
                .unwrap();
            let o1 = x_ctx_o1::CtxO1::build(&contexts, total);
            for (idx, ctx) in contexts.iter().enumerate() {
                for n in ctx.query_offset as usize..ctx.query_offset as usize + ctx.len {
                    let want = (idx, format!("{lengths:?} {strand:?} offset {n}"));
                    assert_eq!(index.context_for_offset(n), want.0, "{}", want.1);
                    // (With both strands NCBI's bounds of `min_length` need every query to
                    // have letters; the one-strand layouts have `min_length` 0.)
                    if strand != QueryStrand::Both || lengths.iter().all(|&len| len > 0) {
                        assert_eq!(bsearch_context_info(n, &index), want.0, "{}", want.1);
                    }
                    if !index.direct_map.is_empty() {
                        assert_eq!(index.direct_map[n] as usize, want.0, "{}", want.1);
                    }
                    if index.x_total_len > n {
                        assert_eq!(index.x_direct(n), want.0, "{}", want.1);
                    }
                    if let Some(o1) = &o1 {
                        if let Some(context) = o1.get(n) {
                            assert_eq!(context, want.0, "{}", want.1);
                        }
                    }
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:225-229
// The diagonal hash starts above 8000 on the LAST context's end: with `-strand plus` that
// is the empty minus context, one past the plus context, so a query of exactly 8000 letters
// takes the hash; with `minus` the last context is the searched one; with both strands the
// value is odd (two contexts and the gap between them).
#[test]
fn the_diagonal_hash_starts_at_the_end_of_the_last_context() {
    let hash = |lengths: &[usize], strand| {
        diag_container_query_length(&strand_contexts(lengths, strand)) > 8000
    };
    // plus
    assert!(!hash(&[7999], QueryStrand::Plus));
    assert!(hash(&[8000], QueryStrand::Plus));
    assert!(!hash(&[4000, 3998], QueryStrand::Plus));
    assert!(hash(&[4000, 3999], QueryStrand::Plus));
    // minus
    assert!(!hash(&[8000], QueryStrand::Minus));
    assert!(hash(&[8001], QueryStrand::Minus));
    assert!(!hash(&[4000, 3999], QueryStrand::Minus));
    assert!(hash(&[4000, 4000], QueryStrand::Minus));
    // both
    assert!(!hash(&[3999], QueryStrand::Both));
    assert!(hash(&[4000], QueryStrand::Both));
    // The value is the searched block's length (the query buffer without its sentinels)
    // plus 1 for `plus`, and that length otherwise.
    for (lengths, strand, extra) in [
        (vec![10usize, 3, 10], QueryStrand::Plus, 1),
        (vec![10, 3, 10], QueryStrand::Minus, 0),
        (vec![10, 3, 10], QueryStrand::Both, 0),
    ] {
        let contexts = strand_contexts(&lengths, strand);
        let last = contexts.last().unwrap();
        let concat_len =
            (last.query_offset as usize + last.len + if last.len == 0 { 1 } else { 2 })
                .saturating_sub(2);
        assert_eq!(
            diag_container_query_length(&contexts),
            concat_len + extra,
            "{strand:?}"
        );
    }
    assert_eq!(diag_container_query_length(&[]), 0);
}

#[test]
fn evalue_reap_keeps_unless_greater_than_cutoff() {
    assert!(hsp_survives_evalue_reap(10.0, 10.0));
    assert!(!hsp_survives_evalue_reap(1e-180, 0.0));
    assert!(hsp_survives_evalue_reap(0.0, 0.0));
    assert!(!hsp_survives_evalue_reap(10.000001, 10.0));
    assert!(hsp_survives_evalue_reap(f64::NAN, 10.0));
    assert!(hsp_survives_evalue_reap(10.0, f64::NAN));
    assert!(hsp_survives_evalue_reap(f64::INFINITY, f64::INFINITY));
    assert!(!hsp_survives_evalue_reap(f64::INFINITY, f64::MAX));
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
fn make_equal_prelim_hit(query_idx: u32, seed_qs: i32, seed_ss: usize) -> PrelimHit {
    PrelimHit {
        context_idx: 0,
        query_idx,
        query_frame: 1,
        query_context_offset: 0,
        prelim_qs: 10,
        prelim_qe: 30,
        prelim_ss: 20,
        prelim_se: 40,
        prelim_score: 50,
        seed_qs,
        seed_ss,
        prelim_evalue: 0.0,
    }
}

#[allow(clippy::too_many_arguments)]
fn prelim_hit(
    context_idx: u32,
    query_frame: i32,
    prelim_qs: i32,
    prelim_qe: i32,
    prelim_ss: usize,
    prelim_se: usize,
    prelim_score: i32,
) -> PrelimHit {
    PrelimHit {
        context_idx,
        query_idx: 0,
        query_frame,
        query_context_offset: 0,
        prelim_qs,
        prelim_qe,
        prelim_ss,
        prelim_se,
        prelim_score,
        seed_qs: prelim_qs,
        seed_ss: prelim_ss,
        prelim_evalue: 0.0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2912-2945
#[test]
fn a_query_split_merges_hsps_across_the_split_on_either_strand() {
    let split = HspListSplit::Query {
        offsets: [1000, 5000],
    };
    // Plus strand: the earlier chunk's HSP ends after the split point, the later
    // chunk's starts before its overlap ends.
    let mut combined = vec![prelim_hit(0, 1, 900, 1050, 100, 250, 150)];
    let incoming = vec![prelim_hit(0, 1, 1000, 1200, 200, 400, 200)];
    merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
    assert_eq!(combined.len(), 1);
    assert_eq!((combined[0].prelim_qs, combined[0].prelim_qe), (900, 1200));
    // The score density of the two (350 / 350) over the merged 300 residues.
    assert_eq!(combined[0].prelim_score, 300);
    // Minus strand: the later chunk comes first on the strand, so the tests swap.
    let mut combined = vec![prelim_hit(1, -1, 4950, 5200, 1000, 1250, 150)];
    let incoming = vec![prelim_hit(1, -1, 4800, 5050, 850, 1100, 150)];
    merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
    assert_eq!(combined.len(), 1);
    assert_eq!((combined[0].prelim_qs, combined[0].prelim_qe), (4800, 5200));
    // The same pair on the plus strand is not in the strips (the earlier HSP does not
    // end after the split point): both stay.
    let mut combined = vec![prelim_hit(0, 1, 800, 950, 1000, 1150, 150)];
    let incoming = vec![prelim_hit(0, 1, 700, 850, 900, 1050, 150)];
    merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
    assert_eq!(combined.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1355
#[test]
fn stable_prelim_sort_preserves_equal_complete_record_order() {
    let mut hsps = vec![
        make_equal_prelim_hit(7, 101, 201),
        make_equal_prelim_hit(3, 102, 202),
        make_equal_prelim_hit(9, 103, 203),
    ];

    qsort_prelim_hits_by(&mut hsps, score_compare_prelim_hits);

    assert_eq!(
        hsps.iter().map(|hsp| hsp.query_idx).collect::<Vec<_>>(),
        vec![7, 3, 9]
    );
    assert_eq!(
        hsps.iter()
            .map(|hsp| (hsp.seed_qs, hsp.seed_ss))
            .collect::<Vec<_>>(),
        vec![(101, 201), (102, 202), (103, 203)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2532
#[test]
fn prelim_purge_retains_the_first_complete_equal_record() {
    let first = make_equal_prelim_hit(7, 101, 201);
    let second = make_equal_prelim_hit(3, 102, 202);

    let kept = purge_prelim_hits_with_common_endpoints(vec![first, second]);

    assert_eq!(kept.len(), 1);
    assert_eq!(kept[0].query_idx, 7);
    assert_eq!((kept[0].seed_qs, kept[0].seed_ss), (101, 201));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:689-706
// Each range starts at the last letter of the previous mask, as NCBI's ranges read
// from `TSequenceRanges::get_data` (seqdb.hpp:280-282).
#[test]
fn test_build_subject_seq_ranges_from_masks() {
    let masks = vec![MaskedInterval::new(2, 5), MaskedInterval::new(8, 12)];
    let ranges = build_subject_seq_ranges_from_masks(&masks, 10);
    assert_eq!(ranges, vec![(0, 2), (4, 8), (9, 10)]);
    let masks = vec![MaskedInterval::new(0, 5)];
    assert_eq!(
        build_subject_seq_ranges_from_masks(&masks, 32),
        vec![(0, 0), (4, 32)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
#[test]
fn test_adjust_blastn_offsets_minus_query_keeps_internal_subject_order() {
    let adjusted = adjust_blastn_offsets(10, 20, 30, 40, 100, -1);
    assert_eq!(adjusted, (81, 90, 40, 31));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1887-1890
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1923-1928
#[test]
fn test_blastn_context_statistics_uses_supplied_eff_searchsp() {
    let params = crate::stats::KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 1.0,
        beta: -1.0,
    };
    let eff_searchsp = 44_573_947_800_i64;
    let raw_score = 46;

    let (bit_score, evalue) =
        calculate_blastn_context_statistics(raw_score, &params, eff_searchsp, false);

    let expected_bit_score =
        (params.lambda * raw_score as f64 - params.k.glibc_ln()) / NCBIMATH_LN2;
    let expected_evalue = (eff_searchsp as f64)
        * (-(params.lambda) * raw_score as f64 + params.k.glibc_ln()).glibc_exp();

    assert_eq!(bit_score, expected_bit_score);
    assert_eq!(evalue, expected_evalue);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1864-1870
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1923-1928
#[test]
fn test_blastn_context_statistics_rounds_evalue_score_only() {
    let params = crate::stats::KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.8,
        beta: -2.0,
    };
    let eff_searchsp = 86_850_803_746_i64;
    let raw_score = 45;

    let (bit_score, evalue) =
        calculate_blastn_context_statistics(raw_score, &params, eff_searchsp, true);

    let expected_bit_score =
        (params.lambda * raw_score as f64 - params.k.glibc_ln()) / NCBIMATH_LN2;
    let expected_evalue_score = raw_score & !1;
    let expected_evalue = (eff_searchsp as f64)
        * (-(params.lambda) * expected_evalue_score as f64 + params.k.glibc_ln()).glibc_exp();

    assert_eq!(bit_score, expected_bit_score);
    assert_eq!(evalue, expected_evalue);
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:204-210;
// c++/src/algo/blast/core/blast_parameters.c:925-946
// NCBI BLAST+ 2.17.0+ raw output is frozen for three unequal subjects.
// Their combined search space excludes weak word7 seeds which could grow
// into extra gapped HSPs under an incorrectly lowered per-subject cutoff.
#[test]
fn word7_uses_subject_set_cutoffs() {
    let sequence: String = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/small_test.fasta"
    ))
    .lines()
    .filter(|line| !line.starts_with('>'))
    .map(str::trim)
    .collect();
    let query: String = (0..3)
        .map(|index| format!(">seq{index}\n{}\n", &sequence[..900]))
        .collect();
    let subject: String = [899, 900, 903]
        .into_iter()
        .enumerate()
        .map(|(index, length)| format!(">seq{index}\n{}\n", &sequence[..length]))
        .collect();
    let expected = include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/unit/helpers/ncbi_reference_data/blastn_word7_subject_set.out"
    ));
    let threads: &[usize] = if cfg!(feature = "parallel") {
        &[1, 4, 8, 16]
    } else {
        &[1]
    };
    for n in threads {
        let options = crate::cli::try_parse_from([
            "losat",
            "blastn",
            "-query",
            "memory-query",
            "-subject",
            "memory-subject",
            "-task",
            "blastn",
            "-word_size",
            "7",
            "-outfmt",
            "6",
            "-num_threads",
            &n.to_string(),
        ])
        .unwrap();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        let mut output = Vec::new();
        let mut stderr = std::io::stderr();
        run_local(
            options.blastn(),
            &reader_records(&query, ReaderConfig::query("BLASTN", false, false)),
            &crate::api::QueryEnd::Input,
            &reader_records(&subject, ReaderConfig::subject("BLASTN", false, false)),
            &mut ReportOutputs::single("6", OutputSink::Writer(&mut output), &mut stderr),
        )
        .unwrap();
        assert_eq!(output, expected.as_slice(), "n{n}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-406,509-513,583-612
// Repeated, unequal-length gapped matches generate dependent prelim HSPs
// across batch boundaries. Raw tabular bytes guard accepted HSP order and
// endpoints; the DP check below compares scores and complete edit scripts.
#[cfg(feature = "parallel")]
#[test]
fn speculative_traceback_preserves_ordered_gapped_hsps() {
    use rayon::prelude::*;
    let mut state = 0x12345678_u32;
    let query: Vec<u8> = (0..512)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state & 3) as usize]
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4170-4189
    // Subjects shorter than `MAX_SUBJECT_OFFSET` keep a zero `start_shift` and return early.
    // Exercise a real subject shift and more than one 16-HSP batch.
    let mut subject = vec![b'N'; 91_000];
    for i in 0..20 {
        subject.extend_from_slice(&query[..240]);
        subject.extend_from_slice(b"GATTACA");
        subject.extend_from_slice(&query[240..512 - i * 13]);
        subject.extend_from_slice(&[b'N'; 64]);
    }
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-traceback-regression-{}-{}",
        crate::utils::test_dir::process_id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&directory).unwrap();
    let query_path = directory.join("query.fa");
    let subject_path = directory.join("subject.fa");
    let output_path = directory.join("output.txt");
    std::fs::write(
        &query_path,
        format!(">query\n{}\n", String::from_utf8(query).unwrap()),
    )
    .unwrap();
    std::fs::write(
        &subject_path,
        format!(">subject\n{}\n", String::from_utf8(subject).unwrap()),
    )
    .unwrap();
    let execute = |threads: usize, outfmt: &str| {
        let args = crate::cli::try_parse_from([
            "losat",
            "blastn",
            "-query",
            query_path.to_str().unwrap(),
            "-subject",
            subject_path.to_str().unwrap(),
            "-task",
            "blastn",
            "-word_size",
            "11",
            "-num_threads",
            &threads.to_string(),
            "-outfmt",
            outfmt,
            "-out",
            output_path.to_str().unwrap(),
        ])
        .unwrap()
        .blastn();
        run(
            args,
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap();
        std::fs::read(&output_path).unwrap()
    };
    let sequential = execute(1, "6");
    let text = std::str::from_utf8(&sequential).unwrap();
    assert!(
        text.lines().count() >= 20,
        "must cross the sixteen-HSP batch boundary"
    );
    assert!(text
        .lines()
        .any(|line| line.split('\t').nth(5).unwrap() != "0"));
    for threads in [4, 2, 4] {
        assert_eq!(execute(threads, "6"), sequential, "n{threads}");
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-513
    // Multiple seeds extend to the same gapped HSP. Compute in reverse input
    // order with independent scratch, then replay in the original seed order.
    let q = encode_iupac_to_blastna(
        reader_records(
            &std::fs::read_to_string(&query_path).unwrap(),
            ReaderConfig::query("BLASTN", false, false),
        )[0]
        .seq(),
    );
    let subj = encode_iupac_to_blastna(
        reader_records(
            &std::fs::read_to_string(&subject_path).unwrap(),
            ReaderConfig::subject("BLASTN", false, false),
        )[0]
        .seq(),
    );
    let matrix = build_blastna_matrix(2, -3);
    let seeds = [60, 200, 300, 400, 80, 180, 320, 420];
    let dp = |index: usize, scratch: &mut GapAlignScratch| {
        let qs = seeds[index];
        extend_gapped_heuristic_with_traceback_with_scratch(
            &q,
            &subj,
            qs,
            91_000 + qs + if qs >= 240 { 7 } else { 0 },
            1,
            2,
            -3,
            &matrix,
            5,
            2,
            60,
            scratch,
        )
    };
    let mut scratch = GapAlignScratch::new();
    let ordered: Vec<_> = (0..8).map(|i| dp(i, &mut scratch)).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:177-188
    let reversed: Vec<_> = crate::utils::threading::with_search_pool(4, "blastn", |pool| {
        Ok(pool.install(|| {
            (0..8)
                .into_par_iter()
                .rev()
                .map_init(GapAlignScratch::new, |scratch, i| dp(i, scratch))
                .collect()
        }))
    })
    .unwrap();
    assert_eq!(ordered, reversed.into_iter().rev().collect::<Vec<_>>());
    assert!(ordered.iter().all(|hsp| hsp.7 > 0 && !hsp.9.is_empty()));
    std::fs::remove_dir_all(directory).unwrap();
}

fn query_prelim_hit(query_idx: u32, prelim_score: i32) -> PrelimHit {
    PrelimHit {
        query_idx,
        ..prelim_hit(query_idx * 2, 1, 10, 30, 20, 40, prelim_score)
    }
}

fn prelim_scores(hit_list: &HitList<PrelimHspList>) -> Vec<(u32, Vec<i32>)> {
    hit_list
        .hsplist_array
        .iter()
        .map(|list| {
            (
                list.oid,
                list.hsps.iter().map(|hsp| hsp.prelim_score).collect(),
            )
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:116-161
#[test]
fn the_collector_keeps_prelim_hitlist_size_subjects_per_query() {
    let subject_hits = vec![
        vec![query_prelim_hit(0, 30), query_prelim_hit(1, 90)],
        vec![query_prelim_hit(0, 50)],
        vec![query_prelim_hit(0, 40), query_prelim_hit(0, 45)],
        Vec::new(),
    ];
    let mut hit_lists = collect_prelim_hit_lists(subject_hits, 3, 2);
    assert!(hit_lists[2].is_none());
    let mut query0 = hit_lists[0].take().unwrap();
    // Subject 2 replaces subject 0 (the worst first score); a list that enters the
    // heap is sorted by e-value, then score (`Blast_HSPListSortByEvalue`).
    query0.sort_by_evalue();
    assert_eq!(
        prelim_scores(&query0),
        vec![(1, vec![50]), (2, vec![45, 40])]
    );
    let mut query1 = hit_lists[1].take().unwrap();
    query1.sort_by_evalue();
    assert_eq!(prelim_scores(&query1), vec![(0, vec![90])]);
}

// `LOSAT_X_BNHITGROUP`: the lists made from a subject's HSPs alone equal the lists made with
// one list per query of the batch, in every hit list, with the HSPs of a list in their order
// (each HSP has its own score). Queries in random order, in runs and all in one query.
#[test]
fn grouped_hit_lists_equal_the_lists_per_query() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = move |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    for round in 0..200u64 {
        let num_queries = 1 + next(40) as usize;
        let size = 1 + next(5) as usize;
        let mut per_query = PrelimHitCollector::new(num_queries, size);
        let mut grouped = PrelimHitCollector::new(num_queries, size);
        let mut score = 0;
        for oid in 0..30 {
            let count = next(25);
            let run = 1 + next(4);
            let mut hits = Vec::new();
            for i in 0..count {
                let query = match round % 3 {
                    0 => next(num_queries as u64),
                    1 => (i / run).min(num_queries as u64 - 1),
                    _ => num_queries as u64 - 1,
                };
                score += 1;
                hits.push(query_prelim_hit(query as u32, score));
            }
            if hits.is_empty() {
                continue;
            }
            per_query.add_per_query(oid, hits.clone());
            grouped.x_add_grouped(oid, hits);
        }
        for (a, b) in per_query.hit_lists.iter().zip(&grouped.hit_lists) {
            assert_eq!(a.as_ref().map(prelim_scores), b.as_ref().map(prelim_scores));
            assert_eq!(a.is_some(), b.is_some());
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2163-2210
#[test]
fn merging_prelim_hit_lists_appends_shared_subjects_and_keeps_the_size() {
    let list = |oid: u32, scores: &[i32]| PrelimHspList {
        oid,
        hsps: scores
            .iter()
            .map(|&score| query_prelim_hit(0, score))
            .collect(),
        best_evalue: 0.0,
    };
    let mut hitlist1 = HitList::new(3);
    hitlist1.update(list(4, &[50]));
    hitlist1.update(list(2, &[70]));
    let mut combined = HitList::new(3);
    combined.update(list(9, &[80]));
    combined.update(list(2, &[60]));
    combined.update(list(7, &[40]));
    let mut combined = Some(combined);
    merge_prelim_hit_list(hitlist1, &mut combined, [0, 0], 100);
    let mut merged = combined.unwrap();
    assert_eq!(merged.hsplist_max, 3);
    merged.sort_by_evalue();
    assert_eq!(
        prelim_scores(&merged),
        vec![(9, vec![80]), (2, vec![70, 60]), (4, vec![50])]
    );
    let mut empty = None;
    let mut single = HitList::new(3);
    single.update(list(5, &[10]));
    merge_prelim_hit_list(single, &mut empty, [0, 0], 100);
    assert_eq!(prelim_scores(empty.as_ref().unwrap()), vec![(5, vec![10])]);
}

/// The NCBI errors of BLASTN's `-outfmt` at NCBI's positions (blast_args.cpp:2745-2748,
/// 2874-2886, 2975-2977; blastn_app.cpp:255-258).
#[test]
fn output_formats_raise_ncbis_errors_where_ncbi_does() {
    let args = |words: &[&str]| {
        let mut argv = vec!["losat", "blastn", "-query", "q", "-subject", "s"];
        argv.extend_from_slice(words);
        crate::cli::try_parse_from(argv).unwrap().blastn()
    };
    let native = |error: anyhow::Error| {
        let native = error.downcast::<crate::cli::NativeError>().unwrap();
        (native.exit, String::from_utf8(native.message).unwrap())
    };
    let options = |text: &str| {
        (
            1,
            format!(
                "LOSAT query/options error: {text}\n{}\n",
                crate::wording::USER_MANUAL_FOOTER
            ),
        )
    };
    // The CLI: a format that NCBI rejects later has no report (`cli_output_formats`).
    let file = Some(std::path::Path::new("o.json"));
    for outfmt in ["13", "14", "19", "21"] {
        assert!(
            cli_output_formats(outfmt, None, false).unwrap().is_empty(),
            "{outfmt}"
        );
    }
    assert!(
        cli_output_formats("13", Some(std::path::Path::new("-")), false)
            .unwrap()
            .is_empty()
    );
    // 13 and 14 to a file are written (OF5c part 1): the XML2 / JSON documents as files.
    for (outfmt, kind) in [
        ("13", crate::report::blast_xml2::Xml2Kind::Files13),
        ("14", crate::report::blast_xml2::Xml2Kind::Files14),
    ] {
        assert_eq!(
            cli_output_formats(outfmt, file, false).unwrap(),
            [BlastnOutputFormat::Xml2(kind)],
            "{outfmt}"
        );
    }
    // Line RG T3b: the query-anchored formats 1 to 4, outfmt 0's report with another section.
    assert_eq!(
        cli_output_formats("1", file, false).unwrap(),
        [BlastnOutputFormat::Pairwise(
            crate::report::pairwise::PairwiseView::QueryAnchored {
                identity: true,
                anchored: true
            }
        )]
    );
    // Outfmt 18, the organism report (line RG T2), rides the custom tabular report's hit
    // lists; NCBI keeps no field list for it.
    for outfmt in ["18", "18 qseqid"] {
        let formats = cli_output_formats(outfmt, file, false).unwrap();
        assert!(
            matches!(
                formats.as_slice(),
                [BlastnOutputFormat::CustomTabular(spec)]
                    if spec.kind == crate::report::tabular::TabularKind::TaxReport
                        && spec.fields.is_empty()
            ),
            "{outfmt}: {formats:?}"
        );
    }
    // SAM (OF6) reads its words and ignores the rest and a delimiter.
    assert_eq!(
        cli_output_formats("17", file, false).unwrap(),
        [BlastnOutputFormat::Sam(
            crate::report::sam::SamSpec::default()
        )]
    );
    assert_eq!(
        cli_output_formats(" 17 delim=, SR SQ foo", None, false).unwrap(),
        [BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
            subject_as_reference: true,
            sequence: true
        })]
    );
    assert_eq!(
        cli_output_formats("0 qaccver", None, false).unwrap(),
        [BlastnOutputFormat::Pairwise(
            crate::report::pairwise::PairwiseView::Pairwise
        )]
    );
    // A field list, a delimiter, and outfmt 10 and 20 are the custom tabular report (line CLI
    // C2); the default 6 and 7 (`std` with a tab) keep their own writers.
    let custom = |outfmt: &str| match cli_output_formats(outfmt, None, false).unwrap().as_slice() {
        [BlastnOutputFormat::CustomTabular(spec)] => (*spec).clone(),
        other => panic!("{outfmt}: {other:?}"),
    };
    assert_eq!(
        custom("7 qaccver").fields,
        [crate::report::tabular::TabularField::QueryAccessionVersion]
    );
    assert_eq!(custom("6 delim=, std").delimiter, ",");
    assert_eq!(
        custom("20").kind,
        crate::report::tabular::TabularKind::CsvWithHeader
    );
    assert_eq!(
        cli_output_formats("6 std", None, false).unwrap(),
        [BlastnOutputFormat::Tabular]
    );
    // A tab given as `delim=` keeps the default writer; NCBI takes `delim=\t` as written, a
    // backslash and a `t`.
    assert_eq!(
        cli_output_formats("7 delim=\t std", None, false).unwrap(),
        [BlastnOutputFormat::TabularWithComments]
    );
    assert_eq!(custom(r"7 delim=\t").delimiter, r"\t");
    assert_eq!(
        native(cli_output_formats("22", None, false).unwrap_err()),
        (
            255,
            "Error: Formatting choice is out of range\n".to_string()
        )
    );
    // The formatting handler raises its errors before its warning; 13 passes it.
    for (outfmt, error) in [
        (
            "19",
            Some("AIRR rearrangement format is only applicable to igblastn"),
        ),
        (
            "21",
            Some("FASTA output format is only applicable to magicblast"),
        ),
        ("13", None),
    ] {
        let mut diagnostics = Vec::new();
        let mut sink = Vec::new();
        let mut outputs =
            ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
        let result = process_options(&mut args(&["-max_target_seqs", "1"]), &mut outputs);
        drop(outputs);
        match error {
            Some(error) => {
                assert_eq!(native(result.unwrap_err()), options(error), "{outfmt}");
                assert!(diagnostics.is_empty(), "{outfmt}");
            }
            None => {
                result.unwrap();
                assert_eq!(
                    String::from_utf8(diagnostics).unwrap(),
                    "Warning: [blastn] Examining 5 or more matches is recommended\n"
                );
            }
        }
    }
    // An ungapped search writes the Seq-aligns (8, 9, 11, 12: Dense-diags) and SAM (17, a
    // subject's HSPs joined into one record).
    for outfmt in ["8", "9", "11", "12", "17", "17 SR"] {
        let mut diagnostics = Vec::new();
        let mut sink = Vec::new();
        let mut outputs =
            ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
        let result = process_options(
            &mut args(&["-ungapped", "-max_target_seqs", "1"]),
            &mut outputs,
        );
        drop(outputs);
        result.unwrap_or_else(|error| panic!("{outfmt}: {error}"));
    }
    // `run_local`'s formats that NCBI rejects later are left out, as the command line's
    // (`cli_output_formats`): NCBI's formatting handler (`process_options`, above) and
    // `search` (the XInclude formats without a file) raise their errors after `-out` is made.
    for outfmt in ["13", "14", "19", "21"] {
        assert_eq!(
            parse_output_formats([("6", None), (outfmt, None)], false).unwrap(),
            [BlastnOutputFormat::Tabular],
            "{outfmt}"
        );
    }
    assert_eq!(
        parse_output_formats([("0", None), ("6", None), ("07", None)], false).unwrap(),
        [
            BlastnOutputFormat::Pairwise(crate::report::pairwise::PairwiseView::Pairwise),
            BlastnOutputFormat::Tabular,
            BlastnOutputFormat::TabularWithComments
        ]
    );
    assert_eq!(
        parse_output_formats([("17", None), ("17 SQ SR", None)], false).unwrap(),
        [
            BlastnOutputFormat::Sam(crate::report::sam::SamSpec::default()),
            BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
                subject_as_reference: true,
                sequence: true
            })
        ]
    );
}

// SubjectHits::ByQuery: the lists of a subject's queries, each sorted by score, give the hit
// lists that their join sorted by score (stable, as `x_qgroup_par::sort_stable`) gave as one
// subject list: ties of score and place within and across queries, queries without HSPs, a
// query past the batch, several subjects and hit lists that drop subjects.
#[test]
fn query_lists_give_the_hit_lists_of_the_sorted_subject_list() {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut below = |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n.max(1)
    };
    type Row = (u32, u32, i32, i32, i32, usize, usize, u64, u64);
    let rows = |hit_lists: &[Option<BlastnHitList>]| -> Vec<Option<Vec<(u32, u32, Vec<Row>)>>> {
        hit_lists
            .iter()
            .map(|list| {
                list.as_ref().map(|list| {
                    list.hsplist_array
                        .iter()
                        .map(|hsps| {
                            let rows = hsps
                                .hsps
                                .iter()
                                .map(|h| {
                                    (
                                        h.q_idx,
                                        h.s_idx,
                                        h.raw_score,
                                        h.internal_q_offset_0,
                                        h.internal_q_end_0,
                                        h.internal_s_offset_0,
                                        h.internal_s_end_0,
                                        h.e_value.to_bits(),
                                        h.bit_score.to_bits(),
                                    )
                                })
                                .collect();
                            (hsps.oid, hsps.query_index, rows)
                        })
                        .collect()
                })
            })
            .collect()
    };
    for _ in 0..300 {
        let queries = 1 + below(10) as usize;
        let subjects = 1 + below(4) as u32;
        let hitlist_size = 1 + below(3) as usize;
        let mut joined: Vec<Option<BlastnHitList>> = (0..queries).map(|_| None).collect();
        let mut by_query: Vec<Option<BlastnHitList>> = (0..queries).map(|_| None).collect();
        for s_idx in 0..subjects {
            let mut lists = Vec::new();
            // one query past the batch now and then (both paths drop its HSPs)
            for q_idx in 0..queries as u32 + u32::from(below(8) == 0) {
                if below(3) == 0 {
                    continue;
                }
                let mut list: Vec<BlastnHsp> = (0..1 + below(6))
                    .map(|_| {
                        let q = below(4) as usize;
                        let s = below(4) as usize;
                        BlastnHsp {
                            identity: 100.0,
                            length: 10,
                            mismatch: 0,
                            gapopen: 0,
                            q_start: q + 1,
                            q_end: q + 10,
                            s_start: s + 1,
                            s_end: s + 10,
                            e_value: below(5) as f64 * 1e-3,
                            bit_score: below(50) as f64,
                            num_ident: 10,
                            query_frame: 1,
                            query_length: 100,
                            q_idx,
                            s_idx,
                            raw_score: below(4) as i32,
                            internal_q_offset_0: q as i32,
                            internal_q_end_0: (q + 1 + below(2) as usize) as i32,
                            internal_s_offset_0: s,
                            internal_s_end_0: s + 1 + below(2) as usize,
                            internal_query_context_offset: 0,
                            gap_info: None,
                            num_positives: 10,
                        }
                    })
                    .collect();
                sort_hsps_by_score(&mut list);
                lists.push(list);
            }
            let mut sorted: Vec<BlastnHsp> = lists.iter().flatten().cloned().collect();
            sorted.sort_by(crate::algorithm::blastn::hsp::score_compare_hsps);
            update_hitlists_with_subject_hits(
                &mut joined,
                SubjectHits::Mixed(sorted),
                hitlist_size,
            );
            update_hitlists_with_subject_hits(
                &mut by_query,
                SubjectHits::ByQuery(lists),
                hitlist_size,
            );
        }
        assert_eq!(rows(&by_query), rows(&joined));
    }
}

/// One chunk of NCBI's `s_GetNextSubjectChunk` model: offset, length, the seq_ranges and
/// the overlap of `Blast_HSPListsMerge`.
type ModelChunk = (i32, i32, Vec<(i32, i32)>, i32);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310 (cited at
// `SubjectSplitState::next_chunk`), 579-584 (the overlap), 145-179 (`s_BackupSubject`).
/// NCBI's chunks of a nucleotide subject of `length` letters with the hard ranges `hard`
/// (`None`: `&full_range`) and the soft ranges `soft` (`None`: `&full_range`), statement by
/// statement in `Int4` (the mask type is hard with `hard`, soft with `soft`). The overlap
/// after the last hard range reads one range past the list in NCBI (stale memory); the model
/// reads the full range's end there, LOSAT's choice (decision H5).
fn ncbi_chunks(
    length: i32,
    hard: Option<&[(i32, i32)]>,
    soft: Option<&[(i32, i32)]>,
) -> Vec<ModelChunk> {
    ncbi_ranged_chunks(length, hard, soft)
        .into_iter()
        .map(|(chunk, _)| chunk)
        .collect()
}

/// `ncbi_chunks` with each chunk's hard range (`hm_index` when it is cut).
fn ncbi_ranged_chunks(
    length: i32,
    hard: Option<&[(i32, i32)]>,
    soft: Option<&[(i32, i32)]>,
) -> Vec<(ModelChunk, usize)> {
    const MAX: i32 = 5_000_000;
    const OVERLAP: i32 = 100;
    let full = [(0, length)];
    let hard_ranges = hard.unwrap_or(&full);
    let soft_ranges = soft.unwrap_or(&full);
    let soft_masking = soft.is_some();
    let (mut hm_index, mut sm_index) = (0usize, 0usize);
    let mut offset = hard_ranges[0].0;
    let mut next = offset;
    let mut chunks = Vec::new();
    loop {
        if next >= length {
            return chunks;
        }
        let residual = next % 4;
        offset = next - residual;
        let range = hm_index;
        let subject_length;
        if offset + MAX < hard_ranges[hm_index].1 {
            subject_length = MAX;
            next = offset + MAX - OVERLAP;
        } else {
            subject_length = hard_ranges[hm_index].1 - offset;
            hm_index += 1;
            next = if hm_index < hard_ranges.len() {
                hard_ranges[hm_index].0
            } else {
                length
            };
        }
        let seq_ranges: Vec<(i32, i32)> = if offset == 0 && residual == 0 && next == length {
            soft_ranges.to_vec()
        } else if !soft_masking {
            vec![(residual, subject_length)]
        } else {
            let start = offset;
            let end = start + subject_length;
            let mut i = sm_index;
            while soft_ranges[i].1 < start {
                i += 1;
            }
            let first = i;
            while i < soft_ranges.len() && soft_ranges[i].0 < end {
                i += 1;
            }
            let len = i - first;
            sm_index = i - 1;
            if len == 0 {
                continue;
            }
            let mut ranges: Vec<(i32, i32)> = soft_ranges[first..first + len]
                .iter()
                .map(|&(left, right)| (left - offset, right - offset))
                .collect();
            if ranges[0].0 < 0 {
                ranges[0].0 = 0;
            }
            if ranges[len - 1].1 > subject_length {
                ranges[len - 1].1 = subject_length;
            }
            ranges
        };
        let range_left = if hm_index < hard_ranges.len() {
            hard_ranges[hm_index].0
        } else {
            length
        };
        let overlap = if offset == range_left { 0 } else { OVERLAP };
        chunks.push(((offset, subject_length, seq_ranges, overlap), range));
    }
}

/// `SubjectSplitState`'s chunks, as the model gives them. A hard-masked subject has no soft
/// ranges (`read_db_seq_ranges`): its chunks carry their ranges inline (`next_hard_chunk`).
fn losat_chunks(
    length: i32,
    hard: Option<&[(i32, i32)]>,
    soft: Option<&[(i32, i32)]>,
) -> Vec<ModelChunk> {
    let masked = hard.is_some() || soft.is_some();
    let soft_ranges = match soft {
        Some(soft) => soft.to_vec(),
        None => Vec::new(),
    };
    let mut state = SubjectSplitState::new(length as usize, soft_ranges);
    let mut hard_state = HardSplitState::default();
    if let Some(hard) = hard {
        hard_state.start(length as usize, hard.to_vec(), true);
    }
    let mut scratch = Vec::new();
    let mut chunks = Vec::new();
    loop {
        let next = if hard.is_some() {
            next_subject_chunk::<true>(&mut state, &mut hard_state, masked, DBSEQ_CHUNK_OVERLAP)
        } else {
            next_subject_chunk::<false>(&mut state, &mut hard_state, masked, DBSEQ_CHUNK_OVERLAP)
        };
        match next {
            SubjectChunkStatus::Done => return chunks,
            SubjectChunkStatus::NoRange => continue,
            SubjectChunkStatus::Ok(chunk) => {
                let ranges = chunk
                    .seq_ranges
                    .as_slice(
                        &state.soft_ranges,
                        &mut scratch,
                        chunk.offset as i32,
                        chunk.length as i32,
                    )
                    .to_vec();
                chunks.push((
                    chunk.offset as i32,
                    chunk.length as i32,
                    ranges,
                    chunk.overlap as i32,
                ));
            }
        }
    }
}

/// The complement of sorted, disjoint masks (as `core::subject_split::complement_masks`).
fn complement(masks: &[(i32, i32)], length: i32) -> Vec<(i32, i32)> {
    let mut ranges = Vec::new();
    let mut left = 0;
    for &(start, end) in masks {
        ranges.push((left, start));
        left = end;
    }
    ranges.push((left, length));
    ranges
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310, 579-584 (quoted
// at `SubjectSplitState::next_chunk`)
// `-db_hard_mask`'s chunks (series NP, M4a-2b decision H1) against the model of NCBI's
// `s_GetNextSubjectChunk` above: a mask at 0 (a chunk of no letters), adjacent masks, residuals
// 0 to 3 (the chunk starts up to three masked letters before its range), one-letter masks, a
// range longer than MAX_DBSEQ_LEN (full chunks with the overlap 100, the first of a range that
// starts on a multiple of four with 0), masks at the end; then random sorted masks. The soft and
// unmasked chunks of the same subjects stay NCBI's too. One chunk differs from NCBI's by
// design (series NP, M4a-2c, `blastinput::db_mask::TAIL_MASK_SEARCHED_WORD`): where the "no
// chunking" branch gives a chunk shorter than the subject the full range (a mask that runs to
// the end after the chunk from 0), LOSAT scans the chunk, `[0, length)`; every other chunk,
// the "no chunking" branch's of a chunk of the whole subject included, is NCBI's.
#[test]
fn hard_masked_chunks_follow_ncbi_s_get_next_subject_chunk() {
    let check = |length: i32, masks: &[(i32, i32)]| {
        let hard = complement(masks, length);
        let found = losat_chunks(length, Some(&hard), None);
        let ncbi_but_the_exception: Vec<ModelChunk> = ncbi_chunks(length, Some(&hard), None)
            .into_iter()
            .map(|(offset, chunk_length, ranges, overlap)| {
                if offset == 0 && chunk_length < length && ranges == [(0, length)] {
                    (offset, chunk_length, vec![(0, chunk_length)], overlap)
                } else {
                    (offset, chunk_length, ranges, overlap)
                }
            })
            .collect();
        assert_eq!(found, ncbi_but_the_exception, "{length} {masks:?}");
        for chunk in &found {
            assert!(chunk.1 >= 0, "{length} {masks:?}: {chunk:?}");
        }
        let soft = complement(masks, length);
        assert_eq!(
            losat_chunks(length, None, Some(&soft)),
            ncbi_chunks(length, None, Some(&soft)),
            "soft {length} {masks:?}"
        );
        assert_eq!(
            losat_chunks(length, None, None),
            ncbi_chunks(length, None, None)
        );
        found
    };
    // A mask at 0: a chunk of length 0 at 0, then the chunk after the mask (residual 2).
    assert_eq!(
        check(1200, &[(0, 102)]),
        [(0, 0, vec![(0, 0)], 100), (100, 1100, vec![(2, 1100)], 100)]
    );
    // Adjacent masks: the hard range between them has no letters (residual 0, length 0).
    assert_eq!(
        check(1200, &[(100, 200), (200, 300)]),
        [
            (0, 100, vec![(0, 100)], 100),
            (200, 0, vec![(0, 0)], 100),
            (300, 900, vec![(0, 900)], 100),
        ]
    );
    for residual in 0..4 {
        let end = 600 + residual;
        let found = check(1200, &[(500, end)]);
        assert_eq!(found[1].0, 600);
        assert_eq!(found[1].2, vec![(residual, 600)]);
    }
    // One-letter masks whose chunks start before them (their HSP lists can merge).
    check(1200, &[(301, 302), (500, 503), (703, 706)]);
    // A range longer than MAX_DBSEQ_LEN: full chunks of it, the first with the overlap 0.
    let long = check(12_000_000, &[(100, 204), (11_000_000, 11_000_010)]);
    assert_eq!(long[1], (204, 5_000_000, vec![(0, 5_000_000)], 0));
    assert_eq!(long[2].0, 204 + 5_000_000 - 100);
    assert_eq!(long[2].3, 100);
    check(12_000_000, &[(100, 202)]);
    // A last mask at the end (the empty last range is not a chunk), masks at both ends.
    check(1200, &[(10, 20), (1100, 1200)]);
    check(1200, &[(0, 100), (1100, 1200)]);
    // Empty masks at 0 before a mask that runs to the end (review F1 of M4a-2b): the chunk
    // from 0 takes NCBI's "no chunking" branch, whose full range LOSAT does not scan past the
    // chunk (M4a-2c); so does a mask that runs to the end alone, and a subject masked whole
    // (a chunk of no letters).
    let found = check(1200, &[(0, 0), (1000, 1200)]);
    assert_eq!(found[1], (0, 1000, vec![(0, 1000)], 100));
    assert_eq!(
        check(1200, &[(1000, 1200)]),
        [(0, 1000, vec![(0, 1000)], 100)]
    );
    assert_eq!(check(1200, &[(0, 1200)]), [(0, 0, vec![(0, 0)], 100)]);
    // A mask of no letters at the end: the one chunk is the whole subject, NCBI's own.
    assert_eq!(
        check(1200, &[(1200, 1200)]),
        [(0, 1200, vec![(0, 1200)], 100)]
    );
    // Random sorted, disjoint masks (empty ones too, some sets starting with empty masks at
    // 0), every shape: the chunks are NCBI's arithmetic, refused or not.
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % bound) as i32
    };
    for case in 0..2000 {
        let length = if case % 50 == 0 {
            5_000_000 + next(3_000_000)
        } else {
            1 + next(3000)
        };
        let mut masks = vec![(0, 0); (next(4) - 1).max(0) as usize];
        let mut at = 0;
        for _ in 0..next(6) {
            let start = at + next(400) * (1 + length / 3000);
            let end = start + next(60);
            if end > length {
                break;
            }
            masks.push((start, end));
            at = end;
        }
        check(length, &masks);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310, 579-584 (quoted
// at `SubjectSplitState::next_chunk`); na_ungapped.c:1647-1668, masksubj.inl:43-58,
// blast_engine.c:493-500, 565 and blast_hits.c:2873-2911, 2956 (cited at `HardMaskShape::of_scan`)
// The refusals of `-db_hard_mask` against NCBI's chunks (the model above, not the chunker
// that the refusals run; reviews F1-F3 and T1 of M4a-2b, review 2 F1 and T1): a shape is
// refused as `TailMask` exactly when a chunk shorter than the sequence scans the full range
// and the scan reads a word (`length >= word`), at a word below `TAIL_MASK_SEARCHED_WORD`
// (series NP, M4a-2c: at 11 and more LOSAT searches it), and as `ShortLastMask` exactly when NCBI's
// merge after the final chunk can take an HSP of an earlier chunk by the overlap that it reads
// past its list: the final chunk is of the last hard range and scans a word (else NCBI does
// not merge it), and an earlier chunk of another hard range that scans a word ends after the
// final chunk's offset (`hspcnt1`; a chunk of the same range is LOSAT's overlap, decision
// B-d); an accepted shape has neither. A chunk scans a word when one of its `seq_ranges`
// holds `word` letters (`s_DetermineScanningOffsets`). Random sorted masks of sequences of 1
// to 64 letters (a third of the sets with empty masks at 0 first, words of 4 to 16 letters),
// and of sequences of up to 15,000,000 letters with one to four masks of 0 to 5 letters (some
// longer) at random places (words of 4 to 28 letters), where ranges over MAX_DBSEQ_LEN letters
// are cut into several chunks and the final chunk can follow a full chunk of its own range.
#[test]
fn hard_mask_refusals_are_the_merges_that_ncbi_decides_past_its_list() {
    use crate::blastinput::db_mask::{HardMaskShape, HardScan, TAIL_MASK_SEARCHED_WORD};
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    let mut next = |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % bound) as i32
    };
    // tails, shorts, accepted; of the long sequences: shorts, accepted with the final chunk
    // after a chunk of its own range (B-d), and accepted where only chunks that scan no word
    // end after the final chunk's offset (review 2 F1).
    let mut counts = [0usize; 6];
    for case in 0..40_000 {
        let long = case % 2 == 1;
        let (length, word) = if long {
            (1 + next(15_000_000), 4 + next(25))
        } else {
            (1 + next(64), 4 + next(13))
        };
        let mut masks = vec![
            (0, 0);
            if next(3) == 0 {
                1 + next(3) as usize
            } else {
                0
            }
        ];
        let mut at = 0;
        for _ in 0..1 + next(4) {
            let (start, end) = if long {
                let start = at + next((length - at).max(1) as u64);
                let letters = if next(8) == 0 { next(5_000) } else { next(6) };
                (start, (start + letters).min(length))
            } else {
                let start = at + next(8);
                (start, (start + next(12)).min(length))
            };
            if start > length {
                break;
            }
            masks.push((start, end));
            at = end;
        }
        let hard = complement(&masks, length);
        let chunks = ncbi_ranged_chunks(length, Some(&hard), None);
        let scans = |chunk: &ModelChunk| chunk.2.iter().any(|&(left, right)| left + word <= right);
        let scans_past = chunks
            .iter()
            .any(|(chunk, _)| chunk.1 < length && chunk.2 == [(0, length)] && length >= word);
        // The final chunk of the last hard range, after which NCBI reads the overlap past
        // the list (the loop ends there), and the earlier chunks that end after its offset.
        let (merges_past, after_own_range, only_no_word) = match chunks.split_last() {
            Some(((last, range), earlier)) if *range == hard.len() - 1 => {
                let reach = |chunk: &ModelChunk| chunk.0 + chunk.1 > last.0;
                (
                    scans(last)
                        && earlier
                            .iter()
                            .any(|(chunk, other)| other != range && scans(chunk) && reach(chunk)),
                    earlier.last().is_some_and(|(_, other)| other == range),
                    scans(last)
                        && earlier
                            .iter()
                            .any(|(chunk, other)| other != range && reach(chunk))
                        && !earlier
                            .iter()
                            .any(|(chunk, other)| other != range && scans(chunk) && reach(chunk)),
                )
            }
            _ => (false, false, false),
        };
        let stored: Vec<(u32, u32)> = masks
            .iter()
            .map(|&(start, end)| (start as u32, end as u32))
            .collect();
        let scan = HardScan::Nucleotide {
            scanned_word: word as usize,
        };
        let shape = HardMaskShape::of_scan(&stored, None, &|| Some(length as usize), scan);
        // The same shape from the bounds of the length that the index gives (the packed bytes
        // and 0 to 3 letters in the last one; review 2 S4).
        let low = (length - length % 4) as usize;
        assert_eq!(
            HardMaskShape::of_scan(
                &stored,
                Some((low, low + 3)),
                &|| Some(length as usize),
                scan
            ),
            shape,
            "{length} {masks:?} word {word}: bounds"
        );
        match shape {
            Ok(Some(HardMaskShape::TailMask { .. })) => counts[0] += 1,
            Ok(Some(HardMaskShape::ShortLastMask { .. })) => {
                counts[1] += 1;
                counts[3] += usize::from(long);
            }
            Ok(None) => {
                counts[2] += 1;
                counts[4] += usize::from(long && after_own_range);
                counts[5] += usize::from(only_no_word);
            }
            other => panic!("{length} {masks:?}: {other:?}"),
        }
        assert_eq!(
            matches!(shape, Ok(Some(HardMaskShape::TailMask { .. }))),
            scans_past && (word as usize) < TAIL_MASK_SEARCHED_WORD,
            "{length} {masks:?} word {word}: {chunks:?}"
        );
        assert_eq!(
            matches!(shape, Ok(Some(HardMaskShape::ShortLastMask { .. }))),
            merges_past && !scans_past,
            "{length} {masks:?} word {word}: {chunks:?}"
        );
    }
    assert!(
        counts[0] > 300
            && counts[1] > 1000
            && counts[2] > 10_000
            && counts[3] > 100
            && counts[4] > 1000
            && counts[5] > 100,
        "{counts:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2388
// A saved hit of context 2 (offset 1302 of the concatenated query) whose walk took 1082 letters
// of the contexts before it keeps NCBI's context-relative `Int4` start, -1082, its length, and
// the greedy seed `q_start + length/2` (blast_gapalign.c:4012-4017), negative here; the sort
// compares the concatenated query's start (`score_compare_match`, blast_extend.c:271-304),
// which no hit of the context before it (offset 1001) shares.
#[test]
fn a_hit_that_starts_left_of_its_context_keeps_ncbis_signed_start() {
    let ctx = QueryContext {
        query_idx: 1,
        frame: 1,
        query_offset: 1302,
        len: 300,
        masks: Vec::new(),
    };
    let ungapped = UngappedData {
        q_start: 220,
        s_start: 0,
        length: 1382,
        score: i32::MAX - 40,
    };
    let hit = UngappedHit::saved(2, &ctx, 1302 + 250, 1332, &ungapped);
    assert_eq!((hit.qs, hit.qe, hit.length()), (-1082, 300, 1382));
    assert_eq!((hit.seed_q_off, hit.seed_s_off), (250, 1332));
    assert_eq!(hit.greedy_seed(), (-391, 691));
    let in_context = UngappedHit::saved(
        2,
        &ctx,
        1302 + 250,
        1332,
        &UngappedData {
            q_start: 1302 + 200,
            s_start: 1232,
            length: 100,
            score: 90,
        },
    );
    assert_eq!((in_context.qs, in_context.qe), (200, 300));
    assert_eq!(in_context.greedy_seed(), (250, 1282));
    let ctx1 = QueryContext {
        query_idx: 0,
        frame: -1,
        query_offset: 1001,
        len: 300,
        masks: Vec::new(),
    };
    let before = UngappedHit::saved(
        1,
        &ctx1,
        1001 + 10,
        1332,
        &UngappedData {
            q_start: 221,
            ..ungapped
        },
    );
    // The same score, subject start and length: the concatenated query's start decides.
    assert_eq!(
        score_compare_ungapped_hits(&hit, &before),
        std::cmp::Ordering::Less
    );
    assert_eq!(
        score_compare_ungapped_hits(&before, &hit),
        std::cmp::Ordering::Greater
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:699-715 (the re-evaluation
// loop, `Blast_HSPReevaluateWithAmbiguitiesUngapped`) and :453-456 (the kept offsets)
// An `-ungapped` HSP of context 1 whose NCBI start lies 21 letters before it (20 letters of
// context 0 and the sentinel) is read from there: a part before the context below the cutoff
// is dropped at the sentinel (the HSP is the context's), one at the cutoff or more stays the
// HSP unless the context's letters score more (then the HSP is theirs), and where it stays
// NCBI's offsets are negative, kept (`PrelimHit::prelim_qs`, `prelim_qe`) with the identities
// of its letters before the context (`-perc_identity 100` keeps it).
#[test]
fn an_ungapped_hsp_is_reevaluated_from_ncbis_start_before_its_context() {
    let matrix = build_blastna_matrix(2, -3);
    // BLASTNA A, C, G, T = 0..3; `other` is a mismatch of each.
    let other = |b: u8| (b + 1) % 4;
    let ctx0: Vec<u8> = (0..20u8).map(|i| (i * 7 + 3) % 4).collect();
    let ctx1: Vec<u8> = (0..40u8).map(|i| (i * 5 + 1) % 4).collect();
    let concat: Vec<u8> = [ctx0.as_slice(), &[15], ctx1.as_slice()].concat();
    // `before` letters of ctx0's tail match the subject (the rest mismatch), then the letter
    // against the sentinel, then the first `inside` letters of ctx1 match.
    let run = |before: usize, inside: usize, percent_identity: f64| {
        let mut subject: Vec<u8> = ctx0
            .iter()
            .enumerate()
            .map(|(i, &b)| if i >= 20 - before { b } else { other(b) })
            .collect();
        subject.push(0);
        subject.extend(
            ctx1.iter()
                .enumerate()
                .map(|(i, &b)| if i < inside { b } else { other(b) }),
        );
        // NCBI's offsets -21 and 40.
        let mut hits = vec![PrelimHit {
            context_idx: 1,
            query_idx: 0,
            query_frame: -1,
            query_context_offset: 21,
            prelim_qs: -21,
            prelim_qe: 40,
            prelim_ss: 0,
            prelim_se: 61,
            prelim_score: i32::MAX - 5,
            seed_qs: 0,
            seed_ss: 21,
            prelim_evalue: 0.0,
        }];
        ungapped_subject_stage(
            &mut hits,
            &subject,
            &concat,
            &concat,
            &matrix,
            &[10, 10],
            percent_identity,
            |_| 0.0,
            PrelimReap::Evalue(10.0),
        );
        hits.iter()
            .map(|h| {
                (
                    h.prelim_qs,
                    h.prelim_qe,
                    h.prelim_ss,
                    h.prelim_se,
                    h.prelim_score,
                    h.prelim_qe - h.prelim_qs,
                )
            })
            .collect::<Vec<_>>()
    };
    // Below the cutoff before the context (2 matches, 4 < 10): the context's 15 matches.
    assert_eq!(run(2, 15, 0.0), vec![(0, 15, 21, 36, 30, 15)]);
    // At the cutoff or more before it (12 matches, 24) and 10 matches inside (20): the piece
    // before it, at NCBI's offsets -13 and -1, all 12 letters identical.
    assert_eq!(run(12, 10, 100.0), vec![(-13, -1, 8, 20, 24, 12)]);
    // The same piece before it, and 30 matches inside (60 > 24): the context's piece.
    assert_eq!(run(12, 30, 0.0), vec![(0, 30, 21, 51, 60, 30)]);
    // Ties keep the piece before the context (`sum > score`): 12 matches each side.
    assert_eq!(run(12, 12, 0.0), vec![(-13, -1, 8, 20, 24, 12)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:466-469,503-507,599
// A megablast traceback start more than 3000 letters outside its 30-letter context on a
// subject of 100000 letters (remaining-gaps X1b-o review 2 M1). Past the context's end
// `AdjustSubjectRange` gives a range shorter than the start (`s_avail` -1001) and the greedy
// aligns to the left of the start only; before the context it gives a start of -1001 in the
// range and the greedy aligns to the right only, also when the shift passes the subject's end
// (range length -1). Each alignment is the run of letters set to match, and its subject
// offsets after the shift are those letters' (no HSP is dropped and no slice panics).
#[test]
fn a_traceback_start_far_outside_its_context_aligns_on_a_long_subject() {
    // Contexts 0 (5000 letters), 1 (30 A) and 2 (5000), sentinels between them; contexts 0
    // and 2 cycle through the four letters, so no gap meets a run of matches.
    let cycle = |n: usize, a: usize| (0..n).map(|i| ((i * 7 + a) % 4) as u8).collect::<Vec<_>>();
    let concat = [
        cycle(5000, 1),
        vec![15],
        vec![0; 30],
        vec![15],
        cycle(5000, 2),
    ]
    .concat();
    let ctx = QueryContext {
        query_idx: 0,
        frame: 1,
        query_offset: 5001,
        len: 30,
        masks: Vec::new(),
    };
    let query = context_query(&concat, &ctx);
    // A subject of 100000 C (no match of the A context) with `run` (concatenated query
    // places) copied at `at`; the traceback from NCBI's start `(q_start, s_start)`, its
    // range, and the alignment's query offsets, subject offsets after the shift and score.
    let trace = |q_start: i32, s_start: usize, run: std::ops::Range<usize>, at: usize| {
        let mut subject = vec![1u8; 100_000];
        subject[at..at + run.len()].copy_from_slice(&concat[run]);
        let mut s_adj = s_start as i32;
        let mut s_len = subject.len() as i32;
        let shift = adjust_subject_range(&mut s_adj, &mut s_len, q_start, query.len);
        let start: traceback_hsp::TraceStart = (q_start, s_start, shift, s_len, s_adj);
        let (qs, qe, ss, se, score, _, ops) = greedy_gapped_alignment_with_traceback(
            query,
            traceback_hsp::adjusted_subject(&subject, &start),
            q_start,
            s_adj,
            1,
            -2,
            0,
            0,
            54,
            &mut GreedyAlignScratch::new(),
        )
        .expect("NCBI aligns from a start outside its context");
        (
            (shift, s_len, s_adj),
            (qs, qe, ss + shift, se + shift, score),
            ops,
        )
    };
    // Past the end: context 2's first 4000 letters end at offset 4031, copied to end at 4093.
    assert_eq!(
        trace(4031, 4093, 5032..9032, 93),
        (
            (0, 3092, 4093),
            (31, 4031, 93, 4093, 4000),
            vec![GapEditOp::Sub(4000)]
        )
    );
    // Before the start: context 0's last 4000 letters from offset -4001, copied from 50000.
    assert_eq!(
        trace(-4001, 50_000, 1000..5000, 50_000),
        (
            (51_001, 6030, -1001),
            (-4001, -1, 50_000, 54_000, 4000),
            vec![GapEditOp::Sub(4000)]
        )
    );
    // The same start near the subject's end: the shift passes it, and the right extension
    // reads the subject's last 1000 letters.
    assert_eq!(
        trace(-4001, 99_000, 1000..2000, 99_000),
        (
            (100_001, -1, -1001),
            (-4001, -3001, 99_000, 100_000, 1000),
            vec![GapEditOp::Sub(1000)]
        )
    );
}

// An HSP outside its query context (`BlastnHsp::internal_q_offset_0`: a start below 0 or an
// end past the context's length) reaches the reports; the search fails or refuses a batch
// where one has a Seq-align query start below 0 (`BlastnHsp::seqalign_query_start`,
// `left_of_query_refusal`: a plus-strand start left of the context, or a minus-strand
// end past it), and the reports show one past the query's end (remaining-gaps X1b-p). The
// drops before it compare NCBI's signed offsets: an HSP left of its context sorts before one
// at offset 0 with the same score and subject range.
#[test]
fn an_hsp_outside_its_context_sorts_by_ncbis_offsets() {
    let hsp = |offset: i32, end: i32| BlastnHsp {
        identity: 100.0,
        length: (end - offset) as usize,
        mismatch: 0,
        gapopen: 0,
        q_start: 0,
        q_end: 0,
        s_start: 0,
        s_end: 0,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: (end - offset) as usize,
        query_frame: 1,
        query_length: 100,
        q_idx: 0,
        s_idx: 0,
        raw_score: 24,
        internal_q_offset_0: offset,
        internal_q_end_0: end,
        internal_s_offset_0: 8,
        internal_s_end_0: 20,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: (end - offset) as usize,
    };
    let left = hsp(-13, -1);
    assert_eq!(left.seqalign_query_start(), -13);
    assert_eq!(hsp(100, 112).seqalign_query_start(), 100);
    let mut minus = hsp(-12, -2);
    minus.query_frame = -1;
    assert_eq!(minus.seqalign_query_start(), 102);
    // `starts_left_of_record` is `record_query_start(..) < 0`, tested first as X1b-o tested an
    // HSP outside its context: inside it (and on its right) no start is below 0.
    assert!(left.outside_query_context() && left.starts_left_of_record(0, false));
    assert!(left.starts_left_of_record(12, false) && left.starts_left_of_record(12, true));
    assert!(!left.starts_left_of_record(13, true));
    assert!(minus.outside_query_context() && !minus.starts_left_of_record(0, true));
    for inside in [hsp(0, 12), hsp(88, 100)] {
        assert!(!inside.outside_query_context());
        assert!(!inside.starts_left_of_record(0, false) && !inside.starts_left_of_record(0, true));
    }
    assert!(hsp(-1, 11).outside_query_context() && hsp(100, 112).outside_query_context());
    assert!(!hsp(100, 112).starts_left_of_record(0, false));
    (minus.internal_q_offset_0, minus.internal_q_end_0) = (95, 104);
    assert_eq!(minus.seqalign_query_start(), -4);
    assert!(minus.starts_left_of_record(0, false) && minus.starts_left_of_record(3, true));
    assert!(!minus.starts_left_of_record(4, true) && minus.starts_left_of_record(4, false));
    assert_eq!(
        crate::algorithm::blastn::hsp::score_compare_hsps(&left, &hsp(0, 12)),
        std::cmp::Ordering::Less
    );
}

// No NCBI counterpart (caller contract of LOSAT's `run_local` entry): a `-db` search must call
// `run_local_database`, and callers of `run_local` with WindowMasker or `-filtering_db` must
// populate the counts source, registry, or resolved repeat filter before calling `run_local`.
#[test]
fn run_local_refuses_database_and_unresolved_filter_inputs() {
    let parse = |words: &[&str]| {
        let mut argv = vec!["losat", "blastn", "-query", "q"];
        argv.extend_from_slice(words);
        crate::cli::try_parse_from(argv).unwrap().blastn()
    };
    let call = |args: BlastnArgs, outfmt: &str| {
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let result = {
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
            run_local(args, &[], &crate::api::QueryEnd::Input, &[], &mut outputs)
        };
        (result, sink, diagnostics)
    };

    let assert_refused = |args: BlastnArgs, outfmt: &str, expected: &str, label: &str| {
        let (result, sink, diagnostics) = call(args, outfmt);
        assert_eq!(format!("{:#}", result.unwrap_err()), expected, "{label}");
        assert!(sink.is_empty() && diagnostics.is_empty(), "{label}");
    };

    let db_msg = "-db (a database search) is run by run_local_database, not run_local";
    for outfmt in ["6", "11"] {
        assert_refused(parse(&["-db", "nt"]), outfmt, db_msg, outfmt);
    }

    let wm_db_msg = "-window_masker_db (a WindowMasker counts file) needs its counts source set before calling run_local (BlastnArgs::window_masker_counts)";
    assert_refused(
        parse(&["-subject", "s", "-window_masker_db", "counts.obinary"]),
        "6",
        wm_db_msg,
        "window_masker_db",
    );

    let wm_taxid_msg = "-window_masker_taxid (a directory of WindowMasker counts files) needs its registry and counts source set before calling run_local (BlastnArgs::window_masker_registry, BlastnArgs::window_masker_counts)";
    let taxid_cases: [(&str, fn(&mut BlastnArgs)); 3] = [
        ("missing both", |_| {}),
        ("missing counts", |a| {
            a.window_masker_registry = Some(Ok(
                crate::blastinput::ncbi_environment::WindowMaskerRegistry::default(),
            ));
        }),
        ("missing registry", |a| {
            a.window_masker_counts = Some(Arc::new(crate::utils::winmask::FileCounts));
        }),
    ];
    for (label, configure) in taxid_cases {
        let mut args = parse(&["-subject", "s", "-window_masker_taxid", "9606"]);
        configure(&mut args);
        assert_refused(args, "6", wm_taxid_msg, label);
    }

    // `-window_masker_taxid 0` disables WindowMasker, so it passes the guard and reaches the
    // empty-subject check (`check_subjects_not_empty`).
    let (result, sink, diagnostics) =
        call(parse(&["-subject", "s", "-window_masker_taxid", "0"]), "6");
    let native = result
        .unwrap_err()
        .downcast::<crate::cli::NativeError>()
        .expect("taxid 0 reaches check_subjects_not_empty");
    assert_eq!(native.exit, 3);
    assert!(String::from_utf8(native.message)
        .unwrap()
        .contains("Empty CBlastQueryVector"));
    assert!(sink.is_empty() && diagnostics.is_empty());

    let rpt_msg = "-filtering_db (a BLAST database of repeats) needs its repeat filter resolved with the run's settings (BlastnArgs::resolve_repeat_filter)";
    assert_refused(
        parse(&["-subject", "s", "-filtering_db", "repeats"]),
        "6",
        rpt_msg,
        "filtering_db",
    );
}

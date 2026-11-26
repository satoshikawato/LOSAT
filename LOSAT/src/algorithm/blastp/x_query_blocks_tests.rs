use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803
// The lookup query with the contexts' Karlin blocks computed on a pool of 1, 3 and 7 threads is
// the reference's (`prepare_blosum62_lookup_query_for_word_size`): the query, the segments and
// every context (Karlin blocks bit for bit, validity). Random queries with low-complexity runs
// (SEG-masked, several segments), a query of X only (no Karlin block: invalid), a short one.
// The same with the blocks computed as extra tasks of the flat table's build
// (LOSAT_X_BPLKPFLAT, `XLookupFlat`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn lookup_query_with_karlin_blocks_on_a_pool_is_the_references() {
    use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_build_flat_in_pool, XFlatBuild,
    };
    use crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts;
    use crate::algorithm::tblastx::lookup::{
        prepare_blosum62_lookup_query_for_word_size, x_prepare_blosum62_lookup_query_in_pool,
    };
    let letters = b"ACDEFGHIKLMNPQRSTVWYBZXUO*";
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut below = |n: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as usize
    };
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    for _ in 0..40 {
        let length = 1 + below(400);
        let mut sequence: Vec<u8> = (0..length).map(|_| letters[below(20)]).collect();
        if length > 60 && below(2) == 0 {
            let at = below(length - 40);
            let run = letters[below(letters.len())];
            sequence[at..at + 30].fill(run);
        }
        if below(10) == 0 {
            sequence[below(length)] = letters[20 + below(6)];
        }
        sequences.push(sequence);
    }
    sequences.push(b"XXXXXXXXXXXX".to_vec());
    sequences.push(b"MK".to_vec());
    let seg = crate::utils::seg::SegParams::default();
    for seg in [None, Some(&seg)] {
        let frames: Vec<Vec<QueryFrame>> = sequences
            .iter()
            .map(|sequence| vec![encode_protein_query_frame_with_seg(sequence, seg)])
            .collect();
        let (query, locations, contexts) =
            prepare_blosum62_lookup_query_for_word_size(&frames, 5, false);
        assert!(contexts.iter().any(|context| !context.is_valid));
        for threads in [1, 3, 7] {
            let (pool_query, pool_locations, pool_contexts) =
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    Ok(x_prepare_blosum62_lookup_query_in_pool(
                        &frames, 5, false, pool,
                    ))
                })
                .expect("pool");
            assert_eq!(pool_query, query);
            assert_eq!(pool_locations, locations);
            x_assert_same_contexts(&pool_contexts, &contexts);
            // LOSAT_X_BPLKPFLAT: the Karlin blocks as extra tasks of the pool's build of the
            // flat table (every one of them at 3 and 7 threads; at one thread the serial builder
            // runs none and `finish` computes them all).
            let from_tasks = crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                let (flat, flat_query, flat_locations, mut flat_contexts) =
                    XLookupFlat::with_karlin_tasks(&frames, 5, XLkpFlat::On, None);
                assert_eq!(flat_query, query);
                assert_eq!(flat_locations, locations);
                let built = x_build_flat_in_pool(
                    5,
                    1100,
                    &flat_query,
                    &flat_locations,
                    &[0..flat_locations.len()],
                    pool,
                    flat.build_args(0),
                );
                assert!(matches!(built, Some(XFlatBuild::Table(_))));
                let from_tasks = flat.finish(&mut flat_contexts);
                x_assert_same_contexts(&flat_contexts, &contexts);
                Ok(from_tasks)
            })
            .expect("pool");
            let expected = if threads > 1 { contexts.len() } else { 0 };
            assert_eq!(from_tasks, expected, "threads {threads}");
        }
    }
}

fn batch(input: Range<usize>, searched: Range<usize>) -> RangedBatch {
    RangedBatch { input, searched }
}

#[test]
fn blocks_close_at_the_first_batch_end_past_the_size() {
    let lengths = [6, 6, 3, 9, 2, 2, 7];
    let batches = [
        batch(0..2, 0..2),
        batch(2..4, 2..4),
        batch(4..6, 4..6),
        batch(6..7, 6..7),
    ];
    // Letters per batch: 12, 12, 4, 7.
    assert_eq!(
        x_query_blocks(&batches, &lengths, 20, None),
        Some(vec![0..4, 4..7])
    );
    assert_eq!(
        x_query_blocks(&batches, &lengths, 12, None),
        Some(vec![0..2, 2..4, 4..7])
    );
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1, None),
        Some(vec![0..2, 2..4, 4..6, 6..7])
    );
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1000, None),
        Some(vec![0..7])
    );
}

#[test]
fn blocks_close_before_a_batch_that_takes_their_words_past_the_most() {
    let lengths = [6, 6, 3, 9, 2, 2, 7];
    let batches = [
        batch(0..2, 0..2),
        batch(2..4, 2..4),
        batch(4..6, 4..6),
        batch(6..7, 6..7),
    ];
    let words = [40u64, 30, 10, 50];
    // The words close a block before the batch that takes it past the most (a block may
    // reach the most).
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1000, Some((&words, 80))),
        Some(vec![0..6, 6..7])
    );
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1000, Some((&words, 70))),
        Some(vec![0..4, 4..7])
    );
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1000, Some((&words, 69))),
        Some(vec![0..2, 2..6, 6..7])
    );
    // A batch above the most is a block of its own; it is never split.
    assert_eq!(
        x_query_blocks(&batches, &lengths, 1000, Some((&words, 1))),
        Some(vec![0..2, 2..4, 4..6, 6..7])
    );
    // Words that never reach the most leave the blocks of the letters.
    assert_eq!(
        x_query_blocks(&batches, &lengths, 20, Some((&words, 1000))),
        Some(vec![0..4, 4..7])
    );
    // A block closed by its letters counts its words anew: the letters close 0..2 and 2..4
    // (12 letters each), and 4..7 holds 10 + 50 words, within 75 (with the words of the
    // blocks before kept, 40 + 30 + 10 would pass 75 and close 4..6).
    assert_eq!(
        x_query_blocks(&batches, &lengths, 12, Some((&words, 75))),
        Some(vec![0..2, 2..4, 4..7])
    );
    // A block closed by its words counts its letters anew: 40 + 30 > 55 closes 0..2, and
    // 2..6 has 16 letters, below 20, when 10 + 50 > 55 closes it.
    assert_eq!(
        x_query_blocks(&batches, &lengths, 20, Some((&words, 55))),
        Some(vec![0..2, 2..6, 6..7])
    );
    // One count per batch, or no blocks.
    assert_eq!(
        x_query_blocks(&batches, &lengths, 20, Some((&words[..3], 1000))),
        None
    );
}

#[test]
fn batches_that_do_not_cover_the_queries_give_no_blocks() {
    let lengths = [5, 5, 5];
    assert_eq!(
        x_query_blocks(&[batch(0..2, 0..2), batch(2..3, 2..2)], &lengths, 1, None),
        None
    );
    assert_eq!(
        x_query_blocks(&[batch(0..2, 0..2)], &lengths, 1, None),
        None
    );
    assert_eq!(
        x_query_blocks(&[batch(0..1, 0..1), batch(1..3, 2..3)], &lengths, 1, None),
        None
    );
    assert_eq!(x_query_blocks(&[], &[], 1, None), Some(vec![]));
}

// The XML checks of the block search agree with the writer's own functions on the probe titles
// of `report/blast_xml/names.rs`: names that pass (periods, `TPA:`, doubled spaces, an empty
// title, a byte >= 0x80 in UTF-8 and Latin-1, a control byte: the stream writes them), and
// punctuation that `x_CleanAndCompress` reads past the end of (refused). A one-record check and the same record among
// others, on a pool of 1 and 3 threads.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn the_xml_checks_agree_with_the_writers_functions() {
    use crate::report::blast_xml::{names, XmlSubjects};
    let titles: [&[u8]; 12] = [
        b"q1 title ending with periods...",
        b"TPA: q3 prefixed",
        b"q4  double  space ,  comma ;",
        b"",
        b"s2 &amp; <y> \"q\"",
        b"tab\there",
        "caf\u{e9}".as_bytes(),
        b"caf\xe9",
        b"a\x01b",
        b", ,",
        b",",
        b"x \xc3\xa9 y",
    ];
    let record = |index: usize, title: &[u8]| {
        FastaRecord::new(
            format!("Query_{}", index + 1),
            title,
            b"MKVLAAGIVALLLAAGCSS",
        )
    };
    let ok_query = |record: &FastaRecord| {
        names::query_texts(
            "BLASTP",
            std::slice::from_ref(record),
            true,
            &[Vec::new()],
            |_| false,
            |_| false,
        )
        .is_ok()
    };
    let ok_subject = |record: &FastaRecord| {
        XmlSubjects::Records {
            records: std::slice::from_ref(record),
            length: &|_| 0,
            protein: true,
        }
        .names(0, "BLASTP")
        .is_ok()
    };
    let records: Vec<FastaRecord> = titles
        .iter()
        .enumerate()
        .map(|(index, title)| record(index, title))
        .collect();
    let passing = records.iter().filter(|record| ok_query(record)).count();
    assert!(
        passing >= 9 && passing < records.len(),
        "the probe titles pass and fail: {passing}"
    );
    for threads in [1, 3] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            for record in &records {
                let one = std::slice::from_ref(record);
                assert_eq!(
                    x_block_xml_run_allowed(one, pool),
                    ok_query(record),
                    "query title {:?}, threads {threads}",
                    String::from_utf8_lossy(&record.title)
                );
                assert_eq!(
                    x_block_xml_subjects_allowed(one, pool),
                    ok_subject(record),
                    "subject title {:?}, threads {threads}",
                    String::from_utf8_lossy(&record.title)
                );
            }
            // One refusable record among many fails the whole set, wherever it lies.
            let many = |bad: Option<usize>| -> Vec<FastaRecord> {
                (0..50)
                    .map(|index| {
                        let title: &[u8] = if Some(index) == bad {
                            b", ,"
                        } else {
                            b"ok title"
                        };
                        record(index, title)
                    })
                    .collect()
            };
            assert!(x_block_xml_run_allowed(&many(None), pool));
            assert!(x_block_xml_subjects_allowed(&many(None), pool));
            for bad in [0, 17, 49] {
                assert!(!x_block_xml_run_allowed(&many(Some(bad)), pool));
                assert!(!x_block_xml_subjects_allowed(&many(Some(bad)), pool));
            }
            Ok(())
        })
        .expect("pool");
    }
}

// The redo's store of SEG-masked subject ranges is the blocks' only for the subjects encoded
// for all blocks (`x_block_seg_store`): its key is a subject's index in the redo's subjects,
// which in another set (another encoding, a `-db` search's listed subjects, whose ranks differ
// between blocks) may name another subject. The `-db` blocks: `tests/psiblast_taxdb_blocks.rs`.
#[test]
fn the_seg_store_of_the_blocks_serves_only_their_subject_encoding() {
    use super::super::encoding::encode_protein_sequence;
    use super::super::x_redo_subjects::XRedoSubjects;
    let records: Vec<FastaRecord> = [&b"MKVLAAGIVALLLAAGCSS"[..], b"MSTNPKPQRKTKRNTNRR"]
        .iter()
        .enumerate()
        .map(|(index, seq)| FastaRecord::new(format!("Subject_{}", index + 1), b"", seq))
        .collect();
    let encode = || {
        records
            .iter()
            .map(|record| encode_protein_sequence(record.seq()))
            .collect::<Vec<_>>()
    };
    let other = encode();
    let subjects = XBlockSubjects::new(&records, encode());
    crate::utils::threading::with_search_pool(1, "blastp", |pool| {
        assert!(x_block_seg_store(pool, XRedoSubjects::Encoded(&subjects.encoded)).is_none());
        let pool = pool.x_with_ext(&subjects);
        let store = x_block_seg_store(&pool, XRedoSubjects::Encoded(&subjects.encoded))
            .expect("the blocks' store");
        assert!(std::sync::Arc::ptr_eq(&store, &subjects.x_seg_store));
        assert!(x_block_seg_store(&pool, XRedoSubjects::Encoded(&other)).is_none());
        assert!(x_block_seg_store(&pool, XRedoSubjects::Encoded(&subjects.encoded[..1])).is_none());
        Ok(())
    })
    .unwrap();
}

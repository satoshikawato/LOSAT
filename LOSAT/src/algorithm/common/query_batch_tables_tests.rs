use super::*;
use crate::algorithm::common::diag_entry::batch_context_distance;

/// Two batches laid out as `q1q.faa` lays them out (fixture `th10` of M/m4b3b/review/probe_r3): `q1` of 60
/// letters, the sentinel, then a second query at 61, here of 240 letters (the fixture's `q` has 48) so
/// that the two batches get tables of different sizes.
const LAYOUT: [(i32, i32); 2] = [(0, 60), (61, 240)];

#[test]
fn each_batch_has_ncbis_ring() {
    // NCBI's table of a batch: the power of two just at least its length plus the window.
    let tables = QueryBatchTables::new(&LAYOUT, 40);
    let masks: Vec<i32> = tables.rings.iter().map(|ring| ring.mask).collect();
    assert_eq!(masks, [127, 511]);
    assert_eq!(QueryBatchTables::new(&[(0, 88)], 40).rings[0].mask, 127);
    assert_eq!(QueryBatchTables::new(&[(0, 89)], 40).rings[0].mask, 255);
    // The one-hit search (window 0) and an empty layout (one batch of nothing).
    assert_eq!(QueryBatchTables::new(&[(0, 64)], 0).rings[0].mask, 63);
    let empty = QueryBatchTables::new(&[], 40);
    assert_eq!(empty.rings[0].mask, 63);
    assert!(empty.is_of(&[], 40));
    assert!(!empty.is_of(&LAYOUT, 40));
    assert!(QueryBatchTables::new(&LAYOUT, 40).is_of(&LAYOUT, 40));
    assert!(!QueryBatchTables::new(&LAYOUT, 40).is_of(&LAYOUT, 0));
}

#[test]
fn a_hit_reads_its_own_batchs_ring_in_its_batchs_offsets() {
    let mut tables = QueryBatchTables::new(&LAYOUT, 40);
    tables.ncbi_table(
        NcbiStart {
            offset: 10,
            cleared: false,
        },
        40,
    );
    assert_eq!(tables.offset, 10);
    // No ring is made before a hit reads it.
    assert_eq!(tables.made(), [0, 0]);
    // Query offset 66 is offset 5 of the second batch: diagonal 5 - 2 of its ring, and the
    // one-hit diagonal 2 - 5.
    let (entry, start) = tables.two_hit_entry(66, 2);
    assert_eq!(start, 61);
    assert_eq!(*entry, DiagStruct::default());
    entry.set_last_hit(1234);
    assert_eq!(tables.made(), [0, 512]);
    assert_eq!(tables.rings[1].entries[3].last_hit(), 1234);
    tables.one_hit_entry(66, 2).set_last_hit(77);
    assert_eq!(tables.rings[1].entries[509].last_hit(), 77);
    // The first batch's hits at query offset 3: its own ring, its start 0.
    let (entry, start) = tables.two_hit_entry(3, 0);
    assert_eq!(start, 0);
    assert_eq!(*entry, DiagStruct::default());
    assert_eq!(tables.made(), [128, 512]);
    // A diagonal of the first query that runs into the second batch's offsets is the second
    // batch's (as in NCBI's search of that batch, where it is a diagonal of its own queries).
    let (_, start) = tables.two_hit_entry(61, 0);
    assert_eq!(start, 61);
}

#[test]
fn every_sequence_starts_on_the_table_of_its_start() {
    let mut tables = QueryBatchTables::new(&LAYOUT, 40);
    tables.ncbi_table(
        NcbiStart {
            offset: 10,
            cleared: false,
        },
        40,
    );
    tables.two_hit_entry(66, 2).0.set_last_hit(1234);
    // The next sequence on the tables: the ring is made again where its first hit reads it, a
    // cleared table's entries the window less than 0.
    tables.ncbi_table(
        NcbiStart {
            offset: 500,
            cleared: true,
        },
        40,
    );
    assert_eq!(tables.offset, 500);
    assert_eq!(tables.made(), [0, 0]);
    let (entry, _) = tables.two_hit_entry(66, 2);
    assert_eq!(*entry, DiagStruct::clear(40));
    assert!(tables.rings[1]
        .entries
        .iter()
        .all(|e| *e == DiagStruct::clear(40)));
}

#[test]
fn the_exit_of_a_chunk_moves_or_clears_every_ring() {
    let mut tables = QueryBatchTables::new(&LAYOUT, 40);
    tables.ncbi_table(
        NcbiStart {
            offset: 10,
            cleared: false,
        },
        40,
    );
    tables.finish_chunk(-25, 40);
    assert_eq!(tables.offset, 25);
    tables.two_hit_entry(66, 2).0.set_last_hit(1234);
    tables.offset = i32::MAX / 4;
    tables.finish_chunk(100, 40);
    assert_eq!(tables.offset, 40);
    assert_eq!(*tables.two_hit_entry(66, 2).0, DiagStruct::clear(40));
    assert_eq!(*tables.two_hit_entry(3, 0).0, DiagStruct::clear(40));
    let Ok(()) = tables.move_offset(7);
    assert_eq!(tables.offset, 47);
}

#[test]
fn fresh_tables_allocate_on_first_hit_and_clear_for_negative_window() {
    let mut tables = QueryBatchTables::new(&LAYOUT, -40);
    assert_eq!(tables.offset, -40);
    assert_eq!(tables.made(), [0, 0]);
    let entry = tables.one_hit_entry(66, 2);
    assert_eq!(*entry, DiagStruct::default());
    entry.set_last_hit(99);
    assert_eq!(tables.made(), [0, 256]);
    tables.clear_table(-40);
    assert_eq!(tables.offset, -40);
    assert_eq!(tables.made(), [0, 0]);
    assert_eq!(*tables.one_hit_entry(66, 2), DiagStruct::clear(-40));
    assert_eq!(*tables.one_hit_entry(3, 0), DiagStruct::clear(-40));
}

#[test]
fn the_context_test_reads_the_batchs_offsets() {
    // A distance whose offset less it falls before the batch's start passes NCBI's test of the
    // batch (its offset less the distance lies below 0): it becomes one past the offset.
    assert_eq!(batch_context_distance(77, 30, 61), 78);
    assert_eq!(batch_context_distance(77, 17, 61), 78);
    // At or after the batch's start, and below 0, the distance is as it was.
    assert_eq!(batch_context_distance(77, 16, 61), 16);
    assert_eq!(batch_context_distance(77, 3, 61), 3);
    assert_eq!(batch_context_distance(20, 30, 61), 30);
    assert_eq!(batch_context_distance(20, 30, 0), 30);
    assert_eq!(batch_context_distance(20, 10, 0), 10);
    for (q, diff, start) in [
        (77_i32, 30_i32, 61_i32),
        (77, 16, 61),
        (20, 30, 61),
        (5, 3, 0),
        (90, 40, 61),
    ] {
        // The unsigned test of the search of all queries on the result equals NCBI's on the
        // batch's offsets, for every context start at or after the batch's.
        for context in [start, start + 10, 85] {
            let all =
                (q.wrapping_sub(batch_context_distance(q, diff, start)) as u32) < context as u32;
            let ncbi = ((q - start).wrapping_sub(diff) as u32) < (context - start) as u32;
            assert_eq!(
                all, ncbi,
                "q {q} diff {diff} batch {start} context {context}"
            );
        }
    }
}

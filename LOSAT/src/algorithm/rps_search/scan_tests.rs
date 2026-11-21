use super::*;
use crate::seqdb::rps::RPS_BACKBONE_SIZE;

fn word(letters: [u8; 3]) -> usize {
    (usize::from(letters[0]) << 10) | (usize::from(letters[1]) << 5) | usize::from(letters[2])
}

struct Table {
    cells: Vec<[i32; 4]>,
    overflow: Vec<i32>,
    pv: Vec<u64>,
}

impl Table {
    fn new() -> Self {
        Self {
            cells: vec![[0; 4]; RPS_BACKBONE_SIZE + 1],
            overflow: Vec::new(),
            pv: vec![0; RPS_BACKBONE_SIZE / 64],
        }
    }

    fn set(&mut self, w: usize, hits: &[i32]) {
        let mut cell = [hits.len() as i32, 0, 0, 0];
        if hits.len() <= 3 {
            cell[1..=hits.len()].copy_from_slice(hits);
        } else {
            cell[1] = hits[0];
            cell[2] = (self.overflow.len() * 4) as i32;
            self.overflow.extend_from_slice(&hits[1..]);
        }
        self.cells[w] = cell;
        self.pv[w / 64] |= 1 << (w % 64);
    }

    fn view(&self) -> RpsTable<'_> {
        RpsTable {
            cells: &self.cells,
            overflow: &self.overflow,
            pv: &self.pv,
        }
    }
}

// Hits go to the bucket of their row in the order of the scan; a cell's overflow run keeps
// its order.
#[test]
fn hits_fill_the_buckets_in_scan_order() {
    let mut table = Table::new();
    table.set(word([1, 2, 3]), &[4100, 10]);
    table.set(word([2, 3, 1]), &[20, 2060, 30, 2050, 5000]);
    let subject = [1u8, 2, 3, 1, 2, 3];
    let mut buckets = RpsBuckets::new(6000, Default::default());
    let mut offset = 0;
    let hits = scan_subject(
        &table.view(),
        &subject,
        &mut offset,
        &mut buckets,
        MAX_SCAN_HITS,
    );
    assert_eq!(hits, 2 + 5 + 2);
    assert_eq!(offset, 4);
    let pairs = |b: usize| -> Vec<(u32, u32)> { buckets.bucket_hits(b) };
    // Word 1 2 3 at s 0 and 3 (rows of the last letter less 2); word 2 3 1 at s 1.
    assert_eq!(pairs(0), vec![(8, 0), (18, 1), (28, 1), (8, 3)]);
    assert_eq!(pairs(1), vec![(2058, 1), (2048, 1)]);
    assert_eq!(pairs(2), vec![(4098, 0), (4998, 1), (4098, 3)]);
}

// A scan stops before the word whose hits do not fit, and the next one starts there.
#[test]
fn a_full_scan_restarts_at_the_word_that_did_not_fit() {
    let mut table = Table::new();
    table.set(word([1, 1, 1]), &[2, 3, 4]);
    let subject = [1u8; 6];
    let mut buckets = RpsBuckets::new(10, Default::default());
    let mut offset = 0;
    // Room for seven hits: two words (six hits), then the third does not fit.
    let first = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
    assert_eq!((first, offset), (6, 2));
    assert_eq!(buckets.bucket_hits(0).len(), 6);
    let second = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
    assert_eq!((second, offset), (6, 4));
    assert_eq!(buckets.bucket_hits(0)[0], (0, 2));
    let last = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
    assert_eq!((last, offset), (0, 4));
}

/// A deterministic sequence of pseudo-random numbers (no NCBI counterpart: test data).
fn lcg(state: &mut u64) -> u32 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    (*state >> 33) as u32
}

// D4 and M1 of PS12s: the blocks of a bucket hold its hits in the order of the scan, and
// `next_block` reads them bucket by bucket, as one `Vec` of NCBI's pairs per bucket holds and
// reads them; a drained bucket's closed blocks go back to the pool, which keeps its blocks
// from scan to scan and volume to volume (fewer or more buckets): what the largest scan's
// buckets held.
// Query offsets cross segments of 2^21 (also inside a bucket's block) and read back whole.
#[test]
fn blocks_read_the_hits_in_the_order_of_vectors_per_bucket() {
    let mut state = 12_345u64;
    // One group of every bucket (the hits are added to it directly).
    let mut buckets = RpsBuckets::with_group(0, Default::default(), usize::MAX / 2);
    let mut largest = 0usize;
    let segment = 1u32 << SEGMENT_SHIFT;
    for (scan, &(rows, hits, first_s)) in [
        (5_000usize, 10usize, 0u32),
        (300_000, 9_000, 0),
        (300_000, 0, 0),
        (40_000, 3 * BLOCK_ENTRIES + 1, segment - 2_000),
        (600_000, 25_000, 0),
        (2_047, 2 * BLOCK_ENTRIES, 3 * segment - 700),
        (2_047, 5, u32::MAX - 4),
        (70_000, 6_000, 5 * segment + 12),
    ]
    .iter()
    .enumerate()
    {
        buckets.reset_for(rows);
        let cursor = buckets.open_group(0);
        let count = rows / RPS_BUCKET_SIZE + 1;
        assert_eq!(buckets.heads.len(), count, "scan {scan}");
        let mut expected: Vec<Vec<(u32, u32)>> = vec![Vec::new(); count];
        for i in 0..hits as u32 {
            // Words of one to three hits at increasing query offsets.
            let s_off = first_s.wrapping_add(i / 2);
            let q_off = lcg(&mut state) % rows.max(1) as u32;
            let s_bits = buckets.entry_bits(s_off);
            buckets.add(q_off, s_bits);
            expected[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s_off));
        }
        for (b, hits) in expected.iter().enumerate() {
            assert_eq!(&buckets.bucket_hits(b), hits, "scan {scan} bucket {b}");
        }
        // Every closed block is full but where a segment started after it.
        let crossed = first_s >> SEGMENT_SHIFT
            != first_s.wrapping_add((hits as u32).saturating_sub(1) / 2) >> SEGMENT_SHIFT;
        for chain in &buckets.chains {
            for closed in chain {
                assert!(
                    crossed || closed.len as usize == BLOCK_ENTRIES,
                    "scan {scan}"
                );
            }
        }
        let held = buckets.heads.len() + buckets.chains.iter().map(Vec::len).sum::<usize>();
        largest = largest.max(held);
        let read = buckets.read_group(cursor);
        assert_eq!(read, expected.concat(), "scan {scan}");
        assert!(buckets.chains.iter().all(Vec::is_empty), "scan {scan}");
        assert!(
            buckets.heads.iter().all(|head| head.fill == 0),
            "scan {scan}"
        );
        // The pool holds what the buckets held at most at one time.
        assert_eq!(buckets.blocks(), largest, "scan {scan}");
        assert_eq!(
            buckets.kept_bytes(),
            buckets.blocks() * BLOCK_ENTRIES * std::mem::size_of::<u32>()
        );
    }
}

// M1 of PS12s: a scan whose query offsets cross a segment of 2^21 gives every hit back whole,
// bucket by bucket in the order of the scan, as NCBI's pairs, whatever the groups (PS12m; the
// cells list their hits by bucket, as a table made by makeprofiledb does). Ten buckets and
// groups of 1 to 4: a word whose hits lie in a later group is skipped in the earlier ones,
// also after the segment changed.
#[test]
fn a_scan_across_a_segment_reads_the_query_offsets_whole() {
    let cells: [([u8; 3], &[i32]); 3] = [
        ([1, 2, 3], &[10, 4100, 12_000]),
        ([2, 3, 1], &[20, 30, 2050, 2060, 5000]),
        ([3, 1, 2], &[18_000, 18_100]),
    ];
    let mut table = Table::new();
    for (letters, entries) in cells {
        table.set(word(letters), entries);
    }
    let boundary = 1usize << SEGMENT_SHIFT;
    let mut subject = vec![20u8; boundary + 12];
    for start in [
        boundary - 7,
        boundary - 4,
        boundary - 1,
        boundary + 2,
        boundary + 6,
    ] {
        subject[start..start + 3].copy_from_slice(&[1, 2, 3]);
        subject[start + 3] = 1;
    }
    subject[boundary - 6..boundary - 3].copy_from_slice(&[3, 1, 2]);
    let mut expected: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 10];
    for s in 0..=subject.len() - 3 {
        let w = (usize::from(subject[s]) << 10)
            | (usize::from(subject[s + 1]) << 5)
            | usize::from(subject[s + 2]);
        let entries = cells
            .iter()
            .find(|(letters, _)| word(*letters) == w)
            .map_or(&[][..], |(_, entries)| entries);
        for &entry in entries {
            let q_off = (entry - 2) as u32;
            expected[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s as u32));
        }
    }
    assert!(expected[0].iter().any(|&(_, s)| s as usize >= boundary));
    assert!(expected[8].iter().any(|&(_, s)| s as usize >= boundary - 7));
    for group_buckets in [RPS_GROUP_BUCKETS, 1, 2, 3, 4] {
        let mut buckets = RpsBuckets::with_group(20_000, Default::default(), group_buckets);
        let mut offset = 0;
        let hits = scan_subject(
            &table.view(),
            &subject,
            &mut offset,
            &mut buckets,
            MAX_SCAN_HITS,
        );
        assert_eq!(hits as usize, expected.iter().map(Vec::len).sum::<usize>());
        let read = buckets.drain(&table.view(), &subject);
        assert_eq!(read, expected.concat(), "{group_buckets} buckets a group");
    }
}

// M2 of PS12s: the buckets of two workers share one pool: the blocks one worker's buckets
// give back once read serve the other's scan, so the run holds the blocks its buckets hold
// at the same time.
#[test]
fn workers_share_the_blocks_given_back_once_read() {
    let pool: std::sync::Arc<RpsBlockPool> = Default::default();
    let mut a = RpsBuckets::new(5_000, pool.clone());
    let mut b = RpsBuckets::new(5_000, pool.clone());
    let fill = |buckets: &mut RpsBuckets| {
        let cursor = buckets.open_group(0);
        for s_off in 0..5 * BLOCK_ENTRIES as u32 {
            let s_bits = buckets.entry_bits(s_off);
            buckets.add(s_off % 3_000, s_bits);
        }
        cursor
    };
    let closed = |buckets: &RpsBuckets| buckets.chains.iter().map(Vec::len).sum::<usize>();
    let cursor = fill(&mut a);
    let closed_a = closed(&a);
    assert!(closed_a > 0);
    assert_eq!(pool.len(), 0);
    assert_eq!(a.read_group(cursor).len(), 5 * BLOCK_ENTRIES);
    assert_eq!(pool.len(), closed_a);
    // The other worker's scan takes those blocks: no block more.
    let cursor = fill(&mut b);
    assert_eq!(closed(&b), closed_a);
    assert_eq!(pool.len(), 0);
    assert_eq!(b.read_group(cursor).len(), 5 * BLOCK_ENTRIES);
    assert_eq!(pool.len(), closed_a);
}

/// NCBI's model of one scan (`BlastRPSScanSubject` with `s_AddToRPSBucket`): for each bucket
/// the pairs appended in order, the scan's total and the offset the next scan starts at. `lists`
/// are the cells' hits by word (empty: no hit), in the order of the cell.
fn reference_scan(
    lists: &[Vec<i32>],
    subject: &[u8],
    offset: &mut i32,
    max_hits: i32,
    num_buckets: usize,
) -> (i32, Vec<Vec<(u32, u32)>>) {
    let mut buckets = vec![Vec::new(); num_buckets];
    let mut total = 0i32;
    let mut s = *offset;
    while s as usize + 3 <= subject.len() {
        let at = s as usize;
        let w = word([subject[at], subject[at + 1], subject[at + 2]]);
        let list = &lists[w];
        if !list.is_empty() {
            if list.len() as i32 > max_hits - total {
                break;
            }
            for &entry in list {
                let q_off = (entry - 2) as u32;
                buckets[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s as u32));
            }
            total += list.len() as i32;
        }
        s += 1;
    }
    *offset = s;
    (total, buckets)
}

/// A random table of `rows` rows over the letters 1 to 4: a hit list for about half of the
/// words, of 1, 2, 3, 4 (the form boundary), 50 or (twice at most) 6,000 entries in ascending
/// order, with equal buckets and repeated entries too.
fn random_lists(state: &mut u64, rows: usize) -> Vec<Vec<i32>> {
    let mut lists = vec![Vec::new(); crate::seqdb::rps::RPS_BACKBONE_SIZE];
    let mut huge = 0;
    for a in 1..=4u8 {
        for b in 1..=4u8 {
            for c in 1..=4u8 {
                if lcg(state).is_multiple_of(2) {
                    continue;
                }
                let mut len = [1usize, 2, 3, 4, 50, 3, 4, 2, 6_000][lcg(state) as usize % 9];
                if len == 6_000 {
                    huge += 1;
                    if huge > 2 {
                        len = 50;
                    }
                }
                let mut list: Vec<i32> = (0..len)
                    .map(|_| 2 + (lcg(state) as usize % (rows - 2)) as i32)
                    .collect();
                list.sort_unstable();
                lists[word([a, b, c])] = list;
            }
        }
    }
    lists
}

fn table_of(lists: &[Vec<i32>]) -> Table {
    let mut table = Table::new();
    for (w, list) in lists.iter().enumerate() {
        if !list.is_empty() {
            table.set(w, list);
        }
    }
    table
}

/// Scans `subject` completely (scan after scan) into buckets of groups of `group_buckets` and
/// asserts, for each scan, NCBI's total, next offset and hits per bucket in order.
fn assert_scans_equal_ncbi(
    table: &RpsTable<'_>,
    lists: &[Vec<i32>],
    subject: &[u8],
    rows: usize,
    max_hits: i32,
    group_buckets: usize,
) -> usize {
    let num_buckets = rows / RPS_BUCKET_SIZE + 1;
    let mut buckets = RpsBuckets::with_group(rows, Default::default(), group_buckets);
    let (mut offset, mut ref_offset, mut scans) = (0i32, 0i32, 0);
    while offset as usize + 3 <= subject.len() {
        let total = scan_subject(table, subject, &mut offset, &mut buckets, max_hits);
        let (ref_total, expected) =
            reference_scan(lists, subject, &mut ref_offset, max_hits, num_buckets);
        assert_eq!((total, offset), (ref_total, ref_offset));
        if group_buckets >= num_buckets {
            // One group: the scan made every hit, as the scan of dev did, and no pass follows.
            assert!(!buckets.beyond);
            assert_eq!(buckets.groups(), 1);
        }
        assert_eq!(
            buckets.drain(table, subject),
            expected.concat(),
            "scan {scans} of {group_buckets} buckets a group"
        );
        scans += 1;
    }
    scans
}

// PS12m (second review, S2): a scan whose every word has one hit past the first group keeps
// nothing per word (a cell without an overflow run): the later group is made by a pass over the subject, and
// the hits are the scan's, in its order (the dev scan's state is the same hits).
#[test]
fn a_repeated_word_with_one_hit_past_the_first_group_needs_no_state_per_word() {
    let mut table = Table::new();
    let q_off = 100 * RPS_BUCKET_SIZE + 77; // bucket 100: group 1
    table.set(word([1, 1, 1]), &[q_off as i32 + 2]);
    let subject = vec![1u8; 20_000];
    let mut buckets = RpsBuckets::new(300_000, Default::default());
    let mut offset = 0;
    let hits = scan_subject(
        &table.view(),
        &subject,
        &mut offset,
        &mut buckets,
        MAX_SCAN_HITS,
    );
    assert_eq!(hits, 20_000 - 2);
    assert_eq!(buckets.groups(), 3);
    // The word's one entry lies in group 1: group 2 has no pass.
    assert!(buckets.needs_pass(0) && buckets.needs_pass(1) && !buckets.needs_pass(2));
    let read = buckets.drain(&table.view(), &subject);
    let expected: Vec<(u32, u32)> = (0..20_000 - 2).map(|s| (q_off as u32, s)).collect();
    assert_eq!(read, expected);
}

// PS12m: whatever the group of buckets (1, 2, 3, 7, 64, or all), every scan gives the word
// finder NCBI's hits in NCBI's order: the same total and next offset, and the concatenation
// of the buckets' pairs. Random tables (cells of 1, 2, 3, 4, 50 and 6,000 hits, up to 342
// buckets), random subjects, and `max_hits` from the longest cell to a few times it, so
// that a subject takes several scans and a scan ends inside a run of heavy words.
#[test]
fn groups_of_buckets_read_the_hits_in_ncbis_order() {
    let mut state = 7u64;
    let mut cut_scans = 0;
    // 131,071 and 131,074 rows: volumes of exactly 64 and of 65 buckets, the production
    // group's edge (one group, and a second group of one bucket).
    for case in 0..28 {
        let rows = [
            1_000usize, 5_000, 40_000, 140_000, 700_000, 131_071, 131_074,
        ][case % 7];
        let lists = random_lists(&mut state, rows);
        let table = table_of(&lists);
        let longest = lists.iter().map(Vec::len).max().unwrap() as i32;
        let length = 20 + lcg(&mut state) as usize % 1_480;
        let subject: Vec<u8> = (0..length)
            .map(|_| 1 + (lcg(&mut state) % 4) as u8)
            .collect();
        let max_hits = longest + (lcg(&mut state) % 40_000) as i32;
        let buckets = rows / RPS_BUCKET_SIZE + 1;
        for group_buckets in [1, 2, 3, 7, RPS_GROUP_BUCKETS, buckets, buckets + 5] {
            let scans = assert_scans_equal_ncbi(
                &table.view(),
                &lists,
                &subject,
                rows,
                max_hits,
                group_buckets,
            );
            cut_scans += scans.saturating_sub(1);
        }
    }
    assert!(
        cut_scans > 20,
        "{cut_scans} scans ended before the subject did"
    );
}

// PS12m: cells whose hits are not in groups of buckets in order are regrouped when the volume
// is opened (`seqdb::rps`); read from the regrouped cells, every group gives the hits that
// NCBI's scan of the cells as they are gives (the model is built from the original, shuffled
// lists). The regrouping is by the group size read (production: `RPS_GROUP_BUCKETS`).
#[test]
fn regrouped_cells_read_the_hits_in_the_order_of_the_original_cells() {
    let mut state = 99u64;
    let mut regrouped_cases = 0;
    for case in 0..10 {
        let rows = [5_000usize, 40_000, 140_000, 700_000, 300_000][case % 5];
        let mut lists = random_lists(&mut state, rows);
        for list in lists.iter_mut().filter(|list| list.len() > 1) {
            // A shuffle (Fisher-Yates), so that most long lists lose their order.
            for i in (1..list.len()).rev() {
                list.swap(i, lcg(&mut state) as usize % (i + 1));
            }
        }
        let table = table_of(&lists);
        let length = 20 + lcg(&mut state) as usize % 1_000;
        let subject: Vec<u8> = (0..length)
            .map(|_| 1 + (lcg(&mut state) % 4) as u8)
            .collect();
        let longest = lists.iter().map(Vec::len).max().unwrap() as i32;
        let buckets = rows / RPS_BUCKET_SIZE + 1;
        for group_buckets in [1, 3, RPS_GROUP_BUCKETS, buckets] {
            let regrouped =
                crate::seqdb::rps::regroup_cells(&table.cells, &table.overflow, group_buckets)
                    .unwrap();
            if group_buckets < buckets {
                assert!(
                    regrouped.is_some(),
                    "shuffled lists in {group_buckets} of {buckets} buckets"
                );
                regrouped_cases += 1;
            }
            let view = match &regrouped {
                Some((cells, overflow)) => RpsTable {
                    cells,
                    overflow,
                    pv: &table.pv,
                },
                None => table.view(),
            };
            assert_scans_equal_ncbi(
                &view,
                &lists,
                &subject,
                rows,
                longest + 3_000,
                group_buckets,
            );
        }
    }
    assert!(regrouped_cases >= 10, "{regrouped_cases} regrouped cases");
}

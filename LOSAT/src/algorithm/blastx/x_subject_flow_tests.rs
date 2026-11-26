use super::{offsets_stay_below_2_30, XTasks};
use crate::algorithm::blastx::{input::FastaRecord, subjects::Subjects};

fn records(lengths: &[usize]) -> Vec<FastaRecord> {
    lengths
        .iter()
        .enumerate()
        .map(|(i, &length)| FastaRecord {
            internal_id: format!("s{i}"),
            title: Vec::new(),
            sequence: vec![b'A'; length],
            lowercase_masks: Vec::new(),
            source: None,
            around: None,
            skipped: None,
        })
        .collect()
}

// The cut covers the subjects in order, no task holds more than `k` subjects, and a task of
// two or more subjects holds at most `k` times the mean length; equal lengths keep the grid.
#[test]
fn the_tasks_cut_by_residues_cover_the_subjects_in_order() {
    let mut lengths = vec![300; 1_000];
    for length in &mut lengths[400..420] {
        *length = 20_000;
    }
    lengths[999] = 7;
    for k in [1, 2, 7, 64, 2_000] {
        let set = records(&lengths);
        let tasks = XTasks::cut(k, Subjects::Fasta(&set)).unwrap();
        let budget = (k * lengths.iter().sum::<usize>()).div_ceil(lengths.len());
        let mut next = 0;
        for task in 0..tasks.count() {
            let range = tasks.range(task);
            assert_eq!(range.start, next, "k {k} task {task}");
            assert!(!range.is_empty() && range.len() <= k, "k {k} task {task}");
            let residues: usize = lengths[range.clone()].iter().sum();
            assert!(range.len() == 1 || residues <= budget, "k {k} task {task}");
            next = range.end;
        }
        assert_eq!(next, lengths.len(), "k {k}");
        if k == 64 {
            // The twenty long subjects (subjects 400-419 of the grid's task 384..448): two per
            // task under the budget of 64 x the mean length (693.707), 26 tasks for the grid's
            // 16.
            assert_eq!(budget, 44_398);
            assert_eq!(tasks.count(), 26);
            assert_eq!(tasks.range(6), 384..401);
            assert_eq!(tasks.range(7), 401..403);
            assert_eq!(tasks.range(16), 419..483);
        }
    }
    let set = records(&[300; 1_000]);
    let tasks = XTasks::cut(64, Subjects::Fasta(&set)).unwrap();
    assert_eq!(tasks.count(), 16);
    assert_eq!(tasks.range(15), 960..1_000);
    assert!(XTasks::cut(64, Subjects::Fasta(&[])).is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-175
// The offset a subject is scanned at is below INT4_MAX / 4 plus the previous subject and a
// window; its entries add at most the subject's own length.
#[test]
fn the_flow_is_kept_only_where_no_offset_passes_2_30() {
    assert!(offsets_stay_below_2_30(40, 5_000_000));
    assert!(offsets_stay_below_2_30((1 << 28) - 2 * 35_000, 35_000));
    assert!(!offsets_stay_below_2_30((1 << 28) - 2 * 35_000 + 1, 35_000));
    assert!(!offsets_stay_below_2_30(536_870_000, 4_910));
    assert!(!offsets_stay_below_2_30(-1, 10));
    // The bound holds for the worst history: INT4_MAX / 4 - 1 plus two subjects and a
    // window stays below 2^30.
    let worst = i64::from(i32::MAX / 4) - 1 + (1 << 28);
    assert!(worst < 1 << 30);
}

use super::{DiagStruct, Diagonals};

/// An entry with `last_hit` and the flag set.
pub(super) fn flagged(last_hit: i32) -> DiagStruct {
    let mut entry = DiagStruct::default();
    entry.set_last_hit(last_hit);
    entry.set_flag(1);
    entry
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
#[test]
fn allocation_power_of_two_and_window_boundaries() {
    for (query, window, length) in [
        (63, 0, 64),
        (64, 0, 64),
        (65, 0, 128),
        (24, 40, 64),
        (25, 40, 128),
    ] {
        let d = Diagonals::new(query, window).unwrap();
        assert_eq!(d.entries.len(), length);
        assert_eq!(d.offset, window);
    }
    assert!(Diagonals::new(1, -1).is_err());
}

// No NCBI counterpart (buffer reuse, LOSAT_X_BXDIAGKEEP): a used table of another length or
// window, made again by `x_remake`, is the table `new` makes.
#[cfg(feature = "parallel")]
#[test]
fn a_remade_table_is_the_new_one() {
    for (from, to) in [
        ((100, 40), (30, 40)),
        ((30, 40), (100, 40)),
        ((70, 40), (70, 8)),
    ] {
        let mut table = Diagonals::new(from.0, from.1).unwrap();
        table.entries.fill(flagged(12_345));
        table.offset = 99;
        assert!(!table.x_fits(to.0, to.1));
        table.x_remake(to.0, to.1).unwrap();
        let new = Diagonals::new(to.0, to.1).unwrap();
        assert!(table.x_fits(to.0, to.1));
        assert_eq!(table.entries, new.entries);
        assert_eq!((table.offset, table.window), (new.offset, new.window));
    }
    assert!(Diagonals::new(8, 8).unwrap().x_remake(1, -1).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
#[test]
fn reset_occurs_at_entry_threshold_and_clears_stale_flags() {
    for window in [0, 1, 39, 40, 41] {
        let mut d = Diagonals::new(65, window).unwrap();
        d.entries.fill(flagged(123));
        d.offset = 536_870_910; // INT4_MAX / 4 - 1, independent literal boundary.
        d.finish_subject(1).unwrap();
        assert_eq!(d.offset, 536_870_911 + window);
        assert!(d.entries.iter().all(|x| *x == flagged(123)));
        d.finish_subject(5_000_000).unwrap();
        assert_eq!(d.offset, window);
        assert!(d.entries.iter().all(|x| *x == DiagStruct::clear(window)));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// Independent worker slots advance only for their assigned subjects, including
// a partial final wave; no query state is shared between slots or batches.
#[test]
fn uneven_wave_slots_and_new_batch_keep_separate_offsets() {
    let mut slots = (0..4)
        .map(|_| Diagonals::new(128, 40).unwrap())
        .collect::<Vec<_>>();
    for (oid, length) in [10, 20, 30, 40, 50, 60, 70, 80, 90].into_iter().enumerate() {
        slots[oid % 4].finish_subject(length).unwrap();
    }
    assert_eq!(
        slots.iter().map(|x| x.offset).collect::<Vec<_>>(),
        [310, 200, 220, 240]
    );
    assert_eq!(Diagonals::new(128, 40).unwrap().offset, 40);
}
